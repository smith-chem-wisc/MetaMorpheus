using EngineLayer;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Modifications;
using Readers;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class CalibrationTask : MetaMorpheusTask
    {
        public CalibrationTask() : base(MyTask.Calibrate)
        {
            CommonParameters = new CommonParameters(
                trimMsMsPeaks: false,
                doPrecursorDeconvolution: false
                );

            CalibrationParameters = new CalibrationParameters();
        }

        public CalibrationParameters CalibrationParameters { get; set; }
        private readonly int NumRequiredPsms = 20;
        private readonly int NumRequiredMs1Datapoints = 50;
        private readonly int NumRequiredMs2Datapoints = 100;
        double MaxPrecursorTolerance = 40;
        double MaxProductTolerance = 150;
        public const string CalibSuffix = "-calib";
        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            LoadModifications(taskId, out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, true, DecoyType.Reverse, localizeableModificationTypes, CommonParameters);

            // write prose settings
            _ = ProseCreatedWhileRunning.Append("The following calibration settings were used: ");
            _ = ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            _ = ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            _ = ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            _ = ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            _ = ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            _ = ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif)) + "; ");
            _ = ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");
            _ = ProseCreatedWhileRunning.Append("max mods per peptide = " + CommonParameters.DigestionParams.MaxModsForPeptide + "; ");
            _ = ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            _ = ProseCreatedWhileRunning.Append("precursor mass tolerance = " + CommonParameters.PrecursorMassTolerance + "; ");
            _ = ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + ". ");
            _ = ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count(p => !p.IsDecoy) + " non-decoy protein entries including " + proteinList.Count(p => p.IsContaminant) + " contaminant sequences. ");

            // start the calibration task
            Status("Calibrating...", new List<string> { taskId });
            MyTaskResults = new MyTaskResults(this)
            {
                NewSpectra = new List<string>(),
                NewFileSpecificTomls = new List<string>()
            };

            MyFileManager myFileManager = new MyFileManager(true);
            List<string> unsuccessfullyCalibratedFilePaths = new List<string>();

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                // get filename stuff
                string originalUncalibratedFilePath = currentRawFileList[spectraFileIndex];
                string originalUncalibratedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUncalibratedFilePath);
                string calibratedFilePath = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".mzML");

                // mark the file as in-progress
                StartingDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                // load the file
                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files" });

                MsDataFile myMsDataFile = myFileManager.LoadFile(originalUncalibratedFilePath, combinedParams);

                // get datapoints to fit calibration function to
                Status("Acquiring calibration data points...", new List<string> { taskId, "Individual Spectra Files" });

                // We first run get acquisition results with a high q-value cutoff to get a rough estimate of the precursor and product mass errors and the q-value cuttoff for calibrating PSMs
                // We then use these errors to set the precursor and product tolerances for the next round of acquisition results
                // We then run get acquisition results again with the new tolerances to get the final calibration data points
                DataPointAquisitionResults acquisitionResults = null;
                double precursorTolerance = MaxPrecursorTolerance;
                double productTolerance = MaxProductTolerance;
                CalibrationParameters.QValueCutoffForCalibratingPSMs = 0.05;
                // get acquisition results with a subset of the data if possible to set the initital tolerances
                acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, new PpmTolerance(precursorTolerance), new PpmTolerance(150), 10);
                if (acquisitionResults.Psms.Count(p=>!p.IsDecoy) > NumRequiredPsms && acquisitionResults.Ms1List.Count > NumRequiredMs1Datapoints && acquisitionResults.Ms2List.Count > NumRequiredMs2Datapoints)
                {
                    //set initial tolerances and q-value cutoff for GetAcquisitionResults
                    var qValueNotchForMaxRequiredPsm = acquisitionResults.Psms.Where(p => !p.IsDecoy).Select(p=>p.FdrInfo.QValueNotch).ToList()[NumRequiredPsms];
                    CalibrationParameters.QValueCutoffForCalibratingPSMs = (0 < qValueNotchForMaxRequiredPsm && qValueNotchForMaxRequiredPsm < MaxPrecursorTolerance) ? qValueNotchForMaxRequiredPsm : 0.05;
                    var potentialNewPrecursorTolerance = acquisitionResults.PsmPrecursorIqrPpmError * 5 + Math.Abs(acquisitionResults.PsmPrecursorMedianPpmError);
                    var potentialNewProductTolerance = acquisitionResults.PsmProductIqrPpmError * 5 + Math.Abs(acquisitionResults.PsmProductMedianPpmError);
                    precursorTolerance = (0 < potentialNewPrecursorTolerance && potentialNewPrecursorTolerance < MaxPrecursorTolerance) ? potentialNewPrecursorTolerance : MaxPrecursorTolerance;
                    productTolerance = (0 < potentialNewProductTolerance && potentialNewProductTolerance < MaxProductTolerance) ? potentialNewProductTolerance : MaxProductTolerance;
                }

                // This is the first set of acquisition results that will be used to calibrate the file
                acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, new PpmTolerance(precursorTolerance), new PpmTolerance(productTolerance));
                if (acquisitionResults.Psms.Count(p => !p.IsDecoy) > NumRequiredPsms && acquisitionResults.Ms1List.Count > NumRequiredMs1Datapoints && acquisitionResults.Ms2List.Count > NumRequiredMs2Datapoints)
                {
                    //set next tolerances and q-value cutoff for GetAcquisitionResults
                    if(acquisitionResults.Psms.Count(p => !p.IsDecoy && p.FdrInfo.QValueNotch > 0.005) > 1000) // it is okay to reduce the threshold for the q-value cutoff if there are a lot of PSMs
                    {
                        CalibrationParameters.QValueCutoffForCalibratingPSMs = 0.005;
                    }
                    else if (acquisitionResults.Psms.Count(p => !p.IsDecoy && p.FdrInfo.QValueNotch > 0.01) > 1000) // we didn't have a ton of PSMs at 0.005, so let's try 0.01
                    {
                        CalibrationParameters.QValueCutoffForCalibratingPSMs = 0.01;
                    }
                    var potentialNewPrecursorTolerance = acquisitionResults.PsmPrecursorIqrPpmError * 5 + Math.Abs(acquisitionResults.PsmPrecursorMedianPpmError);
                    var potentialNewProductTolerance = acquisitionResults.PsmProductIqrPpmError * 5 + Math.Abs(acquisitionResults.PsmProductMedianPpmError);
                    precursorTolerance = (0 < potentialNewPrecursorTolerance && potentialNewPrecursorTolerance < MaxPrecursorTolerance) ? potentialNewPrecursorTolerance : MaxPrecursorTolerance;
                    productTolerance = (0 < potentialNewProductTolerance && potentialNewProductTolerance < MaxProductTolerance) ? potentialNewProductTolerance : MaxProductTolerance;

                    // generate calibration function and shift data points
                    Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                    CalibrationEngine engine = new(myMsDataFile, acquisitionResults, combinedParams, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                    _ = engine.Run();

                    //update file
                    myMsDataFile = engine.CalibratedDataFile;

                    // use the tolereances and cuttoffs defined above to get acquisition results that will be used to calibrate the file
                    // using 5x the interquartile range of the ppm errors as the tolerance to have more data points for calibration
                    acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, new PpmTolerance(precursorTolerance), new PpmTolerance(productTolerance));

                    if (acquisitionResults.Psms.Count(p => !p.IsDecoy) > NumRequiredPsms && acquisitionResults.Ms1List.Count > NumRequiredMs1Datapoints && acquisitionResults.Ms2List.Count > NumRequiredMs2Datapoints)
                    {

                        //generate calibration function and shift data points AGAIN because it's fast and contributes new data
                        Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                        engine = new CalibrationEngine(myMsDataFile, acquisitionResults, combinedParams, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                        _ = engine.Run();

                        //update file
                        // repeat GetAcquisitionResults which will provide the ultimate precursor and product tolerances for the toml
                        myMsDataFile = engine.CalibratedDataFile;

                        // write the calibrated mzML file
                        myMsDataFile.ExportAsMzML(calibratedFilePath, CalibrationParameters.WriteIndexedMzml);
                        myFileManager.DoneWithFile(originalUncalibratedFilePath);

                        // write toml settings for the calibrated file
                        string newTomlFileName = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".toml");

                        FileSpecificParameters fileSpecificParams = new();

                        // carry over file-specific parameters from the uncalibrated file to the calibrated one
                        if (fileSettingsList[spectraFileIndex] != null)
                        {
                            fileSpecificParams = fileSettingsList[spectraFileIndex].Clone();
                        }

                        //suggest 4 * interquartile range as the ppm tolerance
                        fileSpecificParams.PrecursorMassTolerance = new PpmTolerance((4.0 * acquisitionResults.PsmPrecursorIqrPpmError) + Math.Abs(acquisitionResults.PsmPrecursorMedianPpmError));
                        fileSpecificParams.ProductMassTolerance = new PpmTolerance((4.0 * acquisitionResults.PsmProductIqrPpmError) + Math.Abs(acquisitionResults.PsmProductMedianPpmError));

                        Toml.WriteFile(fileSpecificParams, newTomlFileName, tomlConfig);

                        FinishedWritingFile(newTomlFileName, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });

                        // finished calibrating this file
                        FinishedWritingFile(calibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                        MyTaskResults.NewSpectra.Add(calibratedFilePath);
                        MyTaskResults.NewFileSpecificTomls.Add(newTomlFileName);
                        FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                        ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
                    }

                }
                else // give up
                {
                    // provide a message indicating why we couldn't calibrate
                    CalibrationWarnMessage(acquisitionResults);
                    FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                    ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
                    continue;
                }
            }

            // re-write experimental design (if it has been defined) with new calibrated file names
            string assumedPathToExperDesign = Directory.GetParent(currentRawFileList.First()).FullName;
            assumedPathToExperDesign = Path.Combine(assumedPathToExperDesign, GlobalVariables.ExperimentalDesignFileName);

            if (File.Exists(assumedPathToExperDesign))
            {
                WriteNewExperimentalDesignFile(assumedPathToExperDesign, OutputFolder, currentRawFileList, unsuccessfullyCalibratedFilePaths);
            }

            // finished calibrating all files for the task
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            return MyTaskResults;
        }



        private void CalibrationWarnMessage(DataPointAquisitionResults acquisitionResults)
        {
            // provide a message indicating why we couldn't calibrate
            if (acquisitionResults.Psms.Count < NumRequiredPsms)
            {
                Warn("Calibration failure! Could not find enough high-quality PSMs. Required " + NumRequiredPsms + ", saw " + acquisitionResults.Psms.Count);
            }
            else if (acquisitionResults.Ms1List.Count < NumRequiredMs1Datapoints)
            {
                Warn("Calibration failure! Could not find enough MS1 datapoints. Required " + NumRequiredMs1Datapoints + ", saw " + acquisitionResults.Ms1List.Count);
            }
            else if (acquisitionResults.Ms2List.Count < NumRequiredMs2Datapoints)
            {
                Warn("Calibration failure! Could not find enough MS2 datapoints. Required " + NumRequiredMs2Datapoints + ", saw " + acquisitionResults.Ms2List.Count);
            }
        }

        private DataPointAquisitionResults GetDataAcquisitionResults(MsDataFile myMsDataFile, string currentDataFile, List<Modification> variableModifications, List<Modification> fixedModifications, List<Protein> proteinList, string taskId, CommonParameters combinedParameters, Tolerance initPrecTol, Tolerance initProdTol, int ms2Fraction = 1)
        {
            string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(currentDataFile);
            MassDiffAcceptor searchMode = initPrecTol is PpmTolerance ?
                (MassDiffAcceptor)new SinglePpmAroundZeroSearchMode(initPrecTol.Value) :
                new SingleAbsoluteAroundZeroSearchMode(initPrecTol.Value);

            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = GetMs2Scans(myMsDataFile, currentDataFile, combinedParameters).OrderBy(b => b.PrecursorMass).Where((x, i) => i % ms2Fraction == 0).ToArray();
            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            Log("Searching with searchMode: " + searchMode, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });
            Log("Searching with productMassTolerance: " + initProdTol, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });

            bool writeSpectralLibrary = false;
            _ = new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null, proteinList, searchMode, combinedParameters,
                FileSpecificParameters, null, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }, writeSpectralLibrary).Run();
            List<SpectralMatch> allPsms = allPsmsArray.Where(b => b != null).ToList();

            allPsms = allPsms.OrderByDescending(b => b.Score)
                .ThenBy(b => b.BioPolymerWithSetModsMonoisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.BioPolymerWithSetModsMonoisotopicMass.Value) : double.MaxValue)
                .GroupBy(b => (b.FullFilePath, b.ScanNumber, b.BioPolymerWithSetModsMonoisotopicMass)).Select(b => b.First()).ToList();

            _ = new FdrAnalysisEngine(allPsms, searchMode.NumNotches, CommonParameters, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }, doPEP: false).Run();

            List<SpectralMatch> goodIdentifications = allPsms.Where(b => 
                b.FdrInfo.QValueNotch < CalibrationParameters.QValueCutoffForCalibratingPSMs 
                && b.FullSequence != null
                && !b.IsDecoy).ToList();

            if (!goodIdentifications.Any())
            {
                return new DataPointAquisitionResults(null, new List<SpectralMatch>(), new List<LabeledDataPoint>(), new List<LabeledDataPoint>(), 0, 0, 0, 0);
            }

            //get the deconvoluted ms2scans for the good identifications
            List<Ms2ScanWithSpecificMass> goodScans = new List<Ms2ScanWithSpecificMass>();
            List<SpectralMatch> unfilteredPsms = allPsmsArray.ToList();
            foreach (SpectralMatch psm in goodIdentifications)
            {
                goodScans.Add(listOfSortedms2Scans[unfilteredPsms.IndexOf(psm)]);
            }

            DataPointAquisitionResults currentResult = (DataPointAquisitionResults)new DataPointAcquisitionEngine(
                    goodIdentifications,
                    goodScans,
                    myMsDataFile,
                    initPrecTol,
                    initProdTol,
                    CalibrationParameters.MinMS1IsotopicPeaksNeededForConfirmedIdentification,
                    CommonParameters,
                    FileSpecificParameters,
                    new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            return currentResult;
        }

        private static void WriteNewExperimentalDesignFile(string pathToOldExperDesign, string outputFolder, List<string> originalUncalibratedFileNamesWithExtension,
            List<string> unsuccessfullyCalibratedFilePaths)
        {
            List<SpectraFileInfo> oldExperDesign = ExperimentalDesign.ReadExperimentalDesign(pathToOldExperDesign, originalUncalibratedFileNamesWithExtension, out var errors);

            if (errors.Any())
            {
                foreach (string error in errors)
                {
                    Warn(error);
                }

                return;
            }

            List<SpectraFileInfo> newExperDesign = new();

            foreach (SpectraFileInfo uncalibratedSpectraFile in oldExperDesign)
            {
                string originalUncalibratedFilePath = uncalibratedSpectraFile.FullFilePathWithExtension;
                string originalUncalibratedFilenameWithoutExtension = GlobalVariables.GetFilenameWithoutExtension(originalUncalibratedFilePath);
                string calibratedFilePath = Path.Combine(outputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".mzML");

                SpectraFileInfo calibratedSpectraFile = new(calibratedFilePath,
                    uncalibratedSpectraFile.Condition, uncalibratedSpectraFile.BiologicalReplicate, uncalibratedSpectraFile.TechnicalReplicate, uncalibratedSpectraFile.Fraction);

                if (unsuccessfullyCalibratedFilePaths.Contains(uncalibratedSpectraFile.FullFilePathWithExtension))
                {
                    newExperDesign.Add(uncalibratedSpectraFile);
                }
                else
                {
                    newExperDesign.Add(calibratedSpectraFile);
                }
            }

            _ = ExperimentalDesign.WriteExperimentalDesignToFile(newExperDesign);
        }
    }
}