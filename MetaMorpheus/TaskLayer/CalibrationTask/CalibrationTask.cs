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
                productMassTolerance: new PpmTolerance(25),
                precursorMassTolerance: new PpmTolerance(15),
                trimMsMsPeaks: false,
                doPrecursorDeconvolution: false
                );

            CalibrationParameters = new CalibrationParameters();
        }

        public CalibrationParameters CalibrationParameters { get; set; }

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

                bool couldNotFindEnoughDatapoints = false;

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
                DataPointAquisitionResults aR = null;
                CalibrationParameters.QValueCutoffForCalibratingPSMs = 0.05;
                aR = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, new PpmTolerance(40), new PpmTolerance(150), 10);

                if (aR.Psms.Count(p=>!p.IsDecoy) > NumRequiredPsms && aR.Ms1List.Count < NumRequiredMs1Datapoints && aR.Ms2List.Count < NumRequiredMs2Datapoints)
                {
                    var qValueNotchForMaxRequiredPsm = aR.Psms.Where(p => !p.IsDecoy).Select(p=>p.FdrInfo.QValueNotch).ToList()[NumRequiredPsms];
                    CalibrationParameters.QValueCutoffForCalibratingPSMs = Math.Max(0.005, qValueNotchForMaxRequiredPsm);
                }
                else
                {
                    aR = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, new PpmTolerance(40), new PpmTolerance(150));
                    if (aR.Psms.Count(p => !p.IsDecoy) > NumRequiredPsms && aR.Ms1List.Count < NumRequiredMs1Datapoints && aR.Ms2List.Count < NumRequiredMs2Datapoints)
                    {
                        var qValueNotchForMaxRequiredPsm = aR.Psms.Where(p => !p.IsDecoy).Select(p => p.FdrInfo.QValueNotch).ToList()[NumRequiredPsms];
                        CalibrationParameters.QValueCutoffForCalibratingPSMs = Math.Max(0.005, qValueNotchForMaxRequiredPsm);
                    }
                    else if (aR.Psms.Count(p => !p.IsDecoy) > 0 && aR.Ms1List.Count < 0 && aR.Ms2List.Count < 0)
                    {
                        var qValueNotchForMaxRequiredPsm = aR.Psms.Where(p => !p.IsDecoy).Select(p=>p.FdrInfo.QValueNotch).Max();
                        CalibrationParameters.QValueCutoffForCalibratingPSMs = Math.Max(0.005, qValueNotchForMaxRequiredPsm);
                    }
                    else
                    {
                        couldNotFindEnoughDatapoints = true;
                    }
                }

                if(couldNotFindEnoughDatapoints == true)
                {
                    if (aR.Psms.Count < NumRequiredPsms)
                    {
                        Warn("Calibration failure! Could not find enough high-quality PSMs. Required " + NumRequiredPsms + ", saw " + aR.Psms.Count);
                        FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                        continue;
                    }
                    else if (aR.Ms1List.Count < NumRequiredMs1Datapoints)
                    {
                        Warn("Calibration failure! Could not find enough MS1 datapoints. Required " + NumRequiredMs1Datapoints + ", saw " + aR.Ms1List.Count);
                        FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                        continue;
                    }
                    else if (aR.Ms2List.Count < NumRequiredMs2Datapoints)
                    {
                        Warn("Calibration failure! Could not find enough MS2 datapoints. Required " + NumRequiredMs2Datapoints + ", saw " + aR.Ms2List.Count);
                        FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                        continue;
                    }
                }

                DataPointAquisitionResults acquisitionResults = null;
                double precursorTolerance = Math.Min(40, aR.PsmPrecursorIqrPpmError * 5 + Math.Abs(aR.PsmPrecursorMedianPpmError));
                double productTolerance = Math.Min(40, aR.PsmProductIqrPpmError * 5 + Math.Abs(aR.PsmProductMedianPpmError));
                acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, new PpmTolerance(precursorTolerance), new PpmTolerance(productTolerance));

                // stats before calibration
                int prevPsmCount = acquisitionResults.Psms.Count;
                double preCalibrationPrecursorErrorIqr = acquisitionResults.PsmPrecursorIqrPpmError;
                double preCalibrationProductErrorIqr = acquisitionResults.PsmProductIqrPpmError;

                // generate calibration function and shift data points
                Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                CalibrationEngine engine = new(myMsDataFile, acquisitionResults, combinedParams, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                _ = engine.Run();

                //update file
                myMsDataFile = engine.CalibratedDataFile;

                precursorTolerance = Math.Min(40, acquisitionResults.PsmPrecursorIqrPpmError * 5 + Math.Abs(acquisitionResults.PsmPrecursorMedianPpmError));
                productTolerance = Math.Min(40, acquisitionResults.PsmProductIqrPpmError * 5 + Math.Abs(acquisitionResults.PsmProductMedianPpmError));
                acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, new PpmTolerance(precursorTolerance), new PpmTolerance(productTolerance));

                //generate calibration function and shift data points AGAIN because it's fast and contributes new data
                Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                engine = new CalibrationEngine(myMsDataFile, acquisitionResults, combinedParams, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                _ = engine.Run();

                //update file
                myMsDataFile = engine.CalibratedDataFile;

                precursorTolerance = Math.Min(40, acquisitionResults.PsmPrecursorIqrPpmError * 5 + Math.Abs(acquisitionResults.PsmPrecursorMedianPpmError));
                productTolerance = Math.Min(40, acquisitionResults.PsmProductIqrPpmError * 5 + Math.Abs(acquisitionResults.PsmProductMedianPpmError));
                acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, new PpmTolerance(precursorTolerance), new PpmTolerance(productTolerance));

                // write the calibrated mzML file
                myMsDataFile.ExportAsMzML(calibratedFilePath, CalibrationParameters.WriteIndexedMzml);
                myFileManager.DoneWithFile(originalUncalibratedFilePath);

                // stats after calibration
                int postCalibrationPsmCount = acquisitionResults.Psms.Count;
                double postCalibrationPrecursorErrorIqr = acquisitionResults.PsmPrecursorIqrPpmError;
                double postCalibrationProductErrorIqr = acquisitionResults.PsmProductIqrPpmError;

                //// did the data improve? (not used for anything yet...)
                //bool improvement = ImprovGlobal(preCalibrationPrecursorErrorIqr, preCalibrationProductErrorIqr, prevPsmCount, postCalibrationPsmCount, postCalibrationPrecursorErrorIqr, postCalibrationProductErrorIqr);

                // write toml settings for the calibrated file
                string newTomlFileName = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".toml");

                FileSpecificParameters fileSpecificParams = new();

                // carry over file-specific parameters from the uncalibrated file to the calibrated one
                if (fileSettingsList[spectraFileIndex] != null)
                {
                    fileSpecificParams = fileSettingsList[spectraFileIndex].Clone();
                }

                //suggest 4 * interquartile range as the ppm tolerance
                fileSpecificParams.PrecursorMassTolerance = new PpmTolerance((4.0 * postCalibrationPrecursorErrorIqr) + Math.Abs(acquisitionResults.PsmPrecursorMedianPpmError));
                fileSpecificParams.ProductMassTolerance = new PpmTolerance((4.0 * postCalibrationProductErrorIqr) + Math.Abs(acquisitionResults.PsmProductMedianPpmError));

                Toml.WriteFile(fileSpecificParams, newTomlFileName, tomlConfig);

                FinishedWritingFile(newTomlFileName, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });

                // finished calibrating this file
                FinishedWritingFile(calibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                MyTaskResults.NewSpectra.Add(calibratedFilePath);
                MyTaskResults.NewFileSpecificTomls.Add(newTomlFileName);
                FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
                
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

        private readonly int NumRequiredPsms = 20;
        private readonly int NumRequiredMs1Datapoints = 50;
        private readonly int NumRequiredMs2Datapoints = 100;
        public const string CalibSuffix = "-calib";

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