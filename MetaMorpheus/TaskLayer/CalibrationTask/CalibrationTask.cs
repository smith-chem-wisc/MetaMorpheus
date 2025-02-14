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
        private static readonly int NumRequiredPsms = 16;
        private static readonly int NumRequiredMs1Datapoints = 40;
        private static readonly int NumRequiredMs2Datapoints = 80;
        private static readonly double InitialSearchToleranceMultiplier = 1.5;
        private static readonly double PrecursorMultiplierForToml = 3;
        private static readonly double ProductMultiplierForToml = 6;
        
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
            // re-write experimental design (if it has been defined) with new calibrated file names
            string assumedPathToExperDesign = Directory.GetParent(currentRawFileList.First()).FullName;
            assumedPathToExperDesign = Path.Combine(assumedPathToExperDesign, GlobalVariables.ExperimentalDesignFileName);

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                // get filename stuff
                string originalUncalibratedFilePath = currentRawFileList[spectraFileIndex];
                string originalUncalibratedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUncalibratedFilePath);
                string calibratedFilePath = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".mzML");
                string uncalibratedFilePath = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + ".mzML");

                // mark the file as in-progress
                StartingDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                
                FileSpecificParameters fileSpecificParams = new();
                // carry over file-specific parameters from the uncalibrated file to the calibrated one
                if (fileSettingsList[spectraFileIndex] != null)
                {
                    fileSpecificParams = fileSettingsList[spectraFileIndex].Clone();
                }
                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSpecificParams);
                string uncalibratedTomlFilename = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + ".toml");
                // load the file
                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files" });
                MsDataFile myMsDataFile = myFileManager.LoadFile(originalUncalibratedFilePath, combinedParams);

                // get datapoints to fit calibration function to
                Status("Acquiring calibration data points...", new List<string> { taskId, "Individual Spectra Files" });
                DataPointAquisitionResults acquisitionResultsFirst = null;

                acquisitionResultsFirst = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams,
                    new PpmTolerance(combinedParams.ProductMassTolerance.Value * InitialSearchToleranceMultiplier), new PpmTolerance(combinedParams.PrecursorMassTolerance.Value * InitialSearchToleranceMultiplier));
                // check if we have enough data points to calibrate and then calibrate
                if (FileSuitableForCalibration(acquisitionResultsFirst))
                {
                    // write toml settings for the calibrated file
                    string calibratedTomlFilename = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".toml");

                    // set the mass tolerances for the file specific parameters
                    // we use a multiplier of 4 for the tolerance for files that are not calibrated
                    fileSpecificParams.PrecursorMassTolerance = new PpmTolerance(Math.Round((PrecursorMultiplierForToml * acquisitionResultsFirst.PsmPrecursorIqrPpmError) + Math.Abs(acquisitionResultsFirst.PsmPrecursorMedianPpmError),1));
                    fileSpecificParams.ProductMassTolerance = new PpmTolerance(Math.Round((ProductMultiplierForToml * acquisitionResultsFirst.PsmProductIqrPpmError) + Math.Abs(acquisitionResultsFirst.PsmProductMedianPpmError),1));

                    // generate calibration function and shift data points
                    Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                    CalibrationEngine engine = new(myMsDataFile, acquisitionResultsFirst, combinedParams, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                    _ = engine.Run();

                    // get the calibrated data points again to see if there was an increase
                    DataPointAquisitionResults acquisitionResultsSecond = GetDataAcquisitionResults(engine.CalibratedDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, 
                        fileSpecificParams.PrecursorMassTolerance, fileSpecificParams.ProductMassTolerance);

                    myMsDataFile = engine.CalibratedDataFile;
                    acquisitionResultsFirst = acquisitionResultsSecond;
                    fileSpecificParams.PrecursorMassTolerance = new PpmTolerance(Math.Round((PrecursorMultiplierForToml * acquisitionResultsFirst.PsmPrecursorIqrPpmError) + Math.Abs(acquisitionResultsFirst.PsmPrecursorMedianPpmError),1));
                    fileSpecificParams.ProductMassTolerance = new PpmTolerance(Math.Round((ProductMultiplierForToml * acquisitionResultsFirst.PsmProductIqrPpmError) + Math.Abs(acquisitionResultsFirst.PsmProductMedianPpmError), 1));
                    //Try a second round of calibration
                    // generate calibration function and shift data points
                    Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                    engine = new(myMsDataFile, acquisitionResultsFirst, combinedParams, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                    _ = engine.Run();

                    // get the calibrated data points again to see if there was an increase
                    acquisitionResultsSecond = GetDataAcquisitionResults(engine.CalibratedDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, fileSpecificParams.PrecursorMassTolerance, fileSpecificParams.ProductMassTolerance);

                    if (CalibrationHasValue(acquisitionResultsFirst, acquisitionResultsSecond))
                    {
                        myMsDataFile = engine.CalibratedDataFile;
                        // write toml settings for the calibrated file
                        fileSpecificParams.PrecursorMassTolerance = new PpmTolerance(Math.Round((PrecursorMultiplierForToml * acquisitionResultsSecond.PsmPrecursorIqrPpmError) + Math.Abs(acquisitionResultsSecond.PsmPrecursorMedianPpmError), 1));
                        fileSpecificParams.ProductMassTolerance = new PpmTolerance(Math.Round((ProductMultiplierForToml * acquisitionResultsSecond.PsmProductIqrPpmError) + Math.Abs(acquisitionResultsSecond.PsmProductMedianPpmError), 1));
                    }
                    CalibrationOutput(myMsDataFile, calibratedFilePath, fileSpecificParams, calibratedTomlFilename, taskId, originalUncalibratedFilenameWithoutExtension);
                    myFileManager.DoneWithFile(originalUncalibratedFilePath);

                    // finished calibrating this file

                    FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                    ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
                    continue;
                }
                // if we didn't calibrate, write the uncalibrated file to the output folder
                CalibrationOutput(myMsDataFile, uncalibratedFilePath, fileSpecificParams, uncalibratedTomlFilename, taskId, originalUncalibratedFilenameWithoutExtension);
                unsuccessfullyCalibratedFilePaths.Add(originalUncalibratedFilePath);
                // provide a message indicating why we couldn't calibrate
                CalibrationWarnMessage(acquisitionResultsFirst);
                FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
                continue;
            }
            if (File.Exists(assumedPathToExperDesign))
            {
                WriteNewExperimentalDesignFile(assumedPathToExperDesign, OutputFolder, currentRawFileList, unsuccessfullyCalibratedFilePaths);
            }
            // finished calibrating all files for the task
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            return MyTaskResults;
        }
        private bool CalibrationHasValue(DataPointAquisitionResults acquisitionResultsFirst, DataPointAquisitionResults acquisitionResultsSecond)
        {
            bool improvedCounts = acquisitionResultsSecond.Psms.Count > acquisitionResultsFirst.Psms.Count && acquisitionResultsSecond.Psms.Select(p => p.FullSequence).Distinct().Count() > acquisitionResultsFirst.Psms.Select(p => p.FullSequence).Distinct().Count();
            if(improvedCounts)
            {
                return true;
            }
            bool numPsmsIncreased = acquisitionResultsSecond.Psms.Count >= acquisitionResultsFirst.Psms.Count;
            bool numPeptidesIncreased = acquisitionResultsSecond.Psms.Select(p => p.FullSequence).Distinct().Count() >= acquisitionResultsFirst.Psms.Select(p => p.FullSequence).Distinct().Count();
            bool psmPrecursorMedianPpmErrorDecreased = Math.Abs(acquisitionResultsSecond.PsmPrecursorMedianPpmError) <= 1 || Math.Abs(acquisitionResultsSecond.PsmPrecursorMedianPpmError) <= Math.Abs(acquisitionResultsFirst.PsmPrecursorMedianPpmError);
            bool psmProductMedianPpmErrorDecreased = Math.Abs(acquisitionResultsSecond.PsmProductMedianPpmError) <= 1 || Math.Abs(acquisitionResultsSecond.PsmProductMedianPpmError) <= Math.Abs(acquisitionResultsFirst.PsmProductMedianPpmError);

            return (numPsmsIncreased && numPeptidesIncreased && psmPrecursorMedianPpmErrorDecreased && psmProductMedianPpmErrorDecreased);
        }
        private bool FileSuitableForCalibration(DataPointAquisitionResults acquisitionResults)
        {
            return acquisitionResults.Psms.Count >= NumRequiredPsms && acquisitionResults.Ms1List.Count >= NumRequiredMs1Datapoints && acquisitionResults.Ms2List.Count >= NumRequiredMs2Datapoints;
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

        private void CalibrationOutput(MsDataFile msDataFile, string mzFilePath, FileSpecificParameters fileParams, string tomlName, string taskId, string mzFilenameNoExtension)
        {
            msDataFile.ExportAsMzML(mzFilePath, CalibrationParameters.WriteIndexedMzml);
            MyTaskResults.NewSpectra.Add(mzFilePath);
            Toml.WriteFile(fileParams, tomlName, tomlConfig);
            FinishedWritingFile(tomlName, new List<string> { taskId, "Individual Spectra Files", mzFilenameNoExtension });
            FinishedWritingFile(mzFilePath, new List<string> { taskId, "Individual Spectra Files", mzFilenameNoExtension });
            MyTaskResults.NewFileSpecificTomls.Add(tomlName);
        }

        private DataPointAquisitionResults GetDataAcquisitionResults(MsDataFile myMsDataFile, string currentDataFile, List<Modification> variableModifications, List<Modification> fixedModifications, List<Protein> proteinList, string taskId, CommonParameters combinedParameters, Tolerance initPrecTol, Tolerance initProdTol)
        {
            string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(currentDataFile);
            MassDiffAcceptor searchMode = initPrecTol is PpmTolerance ?
                (MassDiffAcceptor)new SinglePpmAroundZeroSearchMode(initPrecTol.Value) :
                new SingleAbsoluteAroundZeroSearchMode(initPrecTol.Value);

            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = GetMs2Scans(myMsDataFile, currentDataFile, combinedParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            Log("Searching with searchMode: " + searchMode, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });
            Log("Searching with productMassTolerance: " + initProdTol, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });

            // overwrite specific properties in the common parameters that are different for this iteration in the calibraiton task
            bool writeSpectralLibrary = false;
            combinedParameters.PrecursorMassTolerance = initPrecTol;
            combinedParameters.ProductMassTolerance = initProdTol;
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