using EngineLayer;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Omics;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class CalibrationTask : MetaMorpheusTask
    {
        public CalibrationTask() : base(MyTask.Calibrate)
        {
            CommonParameters = new CommonParameters(
                trimMsMsPeaks: false,
                doPrecursorDeconvolution: false,
                productMassTolerance: new PpmTolerance(InitialProductTolerance),
                precursorMassTolerance: new PpmTolerance(InitialPrecursorTolerance)
                );

            CalibrationParameters = new CalibrationParameters();
        }

        public CalibrationParameters CalibrationParameters { get; set; }
        private static readonly int NumRequiredPsms = 16;
        private static readonly int NumRequiredMs1Datapoints = 40;
        private static readonly int NumRequiredMs2Datapoints = 80;
        private static readonly double InitialSearchToleranceMultiplier = 2.0;
        private static readonly double PrecursorMultiplierForToml = 3;
        private static readonly double ProductMultiplierForToml = 6;
        private static readonly double InitialPrecursorTolerance = 10;
        private static readonly double InitialProductTolerance = 30;

        public const string CalibSuffix = "-calib";

        private List<string> _unsuccessfullyCalibratedFilePaths;
        private string _taskId;
        private List<IBioPolymer> _proteinList;
        private List<Modification> _variableModifications;
        private List<Modification> _fixedModifications;
        private MyFileManager _myFileManager;

        protected override MyTaskResults RunSpecific(string outputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            Initialize(taskId, dbFilenameList);
            Status("Calibrating...", new List<string> { _taskId });

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                // Can't calibrate .mgf, .d, or .msalign files
                if (!CanFileTypeBeCalibrated(currentRawFileList[spectraFileIndex]))
                {
                    continue;
                }

                // get original file name, and file names for potential outputs
                string originalUncalibratedFilePath = currentRawFileList[spectraFileIndex];
                string uncalibratedNewFullFilePath = Path.Combine(outputFolder, Path.GetFileName(currentRawFileList[spectraFileIndex]));
                string originalUncalibratedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUncalibratedFilePath);
                string calibratedNewFullFilePath = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".mzML")
                    .ToSafeOutputPath(CalibSuffix + ".mzML");

                // mark the file as in-progress
                StartingDataFile(originalUncalibratedFilePath, new List<string> { _taskId, "Individual Spectra Files", originalUncalibratedFilePath });

                // carry over file-specific parameters from the uncalibrated file to the calibrated one and update combined params
                FileSpecificParameters fileSpecificParams = fileSettingsList[spectraFileIndex] == null 
                    ? new() 
                    : fileSettingsList[spectraFileIndex].Clone();
                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSpecificParams);
                
                // load the file
                Status("Loading spectra file...", new List<string> { _taskId, "Individual Spectra Files" });
                MsDataFile myMsDataFile = _myFileManager.LoadFile(originalUncalibratedFilePath, combinedParams).LoadAllStaticData();

                // First round of calibration
                Status("Acquiring calibration data points...", new List<string> { _taskId, "Individual Spectra Files" });
                DataPointAquisitionResults acquisitionResultsFirst = GetDataAcquisitionResults(myMsDataFile, combinedParams, originalUncalibratedFilePath);

                //not enough points on the first go so try again with a little wider tolerance
                if (!SufficientAcquisitionResults(acquisitionResultsFirst))
                {
                    UpdateCombinedParameters(combinedParams,
                        combinedParams.PrecursorMassTolerance.Value * InitialSearchToleranceMultiplier,
                        combinedParams.ProductMassTolerance.Value * InitialSearchToleranceMultiplier);
                    WarnForWiderTolerance(combinedParams.PrecursorMassTolerance.Value, combinedParams.ProductMassTolerance.Value);
                    acquisitionResultsFirst = GetDataAcquisitionResults(myMsDataFile, combinedParams, originalUncalibratedFilePath);
                }
                // If there still aren't enough points, give up
                if(!SufficientAcquisitionResults(acquisitionResultsFirst))
                {
                    WriteUncalibratedFile(originalUncalibratedFilePath, uncalibratedNewFullFilePath, _unsuccessfullyCalibratedFilePaths, acquisitionResultsFirst, _taskId);
                    continue;
                }

                UpdateCombinedParameters(combinedParams, acquisitionResultsFirst);

                // generate calibration function and shift data points
                Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                CalibrationEngine engine = new(myMsDataFile, acquisitionResultsFirst, combinedParams, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                _ = engine.Run();

                // Second round of calibration
                DataPointAquisitionResults acquisitionResultsSecond = GetDataAcquisitionResults(engine.CalibratedDataFile, combinedParams, originalUncalibratedFilePath);

                // If the second acquisition results are worse, then calibration made things worse. So we should give up 
                // and write the uncalibrated file
                if (!SufficientAcquisitionResults(acquisitionResultsSecond))
                {
                    WriteUncalibratedFile(originalUncalibratedFilePath, uncalibratedNewFullFilePath, _unsuccessfullyCalibratedFilePaths, acquisitionResultsFirst, taskId);
                    continue;
                }

                //only continue with the third round if the second round showed improvement
                if (CalibrationHasValue(acquisitionResultsFirst, acquisitionResultsSecond))
                {
                    myMsDataFile = engine.CalibratedDataFile;
                    UpdateCombinedParameters(combinedParams, acquisitionResultsSecond);

                    // generate calibration function and shift data points
                    Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                    engine = new(myMsDataFile, acquisitionResultsSecond, combinedParams, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                    _ = engine.Run();

                    // Third round of calibration
                    DataPointAquisitionResults acquisitionResultsThird = GetDataAcquisitionResults(engine.CalibratedDataFile, combinedParams, originalUncalibratedFilePath);

                    if (CalibrationHasValue(acquisitionResultsSecond, acquisitionResultsThird))
                    {
                        myMsDataFile = engine.CalibratedDataFile;
                        UpdateCombinedParameters(combinedParams, acquisitionResultsThird);
                    }
                }
                // if calibration did not make things better, write uncalibrated file and update the tolerances
                else
                {
                    WriteUncalibratedFile(originalUncalibratedFilePath, uncalibratedNewFullFilePath, _unsuccessfullyCalibratedFilePaths, acquisitionResultsFirst, taskId);
                    continue;
                }

                // Update file specific params to reflect the new tolerances, then write them out
                fileSpecificParams.PrecursorMassTolerance = combinedParams.PrecursorMassTolerance;
                fileSpecificParams.ProductMassTolerance = combinedParams.ProductMassTolerance;

                // write toml settings for the calibrated file
                string calibratedTomlFilename = Path.Combine(outputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".toml");
                CalibrationOutput(myMsDataFile, calibratedNewFullFilePath, fileSpecificParams, calibratedTomlFilename, taskId, originalUncalibratedFilenameWithoutExtension);

                // finished calibrating this file
                _myFileManager.DoneWithFile(originalUncalibratedFilePath);
                FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
            }

            // re-write experimental design (if it has been defined) with new calibrated file names
            UpdateExperimentalDesignFile(currentRawFileList, outputFolder);

            // finished calibrating all files for the task
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            return MyTaskResults;
        }

        private void UpdateExperimentalDesignFile(List<string> currentRawFileList, string outputFolder)
        {
            string assumedPathToExperDesign = Directory.GetParent(currentRawFileList.First()).FullName;
            assumedPathToExperDesign = Path.Combine(assumedPathToExperDesign, GlobalVariables.ExperimentalDesignFileName);
            if (File.Exists(assumedPathToExperDesign))
            {
                WriteNewExperimentalDesignFile(assumedPathToExperDesign, outputFolder, currentRawFileList, _unsuccessfullyCalibratedFilePaths);
            }
        }

        private DataPointAquisitionResults GetDataAcquisitionResults(MsDataFile myMsDataFile, CommonParameters combinedParameters, string originalDataFile)
        {
            string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(originalDataFile);
            MassDiffAcceptor searchMode = combinedParameters.PrecursorMassTolerance is PpmTolerance ?
                new SinglePpmAroundZeroSearchMode(combinedParameters.PrecursorMassTolerance.Value) :
                new SingleAbsoluteAroundZeroSearchMode(combinedParameters.PrecursorMassTolerance.Value);

            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = GetMs2Scans(myMsDataFile, originalDataFile, combinedParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] allPsmsArray = new SpectralMatch[listOfSortedms2Scans.Length];

            Log("Searching with searchMode: " + searchMode, new List<string> { _taskId, "Individual Spectra Files", fileNameWithoutExtension });
            Log("Searching with precursorMassTolerance: " + combinedParameters.PrecursorMassTolerance, new List<string> { _taskId, "Individual Spectra Files", fileNameWithoutExtension });
            Log("Searching with productMassTolerance: " + combinedParameters.ProductMassTolerance, new List<string> { _taskId, "Individual Spectra Files", fileNameWithoutExtension });

            _ = new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, _variableModifications, _fixedModifications, null, null, null, _proteinList, searchMode, combinedParameters,
                FileSpecificParameters, null, new List<string> { _taskId, "Individual Spectra Files", fileNameWithoutExtension }, false).Run();

            List<SpectralMatch> allPsms = allPsmsArray.Where(b => b != null).OrderByDescending(b => b.Score)
                .ThenBy(b => b.BioPolymerWithSetModsMonoisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.BioPolymerWithSetModsMonoisotopicMass.Value) : double.MaxValue)
                .GroupBy(b => (b.FullFilePath, b.ScanNumber, b.BioPolymerWithSetModsMonoisotopicMass)).Select(b => b.First()).ToList();

            _ = new FdrAnalysisEngine(allPsms, searchMode.NumNotches, CommonParameters, FileSpecificParameters, new List<string> { _taskId, "Individual Spectra Files", fileNameWithoutExtension }, doPEP: false).Run();

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
                    combinedParameters.PrecursorMassTolerance,
                    combinedParameters.ProductMassTolerance,
                    CalibrationParameters.MinMS1IsotopicPeaksNeededForConfirmedIdentification,
                    CommonParameters,
                    FileSpecificParameters,
                    new List<string> { _taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            return currentResult;
        }

        /// <summary>
        /// Writes prose settings and initializes the following private fields used by the calibration engine:
        /// _taskId, _variableModifications, _fixedModifications, _proteinList, _myFileManager, _unsuccessfullyCalibratedFilePaths
        /// </summary>
        private void Initialize(string taskId, List<DbForTask> dbFilenameList)
        {
            MyTaskResults = new MyTaskResults(this)
            {
                NewSpectra = new List<string>(),
                NewFileSpecificTomls = new List<string>()
            };
            _taskId = taskId;
            LoadModifications(_taskId, out _variableModifications, out _fixedModifications, out var localizeableModificationTypes);
            // load proteins
            var dbLoader = new DatabaseLoadingEngine(CommonParameters, this.FileSpecificParameters, [taskId], dbFilenameList, taskId, DecoyType.Reverse, true, localizeableModificationTypes);
            var loadingResults = dbLoader.Run() as DatabaseLoadingEngineResults;
            _proteinList = loadingResults!.BioPolymers;
            
            _myFileManager = new MyFileManager(true);
            _unsuccessfullyCalibratedFilePaths = new List<string>();

            // write prose settings
            WriteProse(_fixedModifications, _variableModifications, _proteinList);
        }

        public void WriteProse(List<Modification> fixedModifications, List<Modification> variableModifications, List<IBioPolymer> bioPolymerList)
        {
            // write prose settings
            _ = ProseCreatedWhileRunning.Append("The following calibration settings were used: ");
            _ = ProseCreatedWhileRunning.Append($"{GlobalVariables.AnalyteType.GetDigestionAgentLabel()} = " + CommonParameters.DigestionParams.DigestionAgent + "; ");
            _ = ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            _ = ProseCreatedWhileRunning.Append($"minimum {GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower()} length = " + CommonParameters.DigestionParams.MinLength + "; ");
            _ = ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxLength == int.MaxValue ?
                $"maximum {GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower()} length = unspecified; " :
                $"maximum {GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower()} length = " + CommonParameters.DigestionParams.MaxLength + "; ");
            if (CommonParameters.DigestionParams is DigestionParams digestionParams)
                _ = ProseCreatedWhileRunning.Append("initiator methionine behavior = " + digestionParams.InitiatorMethionineBehavior + "; ");
            _ = ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif)) + "; ");
            _ = ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");
            _ = ProseCreatedWhileRunning.Append($"max mods per {GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower()} = " + CommonParameters.DigestionParams.MaxMods + "; ");
            _ = ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            _ = ProseCreatedWhileRunning.Append("precursor mass tolerance = " + CommonParameters.PrecursorMassTolerance + "; ");
            _ = ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + ". ");
            _ = ProseCreatedWhileRunning.Append("The combined search database contained " + bioPolymerList.Count(p => !p.IsDecoy) + $" non-decoy {GlobalVariables.AnalyteType.GetBioPolymerLabel().ToLower()} entries including " + bioPolymerList.Count(p => p.IsContaminant) + " contaminant sequences. ");
        }

        public bool CanFileTypeBeCalibrated(string originalUncalibratedFilePath)
        {
            string fileExtension = Path.GetExtension(originalUncalibratedFilePath);
            string originalUncalibratedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUncalibratedFilePath);
            if (fileExtension.Equals(".mgf", StringComparison.OrdinalIgnoreCase) || fileExtension.Equals(".d", StringComparison.OrdinalIgnoreCase) || fileExtension.Equals(".msalign", StringComparison.OrdinalIgnoreCase))
            {
                _unsuccessfullyCalibratedFilePaths.Add(originalUncalibratedFilePath);
                // provide a message indicating why we couldn't calibrate
                Warn("Calibration for " + fileExtension + " files is not supported.");
                FinishedDataFile(originalUncalibratedFilePath, new List<string> { _taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { _taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
                return false;
            }
            return true;
        }

        /// <summary>
        /// Update the combined params such that the precursor and product mass tolerances are equal to
        /// (Precursor = 3/Product = 6)MultiplierForToml * acquisitionResults.PsmPrecursorIqrPpmError + |acquisitionResults.PsmPrecursorMedianPpmError|
        /// 
        /// Currently, this coerces the tolerance values to ppm tolerance. In the future, we should add support for absolute tolerance
        /// </summary>
        public static void UpdateCombinedParameters(CommonParameters combinedParams, DataPointAquisitionResults acquisitionResults)
        {
            double newPrecursorPpmTolerance = Math.Round(PrecursorMultiplierForToml * acquisitionResults.PsmPrecursorIqrPpmError + Math.Abs(acquisitionResults.PsmPrecursorMedianPpmError), 1);
            double newProductPpmTolerance = Math.Round(ProductMultiplierForToml * acquisitionResults.PsmProductIqrPpmError + Math.Abs(acquisitionResults.PsmProductMedianPpmError), 1);
            UpdateCombinedParameters(combinedParams, newPrecursorPpmTolerance, newProductPpmTolerance);
        }

        public static void UpdateCombinedParameters(CommonParameters combinedParameters, double newPrecursorTolerance, double newProductTolerance)
        {
            combinedParameters.PrecursorMassTolerance = new PpmTolerance(newPrecursorTolerance);
            combinedParameters.ProductMassTolerance = new PpmTolerance(newProductTolerance);
        }

        private void WriteUncalibratedFile(string originalUncalibratedFilePath, string uncalibratedNewFullFilePath, List<string> unsuccessfullyCalibratedFilePaths,
            DataPointAquisitionResults acquisitionResults, string taskId)
        {
            // if we didn't calibrate, write the uncalibrated file to the output folder as an mzML
            File.Copy(originalUncalibratedFilePath, uncalibratedNewFullFilePath, true);
            // and add it to the list of all unsuccessfully calibrated files
            unsuccessfullyCalibratedFilePaths.Add(uncalibratedNewFullFilePath); 

            // provide a message indicating why we couldn't calibrate
            CalibrationWarnMessage(acquisitionResults);

            // mark the file as done
            FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
            _myFileManager.DoneWithFile(originalUncalibratedFilePath);

            //Update progress
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", Path.GetFileNameWithoutExtension(originalUncalibratedFilePath) }));
        }

        private bool CalibrationHasValue(DataPointAquisitionResults acquisitionResultsFirst, DataPointAquisitionResults acquisitionResultsSecond)
        {
            int psmsCountFirst = acquisitionResultsFirst.Psms.Count;
            int psmsCountSecond = acquisitionResultsSecond.Psms.Count;
            int peptidesCountFirst = acquisitionResultsFirst.Psms.Select(p => p.FullSequence).Distinct().Count();
            int peptidesCountSecond = acquisitionResultsSecond.Psms.Select(p => p.FullSequence).Distinct().Count();
            bool improvedCounts = psmsCountSecond > psmsCountFirst && peptidesCountSecond > peptidesCountFirst;
            if(improvedCounts)
            {
                return true;
            }
            bool numPsmsIncreased = psmsCountSecond >= psmsCountFirst;
            bool numPeptidesIncreased = peptidesCountSecond >= peptidesCountFirst;
            bool psmPrecursorMedianPpmErrorDecreased = Math.Abs(acquisitionResultsSecond.PsmPrecursorMedianPpmError) <= 1 || Math.Abs(acquisitionResultsSecond.PsmPrecursorMedianPpmError) <= Math.Abs(acquisitionResultsFirst.PsmPrecursorMedianPpmError);
            bool psmProductMedianPpmErrorDecreased = Math.Abs(acquisitionResultsSecond.PsmProductMedianPpmError) <= 1 || Math.Abs(acquisitionResultsSecond.PsmProductMedianPpmError) <= Math.Abs(acquisitionResultsFirst.PsmProductMedianPpmError);

            return (numPsmsIncreased && numPeptidesIncreased && psmPrecursorMedianPpmErrorDecreased && psmProductMedianPpmErrorDecreased);
        }
        
        private bool SufficientAcquisitionResults(DataPointAquisitionResults acquisitionResults)
        {
            return acquisitionResults.Psms.Count >= NumRequiredPsms 
                && acquisitionResults.Ms1List.Count >= NumRequiredMs1Datapoints 
                && acquisitionResults.Ms2List.Count >= NumRequiredMs2Datapoints;
        }

        public void WarnForWiderTolerance(double newPrecursorTolerance, double newProductTolerance)
        {
            Warn("Could not find enough PSMs to calibrate with; opening up tolerances to " +
                 Math.Round(newPrecursorTolerance, 2) + " ppm precursor and " +
                 Math.Round(newProductTolerance, 2) + " ppm product");
        }

        public void CalibrationWarnMessage(DataPointAquisitionResults acquisitionResults)
        {
            // provide a message indicating why we couldn't calibrate
            if (acquisitionResults.Psms.Count < NumRequiredPsms)
            {
                Warn($"Calibration failure! Could not find enough high-quality {GlobalVariables.AnalyteType.GetSpectralMatchLabel()}s. Required " + NumRequiredPsms + ", saw " + acquisitionResults.Psms.Count);
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

            foreach (SpectraFileInfo originalSpectraFile in oldExperDesign)
            {
                var originalFileName = originalSpectraFile.FilenameWithoutExtension;
                var k = unsuccessfullyCalibratedFilePaths.FirstOrDefault(fn => fn.Contains(originalFileName));
                if (k != null)
                {
                    newExperDesign.Add(new SpectraFileInfo(k, originalSpectraFile.Condition, originalSpectraFile.BiologicalReplicate, originalSpectraFile.TechnicalReplicate, originalSpectraFile.Fraction));
                }
                else
                {
                    SpectraFileInfo calibratedSpectraFile = new(Path.Combine(outputFolder, originalFileName + CalibSuffix + ".mzML"),
                    originalSpectraFile.Condition, originalSpectraFile.BiologicalReplicate, originalSpectraFile.TechnicalReplicate, originalSpectraFile.Fraction);
                    newExperDesign.Add(calibratedSpectraFile);
                }
            }

            _ = ExperimentalDesign.WriteExperimentalDesignToFile(newExperDesign);
        }
    }
}