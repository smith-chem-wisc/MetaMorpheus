using EngineLayer;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using FlashLFQ;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
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
            ProseCreatedWhileRunning.Append("The following calibration settings were used: ");
            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("max mods per peptide = " + CommonParameters.DigestionParams.MaxModsForPeptide + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            ProseCreatedWhileRunning.Append("precursor mass tolerance = " + CommonParameters.PrecursorMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + ". ");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count(p => !p.IsDecoy) + " non-decoy protein entries including " + proteinList.Count(p => p.IsContaminant) + " contaminant sequences. ");

            // start the calibration task
            Status("Calibrating...", new List<string> { taskId });
            MyTaskResults = new MyTaskResults(this)
            {
                NewSpectra = new List<string>(),
                NewFileSpecificTomls = new List<string>()
            };

            var myFileManager = new MyFileManager(true);
            List<string> unsuccessfullyCalibratedFilePaths = new List<string>();

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                bool couldNotFindEnoughDatapoints = false;

                // get filename stuff
                var originalUncalibratedFilePath = currentRawFileList[spectraFileIndex];
                var originalUncalibratedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUncalibratedFilePath);
                string calibratedFilePath = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".mzML");

                // mark the file as in-progress
                StartingDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                // load the file
                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files" });

                var myMsDataFile = myFileManager.LoadFile(originalUncalibratedFilePath, CommonParameters);

                // get datapoints to fit calibration function to
                Status("Acquiring calibration data points...", new List<string> { taskId, "Individual Spectra Files" });
                DataPointAquisitionResults acquisitionResults = null;

                for (int i = 1; i <= 5; i++)
                {
                    acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, combinedParams.PrecursorMassTolerance, combinedParams.ProductMassTolerance);

                    // enough data points to calibrate?
                    if (acquisitionResults.Psms.Count >= NumRequiredPsms && acquisitionResults.Ms1List.Count > NumRequiredMs1Datapoints && acquisitionResults.Ms2List.Count > NumRequiredMs2Datapoints)
                    {
                        break;
                    }

                    if (i == 1) // failed round 1
                    {
                        CommonParameters.PrecursorMassTolerance = new PpmTolerance(20);
                        CommonParameters.ProductMassTolerance = new PpmTolerance(50);
                    }
                    else if (i == 2) // failed round 2
                    {
                        CommonParameters.PrecursorMassTolerance = new PpmTolerance(30);
                        CommonParameters.ProductMassTolerance = new PpmTolerance(100);
                    }
                    else if (i == 3) // failed round 3
                    {
                        CommonParameters.PrecursorMassTolerance = new PpmTolerance(40);
                        CommonParameters.ProductMassTolerance = new PpmTolerance(150);
                    }
                    else // failed round 4
                    {
                        if (acquisitionResults.Psms.Count < NumRequiredPsms)
                        {
                            Warn("Calibration failure! Could not find enough high-quality PSMs. Required " + NumRequiredPsms + ", saw " + acquisitionResults.Psms.Count);
                        }
                        if (acquisitionResults.Ms1List.Count < NumRequiredMs1Datapoints)
                        {
                            Warn("Calibration failure! Could not find enough MS1 datapoints. Required " + NumRequiredMs1Datapoints + ", saw " + acquisitionResults.Ms1List.Count);
                        }
                        if (acquisitionResults.Ms2List.Count < NumRequiredMs2Datapoints)
                        {
                            Warn("Calibration failure! Could not find enough MS2 datapoints. Required " + NumRequiredMs2Datapoints + ", saw " + acquisitionResults.Ms2List.Count);
                        }

                        couldNotFindEnoughDatapoints = true;
                        FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                        break;
                    }

                    Warn("Could not find enough PSMs to calibrate with; opening up tolerances to " +
                    Math.Round(CommonParameters.PrecursorMassTolerance.Value, 2) + " ppm precursor and " +
                    Math.Round(CommonParameters.ProductMassTolerance.Value, 2) + " ppm product");
                }

                if (couldNotFindEnoughDatapoints)
                {
                    unsuccessfullyCalibratedFilePaths.Add(Path.GetFileNameWithoutExtension(currentRawFileList[spectraFileIndex]));
                    ReportProgress(new ProgressEventArgs(100, "Failed to calibrate!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
                    continue;
                }

                // stats before calibration
                int prevPsmCount = acquisitionResults.Psms.Count;
                double preCalibrationPrecursorErrorIqr = acquisitionResults.PsmPrecursorIqrPpmError;
                double preCalibrationProductErrorIqr = acquisitionResults.PsmProductIqrPpmError;

                // generate calibration function and shift data points
                Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                CalibrationEngine engine = new CalibrationEngine(myMsDataFile, acquisitionResults, CommonParameters, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                engine.Run();

                //update file
                myMsDataFile = engine.CalibratedDataFile;

                // do another search to evaluate calibration results
                Status("Post-calibration search...", new List<string> { taskId, "Individual Spectra Files" });
                acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, combinedParams.PrecursorMassTolerance, combinedParams.ProductMassTolerance);

                //generate calibration function and shift data points AGAIN because it's fast and contributes new data
                Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                engine = new CalibrationEngine(myMsDataFile, acquisitionResults, CommonParameters, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                engine.Run();

                //update file
                myMsDataFile = engine.CalibratedDataFile;

                // write the calibrated mzML file
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, calibratedFilePath, CalibrationParameters.WriteIndexedMzml);
                myFileManager.DoneWithFile(originalUncalibratedFilePath);

                // stats after calibration
                int postCalibrationPsmCount = acquisitionResults.Psms.Count;
                double postCalibrationPrecursorErrorIqr = acquisitionResults.PsmPrecursorIqrPpmError;
                double postCalibrationProductErrorIqr = acquisitionResults.PsmProductIqrPpmError;

                // did the data improve? (not used for anything yet...)
                bool improvement = ImprovGlobal(preCalibrationPrecursorErrorIqr, preCalibrationProductErrorIqr, prevPsmCount, postCalibrationPsmCount, postCalibrationPrecursorErrorIqr, postCalibrationProductErrorIqr);

                // write toml settings for the calibrated file
                var newTomlFileName = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".toml");

                var fileSpecificParams = new FileSpecificParameters();

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

        private bool ImprovGlobal(double prevPrecTol, double prevProdTol, int prevPsmCount, int thisRoundPsmCount, double thisRoundPrecTol, double thisRoundProdTol)
        {
            if (thisRoundPsmCount > prevPsmCount)
            {
                return true;
            }

            var precRatio = thisRoundPrecTol / prevPrecTol;
            var prodRatio = thisRoundProdTol / prevProdTol;

            if (thisRoundPsmCount == prevPsmCount)
            {
                return precRatio + prodRatio < 2; // Take any improvement in ratios
            }

            var countRatio = (double)thisRoundPsmCount / prevPsmCount;
            return countRatio > 0.9 && precRatio + prodRatio < 1.8;
        }

        private DataPointAquisitionResults GetDataAcquisitionResults(MsDataFile myMsDataFile, string currentDataFile, List<Modification> variableModifications, List<Modification> fixedModifications, List<Protein> proteinList, string taskId, CommonParameters combinedParameters, Tolerance initPrecTol, Tolerance initProdTol)
        {
            var fileNameWithoutExtension = Path.GetFileNameWithoutExtension(currentDataFile);
            MassDiffAcceptor searchMode = initPrecTol is PpmTolerance ?
                (MassDiffAcceptor)new SinglePpmAroundZeroSearchMode(initPrecTol.Value) :
                new SingleAbsoluteAroundZeroSearchMode(initPrecTol.Value);

            var listOfSortedms2Scans = GetMs2Scans(myMsDataFile, currentDataFile, combinedParameters).OrderBy(b => b.PrecursorMass).ToArray();
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            Log("Searching with searchMode: " + searchMode, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });
            Log("Searching with productMassTolerance: " + initProdTol, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });

            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null, proteinList, searchMode, combinedParameters,
                FileSpecificParameters, null, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();
            List<PeptideSpectralMatch> allPsms = allPsmsArray.Where(b => b != null).ToList();

            allPsms = allPsms.OrderByDescending(b => b.Score)
                .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                .GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

            new FdrAnalysisEngine(allPsms, searchMode.NumNotches, CommonParameters, FileSpecificParameters, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            List<PeptideSpectralMatch> goodIdentifications = allPsms.Where(b => b.FdrInfo.QValueNotch < 0.001 && !b.IsDecoy && b.FullSequence != null).ToList();

            if (!goodIdentifications.Any())
            {
                return new DataPointAquisitionResults(null, new List<PeptideSpectralMatch>(), new List<LabeledDataPoint>(), new List<LabeledDataPoint>(), 0, 0, 0, 0);
            }

            //get the deconvoluted ms2scans for the good identifications
            List<Ms2ScanWithSpecificMass> goodScans = new List<Ms2ScanWithSpecificMass>();
            List<PeptideSpectralMatch> unfilteredPsms = allPsmsArray.ToList();
            foreach (PeptideSpectralMatch psm in goodIdentifications)
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
            var oldExperDesign = ExperimentalDesign.ReadExperimentalDesign(pathToOldExperDesign, originalUncalibratedFileNamesWithExtension, out var errors);

            if (errors.Any())
            {
                foreach (var error in errors)
                {
                    Warn(error);
                }

                return;
            }

            var newExperDesign = new List<SpectraFileInfo>();

            foreach (var uncalibratedSpectraFile in oldExperDesign)
            {
                var originalUncalibratedFilePath = uncalibratedSpectraFile.FullFilePathWithExtension;
                var originalUncalibratedFilenameWithoutExtension = GlobalVariables.GetFilenameWithoutExtension(originalUncalibratedFilePath);
                string calibratedFilePath = Path.Combine(outputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".mzML");

                var calibratedSpectraFile = new SpectraFileInfo(calibratedFilePath, 
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

            ExperimentalDesign.WriteExperimentalDesignToFile(newExperDesign);
        }
    }
}