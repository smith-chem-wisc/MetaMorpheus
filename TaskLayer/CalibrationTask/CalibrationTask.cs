using EngineLayer;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class CalibrationTask : MetaMorpheusTask
    {
        private int NumRequiredPsms = 20;
        private int NumRequiredMs1Datapoints = 50;
        private int NumRequiredMs2Datapoints = 100;
        private const string CalibSuffix = "-calib";

        public CalibrationTask() : base(TaskType.Calibrate)
        {
            CommonParameters = new CommonParameters(
                productMassTolerance: new PpmTolerance(25),
                precursorMassTolerance: new PpmTolerance(15),
                trimMsMsPeaks: false,
                doPrecursorDeconvolution: false,
                scoreCutoff: 10);

            CalibrationParameters = new CalibrationParameters();
        }

        public CalibrationParameters CalibrationParameters { get; set; }

        protected override MyTaskResults RunSpecific(string outputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            // load modifications
            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown
                .OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown
                .OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
            List<string> localizeableModificationTypes = GlobalVariables.AllModTypesKnown.ToList();

            // what types of fragment ions to search for
            List<ProductType> ionTypes = new List<ProductType>();
            if (CommonParameters.BIons) { ionTypes.Add(ProductType.BnoB1ions); }
            if (CommonParameters.YIons) { ionTypes.Add(ProductType.Y); }
            if (CommonParameters.ZdotIons) { ionTypes.Add(ProductType.Zdot); }
            if (CommonParameters.CIons) { ionTypes.Add(ProductType.C); }

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, true, DecoyType.Reverse, localizeableModificationTypes);

            // write prose settings
            ProseCreatedWhileRunning.Append("The following calibration settings were used: ");
            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");
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

            object lock1 = new object();

            var myFileManager = new MyFileManager(true);

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                // get filename stuff
                var originalUncalibratedFilePath = currentRawFileList[spectraFileIndex];
                var originalUncalibratedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUncalibratedFilePath);
                string calibratedFilePath = Path.Combine(outputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".mzML");

                // mark the file as in-progress
                StartingDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                MsDataFile myMsDataFile;

                // load the file
                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files" });
                lock (lock1)
                {
                    myMsDataFile = myFileManager.LoadFile(originalUncalibratedFilePath, CommonParameters.TopNpeaks, CommonParameters.MinRatio, CommonParameters.TrimMs1Peaks, CommonParameters.TrimMsMsPeaks);
                }

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
                        this.CommonParameters.SetPrecursorMassTolerance(new PpmTolerance(20));
                        this.CommonParameters.SetProductMassTolerance(new PpmTolerance(50));
                    }
                    else if (i == 2) // failed round 2
                    {
                        this.CommonParameters.SetPrecursorMassTolerance(new PpmTolerance(30));
                        this.CommonParameters.SetProductMassTolerance(new PpmTolerance(100));
                    }
                    else if (i == 3) // failed round 3
                    {
                        this.CommonParameters.SetPrecursorMassTolerance(new PpmTolerance(40));
                        this.CommonParameters.SetProductMassTolerance(new PpmTolerance(150));
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
                        FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                        return MyTaskResults;
                    }

                    Warn("Could not find enough PSMs to calibrate with; opening up tolerances to " +
                    Math.Round(CommonParameters.PrecursorMassTolerance.Value, 2) + " ppm precursor and " +
                    Math.Round(CommonParameters.ProductMassTolerance.Value, 2) + " ppm product");
                }

                // stats before calibration
                int prevPsmCount = acquisitionResults.Psms.Count;
                double preCalibrationPrecursorErrorIqr = acquisitionResults.PsmPrecursorIqrPpmError;
                double preCalibrationProductErrorIqr = acquisitionResults.PsmProductIqrPpmError;

                // generate calibration function and shift data points
                Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                new CalibrationEngine(myMsDataFile, acquisitionResults, CommonParameters, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }).Run();

                // do another search to evaluate calibration results
                Status("Post-calibration search...", new List<string> { taskId, "Individual Spectra Files" });
                acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, combinedParams.PrecursorMassTolerance, combinedParams.ProductMassTolerance);

                // stats after calibration
                int postCalibrationPsmCount = acquisitionResults.Psms.Count;
                double postCalibrationPrecursorErrorIqr = acquisitionResults.PsmPrecursorIqrPpmError;
                double postCalibrationProductErrorIqr = acquisitionResults.PsmProductIqrPpmError;

                // did the data improve? (not used for anything yet...)
                bool improvement = ImprovGlobal(preCalibrationPrecursorErrorIqr, preCalibrationProductErrorIqr, prevPsmCount, postCalibrationPsmCount, postCalibrationPrecursorErrorIqr, postCalibrationProductErrorIqr);

                // write toml settings for the calibrated file
                var newTomlFileName = Path.Combine(outputFolder, originalUncalibratedFilenameWithoutExtension + CalibSuffix + ".toml");

                var fileSpecificParams = new FileSpecificParameters();

                // carry over file-specific parameters from the uncalibrated file to the calibrated one
                if (fileSettingsList[spectraFileIndex] != null)
                {
                    fileSpecificParams = fileSettingsList[spectraFileIndex].Clone();
                }

                // don't write over ppm tolerances if they've been specified by the user already in the file-specific settings
                // otherwise, suggest 4 * interquartile range as the ppm tolerance
                if (fileSpecificParams.PrecursorMassTolerance == null)
                {
                    fileSpecificParams.PrecursorMassTolerance = new PpmTolerance((4.0 * postCalibrationPrecursorErrorIqr) + Math.Abs(acquisitionResults.PsmPrecursorMedianPpmError));
                }
                if (fileSpecificParams.ProductMassTolerance == null)
                {
                    fileSpecificParams.ProductMassTolerance = new PpmTolerance((4.0 * postCalibrationProductErrorIqr) + Math.Abs(acquisitionResults.PsmProductMedianPpmError));
                }

                Toml.WriteFile(fileSpecificParams, newTomlFileName, tomlConfig);

                SucessfullyFinishedWritingFile(newTomlFileName, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });

                // write the calibrated mzML file
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, calibratedFilePath, false);
                myFileManager.DoneWithFile(originalUncalibratedFilePath);

                // finished calibrating this file
                SucessfullyFinishedWritingFile(calibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                MyTaskResults.NewSpectra.Add(calibratedFilePath);
                MyTaskResults.NewFileSpecificTomls.Add(newTomlFileName);
                FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
            }

            // re-write experimental design (if it has been defined) with new calibrated file names
            string assumedPathToExperDesign = Directory.GetParent(currentRawFileList.First()).FullName;
            assumedPathToExperDesign = Path.Combine(assumedPathToExperDesign, GlobalVariables.ExperimentalDesignFileName);
            List<string> newExperimentalDesignOutput = new List<string>();
            if (File.Exists(assumedPathToExperDesign))
            {
                var lines = File.ReadAllLines(assumedPathToExperDesign);

                for (int i = 0; i < lines.Length; i++)
                {
                    // header of experimental design file
                    if (i == 0)
                    {
                        newExperimentalDesignOutput.Add(lines[i]);
                    }
                    else
                    {
                        var split = lines[i].Split('\t');
                        string newline = Path.GetFileNameWithoutExtension(split[0]) + CalibSuffix + "\t";
                        for (int j = 1; j < split.Length; j++)
                        {
                            newline += split[j] + "\t";
                        }

                        newExperimentalDesignOutput.Add(newline);
                    }
                }
            }

            File.WriteAllLines(Path.Combine(outputFolder, GlobalVariables.ExperimentalDesignFileName), newExperimentalDesignOutput);

            // finished calibrating all files for the task
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            return MyTaskResults;
        }

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

        private DataPointAquisitionResults GetDataAcquisitionResults(MsDataFile myMsDataFile, string currentDataFile, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, string taskId, CommonParameters combinedParameters, Tolerance initPrecTol, Tolerance initProdTol)
        {
            var fileNameWithoutExtension = Path.GetFileNameWithoutExtension(currentDataFile);
            MassDiffAcceptor searchMode = initPrecTol is PpmTolerance ?
                (MassDiffAcceptor)new SinglePpmAroundZeroSearchMode(initPrecTol.Value) :
                new SingleAbsoluteAroundZeroSearchMode(initPrecTol.Value);

            FragmentTypes fragmentTypesForCalibration = FragmentTypes.None;
            if (combinedParameters.BIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.b;
            if (combinedParameters.YIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.y;
            if (combinedParameters.CIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.c;
            if (combinedParameters.ZdotIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.zdot;

            var listOfSortedms2Scans = GetMs2Scans(myMsDataFile, currentDataFile, combinedParameters.DoPrecursorDeconvolution, combinedParameters.UseProvidedPrecursorInfo, combinedParameters.DeconvolutionIntensityRatio, combinedParameters.DeconvolutionMaxAssumedChargeState, combinedParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            List<ProductType> productTypes = new List<ProductType>();
            if (combinedParameters.BIons)
                productTypes.Add(ProductType.B);
            if (combinedParameters.YIons)
                productTypes.Add(ProductType.Y);
            if (combinedParameters.CIons)
                productTypes.Add(ProductType.C);
            if (combinedParameters.ZdotIons)
                productTypes.Add(ProductType.Zdot);

            Log("Searching with searchMode: " + searchMode, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });
            Log("Searching with productMassTolerance: " + initProdTol, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });

            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productTypes, searchMode, combinedParameters, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            List<PeptideSpectralMatch> allPsms = allPsmsArray.ToList();

            var compactPeptideToProteinPeptideMatching = ((SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine
                (allPsms, proteinList, fixedModifications, variableModifications, productTypes, new List<DigestionParams> { combinedParameters.DigestionParams },
                combinedParameters.ReportAllAmbiguity, combinedParameters, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run()).CompactPeptideToProteinPeptideMatching;

            foreach (var huh in allPsms)
            {
                if (huh != null)
                {
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
                }
            }

            allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

            new FdrAnalysisEngine(allPsms, searchMode.NumNotches, CommonParameters, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            List<PeptideSpectralMatch> goodIdentifications = allPsms.Where(b => b.FdrInfo.QValueNotch < 0.01 && !b.IsDecoy && b.FullSequence != null).ToList();

            if (!goodIdentifications.Any())
            {
                return new DataPointAquisitionResults(null, new List<PeptideSpectralMatch>(), new List<LabeledDataPoint>(), new List<LabeledDataPoint>(), 0, 0, 0, 0);
            }

            DataPointAquisitionResults currentResult = (DataPointAquisitionResults)new DataPointAcquisitionEngine(
                    goodIdentifications,
                    myMsDataFile,
                    initPrecTol,
                    initProdTol,
                    CalibrationParameters.NumFragmentsNeededForEveryIdentification,
                    CalibrationParameters.MinMS1IsotopicPeaksNeededForConfirmedIdentification,
                    CalibrationParameters.MinMS2IsotopicPeaksNeededForConfirmedIdentification,
                    fragmentTypesForCalibration,
                    CommonParameters,
                    new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            return currentResult;
        }
    }
}