using EngineLayer;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using EngineLayer.FdrAnalysis;
using MathNet.Numerics.Statistics;
using EngineLayer.Localization;
using UsefulProteomicsDatabases;

#if NETFRAMEWORK

using IO.Thermo;

#else
#endif

namespace TaskLayer
{
    public class CalibrationTask : MetaMorpheusTask
    {
        #region Public Constructors

        public CalibrationTask() : base(MyTask.Calibrate)
        {
            CommonParameters = new CommonParameters
            {
                ProductMassTolerance = new PpmTolerance(25),
                PrecursorMassTolerance = new PpmTolerance(15),
                TrimMs1Peaks = false,
                TrimMsMsPeaks = false,
                DoPrecursorDeconvolution = false,
                ScoreCutoff = 10
            };

            CalibrationParameters = new CalibrationParameters();
        }

        #endregion Public Constructors

        #region Public Properties

        public CalibrationParameters CalibrationParameters { get; set; }

        #endregion Public Properties

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            if (CommonParameters.ScoreCutoff < 10)
            {
                Warn("The score cutoff for calibration is less than 10; calibration will continue but may give poor results");
            }

            // load modifications
            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
            List<string> localizeableModificationTypes = CommonParameters.LocalizeAll ? GlobalVariables.AllModTypesKnown.ToList() : CommonParameters.ListOfModTypesLocalize.ToList();

            // what types of fragment ions to search for
            List<ProductType> ionTypes = new List<ProductType>();
            if (CommonParameters.BIons)
                ionTypes.Add(ProductType.BnoB1ions);
            if (CommonParameters.YIons)
                ionTypes.Add(ProductType.Y);
            if (CommonParameters.ZdotIons)
                ionTypes.Add(ProductType.Zdot);
            if (CommonParameters.CIons)
                ionTypes.Add(ProductType.C);

            // load proteins
            Status("Loading proteins...", new List<string> { taskId });
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, DecoyType.Reverse, localizeableModificationTypes, b.IsContaminant, out Dictionary<string, Modification> unknownModifications)).ToList();

            // write prose settings
            proseCreatedWhileRunning.Append("The following calibration settings were used: ");
            proseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            proseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            proseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            proseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; "); proseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            proseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            proseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("max mods per peptide = " + CommonParameters.DigestionParams.MaxModsForPeptide + "; ");
            proseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            if (CommonParameters.ProductMassTolerance is PpmTolerance)
            {
                proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " ppm. ");
            }
            else
            {
                proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            }
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Count(p => p.IsContaminant) + " contaminant sequences. ");

            // start the calibration task
            Status("Calibrating...", new List<string> { taskId });
            myTaskResults = new MyTaskResults(this)
            {
                newSpectra = new List<string>(),
                newFileSpecificTomls = new List<string>()
            };

            object lock1 = new object();

            ParallelOptions parallelOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = CommonParameters.MaxParallelFilesToAnalyze
            };

            var myFileManager = new MyFileManager(true);

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                // get filename stuff
                var originalUncalibratedFilePath = currentRawFileList[spectraFileIndex];
                var originalUncalibratedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUncalibratedFilePath);
                string calibratedFilePath = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + "-calib.mzml");

                // mark the file as in-progress
                StartingDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });

                ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;

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
                    if (acquisitionResults.Psms.Count >= numRequiredPsms && acquisitionResults.Ms1List.Count > numRequiredMs1Datapoints && acquisitionResults.Ms2List.Count > numRequiredMs2Datapoints)
                    {
                        break;
                    }

                    // not enough datapoints seen; open the tolerance and try again
                    var newParameters = ((CommonParameters)CommonParameters).Clone();

                    if (i == 1) // failed round 1
                    {
                        newParameters.PrecursorMassTolerance = new PpmTolerance(20);
                        newParameters.ProductMassTolerance = new PpmTolerance(50);
                    }
                    else if (i == 2) // failed round 2
                    {
                        newParameters.PrecursorMassTolerance = new PpmTolerance(30);
                        newParameters.ProductMassTolerance = new PpmTolerance(100);
                    }
                    else if (i == 3) // failed round 3
                    {
                        newParameters.PrecursorMassTolerance = new PpmTolerance(40);
                        newParameters.ProductMassTolerance = new PpmTolerance(150);
                    }
                    else // failed round 4
                    {
                        if (acquisitionResults.Psms.Count < numRequiredPsms)
                        {
                            Warn("Calibration failure! Could not find enough high-quality PSMs. Required " + numRequiredPsms + ", saw " + acquisitionResults.Psms.Count);
                        }
                        if (acquisitionResults.Ms1List.Count < numRequiredMs1Datapoints)
                        {
                            Warn("Calibration failure! Could not find enough MS1 datapoints. Required " + numRequiredMs1Datapoints + ", saw " + acquisitionResults.Ms1List.Count);
                        }
                        if (acquisitionResults.Ms2List.Count < numRequiredMs2Datapoints)
                        {
                            Warn("Calibration failure! Could not find enough MS2 datapoints. Required " + numRequiredMs2Datapoints + ", saw " + acquisitionResults.Ms2List.Count);
                        }
                        FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                        return;
                    }

                    CommonParameters = newParameters;

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
                new CalibrationEngine(myMsDataFile, acquisitionResults, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }).Run();

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
                var newTomlFileName = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + "-calib.toml");

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
                myTaskResults.newSpectra.Add(calibratedFilePath);
                myTaskResults.newFileSpecificTomls.Add(newTomlFileName);
                FinishedDataFile(originalUncalibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilePath });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
            });

            // finished calibrating all files for the task
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Fields

        private int numRequiredPsms = 20;
        private int numRequiredMs1Datapoints = 50;
        private int numRequiredMs2Datapoints = 100;

        #endregion Private Fields

        #region Private Methods

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

        private DataPointAquisitionResults GetDataAcquisitionResults(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, string currentDataFile, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, string taskId, ICommonParameters combinedParameters, Tolerance initPrecTol, Tolerance initProdTol)
        {
            var fileNameWithoutExtension = Path.GetFileNameWithoutExtension(currentDataFile);
            MassDiffAcceptor searchMode;
            if (initPrecTol is PpmTolerance)
                searchMode = new SinglePpmAroundZeroSearchMode(initPrecTol.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode(initPrecTol.Value);

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

            List<ProductType> lp = new List<ProductType>();
            if (combinedParameters.BIons)
                lp.Add(ProductType.B);
            if (combinedParameters.YIons)
                lp.Add(ProductType.Y);
            if (combinedParameters.CIons)
                lp.Add(ProductType.C);
            if (combinedParameters.ZdotIons)
                lp.Add(ProductType.Zdot);

            Log("Searching with searchMode: " + searchMode, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });
            Log("Searching with productMassTolerance: " + initProdTol, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension });

            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, lp, searchMode, false, combinedParameters, initProdTol, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            List<PeptideSpectralMatch> allPsms = allPsmsArray.ToList();

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = ((SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, lp, new List<IDigestionParams> { combinedParameters.DigestionParams }, combinedParameters.ReportAllAmbiguity, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run()).CompactPeptideToProteinPeptideMatching;

            foreach (var huh in allPsms)
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

            allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

            new FdrAnalysisEngine(allPsms, searchMode.NumNotches, CommonParameters, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            List<PeptideSpectralMatch> goodIdentifications = allPsms.Where(b => b.FdrInfo.QValueNotch < 0.01 && !b.IsDecoy && b.FullSequence != null).ToList();

            if (!goodIdentifications.Any())
            {
                return new DataPointAquisitionResults(null, new List<PeptideSpectralMatch>(), new List<LabeledDataPoint>(), new List<LabeledDataPoint>(), 0, 0, 0, 0);
            }

            var dissociationTypes = MetaMorpheusEngine.DetermineDissociationType(lp);
            foreach (var psm in allPsms)
            {
                var theScan = myMsDataFile.GetOneBasedScan(psm.ScanNumber);
                double thePrecursorMass = psm.ScanPrecursorMass;

                foreach (var huh in lp)
                {
                    var ionMasses = psm.CompactPeptides.First().Key.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { huh });
                    Array.Sort(ionMasses);
                    List<double> matchedIonMassesList = new List<double>();
                    List<double> productMassErrorDaList = new List<double>();
                    List<double> productMassErrorPpmList = new List<double>();
                    LocalizationEngine.MatchIons(theScan, initProdTol, ionMasses, matchedIonMassesList, productMassErrorDaList, productMassErrorPpmList, thePrecursorMass, dissociationTypes, false);
                    double[] matchedIonMassesOnlyMatches = matchedIonMassesList.ToArray();
                    psm.MatchedIonDictOnlyMatches.Add(huh, matchedIonMassesOnlyMatches);
                    psm.ProductMassErrorDa.Add(huh, productMassErrorDaList.ToArray());
                    psm.ProductMassErrorPpm.Add(huh, productMassErrorPpmList.ToArray());
                }
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
                    new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            return currentResult;
        }

        #endregion Private Methods
    }
}