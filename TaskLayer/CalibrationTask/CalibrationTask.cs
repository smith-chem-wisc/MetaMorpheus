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
            };

            CalibrationParameters = new CalibrationParameters();
        }

        #endregion Public Constructors

        #region Public Properties

        public CalibrationParameters CalibrationParameters { get; set; }

        #endregion Public Properties

        #region Protected Internal Methods

        protected internal void WriteMs2DataPoints(List<LabeledMs2DataPoint> items, string outputFolder, string fileName, List<string> nestedIDs)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".ms2dptsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(LabeledMs2DataPoint.TabSeparatedHeader + "\t" + PeptideSpectralMatch.GetTabSeparatedHeader());
                foreach (var dp in items)
                {
                    output.Write(dp.Values());
                    output.WriteLine("\t" + dp.identification);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIDs);
        }

        protected internal void WriteMs1DataPoints(List<LabeledMs1DataPoint> items, string outputFolder, string fileName, List<string> nestedIDs)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".ms1dptsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(LabeledMs1DataPoint.TabSeparatedHeader + "\t" + PeptideSpectralMatch.GetTabSeparatedHeader());
                foreach (var dp in items)
                {
                    output.Write(dp.Values());
                    output.WriteLine("\t" + dp.identification);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIDs);
        }

        #endregion Protected Internal Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList)
        {
            myTaskResults = new MyTaskResults(this)
            {
                newSpectra = new List<string>()
            };

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
            List<string> localizeableModificationTypes = CommonParameters.LocalizeAll ? GlobalVariables.AllModTypesKnown.ToList() : CommonParameters.ListOfModTypesLocalize.ToList();

            Status("Loading proteins...", new List<string> { taskId });
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, UsefulProteomicsDatabases.DecoyType.Reverse, localizeableModificationTypes, b.IsContaminant, out Dictionary<string, Modification> um)).ToList();

            proseCreatedWhileRunning.Append("The following calibration settings were used: ");
            proseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            proseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            proseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            if (CommonParameters.DigestionParams.MaxPeptideLength == null)
            {
                proseCreatedWhileRunning.Append("maximum peptide length = unspecified; ");
            }
            else
            {
                proseCreatedWhileRunning.Append("maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            }
            proseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            proseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");

            proseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            object lock1 = new object();
            ParallelOptions parallelOptions = new ParallelOptions();
            if (CommonParameters.MaxParallelFilesToAnalyze.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxParallelFilesToAnalyze.Value;
            var myFileManager = new MyFileManager(true);

            Status("Calibrating...", new List<string> { taskId });
            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var originalUncalibratedFilePath = currentRawFileList[spectraFileIndex];
                var originalUncalibratedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUncalibratedFilePath);

                ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                string calibratedFilePath = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + "-calib.mzml");

                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;

                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files" });
                // only load one file at a time
                lock (lock1)
                {
                    myMsDataFile = myFileManager.LoadFile(originalUncalibratedFilePath, CommonParameters.TopNpeaks, CommonParameters.MinRatio, CommonParameters.TrimMs1Peaks, CommonParameters.TrimMsMsPeaks);
                }

                Status("Aquiring calibration data points...", new List<string> { taskId, "Individual Spectra Files" });
                // get datapoints to fit calibration function to
                var acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, combinedParams.PrecursorMassTolerance, combinedParams.ProductMassTolerance);
                
                // stats before calibration
                int prevPsmCount = acquisitionResults.Item1.Count;

                var preCalibrationPrecursorErrors = acquisitionResults.Item1.Select(p => (p.ScanPrecursorMass - p.PeptideMonisotopicMass.Value) / p.PeptideMonisotopicMass.Value * 1e6).ToList();
                double preCalibrationPrecursorIqr = Statistics.InterquartileRange(preCalibrationPrecursorErrors);

                var preCalibrationProductErrors = acquisitionResults.Item1.SelectMany(p => p.ProductMassErrorPpm.SelectMany(v => v.Value)).ToList();
                double preCalibrationProductIqr = Statistics.InterquartileRange(preCalibrationProductErrors);
                
                // enough data points to calibrate with?
                if (acquisitionResults.Item2 == null)
                {
                    Warn("Could not find any datapoints to calibrate with!");
                    return;
                }
                if (acquisitionResults.Item2.Ms1List.Count < 4 || acquisitionResults.Item2.Ms2List.Count < 4)
                {
                    Warn("Could not find enough MS1 datapoints to calibrate (" + acquisitionResults.Item2.Ms1List.Count + " found)");
                    Warn("Could not find enough MS2 datapoints to calibrate (" + acquisitionResults.Item2.Ms2List.Count + " found)");
                    return;
                }

                // generate calibration function and shift data points
                Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });
                new CalibrationEngine(myMsDataFile, acquisitionResults.Item2, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }).Run();

                // do another search to evaluate calibration results
                Status("Post-calibration search...", new List<string> { taskId, "Individual Spectra Files" });
                acquisitionResults = GetDataAcquisitionResults(myMsDataFile, originalUncalibratedFilePath, variableModifications, fixedModifications, proteinList, taskId, combinedParams, combinedParams.PrecursorMassTolerance, combinedParams.ProductMassTolerance);

                // stats after calibration
                int postCalibrationPsmCount = acquisitionResults.Item1.Count;

                var postCalibrationPrecursorErrors = acquisitionResults.Item1.Select(p => (p.ScanPrecursorMass - p.PeptideMonisotopicMass) / p.PeptideMonisotopicMass * 1e6).ToList();
                double postCalibrationPrecursorIqr = Statistics.InterquartileRange(postCalibrationPrecursorErrors);

                var postCalibrationProductErrors = acquisitionResults.Item1.SelectMany(p => p.ProductMassErrorPpm.SelectMany(v => v.Value)).ToList();
                double postCalibrationProductIqr = Statistics.InterquartileRange(postCalibrationProductErrors);

                // did the data improve? (not used for anything yet...)
                bool improvement = ImprovGlobal(preCalibrationPrecursorIqr, preCalibrationProductIqr, prevPsmCount, postCalibrationPsmCount, postCalibrationPrecursorIqr, postCalibrationProductIqr);

                // write suggested tolerances for this file
                var tomlFileName = Path.Combine(OutputFolder, originalUncalibratedFilenameWithoutExtension + "-calib.toml");
                FileSpecificTolerances f = new FileSpecificTolerances
                {
                    PrecursorMassTolerance = new PpmTolerance(4.0 * postCalibrationPrecursorIqr),
                    ProductMassTolerance = new PpmTolerance(4.0 * postCalibrationProductIqr)
                };
                Toml.WriteFile(f, tomlFileName, tomlConfig);
                SucessfullyFinishedWritingFile(tomlFileName, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });

                // write the calibrated MZML file
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, calibratedFilePath, false);

                // all done
                SucessfullyFinishedWritingFile(calibratedFilePath, new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension });
                myTaskResults.newSpectra.Add(calibratedFilePath);
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUncalibratedFilenameWithoutExtension }));
            });
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private bool ImprovGlobal(double prevPrecTol, double prevProdTol, int prevPsmCount, int thisRoundPsmCount, double thisRoundPrecTol, double thisRoundProdTol)
        {
            if (thisRoundPsmCount > prevPsmCount)
                return true;

            var precRatio = thisRoundPrecTol / prevPrecTol;
            var prodRatio = thisRoundProdTol / prevProdTol;

            if (thisRoundPsmCount == prevPsmCount)
                return precRatio + prodRatio < 2; // Take any improvement in ratios

            var countRatio = (double)thisRoundPsmCount / prevPsmCount;
            return countRatio > 0.9 && precRatio + prodRatio < 1.8;
        }

        private (List<PeptideSpectralMatch>, DataPointAquisitionResults) GetDataAcquisitionResults(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, string currentDataFile, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, string taskId, ICommonParameters combinedParameters, Tolerance initPrecTol, Tolerance initProdTol)
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

            new FdrAnalysisEngine(allPsms, searchMode.NumNotches, false, new List<string> { taskId, "Individual Spectra Files", fileNameWithoutExtension }).Run();

            List<PeptideSpectralMatch> goodIdentifications = allPsms.Where(b => b.FdrInfo.QValueNotch < 0.01 && !b.IsDecoy && b.FullSequence != null).ToList();

            if (!goodIdentifications.Any())
            {
                Warn("No PSMs below 1% FDR observed!");
                return (new List<PeptideSpectralMatch>(), null);
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

            return (goodIdentifications, currentResult);
        }

        #endregion Private Methods

        #region Private Classes

        private class FileSpecificTolerances
        {
            #region Public Properties

            public Tolerance PrecursorMassTolerance { get; set; }
            public Tolerance ProductMassTolerance { get; set; }

            #endregion Public Properties
        }

        #endregion Private Classes
    }
}