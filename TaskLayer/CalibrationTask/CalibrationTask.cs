using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using SharpLearning.Common.Interfaces;
using SharpLearning.GradientBoost.Learners;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace TaskLayer
{
    public class CalibrationTask : MetaMorpheusTask
    {
        #region Public Constructors

        public CalibrationTask() : base(MyTask.Calibrate)
        {
            CommonParameters = new CommonParameters
            {
                ProductMassTolerance = new PpmTolerance(30),
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
                output.WriteLine(LabeledMs2DataPoint.TabSeparatedHeader + "\t" + Psm.GetTabSeparatedHeader());
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
                output.WriteLine(LabeledMs1DataPoint.TabSeparatedHeader + "\t" + Psm.GetTabSeparatedHeader());
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

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettings)
        {
            myTaskResults = new MyTaskResults(this)
            {
                newSpectra = new List<string>()
            };

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> localizeableModifications;
            if (CommonParameters.LocalizeAll)
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            Status("Loading proteins...", new List<string> { taskId });
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> um)).ToList();

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
            //proseCreatedWhileRunning.Append("parent mass tolerance(s) = {" + searchMode.ToProseString() + "}; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            object lock1 = new object();
            ParallelOptions parallelOptions = new ParallelOptions();
            if (CommonParameters.MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxDegreeOfParallelism.Value;
            Status("Calibrating...", new List<string> { taskId });

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    var currentDataFile = currentRawFileList[spectraFileIndex];

                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = Mzml.LoadAllStaticData(currentDataFile);

                    var initLearners = new List<ILearner<double>>()
                    {
                        new MzMultiplier()
                    };

                    var mzSepLearners = new List<ILearner<double>>()
                    {
                        new LinearCalibrationFunctionMathNet(new int[] { }),

                        new LinearCalibrationFunctionMathNet(new [] { 0 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1 }),
                        new LinearCalibrationFunctionMathNet(new [] { 2 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 2 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2 }),

                        new RegressionAbsoluteLossGradientBoostLearner( iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner( maximumTreeDepth:6, iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner( maximumTreeDepth:9, iterations:1000),
                    };

                    mzSepLearners = mzSepLearners.Select(b => new SeparateMzLearner(b) as ILearner<double>).ToList();

                    var learners = new List<ILearner<double>>()
                    {
                        new LinearCalibrationFunctionMathNet(new int[] { }),

                        new LinearCalibrationFunctionMathNet(new [] { 0 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1 }),
                        new LinearCalibrationFunctionMathNet(new [] { 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 4 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 2, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 2, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 3, 4 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 2, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 2, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 3, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 2, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 2, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 3, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 2, 3, 4 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 3, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 2, 3, 4 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 2, 3, 4 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2, 3, 4 }),

                        new RegressionAbsoluteLossGradientBoostLearner( iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner( maximumTreeDepth:6, iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner( maximumTreeDepth:9, iterations:1000),
                    };

                    var indLearners = new List<ILearner<double>>()
                    {
                        new LinearCalibrationFunctionMathNet(new int[] { }),

                        new LinearCalibrationFunctionMathNet(new [] { 0 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1 }),

                        new RegressionAbsoluteLossGradientBoostLearner( iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner( maximumTreeDepth:6, iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner( maximumTreeDepth:9, iterations:1000),
                    };

                    //List<string> learnerType = new List<string> { "init", "mid", "final" };
                    //List<List<ILearner<double>>> alllearnerse = new List<List<ILearner<double>>> { initLearners, mzSepLearners, learners };
                    //List<string> learnerType = new List<string> { "mid", "final" };
                    //List<List<ILearner<double>>> alllearnerse = new List<List<ILearner<double>>> { mzSepLearners, learners };
                    //List<string> learnerType = new List<string> {  "final" };
                    //List<List<ILearner<double>>> alllearnerse = new List<List<ILearner<double>>> { learners };
                    //List<string> learnerType = new List<string> { "init", "final" };
                    ////List<List<ILearner<double>>> alllearnerse = new List<List<ILearner<double>>> { initLearners, learners };
                    //List<string> learnerType = new List<string> { "mid", "onlyIndividual" };
                    //List<List<ILearner<double>>> alllearners = new List<List<ILearner<double>>> { mzSepLearners, indLearners };
                    List<string> learnerType = new List<string> { "init", "init", "init", "init", "init", "init", "mid", "onlyIndividual" };
                    List<List<ILearner<double>>> alllearners = new List<List<ILearner<double>>> { initLearners, initLearners, initLearners, initLearners, initLearners, initLearners, mzSepLearners, indLearners };

                    for (int i = 0; i < alllearners.Count; i++)
                    {
                        Console.WriteLine(learnerType[i]);

                        DoTheLoop(
                            myMsDataFile,
                            currentDataFile,
                            variableModifications,
                            fixedModifications,
                            proteinList,
                            taskId,
                            OutputFolder,
                            alllearners[i],
                            learnerType[i]);

                        Console.WriteLine("  PrecursorMassTolerance: " + CalibrationParameters.PrecursorMassTolerance);
                        Console.WriteLine("  ProductMassTolerance: " + CommonParameters.ProductMassTolerance);
                    }

                    Warn("");
                    Warn("");
                }
                );

            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private void DoTheLoop(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, string currentDataFile, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, string taskId, string outputFolder, List<ILearner<double>> learners, string learnerType)
        {
            MassDiffAcceptor searchMode;
            if (CalibrationParameters.PrecursorMassTolerance is PpmTolerance)
                searchMode = new SinglePpmAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);

            FragmentTypes fragmentTypesForCalibration = FragmentTypes.None;
            if (CommonParameters.BIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.b;
            if (CommonParameters.YIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.y;
            if (CommonParameters.CIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.c;
            if (CommonParameters.ZdotIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.zdot;

            List<Psm> goodIdentifications = GetGoodIdentifications(myMsDataFile, searchMode, currentDataFile, variableModifications, fixedModifications, proteinList, taskId);

            // Store

            Console.WriteLine("  CurrentGoodIDs: " + goodIdentifications.Count);

            DataPointAquisitionResults prevResult = new DataPointAquisitionResults(null, new List<LabeledMs1DataPoint>(), new List<LabeledMs2DataPoint>(), 0, 0, 0, 0);
            Tolerance prevPrecursorMassToleranceForDatapointAcquisition = CalibrationParameters.PrecursorMassTolerance;
            Tolerance prevProductMassToleranceForDatapointAcquisition = CommonParameters.ProductMassTolerance;

            Tolerance precursorMassToleranceForDatapointAcquisition = CalibrationParameters.PrecursorMassTolerance;
            Tolerance productMassToleranceForDatapointAcquisition = CommonParameters.ProductMassTolerance;
            do
            {
                DataPointAquisitionResults currentResult = (DataPointAquisitionResults)new DataPointAcquisitionEngine(
                    goodIdentifications,
                    myMsDataFile,
                    precursorMassToleranceForDatapointAcquisition,
                    productMassToleranceForDatapointAcquisition,
                    CalibrationParameters.numFragmentsNeededForEveryIdentification,
                    CalibrationParameters.minMS1isotopicPeaksNeededForConfirmedIdentification,
                    CalibrationParameters.minMS2isotopicPeaksNeededForConfirmedIdentification,
                    fragmentTypesForCalibration,
                    new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                precursorMassToleranceForDatapointAcquisition = new PpmTolerance(Math.Max(Math.Abs(currentResult.Ms1InfoPPM.Item1 + 6 * currentResult.Ms1InfoPPM.Item2), Math.Abs(currentResult.Ms1InfoPPM.Item1 - 6 * currentResult.Ms1InfoPPM.Item2)));
                productMassToleranceForDatapointAcquisition = new PpmTolerance(Math.Max(Math.Abs(currentResult.Ms2InfoPPM.Item1 + 6 * currentResult.Ms2InfoPPM.Item2), Math.Abs(currentResult.Ms2InfoPPM.Item1 - 6 * currentResult.Ms2InfoPPM.Item2)));

                #region result output and writing

                WriteMs1DataPoints(currentResult.Ms1List, outputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + learnerType, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                WriteMs2DataPoints(currentResult.Ms2List, outputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + learnerType, new List<string> { taskId, "Individual Spectra Files", currentDataFile });

                Console.WriteLine("    MS1th : " + currentResult.Ms1InfoTh.Item1 + " : " + currentResult.Ms1InfoTh.Item2 + " :  MS2th : " + currentResult.Ms2InfoTh.Item1 + " : " + currentResult.Ms2InfoTh.Item2);
                Console.WriteLine("    MS1ppm : " + currentResult.Ms1InfoPPM.Item1 + " : " + currentResult.Ms1InfoPPM.Item2 + " :  MS2ppm : " + currentResult.Ms2InfoPPM.Item1 + " : " + currentResult.Ms2InfoPPM.Item2);
                Console.WriteLine("    currentResult.Ms1List.Count : " + currentResult.Ms1List.Count + " : currentResult.Ms2List.Count : " + currentResult.Ms2List.Count + " :  currentResult.Count : " + currentResult.Count);
                Console.WriteLine("    currentResult.numMs1MassChargeCombinationsConsidered : " + currentResult.NumMs1MassChargeCombinationsConsidered);
                Console.WriteLine("    currentResult.numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks : " + currentResult.NumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
                Console.WriteLine("    currentResult.numMs2MassChargeCombinationsConsidered : " + currentResult.NumMs2MassChargeCombinationsConsidered);
                Console.WriteLine("    currentResult.numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks : " + currentResult.NumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
                Console.WriteLine("    precursorMassToleranceForDatapointAcquisition : " + precursorMassToleranceForDatapointAcquisition);
                Console.WriteLine("    productMassToleranceForDatapointAcquisition : " + productMassToleranceForDatapointAcquisition);

                #endregion result output and writing

                if (ImprovedResults(prevResult, currentResult))
                {
                    prevResult = currentResult;
                    prevPrecursorMassToleranceForDatapointAcquisition = precursorMassToleranceForDatapointAcquisition;
                    prevProductMassToleranceForDatapointAcquisition = productMassToleranceForDatapointAcquisition;
                }
                else
                    break;
            } while (true);

            CalibrationResults res = (CalibrationResults)new CalibrationEngine(myMsDataFile, prevResult, learners, learnerType, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

            // Update for reuse
            CalibrationParameters.PrecursorMassTolerance = prevPrecursorMassToleranceForDatapointAcquisition;
            searchMode = new SinglePpmAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
            CommonParameters.ProductMassTolerance = prevProductMassToleranceForDatapointAcquisition;

            // Store
            //paths.Add(Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "-" + learnerType + "-SearchRound" + searchRound + ".mzML"));
            //MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, paths.Last(), false);
        }

        private bool ImprovedResults(DataPointAquisitionResults prevResult, DataPointAquisitionResults currentResult)
        {
            //if ((prevResult.Ms1List.Count <= currentResult.Ms1List.Count && prevResult.Ms2List.Count < currentResult.Ms2List.Count)
            //    || (prevResult.Ms1List.Count < currentResult.Ms1List.Count && prevResult.Ms2List.Count <= currentResult.Ms2List.Count)
            //    || (prevResult.Ms1List.Count == currentResult.Ms1List.Count && prevResult.Ms2List.Count == currentResult.Ms2List.Count && ((prevResult.Ms1InfoPPM.Item2 > currentResult.Ms1InfoPPM.Item2 && prevResult.Ms2InfoPPM.Item2 >= currentResult.Ms2InfoPPM.Item2) || (prevResult.Ms1InfoPPM.Item2 >= currentResult.Ms1InfoPPM.Item2 && prevResult.Ms2InfoPPM.Item2 > currentResult.Ms2InfoPPM.Item2))))
            if ((prevResult.Ms1List.Count <= currentResult.Ms1List.Count && prevResult.Ms2List.Count < currentResult.Ms2List.Count)
               || (prevResult.Ms1List.Count < currentResult.Ms1List.Count && prevResult.Ms2List.Count <= currentResult.Ms2List.Count)
               || (prevResult.Ms1InfoPPM.Item2 > currentResult.Ms1InfoPPM.Item2 && prevResult.Ms2InfoPPM.Item2 >= currentResult.Ms2InfoPPM.Item2) || (prevResult.Ms1InfoPPM.Item2 >= currentResult.Ms1InfoPPM.Item2 && prevResult.Ms2InfoPPM.Item2 > currentResult.Ms2InfoPPM.Item2))
                return true;
            return false;
        }

        private bool BreakBecauseOfWorseResults(List<int> numberOfGoodIds, List<Tolerance> precursorMassTolerances, List<Tolerance> productMassTolerances)
        {
            var idRatio = (double)numberOfGoodIds.Last() / numberOfGoodIds[numberOfGoodIds.Count - 2];
            var tolRatio1 = precursorMassTolerances.Last().Value / precursorMassTolerances[precursorMassTolerances.Count - 2].Value;
            var tolRatio2 = productMassTolerances.Last().Value / productMassTolerances[productMassTolerances.Count - 2].Value;

            return
                    // Break if too much decrease in number of ids
                    idRatio <= 0.9 // Break if too much decrease in number of ids
                ||
                   // Break if ids decreased a bit, but tols did not decrase enough
                   (idRatio > 0.9 && idRatio < 1) && (tolRatio1 > 1 || tolRatio2 > 1 || (tolRatio1 >= 0.9 && tolRatio2 >= 0.9))
                ||
                   // Break if ids did not decrease but tols did not decrease
                   (idRatio >= 1) && (tolRatio1 >= 1 || tolRatio2 >= 1);
        }

        private List<Psm> GetGoodIdentifications(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, MassDiffAcceptor searchMode, string currentDataFile, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, string taskId)
        {
            List<ProductType> lp = new List<ProductType>();
            if (CommonParameters.BIons)
                lp.Add(ProductType.B);
            if (CommonParameters.YIons)
                lp.Add(ProductType.Y);
            if (CommonParameters.CIons)
                lp.Add(ProductType.C);
            if (CommonParameters.ZdotIons)
                lp.Add(ProductType.Zdot);

            var listOfSortedms2Scans = GetMs2Scans(myMsDataFile, currentDataFile, CommonParameters.DoPrecursorDeconvolution, CommonParameters.UseProvidedPrecursorInfo, CommonParameters.DeconvolutionIntensityRatio, CommonParameters.DeconvolutionMaxAssumedChargeState, CommonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];

            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, lp, searchMode, false, CommonParameters, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

            List<Psm> allPsms = allPsmsArray.ToList();

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = ((SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, lp, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run()).CompactPeptideToProteinPeptideMatching;

            foreach (var huh in allPsms)
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

            allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

            new FdrAnalysisEngine(allPsms, searchMode, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

            return allPsms.Where(b => b.FdrInfo.QValue < 0.01 && !b.IsDecoy).ToList();
        }

        #endregion Private Methods
    }
}