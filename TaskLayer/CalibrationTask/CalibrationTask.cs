using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using SharpLearning.Common.Interfaces;
using SharpLearning.Ensemble.Learners;
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
            MassDiffAcceptor searchMode;
            if (CalibrationParameters.PrecursorMassTolerance is PpmTolerance)
                searchMode = new SinglePpmAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);

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
            proseCreatedWhileRunning.Append("parent mass tolerance(s) = {" + searchMode.ToProseString() + "}; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            object lock1 = new object();
            ParallelOptions parallelOptions = new ParallelOptions();
            if (CommonParameters.MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxDegreeOfParallelism.Value;
            Status("Calibrating...", new List<string> { taskId });

            Action<List<LabeledMs1DataPoint>> ms1Action = (List<LabeledMs1DataPoint> theList) => {; };
            Action<List<LabeledMs2DataPoint>> ms2Action = (List<LabeledMs2DataPoint> theList) => {; };

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    var currentDataFile = currentRawFileList[spectraFileIndex];

                    int searchRound = 1;

                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = Mzml.LoadAllStaticData(currentDataFile);

                    List<ILearner<double>> linearLearners = new List<ILearner<double>>()
                    {
                        new MzMultiplier()
                    };

                    //List<ILearner<double>> linearLearners = new List<ILearner<double>>()
                    //{
                    //    new LinearCalibrationFunctionMathNet(new int[] { }),

                    //    new LinearCalibrationFunctionMathNet(new [] { 0 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 1 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 2 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 3 }),

                    //    new LinearCalibrationFunctionMathNet(new [] { 0, 1 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 0, 2 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 0, 3 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 1, 2 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 1, 3 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 2, 3 }),

                    //    new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 0, 1, 3 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 0, 2, 3 }),
                    //    new LinearCalibrationFunctionMathNet(new [] { 1, 2, 3 }),

                    //    new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2, 3 }),
                    //};

                    int prevGoodIds;
                    int currentGoodIDs = 0;
                    Tolerance prevPrecursorMassTolerance = null;
                    Tolerance prevProductMassTolerance = null;
                    CalibrationResults res;
                    do
                    {
                        Console.WriteLine("SearchRound: " + searchRound);
                        prevGoodIds = currentGoodIDs;

                        List<Psm> goodIdentifications = GetGoodIdentifications(myMsDataFile, searchMode, currentDataFile, variableModifications, fixedModifications, proteinList, taskId);

                        currentGoodIDs = goodIdentifications.Count;
                        Console.WriteLine("currentGoodIDs: " + currentGoodIDs);

                        if (prevGoodIds >= currentGoodIDs)
                        {
                            Console.WriteLine("Done with linear!");
                            break;
                        }

                        if (CalibrationParameters.WriteIntermediateFiles)
                        {
                            ms1Action = (List<LabeledMs1DataPoint> theList) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + searchRound, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                            ms2Action = (List<LabeledMs2DataPoint> theList) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + searchRound, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                        }

                        res = (CalibrationResults)new CalibrationEngine(myMsDataFile, CommonParameters, goodIdentifications, ms1Action, ms2Action, CalibrationParameters, linearLearners, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                        prevPrecursorMassTolerance = CalibrationParameters.PrecursorMassTolerance;
                        prevProductMassTolerance = CommonParameters.ProductMassTolerance;

                        CalibrationParameters.PrecursorMassTolerance = new PpmTolerance(Math.Max(Math.Abs(res.ms1Info.Item1 + 5 * res.ms1Info.Item2), Math.Abs(res.ms1Info.Item1 - 5 * res.ms1Info.Item2)));
                        searchMode = new SinglePpmAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
                        CommonParameters.ProductMassTolerance = new PpmTolerance(Math.Max(Math.Abs(res.ms2Info.Item1 + 5 * res.ms2Info.Item2), Math.Abs(res.ms2Info.Item1 - 5 * res.ms2Info.Item2)));

                        Console.WriteLine("PrecursorMassTolerance: " + CalibrationParameters.PrecursorMassTolerance);
                        Console.WriteLine("ProductMassTolerance: " + CommonParameters.ProductMassTolerance);

                        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "-LinearCalibSearchRound" + searchRound + ".mzML"), false);

                        searchRound++;
                    } while (true);

                    CalibrationParameters.PrecursorMassTolerance = prevPrecursorMassTolerance;
                    searchMode = new SinglePpmAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
                    CommonParameters.ProductMassTolerance = prevProductMassTolerance;

                    Console.WriteLine("PrecursorMassTolerance: " + CalibrationParameters.PrecursorMassTolerance);
                    Console.WriteLine("ProductMassTolerance: " + CommonParameters.ProductMassTolerance);

                    myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "-LinearCalibSearchRound" + (searchRound - 2) + ".mzML"));

                    //IIndexedLearner<double>[] forEnsemble = new IIndexedLearner<double>[]
                    //{
                    //    new RegressionAbsoluteLossGradientBoostLearner(runParallel: false),
                    //    new RegressionAbsoluteLossGradientBoostLearner(runParallel: false,iterations:1000),
                    //    new RegressionAbsoluteLossGradientBoostLearner(runParallel: false,maximumTreeDepth:6),
                    //    new RegressionAbsoluteLossGradientBoostLearner(runParallel: false,maximumTreeDepth:6,iterations:1000),
                    //};

                    List<ILearner<double>> learners = new List<ILearner<double>>()
                    {
                        new LinearCalibrationFunctionMathNet(new int[] { }),

                        new LinearCalibrationFunctionMathNet(new [] { 0 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1 }),
                        new LinearCalibrationFunctionMathNet(new [] { 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 3 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 2, 3 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 2, 3 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 2, 3 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2, 3 }),

                        //new RegressionAbsoluteLossGradientBoostLearner(runParallel: false),
                        //new RegressionAbsoluteLossGradientBoostLearner(runParallel: false,iterations:1000),
                        //new RegressionAbsoluteLossGradientBoostLearner(runParallel: false,maximumTreeDepth:10),
                        new RegressionAbsoluteLossGradientBoostLearner(runParallel: false,maximumTreeDepth:10,iterations:1000),

                        //new RegressionBackwardEliminationModelSelectingEnsembleLearner(forEnsemble,2),
                        //new RegressionRandomModelSelectingEnsembleLearner(forEnsemble,2),
                        //new RegressionForwardSearchModelSelectingEnsembleLearner(forEnsemble,2),
                    };

                    List<Psm> finalIdentifications = GetGoodIdentifications(myMsDataFile, searchMode, currentDataFile, variableModifications, fixedModifications, proteinList, taskId);

                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        ms1Action = (List<LabeledMs1DataPoint> theList) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "beforeNonlinear", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                        ms2Action = (List<LabeledMs2DataPoint> theList) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "beforeNonlinear", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                    }

                    res = (CalibrationResults)new CalibrationEngine(myMsDataFile, CommonParameters, finalIdentifications, ms1Action, ms2Action, CalibrationParameters, learners, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                    CalibrationParameters.PrecursorMassTolerance = new PpmTolerance(Math.Max(Math.Abs(res.ms1Info.Item1 + 5 * res.ms1Info.Item2), Math.Abs(res.ms1Info.Item1 - 5 * res.ms1Info.Item2)));
                    searchMode = new SinglePpmAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
                    CommonParameters.ProductMassTolerance = new PpmTolerance(Math.Max(Math.Abs(res.ms2Info.Item1 + 5 * res.ms2Info.Item2), Math.Abs(res.ms2Info.Item1 - 5 * res.ms2Info.Item2)));

                    Console.WriteLine("PrecursorMassTolerance: " + CalibrationParameters.PrecursorMassTolerance);
                    Console.WriteLine("ProductMassTolerance: " + CommonParameters.ProductMassTolerance);

                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "-NonLinearCalib.mzML"), false);

                    finalIdentifications = GetGoodIdentifications(myMsDataFile, searchMode, currentDataFile, variableModifications, fixedModifications, proteinList, taskId);

                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        ms1Action = (List<LabeledMs1DataPoint> theList) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "beforeNonlinear2", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                        ms2Action = (List<LabeledMs2DataPoint> theList) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "beforeNonlinear2", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                    }
                    res = (CalibrationResults)new CalibrationEngine(myMsDataFile, CommonParameters, finalIdentifications, ms1Action, ms2Action, CalibrationParameters, learners, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                    CalibrationParameters.PrecursorMassTolerance = new PpmTolerance(Math.Max(Math.Abs(res.ms1Info.Item1 + 5 * res.ms1Info.Item2), Math.Abs(res.ms1Info.Item1 - 5 * res.ms1Info.Item2)));
                    searchMode = new SinglePpmAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
                    CommonParameters.ProductMassTolerance = new PpmTolerance(Math.Max(Math.Abs(res.ms2Info.Item1 + 5 * res.ms2Info.Item2), Math.Abs(res.ms2Info.Item1 - 5 * res.ms2Info.Item2)));

                    Console.WriteLine("PrecursorMassTolerance: " + CalibrationParameters.PrecursorMassTolerance);
                    Console.WriteLine("ProductMassTolerance: " + CommonParameters.ProductMassTolerance);

                    finalIdentifications = GetGoodIdentifications(myMsDataFile, searchMode, currentDataFile, variableModifications, fixedModifications, proteinList, taskId);

                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "-NonLinearCalib2.mzML"), false);

                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        ms1Action = (List<LabeledMs1DataPoint> theList) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "beforeNonlinear3", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                        ms2Action = (List<LabeledMs2DataPoint> theList) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "beforeNonlinear3", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                    }
                    res = (CalibrationResults)new CalibrationEngine(myMsDataFile, CommonParameters, finalIdentifications, ms1Action, ms2Action, CalibrationParameters, learners, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                    CalibrationParameters.PrecursorMassTolerance = new PpmTolerance(Math.Max(Math.Abs(res.ms1Info.Item1 + 5 * res.ms1Info.Item2), Math.Abs(res.ms1Info.Item1 - 5 * res.ms1Info.Item2)));
                    searchMode = new SinglePpmAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
                    CommonParameters.ProductMassTolerance = new PpmTolerance(Math.Max(Math.Abs(res.ms2Info.Item1 + 5 * res.ms2Info.Item2), Math.Abs(res.ms2Info.Item1 - 5 * res.ms2Info.Item2)));

                    Console.WriteLine("PrecursorMassTolerance: " + CalibrationParameters.PrecursorMassTolerance);
                    Console.WriteLine("ProductMassTolerance: " + CommonParameters.ProductMassTolerance);

                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "-NonLinearCalib3.mzML"), false);

                    Warn("");
                    Warn("");
                }
                );

            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

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