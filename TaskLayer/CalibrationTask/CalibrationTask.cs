using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using IO.MzML;
using IO.Thermo;
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
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, true, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> um)).ToList();

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
            if (CommonParameters.MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxDegreeOfParallelism.Value;
            Status("Calibrating...", new List<string> { taskId });

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    var currentDataFile = currentRawFileList[spectraFileIndex];

                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                    lock (lock1) // Lock because reading is sequential
                    {
                        if (Path.GetExtension(currentDataFile).Equals(".mzML", StringComparison.InvariantCultureIgnoreCase))
                            myMsDataFile = Mzml.LoadAllStaticData(currentDataFile);
                        else
                            myMsDataFile = ThermoStaticData.LoadAllStaticData(currentDataFile);
                    }

                    var initLearners = new List<ILearner<double>>
                    {
                        new MzMultiplier()
                    };

                    var mzSepLearners = new List<ILearner<double>>
                    {
                        new LinearCalibrationFunctionMathNet(new int[] { }),

                        new LinearCalibrationFunctionMathNet(new [] { 0 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1 }),
                        new LinearCalibrationFunctionMathNet(new [] { 2 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1 }),
                        new LinearCalibrationFunctionMathNet(new [] { 0, 2 }),
                        new LinearCalibrationFunctionMathNet(new [] { 1, 2 }),

                        new LinearCalibrationFunctionMathNet(new [] { 0, 1, 2 }),

                        new RegressionAbsoluteLossGradientBoostLearner(iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner(maximumTreeDepth:6, iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner(maximumTreeDepth:9, iterations:1000),
                    };

                    mzSepLearners = mzSepLearners.Select(b => new SeparateMzLearner(b) as ILearner<double>).ToList();

                    var learners = new List<ILearner<double>>
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

                        new RegressionAbsoluteLossGradientBoostLearner(iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner(maximumTreeDepth:6, iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner(maximumTreeDepth:9, iterations:1000),
                    };

                    var intLearners = new List<ILearner<double>>
                    {
                        new LinearCalibrationFunctionMathNet(new int[] { }),

                        new LinearCalibrationFunctionMathNet(new [] { 0 }),

                        new RegressionAbsoluteLossGradientBoostLearner(iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner(maximumTreeDepth:6, iterations:1000),
                        new RegressionAbsoluteLossGradientBoostLearner(maximumTreeDepth:9, iterations:1000),
                    };

                    (int count, DataPointAquisitionResults datapointAcquisitionResult) = GetDataAcquisitionResultsAndSetTolerances(myMsDataFile, currentDataFile, variableModifications, fixedModifications, proteinList, taskId);

                    if (datapointAcquisitionResult == null)
                    {
                        Warn("datapointAcquisitionResult is null");
                        return;
                    }
                    if (datapointAcquisitionResult.Ms1List.Count < 4 || datapointAcquisitionResult.Ms2List.Count < 4)
                    {
                        Warn("datapointAcquisitionResult.Ms1List.Count: " + datapointAcquisitionResult.Ms1List.Count);
                        Warn("datapointAcquisitionResult.Ms1List.Count: " + datapointAcquisitionResult.Ms1List.Count);
                        return;
                    }
                    WriteMs1DataPoints(datapointAcquisitionResult.Ms1List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "init", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                    WriteMs2DataPoints(datapointAcquisitionResult.Ms2List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "init", new List<string> { taskId, "Individual Spectra Files", currentDataFile });

                    int prevCount;
                    Tolerance prevPrecTol;
                    Tolerance prevProdTol;
                    string bestFilePath = null;
                    var round = 1;
                    do
                    {
                        new CalibrationEngine(myMsDataFile, datapointAcquisitionResult, initLearners, "mz", new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                        prevCount = count;
                        prevPrecTol = CommonParameters.PrecursorMassTolerance;
                        prevProdTol = CommonParameters.ProductMassTolerance;

                        (count, datapointAcquisitionResult) = GetDataAcquisitionResultsAndSetTolerances(myMsDataFile, currentDataFile, variableModifications, fixedModifications, proteinList, taskId);

                        if (datapointAcquisitionResult == null)
                        {
                            Warn("datapointAcquisitionResult is null");
                            return;
                        }
                        if (datapointAcquisitionResult.Ms1List.Count < 4 || datapointAcquisitionResult.Ms2List.Count < 4)
                        {
                            Warn("datapointAcquisitionResult.Ms1List.Count: " + datapointAcquisitionResult.Ms1List.Count);
                            Warn("datapointAcquisitionResult.Ms1List.Count: " + datapointAcquisitionResult.Ms1List.Count);
                            return;
                        }

                        if (round >= 3 && !ImprovGlobal(prevPrecTol, prevProdTol, prevCount, count))
                            break;

                        WriteMs1DataPoints(datapointAcquisitionResult.Ms1List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "alignment", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                        WriteMs2DataPoints(datapointAcquisitionResult.Ms2List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "alignment", new List<string> { taskId, "Individual Spectra Files", currentDataFile });

                        bestFilePath = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "alignment.mzml");

                        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, bestFilePath, false);
                        SucessfullyFinishedWritingFile(bestFilePath, new List<string> { taskId, "Individual Spectra Files", currentDataFile });

                        round++;
                    } while (true);

                    CommonParameters.PrecursorMassTolerance = prevPrecTol;
                    CommonParameters.ProductMassTolerance = prevProdTol;

                    myMsDataFile = Mzml.LoadAllStaticData(bestFilePath);

                    do
                    {
                        new CalibrationEngine(myMsDataFile, datapointAcquisitionResult, mzSepLearners, "mzRtTicInj", new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                        prevCount = count;
                        prevPrecTol = CommonParameters.PrecursorMassTolerance;
                        prevProdTol = CommonParameters.ProductMassTolerance;

                        (count, datapointAcquisitionResult) = GetDataAcquisitionResultsAndSetTolerances(myMsDataFile, currentDataFile, variableModifications, fixedModifications, proteinList, taskId);

                        if (datapointAcquisitionResult == null)
                        {
                            Warn("datapointAcquisitionResult is null");
                            return;
                        }
                        if (datapointAcquisitionResult.Ms1List.Count < 4 || datapointAcquisitionResult.Ms2List.Count < 4)
                        {
                            Warn("datapointAcquisitionResult.Ms1List.Count: " + datapointAcquisitionResult.Ms1List.Count);
                            Warn("datapointAcquisitionResult.Ms1List.Count: " + datapointAcquisitionResult.Ms1List.Count);
                            return;
                        }

                        if (!ImprovGlobal(prevPrecTol, prevProdTol, prevCount, count))
                            break;

                        WriteMs1DataPoints(datapointAcquisitionResult.Ms1List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "inter-scan", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                        WriteMs2DataPoints(datapointAcquisitionResult.Ms2List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "inter-scan", new List<string> { taskId, "Individual Spectra Files", currentDataFile });

                        bestFilePath = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "inter-scan.mzml");

                        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, bestFilePath, false);
                        SucessfullyFinishedWritingFile(bestFilePath, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                        round++;
                    } while (true);

                    //CommonParameters.PrecursorMassTolerance = prevPrecTol;
                    //CommonParameters.ProductMassTolerance = prevProdTol;

                    //myMsDataFile = Mzml.LoadAllStaticData(bestFilePath);

                    //do
                    //{
                    //    new CalibrationEngine(myMsDataFile, datapointAcquisitionResult, intLearners, "Int", new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                    //    prevCount = count;
                    //    prevPrecTol = CommonParameters.PrecursorMassTolerance;
                    //    prevProdTol = CommonParameters.ProductMassTolerance;

                    //    (count, datapointAcquisitionResult) = GetDataAcquisitionResultsAndSetTolerances(myMsDataFile, currentDataFile, variableModifications, fixedModifications, proteinList, taskId);

                    //    if (datapointAcquisitionResult == null)
                    //    {
                    //        Warn("datapointAcquisitionResult is null");
                    //        return;
                    //    }
                    //    if (datapointAcquisitionResult.Ms1List.Count < 4 || datapointAcquisitionResult.Ms2List.Count < 4)
                    //    {
                    //        Warn("datapointAcquisitionResult.Ms1List.Count: " + datapointAcquisitionResult.Ms1List.Count);
                    //        Warn("datapointAcquisitionResult.Ms1List.Count: " + datapointAcquisitionResult.Ms1List.Count);
                    //        return;
                    //    }

                    //    if (!ImprovGlobal(prevPrecTol, prevProdTol, prevCount, count))
                    //        break;

                    //    WriteMs1DataPoints(datapointAcquisitionResult.Ms1List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "final", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                    //    WriteMs2DataPoints(datapointAcquisitionResult.Ms2List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "final", new List<string> { taskId, "Individual Spectra Files", currentDataFile });

                    //    bestFilePath = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "final.mzml");

                    //    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, bestFilePath, false);
                    //    SucessfullyFinishedWritingFile(bestFilePath, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                    //    round++;
                    //} while (true);

                    //CommonParameters.PrecursorMassTolerance = prevPrecTol;
                    //CommonParameters.ProductMassTolerance = prevProdTol;

                    myTaskResults.newSpectra.Add(bestFilePath);
                }

                );

            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private bool ImprovGlobal(Tolerance prevPrecTol, Tolerance prevProdTol, int prevCount, int count)
        {
            if (count > prevCount)
                return true;

            Tolerance currentPrecTol = CommonParameters.PrecursorMassTolerance;
            Tolerance currentProdTol = CommonParameters.ProductMassTolerance;

            var precRatio = currentPrecTol.Value / prevPrecTol.Value;
            var prodRatio = currentProdTol.Value / prevProdTol.Value;

            if (count == prevCount)
                return precRatio + prodRatio < 2; // Take any improvement in ratios

            var countRatio = (double)count / prevCount;
            return countRatio > 0.9 && precRatio + prodRatio < 1.8;
        }

        private (int, DataPointAquisitionResults) GetDataAcquisitionResultsAndSetTolerances(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, string currentDataFile, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, string taskId)
        {
            MassDiffAcceptor searchMode;
            if (CommonParameters.PrecursorMassTolerance is PpmTolerance)
                searchMode = new SinglePpmAroundZeroSearchMode(CommonParameters.PrecursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode(CommonParameters.PrecursorMassTolerance.Value);

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

            if (!goodIdentifications.Any())
                return (0, null);

            // Store

            DataPointAquisitionResults bestResult = new DataPointAquisitionResults(null, new List<LabeledMs1DataPoint>(), new List<LabeledMs2DataPoint>());

            Tolerance bestPrecursorMassToleranceForDatapointAcquisition = CommonParameters.PrecursorMassTolerance;
            Tolerance bestProductMassToleranceForDatapointAcquisition = CommonParameters.ProductMassTolerance;

            Tolerance testPrecursorMassToleranceForDatapointAcquisition = CommonParameters.PrecursorMassTolerance;
            Tolerance testProductMassToleranceForDatapointAcquisition = CommonParameters.ProductMassTolerance;

            var round = 1;
            do
            {
                DataPointAquisitionResults currentResult = (DataPointAquisitionResults)new DataPointAcquisitionEngine(
                    goodIdentifications,
                    myMsDataFile,
                    testPrecursorMassToleranceForDatapointAcquisition,
                    testProductMassToleranceForDatapointAcquisition,
                    CalibrationParameters.NumFragmentsNeededForEveryIdentification,
                    CalibrationParameters.MinMS1IsotopicPeaksNeededForConfirmedIdentification,
                    CalibrationParameters.MinMS2IsotopicPeaksNeededForConfirmedIdentification,
                    fragmentTypesForCalibration,
                    new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                if (currentResult.Ms1List.Count == 0 || currentResult.Ms2List.Count == 0)
                    return (0, null);

                var computedPrecursorMassToleranceForDatapointAcquisition = new PpmTolerance(Math.Max(Math.Abs(currentResult.Ms1InfoPpm.Item1 + 6 * currentResult.Ms1InfoPpm.Item2), Math.Abs(currentResult.Ms1InfoPpm.Item1 - 6 * currentResult.Ms1InfoPpm.Item2)));
                var computedProductMassToleranceForDatapointAcquisition = new PpmTolerance(Math.Max(Math.Abs(currentResult.Ms2InfoPpm.Item1 + 6 * currentResult.Ms2InfoPpm.Item2), Math.Abs(currentResult.Ms2InfoPpm.Item1 - 6 * currentResult.Ms2InfoPpm.Item2)));

                bool ms1Worse = currentResult.Ms1List.Count <= bestResult.Ms1List.Count && ((currentResult.Ms1InfoPpm.Item2 / bestResult.Ms1InfoPpm.Item2) >= ((double)currentResult.Ms1List.Count / bestResult.Ms1List.Count));
                bool ms2Worse = currentResult.Ms2List.Count <= bestResult.Ms2List.Count && ((currentResult.Ms2InfoPpm.Item2 / bestResult.Ms2InfoPpm.Item2) >= ((double)currentResult.Ms2List.Count / bestResult.Ms2List.Count));

                if (ms1Worse || ms2Worse)
                    break;

                bestPrecursorMassToleranceForDatapointAcquisition = testPrecursorMassToleranceForDatapointAcquisition;
                bestProductMassToleranceForDatapointAcquisition = testProductMassToleranceForDatapointAcquisition;

                testPrecursorMassToleranceForDatapointAcquisition = computedPrecursorMassToleranceForDatapointAcquisition;
                testProductMassToleranceForDatapointAcquisition = computedProductMassToleranceForDatapointAcquisition;

                bestResult = currentResult;
                round++;
            } while (true);

            CommonParameters.PrecursorMassTolerance = bestPrecursorMassToleranceForDatapointAcquisition;
            CommonParameters.ProductMassTolerance = bestProductMassToleranceForDatapointAcquisition;

            return (goodIdentifications.Count, bestResult);
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