using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using SharpLearning.Common.Interfaces;
using SharpLearning.GradientBoost.Learners;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

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

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList)
        {
            myTaskResults = new MyTaskResults(this)
            {
                newSpectra = new List<string>()
            };

            Status("Loading modifications...", new List<string> { taskId, "Individual Spectra Files" });
            List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> localizeableModifications;
            if (CommonParameters.LocalizeAll)
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            Status("Loading proteins...", new List<string> { taskId, "Individual Spectra Files" });
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, UsefulProteomicsDatabases.DecoyType.Reverse, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> um)).ToList();

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
            Status("Calibrating...", new List<string> { taskId, "Individual Spectra Files" });

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var currentDataFile = currentRawFileList[spectraFileIndex];
                ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                lock (lock1) // Lock because reading is sequential
                {
                    if (Path.GetExtension(currentDataFile).Equals(".mzML", StringComparison.OrdinalIgnoreCase))
                        myMsDataFile = Mzml.LoadAllStaticData(currentDataFile);
                    else
#if NETFRAMEWORK
                        myMsDataFile = ThermoStaticData.LoadAllStaticData(currentDataFile);
#else
                        {
                            Warn("No capability for reading " + currentDataFile);
                            return;
                        }
#endif
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

                var finalLearners = new List<ILearner<double>>
                    {
                        new RegressionAbsoluteLossGradientBoostLearner(iterations: 1000),
                        new RegressionAbsoluteLossGradientBoostLearner(maximumTreeDepth: 6, iterations: 1000),
                        new RegressionAbsoluteLossGradientBoostLearner(maximumTreeDepth: 9, iterations: 1000),
                    };

                (int count, DataPointAquisitionResults datapointAcquisitionResult, Tolerance precTol, Tolerance prodTol) = GetDataAcquisitionResultsAndAppropriateTolerances(myMsDataFile, currentDataFile, variableModifications, fixedModifications, proteinList, taskId, combinedParams, combinedParams.PrecursorMassTolerance, combinedParams.ProductMassTolerance);

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
                if (CalibrationParameters.WriteIntermediateFiles)
                {
                    WriteMs1DataPoints(datapointAcquisitionResult.Ms1List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "init", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                    WriteMs2DataPoints(datapointAcquisitionResult.Ms2List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "init", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                }

                int prevCount;
                Tolerance prevPrecTol;
                Tolerance prevProdTol;
                string bestFilePath = null;
                var round = 1;
                do
                {
                    new CalibrationEngine(myMsDataFile, datapointAcquisitionResult, initLearners, "mz", new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                    prevCount = count;
                    prevPrecTol = precTol;
                    prevProdTol = prodTol;

                    (count, datapointAcquisitionResult, precTol, prodTol) = GetDataAcquisitionResultsAndAppropriateTolerances(myMsDataFile, currentDataFile, variableModifications, fixedModifications, proteinList, taskId, combinedParams, precTol, prodTol);

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

                    if (round >= 3 && !ImprovGlobal(prevPrecTol, prevProdTol, prevCount, count, precTol, prodTol))
                    {
                        break;
                    }
                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        WriteMs1DataPoints(datapointAcquisitionResult.Ms1List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "alignment", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                        WriteMs2DataPoints(datapointAcquisitionResult.Ms2List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "alignment", new List<string> { taskId, "Individual Spectra Files", currentDataFile });

                        var intermediateFilePath = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "alignment.mzml");
                        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, intermediateFilePath, false);
                        SucessfullyFinishedWritingFile(intermediateFilePath, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                    }

                    bestFilePath = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "-calib.mzml");
                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, bestFilePath, false);

                    round++;
                } while (true);

                myMsDataFile = Mzml.LoadAllStaticData(bestFilePath);

                do
                {
                    new CalibrationEngine(myMsDataFile, datapointAcquisitionResult, mzSepLearners, "mzRtTicInj", new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

                    prevCount = count;
                    prevPrecTol = precTol;
                    prevProdTol = prodTol;

                    (count, datapointAcquisitionResult, precTol, prodTol) = GetDataAcquisitionResultsAndAppropriateTolerances(myMsDataFile, currentDataFile, variableModifications, fixedModifications, proteinList, taskId, combinedParams, precTol, prodTol);

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

                    if (!ImprovGlobal(prevPrecTol, prevProdTol, prevCount, count, precTol, prodTol))
                        break;

                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        WriteMs1DataPoints(datapointAcquisitionResult.Ms1List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "inter-scan", new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                        WriteMs2DataPoints(datapointAcquisitionResult.Ms2List, OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "inter-scan", new List<string> { taskId, "Individual Spectra Files", currentDataFile });

                        var intermediateFilePath = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "round" + round + "inter-scan.mzml");
                        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, intermediateFilePath, false);
                        SucessfullyFinishedWritingFile(intermediateFilePath, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                    }

                    bestFilePath = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "-calib.mzml");
                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, bestFilePath, false);

                    round++;
                } while (true);

                SucessfullyFinishedWritingFile(bestFilePath, new List<string> { taskId, "Individual Spectra Files", currentDataFile });

                {
                    var tomlFileName = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(currentDataFile) + "-calib.toml");
                    FileSpecificTolerances f = new FileSpecificTolerances
                    {
                        PrecursorMassTolerance = precTol,
                        ProductMassTolerance = prodTol
                    };
                    Toml.WriteFile(f, tomlFileName, tomlConfig);
                    SucessfullyFinishedWritingFile(tomlFileName, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
                }

                myTaskResults.newSpectra.Add(bestFilePath);
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", currentDataFile }));
            });
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private bool ImprovGlobal(Tolerance prevPrecTol, Tolerance prevProdTol, int prevCount, int count, Tolerance currentPrecTol, Tolerance currentProdTol)
        {
            if (count > prevCount)
                return true;

            var precRatio = currentPrecTol.Value / prevPrecTol.Value;
            var prodRatio = currentProdTol.Value / prevProdTol.Value;

            if (count == prevCount)
                return precRatio + prodRatio < 2; // Take any improvement in ratios

            var countRatio = (double)count / prevCount;
            return countRatio > 0.9 && precRatio + prodRatio < 1.8;
        }

        private (int, DataPointAquisitionResults, Tolerance, Tolerance) GetDataAcquisitionResultsAndAppropriateTolerances(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, string currentDataFile, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, string taskId, ICommonParameters combinedParameters, Tolerance initPrecTol, Tolerance initProdTol)
        {
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

            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];

            List<ProductType> lp = new List<ProductType>();
            if (combinedParameters.BIons)
                lp.Add(ProductType.B);
            if (combinedParameters.YIons)
                lp.Add(ProductType.Y);
            if (combinedParameters.CIons)
                lp.Add(ProductType.C);
            if (combinedParameters.ZdotIons)
                lp.Add(ProductType.Zdot);

            Log("Searching with searchMode: " + searchMode, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
            Log("Searching with productMassTolerance: " + initProdTol, new List<string> { taskId, "Individual Spectra Files", currentDataFile });

            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, lp, searchMode, false, combinedParameters, initProdTol, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

            List<Psm> allPsms = allPsmsArray.ToList();

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = ((SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, lp, new List<DigestionParams> { combinedParameters.DigestionParams }, combinedParameters.ReportAllAmbiguity, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run()).CompactPeptideToProteinPeptideMatching;

            foreach (var huh in allPsms)
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

            allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

            new FdrAnalysisEngine(allPsms, searchMode.NumNotches, false, new List<string> { taskId, "Individual Spectra Files", currentDataFile }).Run();

            List<Psm> goodIdentifications = allPsms.Where(b => b.FdrInfo.QValue < 0.01 && !b.IsDecoy).ToList();

            if (!goodIdentifications.Any())
            {
                Warn("No good identifications!");
                return (0, null, null, null);
            }

            // Store

            DataPointAquisitionResults bestResult = new DataPointAquisitionResults(null, new List<LabeledMs1DataPoint>(), new List<LabeledMs2DataPoint>());

            Tolerance bestPrecursorMassToleranceForDatapointAcquisition = initPrecTol;
            Tolerance bestProductMassToleranceForDatapointAcquisition = initProdTol;

            Tolerance testPrecursorMassToleranceForDatapointAcquisition = initPrecTol;
            Tolerance testProductMassToleranceForDatapointAcquisition = initProdTol;

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
                {
                    Warn("currentResult.Ms1List.Count = " + currentResult.Ms1List.Count + " currentResult.Ms2List.Count = " + currentResult.Ms2List.Count);
                    return (0, null, null, null);
                }

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

            Log("Returning precTol:" + bestPrecursorMassToleranceForDatapointAcquisition, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
            Log("Returning prodTol:" + bestProductMassToleranceForDatapointAcquisition, new List<string> { taskId, "Individual Spectra Files", currentDataFile });
            return (goodIdentifications.Count, bestResult, bestPrecursorMassToleranceForDatapointAcquisition, bestProductMassToleranceForDatapointAcquisition);
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