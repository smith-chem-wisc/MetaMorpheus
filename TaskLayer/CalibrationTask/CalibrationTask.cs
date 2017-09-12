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
using SharpLearning.RandomForest.Learners;
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
                output.WriteLine(LabeledMs2DataPoint.TabSeparatedHeaderForMs1 + "\t" + Psm.GetTabSeparatedHeader());
                foreach (var dp in items)
                {
                    output.Write(string.Join("\t", dp.Inputs));
                    output.Write("\t" + dp.Label);
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
                output.WriteLine(LabeledMs1DataPoint.TabSeparatedHeaderForMs1 + "\t" + Psm.GetTabSeparatedHeader());
                foreach (var dp in items)
                {
                    output.Write(string.Join("\t", dp.Inputs));
                    output.Write("\t" + dp.Label);
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
            var searchModes = searchMode;

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

            List<ProductType> lp = new List<ProductType>();
            FragmentTypes fragmentTypesForCalibration = FragmentTypes.None;
            if (CommonParameters.BIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.b;
                lp.Add(ProductType.B);
            }
            if (CommonParameters.YIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.y;
                lp.Add(ProductType.Y);
            }
            if (CommonParameters.CIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.c;
                lp.Add(ProductType.C);
            }
            if (CommonParameters.ZdotIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.zdot;
                lp.Add(ProductType.Zdot);
            }

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
            proseCreatedWhileRunning.Append("parent mass tolerance(s) = {" + searchModes.ToProseString() + "}; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            object lock1 = new object();
            ParallelOptions parallelOptions = new ParallelOptions();
            if (CommonParameters.MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxDegreeOfParallelism.Value;
            Status("Calibrating...", new List<string> { taskId });
            TerminusType terminusType = ProductTypeToTerminusType.IdentifyTerminusType(lp);

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    var origDataFile = currentRawFileList[spectraFileIndex];

                    List<ILearner<double>> learners = new List<ILearner<double>>()
                    {
                        new RegressionAbsoluteLossGradientBoostLearner(),
                        new RegressionHuberLossGradientBoostLearner(),
                        new RegressionSquareLossGradientBoostLearner(),
                        new RegressionExtremelyRandomizedTreesLearner(),

                        new LinearCalibrationFunctionMathNet(new TransformFunction(b => new double[] { b[0] }, 1)),
                        new LinearCalibrationFunctionMathNet(new TransformFunction(b => new double[] { b[0], b[1] }, 2)),
                        //new LinearCalibrationFunctionMathNet(new TransformFunction(b => new double[] { b[0], b[1], Math.Log(b[2]), Math.Log(b[3]) }, 4)),

                        //new RandomForestCalibrationFunction(40, 10), //29-32

                        //new RegressionGradientBoostLearner(),
                        //new RegressionRandomForestLearner(),
                        //new IdentityCalibrationFunction(),
                        //new ConstantCalibrationFunction(),
                        //new RegressionAdaBoostLearner(),
                        //new RegressionDecisionTreeLearner(),
                        //new RegressionQuantileLossGradientBoostLearner(),
                        };

                    List<(Func<DataPointAquisitionResults, DataPointAquisitionResults, bool>, string)> ContinueLoops = new List<(Func<DataPointAquisitionResults, DataPointAquisitionResults, bool>, string)>
                    {
                        ((DataPointAquisitionResults a, DataPointAquisitionResults b) => b.Count > a.Count, "countIncrease"),
                        ((DataPointAquisitionResults a, DataPointAquisitionResults b) => false, "single"),
                    };

                    List<bool> DoLinearMzSeparately = new List<bool>
                    {
                        true,
                        false,
                    };

                    List<ILearner<double>> DoFirst = new List<ILearner<double>>()
                    {
                        new LinearCalibrationFunctionMathNet(new TransformFunction(b => new double[] { b[0] }, 1)),
                        new LinearCalibrationFunctionMathNet(new TransformFunction(b => new double[] { b[0], b[1] }, 2)),
                    };

                    int calibIndex = 0;
                    foreach (var CalibrationSetting in GetAllCalibSettings(learners, ContinueLoops, DoLinearMzSeparately, DoFirst))
                    {
                        calibIndex++;
                        //IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);
                        IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = Mzml.LoadAllStaticData(origDataFile);

                        {
                            var listOfSortedms2Scans = GetMs2Scans(myMsDataFile, origDataFile, CommonParameters.DoPrecursorDeconvolution, CommonParameters.UseProvidedPrecursorInfo, CommonParameters.DeconvolutionIntensityRatio, CommonParameters.DeconvolutionMaxAssumedChargeState, CommonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

                            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];

                            var searchEngine = new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, lp, searchModes, false, CommonParameters, new List<string>());
                            searchEngine.Run();
                            List<Psm> allPsms = allPsmsArray.ToList();
                            // Group and order psms

                            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, ProductTypeToTerminusType.IdentifyTerminusType(lp), new List<DigestionParams> { CommonParameters.DigestionParams }, new List<string>());

                            var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;

                            foreach (var huh in allPsms)
                                if (huh != null && huh.MostProbableProteinInfo == null)
                                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

                            allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

                            // Run FDR analysis on single file results
                            new FdrAnalysisEngine(allPsms, searchModes, new List<string>()).Run();

                            var localizationEngine = new LocalizationEngine(allPsms, lp, myMsDataFile, CommonParameters.ProductMassTolerance, new List<string>(), false);
                            localizationEngine.Run();

                            //Console.WriteLine("scan= " + allPsms.First().ScanNumber + " fullSeq= " + allPsms.First().FullSequence + " daError= " + (allPsms.First().ScanPrecursorMass - allPsms.First().PeptideMonisotopicMass) + " ppmError= " + ((allPsms.First().ScanPrecursorMass - allPsms.First().PeptideMonisotopicMass) / allPsms.First().PeptideMonisotopicMass * 1e6));
                            //Console.Write('\t' + "[");
                            //foreach (var kvp in allPsms.First().MatchedIonDictOnlyMatches)
                            //    Console.Write("[" + string.Join(",", kvp.Value.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "];");
                            //Console.WriteLine("]");
                            //Console.Write('\t' + "[");
                            //foreach (var kvp in allPsms.First().ProductMassErrorDa)
                            //    Console.Write("[" + string.Join(",", kvp.Value.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
                            //Console.WriteLine("]");
                            //Console.Write('\t' + "[");
                            //foreach (var kvp in allPsms.First().ProductMassErrorPpm)
                            //    Console.Write("[" + string.Join(",", kvp.Value.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
                            //Console.WriteLine("]");

                            var goodIdentifications = allPsms.Where(b => b.FdrInfo.QValue < 0.01 && !b.IsDecoy).ToList();

                            Action<List<LabeledMs1DataPoint>, string> ms1Action = (List<LabeledMs1DataPoint> theList, string s) => {; };
                            Action<List<LabeledMs2DataPoint>, string> ms2Action = (List<LabeledMs2DataPoint> theList, string s) => {; };
                            if (CalibrationParameters.WriteIntermediateFiles)
                            {
                                ms1Action = (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                                ms2Action = (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                            }

                            int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
                            int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
                            int numFragmentsNeededForEveryIdentification = 10;
                            Random rnd = new Random(spectraFileIndex);

                            MetaMorpheusEngineResults theResult = new CalibrationEngine(myMsDataFile, CommonParameters.ProductMassTolerance, goodIdentifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, CalibrationParameters.PrecursorMassTolerance, fragmentTypesForCalibration, ms1Action, ms2Action, false, rnd, new List<string>(), 10, 40, CalibrationSetting).Run();

                            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, Path.Combine(OutputFolder, "yah" + calibIndex + ".mzML"), false);
                        }

                        //#region aftercalib

                        //{
                        //    var ms2ScansAfterCalib = GetMs2Scans(myMsDataFile, origDataFile, CommonParameters.DoPrecursorDeconvolution, CommonParameters.UseProvidedPrecursorInfo, CommonParameters.DeconvolutionIntensityRatio, CommonParameters.DeconvolutionMaxAssumedChargeState, CommonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                        //    Psm[] allPsmsArray = new Psm[ms2ScansAfterCalib.Length];
                        //    new ClassicSearchEngine(allPsmsArray, ms2ScansAfterCalib, variableModifications, fixedModifications, proteinList, lp, searchModes, false, CommonParameters, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                        //    List<Psm> allPsms = allPsmsArray.ToList();

                        //    SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngineTest = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, terminusType, new List<DigestionParams> { CommonParameters.DigestionParams }, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        //    var resTest = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngineTest.Run();
                        //    Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatchingTest = resTest.CompactPeptideToProteinPeptideMatching;

                        //    foreach (var huh in allPsms)
                        //        if (huh != null && huh.MostProbableProteinInfo == null)
                        //            huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatchingTest);

                        //    // Group and order psms
                        //    allPsms = allPsmsArray.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

                        //    new FdrAnalysisEngine(allPsms, searchModes, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                        //    var localizationEngine = new LocalizationEngine(allPsms, lp, myMsDataFile, CommonParameters.ProductMassTolerance, new List<string>(), false);
                        //    localizationEngine.Run();

                        //    Console.WriteLine("scan= " + allPsms.First().ScanNumber + " fullSeq= " + allPsms.First().FullSequence + " daError= " + (allPsms.First().ScanPrecursorMass - allPsms.First().PeptideMonisotopicMass) + " ppmError= " + ((allPsms.First().ScanPrecursorMass - allPsms.First().PeptideMonisotopicMass) / allPsms.First().PeptideMonisotopicMass * 1e6));
                        //    Console.Write('\t' + "[");
                        //    foreach (var kvp in allPsms.First().MatchedIonDictOnlyMatches)
                        //        Console.Write("[" + string.Join(",", kvp.Value.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "];");
                        //    Console.WriteLine("]");
                        //    Console.Write('\t' + "[");
                        //    foreach (var kvp in allPsms.First().ProductMassErrorDa)
                        //        Console.Write("[" + string.Join(",", kvp.Value.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
                        //    Console.WriteLine("]");
                        //    Console.Write('\t' + "[");
                        //    foreach (var kvp in allPsms.First().ProductMassErrorPpm)
                        //        Console.Write("[" + string.Join(",", kvp.Value.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
                        //    Console.WriteLine("]");

                        //    Warn(CalibrationSetting.ToString());
                        //    Warn("psms: " + allPsms.Count(a => a.FdrInfo.QValue <= .01 && !a.IsDecoy));
                        //}

                        Warn("");
                        Warn("");
                    }
                }
        );

            ReportProgress(new ProgressEventArgs(100, "Done!",

                    new List<string> { taskId, "Individual Spectra Files" }));

            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private IEnumerable<CalibrationSetting> GetAllCalibSettings(
            List<ILearner<double>> learners,
            List<(Func<DataPointAquisitionResults, DataPointAquisitionResults, bool>, string)> continueLoops,
            List<bool> DoLinearMzSeparately,
            List<ILearner<double>> DoFirst)
        {
            foreach (var ok2 in continueLoops)
                foreach (var ok3 in DoLinearMzSeparately)
                {
                    if (ok3)
                    {
                        yield return new CalibrationSetting
                        {
                            Learners = learners.Select(b => new SeparateMzLearner(b) as ILearner<double>).ToList(),
                            ContinueLoop = ok2,
                            DoFirst = null,
                        };
                        yield return new CalibrationSetting
                        {
                            Learners = learners.Select(b => new SeparateMzLearner(b) as ILearner<double>).ToList(),
                            ContinueLoop = ok2,
                            DoFirst = DoFirst,
                        };
                    }
                    else
                    {
                        yield return new CalibrationSetting
                        {
                            Learners = learners,
                            ContinueLoop = ok2,
                            DoFirst = null,
                        };
                        yield return new CalibrationSetting
                        {
                            Learners = learners,
                            ContinueLoop = ok2,
                            DoFirst = DoFirst,
                        };
                    }
                }
        }

        #endregion Private Methods
    }
}