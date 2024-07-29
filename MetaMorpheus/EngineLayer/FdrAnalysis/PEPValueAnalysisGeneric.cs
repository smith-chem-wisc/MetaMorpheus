using Chemistry;
using EngineLayer.CrosslinkSearch;
using EngineLayer.FdrAnalysis;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using Microsoft.ML.Data;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics.Modifications;
using Omics;

namespace EngineLayer
{
    public static class PEP_Analysis_Cross_Validation
    {
        private static readonly double AbsoluteProbabilityThatDistinguishesPeptides = 0.05;
        private static Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();
        private static Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();
        private static Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();
        
        /// <summary>
        /// A dictionary which stores the chimeric ID string in the key and the number of chimeric identifications as the vale
        /// </summary>
        private static Dictionary<string, int> chimeraCountDictionary = new Dictionary<string, int>();
        public static bool UsePeptideLevelQValueForTraining = true;
        public static double QValueCutoff = 0.005;


        /// <summary>
        /// This method is used to compute the PEP values for all PSMs in a dataset. 
        /// </summary>
        /// <param name="psms"></param>
        /// <param name="searchType"></param>
        /// <param name="fileSpecificParameters"></param>
        /// <param name="outputFolder"></param>
        /// <returns></returns>
        public static string ComputePEPValuesForAllPSMsGeneric(List<SpectralMatch> psms, string searchType, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, string outputFolder)
        {
            string[] trainingVariables = PsmData.trainingInfos[searchType];

            //ensure that the order is always stable.
            psms = psms.OrderByDescending(p => p).ToList();
            List<int> allPeptideIndices = new List<int>();
            List<SpectralMatch> peptides = psms
                .GroupBy(b => b.FullSequence)
                .Select(b => b.FirstOrDefault()).ToList();
            List<int> countOfPeptidesInEachFile = peptides.GroupBy(b => b.FullFilePath).Select(b => b.Count()).ToList();
            bool allFilesContainPeptides = (countOfPeptidesInEachFile.Count == fileSpecificParameters.Count); //rare condition where each file has psms but some files don't have peptides. probably only happens in unit tests.
            QValueCutoff = fileSpecificParameters.Select(t => t.fileSpecificParameters.QValueCutoffForPepCalculation).Min();

            int chargeStateMode = 0;
            Dictionary<string, float> fileSpecificMedianFragmentMassErrors = new Dictionary<string, float>();
            if (peptides.Count() > 100 && allFilesContainPeptides)
            {
                foreach (var peptide in peptides)
                {
                    allPeptideIndices.Add(psms.IndexOf(peptide));
                }
                chargeStateMode = GetChargeStateMode(peptides);
                fileSpecificMedianFragmentMassErrors = GetFileSpecificMedianFragmentMassError(peptides);
            }
            else
            {
                //there are too few psms to do any meaningful training if we used only peptides. So, we will train using psms instead.
                UsePeptideLevelQValueForTraining = false;
                allPeptideIndices = Enumerable.Range(0, psms.Count).ToList();
                chargeStateMode = GetChargeStateMode(psms);
                fileSpecificMedianFragmentMassErrors = GetFileSpecificMedianFragmentMassError(psms);
            }
            

            //These two dictionaries contain the average and standard deviations of hydrophobicitys measured in 1 minute increments accross each raw
            //file separately. An individully measured hydrobophicty calculated for a specific PSM sequence is compared to these values by computing
            //the z-score. That z-score is used as a feature for machine learning.
            //Separate dictionaries are created for peptides with modifications because SSRcalc doesn't really do a good job predicting hyrophobicity

            //The first string in the dictionary is the filename
            //The value of the dictionary is another dictionary that profiles the hydrophobicity behavior.
            //Each key is a retention time rounded to the nearest minute.
            //The value Tuple is the average and standard deviation, respectively, of the predicted hydrophobicities of the observed peptides eluting at that rounded retention time.
            
            if (trainingVariables.Contains("HydrophobicityZScore"))
            {
                if (peptides.Count() > 100 && allFilesContainPeptides)
                {
                    fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = ComputeHydrophobicityValues(peptides, fileSpecificParameters, false);
                    fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = ComputeHydrophobicityValues(peptides, fileSpecificParameters, true);
                    fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE = ComputeMobilityValues(peptides, fileSpecificParameters);
                }
                else
                {
                    fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = ComputeHydrophobicityValues(psms, fileSpecificParameters, false);
                    fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = ComputeHydrophobicityValues(psms, fileSpecificParameters, true);
                    fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE = ComputeMobilityValues(psms, fileSpecificParameters);
                }
            }

            if (trainingVariables.Contains("ChimeraCount"))
                chimeraCountDictionary = psms.GroupBy(p => p.ChimeraIdString)
                    .ToDictionary(p => p.Key, p => p.Count());

            MLContext mlContext = new MLContext();

            //the number of groups used for cross-validation is hard-coded at four. Do not change this number without changes other areas of effected code.
            int numGroups = 4;
            if (psms.Count < 1000 || allPeptideIndices.Count < 500)
            {
                numGroups = 2;
            }
            List<int>[] psmGroupIndices = Get_PSM_Group_Indices(psms, numGroups);

            //the psms will be randomly divided. but then we want to make another array that just contains the subset of peptides that are in those psms. that way we don't compute pep using any peptides that were used in training.
            List<int>[] peptideGroupIndices = Get_Peptide_Group_Indices(psmGroupIndices, allPeptideIndices);
            IEnumerable<PsmData>[] PSMDataGroups = new IEnumerable<PsmData>[numGroups];
            
            for (int i = 0; i < numGroups; i++)
            {
                PSMDataGroups[i] = CreatePsmData(searchType, fileSpecificParameters, psms, peptideGroupIndices[i], fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode);
            }

            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features", numberOfTrees: 400);
            var pipeline = mlContext.Transforms.Concatenate("Features", trainingVariables)
                .Append(trainer);

            List<CalibratedBinaryClassificationMetrics> allMetrics = new List<CalibratedBinaryClassificationMetrics>();
            int sumOfAllAmbiguousPeptidesResolved = 0;

            bool allSetsContainPositiveAndNegativeTrainingExamples = true;
            int groupNumber = 0;
            while (allSetsContainPositiveAndNegativeTrainingExamples == true && groupNumber < numGroups)
            {
                if (PSMDataGroups[groupNumber].Where(p => p.Label == true).Count() == 0 || PSMDataGroups[groupNumber].Where(p => p.Label == false).Count() == 0)
                {
                    allSetsContainPositiveAndNegativeTrainingExamples = false;
                }
                groupNumber++;
            }

            if (allSetsContainPositiveAndNegativeTrainingExamples)
            {
                for (int groupIndexNumber = 0; groupIndexNumber < numGroups; groupIndexNumber++)
                {
                    List<int> allGroupIndexes = Enumerable.Range(0, numGroups).ToList();
                    allGroupIndexes.RemoveAt(groupIndexNumber);

                    //concat doesn't work in a loop, therefore I had to hard code the concat to group 3 out of 4 lists. if the const int numGroups value is changed, then the concat has to be changed accordingly.
                    IDataView dataView = mlContext.Data.LoadFromEnumerable(PSMDataGroups[allGroupIndexes[0]]);
                    if (numGroups > 2)
                    {
                        dataView = mlContext.Data.LoadFromEnumerable(PSMDataGroups[allGroupIndexes[0]].Concat(PSMDataGroups[allGroupIndexes[1]].Concat(PSMDataGroups[allGroupIndexes[2]])));
                    }
                    trainedModels[groupIndexNumber] = pipeline.Fit(dataView);
                    var myPredictions = trainedModels[groupIndexNumber].Transform(mlContext.Data.LoadFromEnumerable(PSMDataGroups[groupIndexNumber]));
                    CalibratedBinaryClassificationMetrics metrics = mlContext.BinaryClassification.Evaluate(data: myPredictions, labelColumnName: "Label", scoreColumnName: "Score");

                    //Parallel operation of the following code requires the method to be stored and then read, once for each thread
                    //if not output directory is specified, the model cannot be stored, and we must force single-threaded operation
                    if (outputFolder != null)
                    {
                        mlContext.Model.Save(trainedModels[groupIndexNumber], dataView.Schema, Path.Combine(outputFolder, "model.zip"));
                    }

                    //model is trained on peptides but here we can use that to compute PEP for all PSMs
                    int ambiguousPeptidesResolved = Compute_PSM_PEP(psms, psmGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], searchType, fileSpecificParameters, fileSpecificMedianFragmentMassErrors, chargeStateMode, outputFolder);

                    allMetrics.Add(metrics);
                    sumOfAllAmbiguousPeptidesResolved += ambiguousPeptidesResolved;
                }

                return AggregateMetricsForOutput(allMetrics, sumOfAllAmbiguousPeptidesResolved);
            }
            else
            {
                return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
            }
        }

        private static List<int>[] Get_Peptide_Group_Indices(List<int>[] psmGroupIndices, List<int> allPeptideIndices)
        {
            List<int>[] peptideGroupIndices = new List<int>[psmGroupIndices.Length];
            for (int i = 0; i < psmGroupIndices.Length; i++)
            {
                peptideGroupIndices[i] = psmGroupIndices[i].Intersect(allPeptideIndices).ToList();
            }
            return peptideGroupIndices;
        }

        public static string AggregateMetricsForOutput(List<CalibratedBinaryClassificationMetrics> allMetrics, int sumOfAllAmbiguousPeptidesResolved)
        {
            List<double> accuracy = allMetrics.Select(m => m.Accuracy).ToList();
            List<double> areaUnderRocCurve = allMetrics.Select(m => m.AreaUnderRocCurve).ToList();
            List<double> areaUnderPrecisionRecallCurve = allMetrics.Select(m => m.AreaUnderPrecisionRecallCurve).ToList();
            List<double> F1Score = allMetrics.Select(m => m.F1Score).ToList();
            List<double> logLoss = allMetrics.Select(m => m.LogLoss).ToList();
            List<double> logLossReduction = allMetrics.Select(m => m.LogLossReduction).ToList();
            List<double> positivePrecision = allMetrics.Select(m => m.PositivePrecision).ToList();
            List<double> positiveRecall = allMetrics.Select(m => m.PositiveRecall).ToList();
            List<double> negativePrecision = allMetrics.Select(m => m.NegativePrecision).ToList();
            List<double> negativeRecall = allMetrics.Select(m => m.NegativeRecall).ToList();

            // log-loss can stochastically take on a value of infinity.
            // correspondingly, log-loss reduction can be negative infinity.
            // when this happens for one or more of the metrics, it can lead to uninformative numbers.
            // so, unless they are all infinite, we remove them from the average. If they are all infinite, we report that.

            logLoss.RemoveAll(x => x == Double.PositiveInfinity);
            logLossReduction.RemoveAll(x => x == Double.NegativeInfinity);

            double logLossAverage = Double.PositiveInfinity;
            double logLossReductionAverage = Double.NegativeInfinity;

            if ((logLoss != null) && (logLoss.Any()))
            {
                logLossAverage = logLoss.Average();
            }

            if ((logLossReduction != null) && (logLossReduction.Any()))
            {
                logLossReductionAverage = logLossReduction.Average();
            }

            StringBuilder s = new StringBuilder();
            s.AppendLine();
            s.AppendLine("************************************************************");
            s.AppendLine("*       Metrics for Determination of PEP Using Binary Classification      ");
            s.AppendLine("*-----------------------------------------------------------");
            s.AppendLine("*       Accuracy:  " + accuracy.Average().ToString());
            s.AppendLine("*       Area Under Curve:  " + areaUnderRocCurve.Average().ToString());
            s.AppendLine("*       Area under Precision recall Curve:  " + areaUnderPrecisionRecallCurve.Average().ToString());
            s.AppendLine("*       F1Score:  " + F1Score.Average().ToString());
            s.AppendLine("*       LogLoss:  " + logLossAverage.ToString());
            s.AppendLine("*       LogLossReduction:  " + logLossReductionAverage.ToString());
            s.AppendLine("*       PositivePrecision:  " + positivePrecision.Average().ToString());
            s.AppendLine("*       PositiveRecall:  " + positiveRecall.Average().ToString());
            s.AppendLine("*       NegativePrecision:  " + negativePrecision.Average().ToString());
            s.AppendLine("*       NegativeRecall:  " + negativeRecall.Average().ToString());
            s.AppendLine("*       Count of Ambiguous Peptides Removed:  " + sumOfAllAmbiguousPeptidesResolved.ToString());
            s.AppendLine("************************************************************");
            return s.ToString();
        }

        public static int Compute_PSM_PEP(List<SpectralMatch> psms, List<int> psmIndices, MLContext mLContext, TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>> trainedModel, string searchType, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, Dictionary<string, float> fileSpecificMedianFragmentMassErrors, int chargeStateMode, string outputFolder)
        {
            int maxThreads = fileSpecificParameters.FirstOrDefault().fileSpecificParameters.MaxThreadsToUsePerFile;
            object lockObject = new object();
            int ambiguousPeptidesResolved = 0;

            //the trained model is not threadsafe. Therefore, to use the same model for each thread saved the model to disk. Then each thread reads its own copy of the model back from disk.
            //If there is no output folder specified, then this can't happen. We set maxthreads eqaul to one and use the model that gets passed into the method.
            if (String.IsNullOrEmpty(outputFolder))
            {
                maxThreads = 1;
            }

            Parallel.ForEach(Partitioner.Create(0, psmIndices.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    ITransformer threadSpecificTrainedModel;
                    if (maxThreads == 1)
                    {
                        threadSpecificTrainedModel = trainedModel;
                    }
                    else
                    {
                        threadSpecificTrainedModel = mLContext.Model.Load(Path.Combine(outputFolder, "model.zip"), out DataViewSchema savedModelSchema);
                    }

                    // one prediction engine per thread, because the prediction engine is not thread-safe
                    var threadPredictionEngine = mLContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(threadSpecificTrainedModel);

                    int ambigousPeptidesRemovedinThread = 0;

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        SpectralMatch psm = psms[psmIndices[i]];

                        if (psm != null)
                        {
                            List<int> indiciesOfPeptidesToRemove = new List<int>();
                            List<double> pepValuePredictions = new List<double>();

                            //Here we compute the pepvalue predection for each ambiguous peptide in a PSM. Ambiguous peptides with lower pepvalue predictions are removed from the PSM.

                            List<int> allBmpNotches = new List<int>();
                            List<IBioPolymerWithSetMods> allBmpPeptides = new List<IBioPolymerWithSetMods>();

                            foreach (var (Notch, Peptide) in psm.BestMatchingBioPolymersWithSetMods)
                            {
                                allBmpNotches.Add(Notch);
                                allBmpPeptides.Add(Peptide);
                                PsmData pd = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, Peptide, Notch, !Peptide.Parent.IsDecoy);
                                var pepValuePrediction = threadPredictionEngine.Predict(pd);
                                pepValuePredictions.Add(pepValuePrediction.Probability);
                                //A score is available using the variable pepvaluePrediction.Score
                            }

                            GetIndiciesOfPeptidesToRemove(indiciesOfPeptidesToRemove, pepValuePredictions);
                            int peptidesRemoved = 0;
                            RemoveBestMatchingPeptidesWithLowPEP(psm, indiciesOfPeptidesToRemove, allBmpNotches, allBmpPeptides, pepValuePredictions, ref peptidesRemoved);
                            ambigousPeptidesRemovedinThread += peptidesRemoved;
                        }
                    }

                    lock (lockObject)
                    {
                        ambiguousPeptidesResolved += ambigousPeptidesRemovedinThread;
                    }
                });
            return ambiguousPeptidesResolved;
        }

        public static List<int>[] Get_PSM_Group_Indices(List<SpectralMatch> psms, int numGroups)
        {
            List<int>[] groupsOfIndicies = new List<int>[numGroups];
            var targetIndexes = psms.Select((item, index) => new { Item = item, Index = index })
                .Where(x => !x.Item.IsDecoy)
                .Select(x => x.Index)
                .ToList();
            RandomizeListInPlace(targetIndexes);
            var decoyIndexes = psms.Select((item, index) => new { Item = item, Index = index })
                .Where(x => x.Item.IsDecoy)
                .Select(x => x.Index)
                .ToList();
            RandomizeListInPlace(decoyIndexes);

            var targetGroups = DivideListIntoGroups(targetIndexes, numGroups);
            var decoyGroups = DivideListIntoGroups(decoyIndexes, numGroups);

            for (int i = 0; i < numGroups; i++)
            {
                groupsOfIndicies[i] = targetGroups[i].Concat(decoyGroups[i]).ToList();
            }

            return groupsOfIndicies;
        }

        static void RandomizeListInPlace<T>(List<T> list)
        {
            Random rng = new Random(42);
            int n = list.Count;
            while (n > 1)
            {
                n--;
                int k = rng.Next(n + 1);
                T value = list[k];
                list[k] = list[n];
                list[n] = value;
            }

        }

        static List<List<T>> DivideListIntoGroups<T>(List<T> list, int n)
        {
            var groups = new List<List<T>>();
            int groupSize = (int)Math.Ceiling(list.Count / (double)n);

            for (int i = 0; i < n; i++)
            {
                groups.Add(list.Skip(i * groupSize).Take(groupSize).ToList());
            }

            return groups;
        }

        public static void RemoveBestMatchingPeptidesWithLowPEP(SpectralMatch psm, List<int> indiciesOfPeptidesToRemove, List<int> notches, List<IBioPolymerWithSetMods> pwsmList, List<double> pepValuePredictions, ref int ambiguousPeptidesRemovedCount)
        {
            foreach (int i in indiciesOfPeptidesToRemove)
            {
                psm.RemoveThisAmbiguousPeptide(notches[i], pwsmList[i]);
                ambiguousPeptidesRemovedCount++;
            }
            psm.PsmFdrInfo.PEP = 1 - pepValuePredictions.Max();
            psm.PeptideFdrInfo.PEP = 1 - pepValuePredictions.Max();
        }

        /// <summary>
        /// Given a set of PEP values, this method will find the indicies of BestMatchingBioPolymersWithSetMods that are not within the required tolerance
        /// This method will also remove the low scoring predictions from the set.
        /// </summary>
        public static void GetIndiciesOfPeptidesToRemove(List<int> indiciesOfPeptidesToRemove, List<double> pepValuePredictions)
        {
            double highestPredictedPEPValue = pepValuePredictions.Max();
            for (int i = 0; i < pepValuePredictions.Count; i++)
            {
                if ((highestPredictedPEPValue - pepValuePredictions[i]) > AbsoluteProbabilityThatDistinguishesPeptides)
                {
                    indiciesOfPeptidesToRemove.Add(i);
                }
            }

            foreach (int i in indiciesOfPeptidesToRemove.OrderByDescending(p => p))
            {
                pepValuePredictions.RemoveAt(i);
            }
        }

        /// <summary>
        /// Here we're getting the most common charge state for precursors that are Targets with q<=0.01.

        public static int GetChargeStateMode(List<SpectralMatch> psms)
        {
            return psms.Where(p => p.IsDecoy != true && p.FdrInfo.QValue <= 0.01).Select(p => p.ScanPrecursorCharge).GroupBy(n => n).OrderByDescending(g => g.Count()).Select(g => g.Key).FirstOrDefault();
        }

        public static Dictionary<string, Dictionary<int, Tuple<double, double>>> ComputeHydrophobicityValues(List<SpectralMatch> psms, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, bool computeHydrophobicitiesforModifiedPeptides)
        {
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);

            //TODO change the tuple so the values have names
            Dictionary<string, Dictionary<int, Tuple<double, double>>> rtHydrophobicityAvgDev = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            List<string> filenames = fileSpecificParameters.Where(s => s.fileSpecificParameters.SeparationType == "HPLC").Select(s => Path.GetFileName(s.fileName)).ToList();

            filenames = filenames.Distinct().ToList();

            foreach (string filename in filenames)
            {
                Dictionary<int, List<double>> hydrophobicities = new Dictionary<int, List<double>>();
                Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations = new Dictionary<int, Tuple<double, double>>();

                foreach (SpectralMatch psm in psms.Where(f => (f.FullFilePath == null || Path.GetFileName(f.FullFilePath) == filename) && f.FdrInfo.QValue <= 0.01 && !f.IsDecoy))
                {
                    List<string> fullSequences = new List<string>();
                    foreach ((int notch, IBioPolymerWithSetMods pwsm) in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        if (fullSequences.Contains(pwsm.FullSequence))
                        {
                            continue;
                        }
                        fullSequences.Add(pwsm.FullSequence);

                        double predictedHydrophobicity = pwsm is PeptideWithSetModifications pep ?  calc.ScoreSequence(pep) : 0;

                        //here i'm grouping this in 2 minute increments becuase there are cases where you get too few data points to get a good standard deviation an average. This is for stability.
                        int possibleKey = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));

                        //First block of if statement is for modified peptides.
                        if (pwsm.AllModsOneIsNterminus.Any() && computeHydrophobicitiesforModifiedPeptides)
                        {
                            if (hydrophobicities.ContainsKey(possibleKey))
                            {
                                hydrophobicities[possibleKey].Add(predictedHydrophobicity);
                            }
                            else
                            {
                                hydrophobicities.Add(possibleKey, new List<double>() { predictedHydrophobicity });
                            }
                        }
                        //this second block of if statment is for unmodified peptides.
                        else if (!pwsm.AllModsOneIsNterminus.Any() && !computeHydrophobicitiesforModifiedPeptides)
                        {
                            if (hydrophobicities.ContainsKey(possibleKey))
                            {
                                hydrophobicities[possibleKey].Add(predictedHydrophobicity);
                            }
                            else
                            {
                                hydrophobicities.Add(possibleKey, new List<double>() { predictedHydrophobicity });
                            }
                        }
                    }
                }

                List<double> allSquaredHyrophobicityDifferences = new List<double>();

                foreach (int retentionTimeBin in hydrophobicities.Keys)
                {
                    //TODO consider using inner-quartile range instead of standard deviation
                    double averageHydrophobicity = hydrophobicities[retentionTimeBin].Average();
                    averagesCommaStandardDeviations.Add(retentionTimeBin, new Tuple<double, double>(averageHydrophobicity, hydrophobicities[retentionTimeBin].StandardDeviation()));
                    foreach (double hydrophobicity in hydrophobicities[retentionTimeBin])
                    {
                        double difference = Math.Abs(hydrophobicity - averageHydrophobicity);
                        if (!double.IsNaN(difference) && difference > 0)
                        {
                            allSquaredHyrophobicityDifferences.Add(Math.Pow(difference, 2));
                        }
                    }
                }

                //some standard deviations are too small or too large because of random reasons, so we replace those small numbers of oddballs with reasonable numbers.
                double globalStDev = 1;
                if (allSquaredHyrophobicityDifferences.Count() > 1)
                {
                    globalStDev = Math.Sqrt(allSquaredHyrophobicityDifferences.Sum() / (allSquaredHyrophobicityDifferences.Count() - 1));
                }

                Dictionary<int, Tuple<double, double>> stDevsToChange = new Dictionary<int, Tuple<double, double>>();
                foreach (KeyValuePair<int, Tuple<double, double>> item in averagesCommaStandardDeviations)
                {
                    //add stability. not allowing stdevs that are too small or too large at one position relative to the global stdev
                    //here we are finding which stdevs are out of whack.
                    if (Double.IsNaN(item.Value.Item2) || item.Value.Item2 < 0.5 || (item.Value.Item2 / globalStDev) > 3)
                    {
                        Tuple<double, double> pair = new Tuple<double, double>(averagesCommaStandardDeviations[item.Key].Item1, globalStDev);
                        stDevsToChange.Add(item.Key, pair);
                    }
                }
                //here we are replacing the stdevs that are out of whack.
                foreach (int key in stDevsToChange.Keys)
                {
                    averagesCommaStandardDeviations[key] = stDevsToChange[key];
                }

                rtHydrophobicityAvgDev.Add(filename, averagesCommaStandardDeviations);
            }
            return rtHydrophobicityAvgDev;
        }

        public static Dictionary<string, Dictionary<int, Tuple<double, double>>> ComputeMobilityValues(List<SpectralMatch> psms, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters)
        {
            Dictionary<string, Dictionary<int, Tuple<double, double>>> rtMobilityAvgDev = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            List<string> filenames = fileSpecificParameters.Where(s => s.fileSpecificParameters.SeparationType == "CZE").Select(s => Path.GetFileName(s.fileName)).ToList();

            filenames = filenames.Distinct().ToList();

            foreach (string filename in filenames)
            {
                Dictionary<int, List<double>> mobilities = new Dictionary<int, List<double>>();
                Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations = new Dictionary<int, Tuple<double, double>>();

                foreach (SpectralMatch psm in psms.Where(f => (f.FullFilePath == null || Path.GetFileName(f.FullFilePath) == filename) && f.FdrInfo.QValue <= 0.01 && !f.IsDecoy))
                {
                    List<string> fullSequences = new List<string>();
                    foreach ((int notch, IBioPolymerWithSetMods pwsm) in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        if (fullSequences.Contains(pwsm.FullSequence))
                        {
                            continue;
                        }
                        fullSequences.Add(pwsm.FullSequence);

                        double predictedMobility = pwsm is PeptideWithSetModifications pep ? 100.0 * GetCifuentesMobility(pep) : 0;

                        //here i'm grouping this in 2 minute increments becuase there are cases where you get too few data points to get a good standard deviation an average. This is for stability.
                        int possibleKey = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));

                        if (mobilities.ContainsKey(possibleKey))
                        {
                            mobilities[possibleKey].Add(predictedMobility);
                        }
                        else
                        {
                            mobilities.Add(possibleKey, new List<double> { predictedMobility });
                        }
                    }
                }

                List<double> allSquaredMobilityDifferences = new List<double>();

                foreach (int retentionTimeBin in mobilities.Keys)
                {
                    //TODO consider using inner-quartile range instead of standard deviation
                    double averageMobility = mobilities[retentionTimeBin].Average();
                    averagesCommaStandardDeviations.Add(retentionTimeBin, new Tuple<double, double>(averageMobility, mobilities[retentionTimeBin].StandardDeviation()));
                    foreach (double hydrophobicity in mobilities[retentionTimeBin])
                    {
                        double difference = Math.Abs(hydrophobicity - averageMobility);
                        if (!double.IsNaN(difference) && difference > 0)
                        {
                            allSquaredMobilityDifferences.Add(Math.Pow(difference, 2));
                        }
                    }
                }

                //some standard deviations are too small or too large because of random reasons, so we replace those small numbers of oddballs with reasonable numbers.
                double globalStDev = 1;
                if (allSquaredMobilityDifferences.Count() > 1)
                {
                    globalStDev = Math.Sqrt(allSquaredMobilityDifferences.Sum() / (allSquaredMobilityDifferences.Count() - 1));
                }

                Dictionary<int, Tuple<double, double>> stDevsToChange = new Dictionary<int, Tuple<double, double>>();

                GetStDevsToChange(stDevsToChange, averagesCommaStandardDeviations, globalStDev);
                UpdateOutOfRangeStDevsWithGlobalAverage(stDevsToChange, averagesCommaStandardDeviations);

                rtMobilityAvgDev.Add(filename, averagesCommaStandardDeviations);
            }
            return rtMobilityAvgDev;
        }

        /// <summary>
        /// This gathers a set of standard deviations that are outside the range of acceptable.
        /// </summary>
        public static void GetStDevsToChange(Dictionary<int, Tuple<double, double>> stDevsToChange, Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations, double globalStDev)
        {
            foreach (KeyValuePair<int, Tuple<double, double>> item in averagesCommaStandardDeviations)
            {
                //add stability. not allowing stdevs that are too small or too large at one position relative to the global stdev
                //here we are finding which stdevs are out of whack.
                if (Double.IsNaN(item.Value.Item2) || item.Value.Item2 < 0.05 || (item.Value.Item2 / globalStDev) > 3)
                {
                    Tuple<double, double> pair = new Tuple<double, double>(averagesCommaStandardDeviations[item.Key].Item1, globalStDev);
                    stDevsToChange.Add(item.Key, pair);
                }
            }
        }

        /// <summary>
        /// here we are replacing the stdevs that are out of whack.
        /// </summary>
        public static void UpdateOutOfRangeStDevsWithGlobalAverage(Dictionary<int, Tuple<double, double>> stDevsToChange, Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations)
        {
            foreach (int key in stDevsToChange.Keys)
            {
                averagesCommaStandardDeviations[key] = stDevsToChange[key];
            }
        }

        private static double GetCifuentesMobility(IBioPolymerWithSetMods pwsm)
        {
            int charge = 1 + pwsm.BaseSequence.Count(f => f == 'K') + pwsm.BaseSequence.Count(f => f == 'R') + pwsm.BaseSequence.Count(f => f == 'H') - CountModificationsThatShiftMobility(pwsm.AllModsOneIsNterminus.Values.AsEnumerable());// the 1 + is for N-terminal

            double mobility = (Math.Log(1 + 0.35 * (double)charge)) / Math.Pow(pwsm.MonoisotopicMass, 0.411);

            return mobility;
        }

        private static float GetSSRCalcHydrophobicityZScore(SpectralMatch psm, IBioPolymerWithSetMods Peptide, Dictionary<string, Dictionary<int, Tuple<double, double>>> d)
        {
            //Using SSRCalc3 but probably any number of different calculators could be used instead. One could also use the CE mobility.
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);
            double hydrophobicityZscore = double.NaN;

            if (d.ContainsKey(Path.GetFileName(psm.FullFilePath)))
            {
                int time = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));
                if (d[Path.GetFileName(psm.FullFilePath)].Keys.Contains(time))
                {
                    double predictedHydrophobicity = Peptide is PeptideWithSetModifications pep ? calc.ScoreSequence(pep) : 0;

                    hydrophobicityZscore = Math.Abs(d[Path.GetFileName(psm.FullFilePath)][time].Item1 - predictedHydrophobicity) / d[Path.GetFileName(psm.FullFilePath)][time].Item2;
                }
            }

            double maxHydrophobicityZscore = 10; // each "Z" is one standard deviation. so, maxHydrophobicityZscore 10 is quite large
            if (double.IsNaN(hydrophobicityZscore) || double.IsInfinity(hydrophobicityZscore) || hydrophobicityZscore > maxHydrophobicityZscore)
            {
                hydrophobicityZscore = maxHydrophobicityZscore;
            }

            return (float)hydrophobicityZscore;
        }

        private static float GetMobilityZScore(SpectralMatch psm, IBioPolymerWithSetMods selectedPeptide)
        {
            double mobilityZScore = double.NaN;

            if (fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE.ContainsKey(Path.GetFileName(psm.FullFilePath)))
            {
                int time = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));
                if (fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[Path.GetFileName(psm.FullFilePath)].Keys.Contains(time))
                {
                    double predictedMobility = 100.0 * GetCifuentesMobility(selectedPeptide);

                    mobilityZScore = Math.Abs(fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[Path.GetFileName(psm.FullFilePath)][time].Item1 - predictedMobility) / fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[Path.GetFileName(psm.FullFilePath)][time].Item2;
                }
            }

            double maxMobilityZscore = 10; // each "Z" is one standard deviation. so, maxHydrophobicityZscore 10 is quite large
            if (double.IsNaN(mobilityZScore) || double.IsInfinity(mobilityZScore) || mobilityZScore > maxMobilityZscore)
            {
                mobilityZScore = maxMobilityZscore;
            }

            return (float)mobilityZScore;
        }

        public static IEnumerable<PsmData> CreatePsmData(string searchType, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
            List<SpectralMatch> psms, List<int> psmIndicies,
            Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_unmodified,
            Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_modified,
            Dictionary<string, float> fileSpecificMedianFragmentMassErrors, int chargeStateMode)
        {
            object psmDataListLock = new object();
            List<PsmData> psmDataList = new List<PsmData>();
            List<double> psmOrder = new List<double>();
            int maxThreads = fileSpecificParameters.FirstOrDefault().fileSpecificParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            Parallel.ForEach(Partitioner.Create(0, psmIndicies.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<PsmData> localPsmDataList = new List<PsmData>();
                    List<double> localPsmOrder = new List<double>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        SpectralMatch psm = psms[psmIndicies[i]];

                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops) { return; }

                        PsmData newPsmData = new PsmData();
                        if (searchType == "crosslink")
                        {
                            CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)psms[i];

                            bool label;
                            if (csm.IsDecoy || csm.BetaPeptide.IsDecoy)
                            {
                                label = false;
                                newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, csm.BestMatchingBioPolymersWithSetMods.First().Peptide, 0, label);
                            }
                            else if (!csm.IsDecoy && !csm.BetaPeptide.IsDecoy && psm.GetFdrInfo(UsePeptideLevelQValueForTraining).QValue <= QValueCutoff)
                            {
                                label = true;
                                newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, csm.BestMatchingBioPolymersWithSetMods.First().Peptide, 0, label);
                            }
                            else
                            {
                                continue;
                            }
                            localPsmDataList.Add(newPsmData);
                            localPsmOrder.Add(i);
                        }
                        else
                        {
                            double bmp = 0;
                            foreach (var (notch, peptideWithSetMods) in psm.BestMatchingBioPolymersWithSetMods)
                            {
                                bool label;
                                double bmpc = psm.BestMatchingBioPolymersWithSetMods.Count();
                                if (peptideWithSetMods.Parent.IsDecoy)
                                {
                                    label = false;
                                    newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, peptideWithSetMods, notch, label);
                                }
                                else if (!peptideWithSetMods.Parent.IsDecoy && psm.GetFdrInfo(UsePeptideLevelQValueForTraining).QValue <= QValueCutoff)
                                {
                                    label = true;
                                    newPsmData = CreateOnePsmDataEntry(searchType, fileSpecificParameters, psm,  timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, peptideWithSetMods, notch, label);
                                }
                                else
                                {
                                    continue;
                                }
                                localPsmDataList.Add(newPsmData);
                                localPsmOrder.Add(i + (bmp / bmpc / 2.0));
                                bmp += 1.0;
                            }
                        }
                    }
                    lock (psmDataListLock)
                    {
                        psmDataList.AddRange(localPsmDataList);
                        psmOrder.AddRange(localPsmOrder);
                    }
                });
            PsmData[] pda = psmDataList.ToArray();
            double[] order = psmOrder.ToArray();

            Array.Sort(order, pda);//this sorts both arrays thru sorting the array in position one. The order array, keeps track of the positon in the original psms list and returns the PsmData array in that same order.

            return pda.AsEnumerable();
        }

        public static PsmData CreateOnePsmDataEntry(string searchType, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, SpectralMatch psm, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_unmodified, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_modified, Dictionary<string, float> fileSpecificMedianFragmentMassErrors, int chargeStateMode, IBioPolymerWithSetMods selectedPeptide, int notchToUse, bool label)
        {
            double normalizationFactor = selectedPeptide.BaseSequence.Length;
            float totalMatchingFragmentCount = 0;
            float internalMatchingFragmentCount = 0;
            float intensity = 0;
            float chargeDifference = 0;
            float deltaScore = 0;
            int notch = 0;
            float ambiguity = 0;
            float modCount = 0;
            float absoluteFragmentMassError = 0;
            float spectralAngle = 0;
            float hasSpectralAngle = 0;
            float chimeraCount = 0;
            float peaksInPrecursorEnvelope = 0;
            float mostAbundantPrecursorPeakIntensity = 0;
            float fractionalIntensity = 0;

            float missedCleavages = 0;
            float longestSeq = 0;
            float complementaryIonCount = 0;
            float hydrophobicityZscore = float.NaN;
            bool isVariantPeptide = false;

            //crosslink specific features
            float alphaIntensity = 0;
            float betaIntensity = 0;
            float longestFragmentIonSeries_Alpha = 0;
            float longestFragmentIonSeries_Beta = 0;
            float isDeadEnd = 0;
            float isLoop = 0;
            float isInter = 0;
            float isIntra = 0;

            double multiplier = 10;
            if (searchType != "crosslink")
            {
                if (searchType == "top-down")
                {
                    normalizationFactor = 1.0;
                }
                // count only terminal fragment ions
                totalMatchingFragmentCount = (float)(Math.Round(psm.BioPolymersWithSetModsToMatchingFragments[selectedPeptide].Count(p => p.NeutralTheoreticalProduct.SecondaryProductType == null) / normalizationFactor * multiplier, 0));
                internalMatchingFragmentCount = (float)(Math.Round(psm.BioPolymersWithSetModsToMatchingFragments[selectedPeptide].Count(p => p.NeutralTheoreticalProduct.SecondaryProductType != null) / normalizationFactor * multiplier, 0));
                intensity = (float)Math.Min(50, Math.Round((psm.Score - (int)psm.Score) / normalizationFactor * Math.Pow(multiplier, 2), 0));
                chargeDifference = -Math.Abs(chargeStateMode - psm.ScanPrecursorCharge);
                deltaScore = (float)Math.Round(psm.DeltaScore / normalizationFactor * multiplier, 0);
                notch = notchToUse;
                modCount = Math.Min((float)selectedPeptide.AllModsOneIsNterminus.Keys.Count(), 10);
                if (psm.BioPolymersWithSetModsToMatchingFragments[selectedPeptide]?.Count() > 0)
                {
                    absoluteFragmentMassError = (float)Math.Min(100.0, Math.Round(10.0 * Math.Abs(GetAverageFragmentMassError(psm.BioPolymersWithSetModsToMatchingFragments[selectedPeptide]) - fileSpecificMedianFragmentMassErrors[Path.GetFileName(psm.FullFilePath)])));
                }

                ambiguity = Math.Min((float)(psm.BioPolymersWithSetModsToMatchingFragments.Keys.Count - 1), 10);
                longestSeq = (float)Math.Round(SpectralMatch.GetLongestIonSeriesBidirectional(psm.BioPolymersWithSetModsToMatchingFragments, selectedPeptide) / normalizationFactor * multiplier, 0);
                complementaryIonCount = (float)Math.Round(SpectralMatch.GetCountComplementaryIons(psm.BioPolymersWithSetModsToMatchingFragments, selectedPeptide) / normalizationFactor * multiplier, 0);
                isVariantPeptide = PeptideIsVariant(selectedPeptide);
                spectralAngle = (float)psm.SpectralAngle;
                if (chimeraCountDictionary.TryGetValue(psm.ChimeraIdString, out int val))
                    chimeraCount = val;
                peaksInPrecursorEnvelope = psm.PrecursorScanEnvelopePeakCount;
                mostAbundantPrecursorPeakIntensity = (float)Math.Round((float)psm.PrecursorScanIntensity / normalizationFactor * multiplier, 0);
                fractionalIntensity = (float)psm.PrecursorFractionalIntensity;

                if (PsmHasSpectralAngle(psm))
                {
                    hasSpectralAngle = 1;
                }

                if (psm.DigestionParams.Protease.Name != "top-down")
                {
                    missedCleavages = selectedPeptide.MissedCleavages;
                    bool fileIsCzeSeparationType = fileSpecificParameters.Any(p => Path.GetFileName(p.fileName) == Path.GetFileName(psm.FullFilePath) && p.fileSpecificParameters.SeparationType == "CZE");

                    if (!fileIsCzeSeparationType)
                    {
                        if (selectedPeptide.BaseSequence.Equals(selectedPeptide.FullSequence))
                        {
                            hydrophobicityZscore = (float)Math.Round(GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_unmodified) * 10.0, 0);
                        }
                        else
                        {
                            hydrophobicityZscore = (float)Math.Round(GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_modified) * 10.0, 0);
                        }
                    }
                    else
                    {
                        hydrophobicityZscore = (float)Math.Round(GetMobilityZScore(psm, selectedPeptide) * 10.0, 0);
                    }
                }
                //this is not for actual crosslinks but for the byproducts of crosslink loop links, deadends, etc.
                if (psm is CrosslinkSpectralMatch)
                {
                    CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)psm;
                    isDeadEnd = Convert.ToSingle((csm.CrossType == PsmCrossType.DeadEnd) || (csm.CrossType == PsmCrossType.DeadEndH2O) || (csm.CrossType == PsmCrossType.DeadEndNH2) || (csm.CrossType == PsmCrossType.DeadEndTris));
                    isLoop = Convert.ToSingle(csm.CrossType == PsmCrossType.Loop);
                }
            }
            else
            {
                CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)psm;
                PeptideWithSetModifications selectedAlphaPeptide = csm.BestMatchingBioPolymersWithSetMods.Select(p => p.Peptide as PeptideWithSetModifications).First();
                PeptideWithSetModifications selectedBetaPeptide = csm.BetaPeptide?.BestMatchingBioPolymersWithSetMods.Select(p => p.Peptide as PeptideWithSetModifications).First();

                float alphaNormalizationFactor = selectedAlphaPeptide.BaseSequence.Length;
                float betaNormalizationFactor = selectedBetaPeptide == null ? (float)0 : selectedBetaPeptide.BaseSequence.Length;
                float totalNormalizationFactor = alphaNormalizationFactor + betaNormalizationFactor;

                totalMatchingFragmentCount = (float)Math.Round(csm.XLTotalScore / totalNormalizationFactor * 10, 0);

                //Compute fragment mass error
                int alphaCount = 0;
                float alphaError = 0;
                if (csm.BioPolymersWithSetModsToMatchingFragments[selectedAlphaPeptide]?.Count > 0)
                {
                    alphaCount = csm.BioPolymersWithSetModsToMatchingFragments[selectedAlphaPeptide].Count;
                    alphaError = Math.Abs(GetAverageFragmentMassError(csm.BioPolymersWithSetModsToMatchingFragments[selectedAlphaPeptide]));
                }
                int betaCount = 0;
                float betaError = 0;
                if (csm.BetaPeptide.BioPolymersWithSetModsToMatchingFragments[selectedBetaPeptide]?.Count > 0)
                {
                    betaCount = csm.BetaPeptide.BioPolymersWithSetModsToMatchingFragments[selectedBetaPeptide].Count;
                    betaError = Math.Abs(GetAverageFragmentMassError(csm.BetaPeptide.BioPolymersWithSetModsToMatchingFragments[selectedBetaPeptide]));
                }

                float averageError = 0;
                if ((alphaCount + betaCount) > 0)
                {
                    averageError = (alphaCount * alphaError + betaCount * betaError) / (alphaCount + betaCount);
                }

                absoluteFragmentMassError = (float)Math.Min(100, Math.Round(averageError - fileSpecificMedianFragmentMassErrors[Path.GetFileName(csm.FullFilePath)] * 10.0, 0));
                //End compute fragment mass error

                deltaScore = (float)Math.Round(csm.DeltaScore / totalNormalizationFactor * 10.0, 0);
                chargeDifference = -Math.Abs(chargeStateMode - psm.ScanPrecursorCharge);
                alphaIntensity = (float)Math.Min(100, Math.Round((csm.Score - (int)csm.Score) / alphaNormalizationFactor * 100.0, 0));
                betaIntensity = csm.BetaPeptide == null ? (float)0 : (float)Math.Min(100.0, Math.Round((csm.BetaPeptide.Score - (int)csm.BetaPeptide.Score) / betaNormalizationFactor * 100.0, 0));
                longestFragmentIonSeries_Alpha = (float)Math.Round(SpectralMatch.GetLongestIonSeriesBidirectional(csm.BioPolymersWithSetModsToMatchingFragments, selectedAlphaPeptide) / alphaNormalizationFactor * 10.0, 0);
                longestFragmentIonSeries_Beta = selectedBetaPeptide == null ? (float)0 : SpectralMatch.GetLongestIonSeriesBidirectional(csm.BetaPeptide.BioPolymersWithSetModsToMatchingFragments, selectedBetaPeptide) / betaNormalizationFactor;
                longestFragmentIonSeries_Beta = (float)Math.Round(longestFragmentIonSeries_Beta * 10.0, 0);
                isInter = Convert.ToSingle(csm.CrossType == PsmCrossType.Inter);
                isIntra = Convert.ToSingle(csm.CrossType == PsmCrossType.Intra);
            }

            psm.PsmData_forPEPandPercolator = new PsmData
            {
                TotalMatchingFragmentCount = totalMatchingFragmentCount,
                Intensity = intensity,
                PrecursorChargeDiffToMode = chargeDifference,
                DeltaScore = deltaScore,
                Notch = notch,
                ModsCount = modCount,
                AbsoluteAverageFragmentMassErrorFromMedian = absoluteFragmentMassError,
                MissedCleavagesCount = missedCleavages,
                Ambiguity = ambiguity,
                LongestFragmentIonSeries = longestSeq,
                ComplementaryIonCount = complementaryIonCount,
                HydrophobicityZScore = hydrophobicityZscore,
                IsVariantPeptide = Convert.ToSingle(isVariantPeptide),

                AlphaIntensity = alphaIntensity,
                BetaIntensity = betaIntensity,
                LongestFragmentIonSeries_Alpha = longestFragmentIonSeries_Alpha,
                LongestFragmentIonSeries_Beta = longestFragmentIonSeries_Beta,
                IsDeadEnd = isDeadEnd,
                IsLoop = isLoop,
                IsInter = isInter,
                IsIntra = isIntra,

                Label = label,

                SpectralAngle = spectralAngle,
                HasSpectralAngle = hasSpectralAngle,
                PeaksInPrecursorEnvelope = peaksInPrecursorEnvelope,
                ChimeraCount = chimeraCount,
                MostAbundantPrecursorPeakIntensity = mostAbundantPrecursorPeakIntensity,
                PrecursorFractionalIntensity = fractionalIntensity,
                InternalIonCount = internalMatchingFragmentCount,
            };

            return psm.PsmData_forPEPandPercolator;
        }

        private static bool PeptideIsVariant(IBioPolymerWithSetMods bpwsm)
        {
            if (bpwsm is not  PeptideWithSetModifications pwsm)
                return false;
            
            bool identifiedVariant = false;
            if (pwsm.Protein.AppliedSequenceVariations.Count() > 0)
            {
                foreach (var variant in pwsm.Protein.AppliedSequenceVariations)
                {
                    if (pwsm.IntersectsAndIdentifiesVariation(variant).identifies)
                    {
                        identifiedVariant = true;
                        break;
                    }
                }
            }
            return identifiedVariant;
        }

        private static bool PsmHasSpectralAngle(SpectralMatch psm)
        {
            return psm.SpectralAngle >= 0;
        }

            public static bool ContainsModificationsThatShiftMobility(IEnumerable<Modification> modifications)
        {
            List<string> shiftingModifications = new List<string> { "Acetylation", "Ammonia loss", "Carbamyl", "Deamidation", "Formylation",
                "N2-acetylarginine", "N6-acetyllysine", "N-acetylalanine", "N-acetylaspartate", "N-acetylcysteine", "N-acetylglutamate", "N-acetylglycine",
                "N-acetylisoleucine", "N-acetylmethionine", "N-acetylproline", "N-acetylserine", "N-acetylthreonine", "N-acetyltyrosine", "N-acetylvaline",
                "Phosphorylation", "Phosphoserine", "Phosphothreonine", "Phosphotyrosine", "Sulfonation" };

            return shiftingModifications.Concat(modifications.Select(m => m.OriginalId).Distinct()).GroupBy(s => s).Where(s => s.Count() > 1).Any();
        }

        public static int CountModificationsThatShiftMobility(IEnumerable<Modification> modifications)
        {
            List<string> shiftingModifications = new List<string> { "Acetylation", "Ammonia loss", "Carbamyl", "Deamidation", "Formylation",
                "N2-acetylarginine", "N6-acetyllysine", "N-acetylalanine", "N-acetylaspartate", "N-acetylcysteine", "N-acetylglutamate", "N-acetylglycine",
                "N-acetylisoleucine", "N-acetylmethionine", "N-acetylproline", "N-acetylserine", "N-acetylthreonine", "N-acetyltyrosine", "N-acetylvaline",
                "Phosphorylation", "Phosphoserine", "Phosphothreonine", "Phosphotyrosine", "Sulfonation" };

            return modifications.Select(n => n.OriginalId).Intersect(shiftingModifications).Count();
        }

        public static Dictionary<string, float> GetFileSpecificMedianFragmentMassError(List<SpectralMatch> psms)
        {
            Dictionary<string, float> fileSpecificMassErrors = new Dictionary<string, float>();
            foreach (string filename in psms.Select(p => Path.GetFileName(p.FullFilePath)).Distinct())
            {
                fileSpecificMassErrors.Add(filename, GetMedianAverageMassError(psms.Where(p => Path.GetFileName(p.FullFilePath) == filename)));
            }
            return fileSpecificMassErrors;
        }

        public static float GetMedianAverageMassError(IEnumerable<SpectralMatch> psms)
        {
            List<float> averageMassErrors = new List<float>();
            foreach (SpectralMatch psm in psms)
            {
                {
                    foreach (KeyValuePair<IBioPolymerWithSetMods, List<MatchedFragmentIon>> peptide_MFI in psm.BioPolymersWithSetModsToMatchingFragments)
                    {
                        if (peptide_MFI.Value != null && peptide_MFI.Value.Count > 0)
                        {
                            averageMassErrors.Add(GetAverageFragmentMassError(peptide_MFI.Value));
                        }
                    }
                }
            }
            return averageMassErrors.Median();
        }

        public static float GetAverageFragmentMassError(IEnumerable<MatchedFragmentIon> matchedIons)
        {
            var matchedIonsGroupedByProductType = matchedIons.GroupBy(i => i.NeutralTheoreticalProduct.ProductType).OrderBy(i => i.Key).ToList();
            List<float> massErrors = new List<float>();
            foreach (var productType in matchedIonsGroupedByProductType)
            {
                var products = productType.OrderBy(p => p.NeutralTheoreticalProduct.FragmentNumber)
                    .ToList();

                for (int i = 0; i < products.Count; i++)
                {
                    MatchedFragmentIon ion = products[i];

                    float massError = (float)(ion.Mz.ToMass(ion.Charge) - ion.NeutralTheoreticalProduct.NeutralMass);
                    float ppmMassError = (float)(massError / ion.NeutralTheoreticalProduct.NeutralMass * 1e6);
                    massErrors.Add(ppmMassError);
                }
            }

            return massErrors.Average();
        }
    }
}