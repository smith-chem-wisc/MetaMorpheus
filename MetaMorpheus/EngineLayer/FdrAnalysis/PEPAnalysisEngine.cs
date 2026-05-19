using Chemistry;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.SSRCalc;
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
using Easy.Common.Extensions;
using System.Threading;
using EngineLayer.SpectrumMatch;

namespace EngineLayer
{
    public class PepAnalysisEngine
    {
        private int _randomSeed = 42;

        /// <summary>
        /// This method contains the hyper-parameters that will be used when training the machine learning model
        /// </summary>
        /// <returns> Options object to be passed in to the FastTree constructor </returns>
        public Microsoft.ML.Trainers.FastTree.FastTreeBinaryTrainer.Options BGDTreeOptions =>
            new Microsoft.ML.Trainers.FastTree.FastTreeBinaryTrainer.Options
            {
                NumberOfThreads = 1,
                NumberOfTrees = 400,
                MinimumExampleCountPerLeaf = 10,
                NumberOfLeaves = 20,
                LearningRate = 0.2,
                LabelColumnName = "Label",
                FeatureColumnName = "Features",
                Seed = _randomSeed,
                FeatureSelectionSeed = _randomSeed,
                RandomStart = false
            };

        private static readonly double AbsoluteProbabilityThatDistinguishesPeptides = 0.05;

        //These two dictionaries contain the average and standard deviations of hydrophobicitys measured in 1 minute increments accross each raw
        //file separately. An individully measured hydrobophicty calculated for a specific PSM sequence is compared to these values by computing
        //the z-score. That z-score is used as a feature for machine learning.
        //Separate dictionaries are created for peptides with modifications because SSRcalc doesn't really do a good job predicting hyrophobicity

        //The first string in the dictionary is the filename
        //The value of the dictionary is another dictionary that profiles the hydrophobicity behavior.
        //Each key is a retention time rounded to the nearest minute.
        //The value Tuple is the average and standard deviation, respectively, of the predicted hydrophobicities of the observed peptides eluting at that rounded retention time.
        public Dictionary<string, Dictionary<int, Tuple<double, double>>> FileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified { get; private set; }
        public Dictionary<string, Dictionary<int, Tuple<double, double>>> FileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified { get; private set; }
        public Dictionary<string, Dictionary<int, Tuple<double, double>>> FileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE { get; private set; }

        /// <summary>
        /// A dictionary which stores the chimeric ID string in the key and the number of chimeric identifications as the vale
        /// </summary>
        private Dictionary<string, int> chimeraCountDictionary = new Dictionary<string, int>();
        public Dictionary<string, float> FileSpecificMedianFragmentMassErrors { get; private set; }
        public Dictionary<string, CommonParameters> FileSpecificParametersDictionary { get; private set; }
        public int ChargeStateMode { get; private set; }

        public double QValueCutoff { get; }
        public bool UsePeptideLevelQValueForTraining = true;
        public string[] TrainingVariables { get; }
        public string OutputFolder { get; }
        public List<SpectralMatch> AllPsms { get; }
        public string SearchType { get; }
        public IRetentionTimePredictor RetentionTimePredictor { get; }

        /// <summary>
        /// When true, after the first cross-fit fit the positive training pool is re-selected
        /// using PEP-derived q-values and the model is retrained, repeating until the positive
        /// pool stabilizes or <see cref="MaxTrainingIterations"/> is reached. When false, PEP
        /// training is single-pass (the historical behavior).
        /// </summary>
        public bool IterativeTraining { get; }

        /// <summary>
        /// Hard upper bound on the number of training iterations when <see cref="IterativeTraining"/>
        /// is true. Ignored when iterative training is off (a single pass always runs).
        /// </summary>
        public int MaxTrainingIterations { get; }

        /// <summary>
        /// Positive (target) training-pool size recorded once per training iteration, in order.
        /// Always has at least one entry after <see cref="ComputePEPValuesForAllPSMs"/> runs.
        /// </summary>
        public List<int> PositivePoolCountPerIteration { get; } = new List<int>();

        /// <summary>
        /// Per-fold PEP-derived q-value maps used to re-select positive training examples on the
        /// most recent relabeling step. Each fold's map is keyed exclusively by that fold's PSMs,
        /// and the maps are pairwise disjoint — fold k is relabeled only from predictions of the
        /// model trained on folds != k, never from a single global PEP. Null until a relabel runs.
        /// </summary>
        internal List<Dictionary<SpectralMatch, double>> PerFoldRelabelPepQValues { get; private set; }

        /// <summary>
        /// Fractional change in positive-pool size below which iterative training is considered
        /// converged and stops early.
        /// </summary>
        internal const double ConvergenceThreshold = 0.01;

        /// <summary>
        /// This method is used to compute the PEP values for all PSMs in a dataset. 
        /// </summary>
        /// <param name="psms"></param>
        /// <param name="searchType"></param>
        /// <param name="fileSpecificParameters"></param>
        /// <param name="outputFolder"></param>
        /// <returns></returns>
        public void SetFileSpecificParameters(List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters)
        {
            FileSpecificParametersDictionary = fileSpecificParameters.ToDictionary(p => Path.GetFileName(p.fileName), p => p.fileSpecificParameters);
        }

        public PepAnalysisEngine(List<SpectralMatch> psms, string searchType, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, string outputFolder, IRetentionTimePredictor? rtPredictor = null, bool iterativeTraining = false, int maxTrainingIterations = 3)
        {
            // This creates a new list of PSMs, but does not clone the Psms themselves.
            // This allows the PSMs to be modified and the order to be preserved
            AllPsms = psms.OrderByDescending(p => p).ToList();
            TrainingVariables = PsmData.trainingInfos[searchType];
            RetentionTimePredictor = rtPredictor ?? new SSRCalc3RetentionTimePredictor();
            OutputFolder = outputFolder;
            SearchType = searchType;
            SetFileSpecificParameters(fileSpecificParameters);
            BuildFileSpecificDictionaries(psms, TrainingVariables);
            double minQ = searchType == "top-down" ? 0.025 : 0.005; // Less stringent FDR cut-off for top-down
            QValueCutoff = Math.Max(fileSpecificParameters.Select(t => t.fileSpecificParameters.QValueCutoffForPepCalculation).Min(), minQ);
            // If we have more than 100 peptides, we will train on the peptide level. Otherwise, we will train on the PSM level
            UsePeptideLevelQValueForTraining = psms.Select(psm => psm.FullSequence).Distinct().Count(seq => seq.IsNotNullOrEmpty()) >= 100;
            IterativeTraining = iterativeTraining;
            MaxTrainingIterations = Math.Max(1, maxTrainingIterations);
        }

        public string ComputePEPValuesForAllPSMs()
        {
            List<SpectralMatchGroup> peptideGroups = UsePeptideLevelQValueForTraining
                ? SpectralMatchGroup.GroupByBaseSequence(AllPsms)
                : SpectralMatchGroup.GroupByIndividualPsm(AllPsms);

            if (UsePeptideLevelQValueForTraining && (peptideGroups.Count(g => g.BestMatch.IsDecoy) < 4 || peptideGroups.Count(g => !g.BestMatch.IsDecoy) < 4))
            {
                // If we don't have enough peptides to train at the peptide level, we will train at the PSM level
                peptideGroups = SpectralMatchGroup.GroupByIndividualPsm(AllPsms);
                UsePeptideLevelQValueForTraining = false;
            }

            int numGroups = 4;
            List<int>[] peptideGroupIndices = GetPeptideGroupIndices(peptideGroups, numGroups);
            int maxThreads = FileSpecificParametersDictionary.Values.FirstOrDefault().MaxThreadsToUsePerFile;

            // When iterative training is off, a single pass runs and the result is identical to the
            // historical single-pass pipeline. When on, the positive training pool is re-selected from
            // PEP-derived q-values after each fit and the model retrained, up to MaxTrainingIterations
            // or until the positive pool stabilizes (see ConvergenceThreshold).
            int maxIterations = IterativeTraining ? MaxTrainingIterations : 1;

            MLContext mlContext = null;
            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels = null;
            List<CalibratedBinaryClassificationMetrics> allMetrics = null;
            int sumOfAllAmbiguousPeptidesResolved = 0;
            int positiveTrainingCount = 0;
            int negativeTrainingCount = 0;

            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                // Build the labeled training data for this iteration. Iteration 0 labels positives by
                // raw-score q-value (the historical rule); later iterations relabel from the previous
                // iteration's per-fold, cross-fit-isolated PEP q-values held in PerFoldRelabelPepQValues.
                IEnumerable<PsmData>[] PSMDataGroups = new IEnumerable<PsmData>[numGroups];
                bool allGroupsHavePositiveAndNegativeTrainingExamples = true;
                Parallel.ForEach(
                    Enumerable.Range(0, numGroups),
                    new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                    group => {
                        PSMDataGroups[group] = CreatePsmData(SearchType, peptideGroups, peptideGroupIndices[group], PerFoldRelabelPepQValues?[group]);
                        if (!PSMDataGroups[group].Any(p => p.Label) || !PSMDataGroups[group].Any(p => !p.Label))
                        {
                            allGroupsHavePositiveAndNegativeTrainingExamples = false;
                        }

                    });

                if (!allGroupsHavePositiveAndNegativeTrainingExamples)
                {
                    if (iteration == 0)
                    {
                        return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
                    }
                    // A relabeling step left a fold without both classes. Keep the previous iteration's
                    // models, finish by removing ambiguous peptides with them, and stop iterating.
                    sumOfAllAmbiguousPeptidesResolved = 0;
                    for (int k = 0; k < numGroups; k++)
                    {
                        sumOfAllAmbiguousPeptidesResolved += Compute_PSM_PEP(peptideGroups, peptideGroupIndices[k], mlContext, trainedModels[k], SearchType, OutputFolder, removeAmbiguous: true);
                    }
                    break;
                }

                positiveTrainingCount = PSMDataGroups.SelectMany(p => p).Count(p => p.Label);
                negativeTrainingCount = PSMDataGroups.SelectMany(p => p).Count(p => !p.Label);
                PositivePoolCountPerIteration.Add(positiveTrainingCount);

                // Converged once the positive pool barely moves between iterations.
                bool converged = iteration > 0
                    && HasConverged(PositivePoolCountPerIteration[iteration - 1], positiveTrainingCount);
                bool isFinalPass = converged || iteration == maxIterations - 1;

                mlContext = new MLContext(seed: _randomSeed);
                trainedModels = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

                var trainer = mlContext.BinaryClassification.Trainers.FastTree(BGDTreeOptions);
                var pipeline = mlContext.Transforms.Concatenate("Features", TrainingVariables)
                    .Append(trainer);

                allMetrics = new List<CalibratedBinaryClassificationMetrics>();
                sumOfAllAmbiguousPeptidesResolved = 0;

                for (int groupIndexNumber = 0; groupIndexNumber < numGroups; groupIndexNumber++)
                {
                    // Train on every fold except the held-out one, then predict the held-out fold. This
                    // cross-fit structure means each PSM's PEP comes from a model that never saw it.
                    IEnumerable<int> trainingFolds = Enumerable.Range(0, numGroups).Where(g => g != groupIndexNumber);
                    IDataView dataView = mlContext.Data.LoadFromEnumerable(trainingFolds.SelectMany(g => PSMDataGroups[g]));
                    trainedModels[groupIndexNumber] = pipeline.Fit(dataView);
                    var myPredictions = trainedModels[groupIndexNumber].Transform(mlContext.Data.LoadFromEnumerable(PSMDataGroups[groupIndexNumber]));
                    CalibratedBinaryClassificationMetrics metrics = mlContext.BinaryClassification.Evaluate(data: myPredictions, labelColumnName: "Label", scoreColumnName: "Score");

                    // Ambiguous-peptide removal mutates the PSM, so it is deferred to the final pass.
                    int ambiguousPeptidesResolved = Compute_PSM_PEP(peptideGroups, peptideGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], SearchType, OutputFolder, removeAmbiguous: isFinalPass);

                    allMetrics.Add(metrics);
                    sumOfAllAmbiguousPeptidesResolved += ambiguousPeptidesResolved;
                }

                if (isFinalPass)
                {
                    break;
                }

                // Prepare the next iteration's relabeling: one PEP q-value map per fold, each computed
                // only from that fold's held-out cross-fit PEPs. Fold k is therefore relabeled solely
                // from the model trained on folds != k, never from a single global PEP.
                PerFoldRelabelPepQValues = ComputePerFoldPepQValues(peptideGroups, peptideGroupIndices);
            }

            return AggregateMetricsForOutput(allMetrics, sumOfAllAmbiguousPeptidesResolved, positiveTrainingCount, negativeTrainingCount, QValueCutoff, PositivePoolCountPerIteration);
        }

        /// <summary>
        /// True when the positive training pool changed by less than <see cref="ConvergenceThreshold"/>
        /// (as a fraction of the previous pool) between two iterations. Iterative training stops once
        /// this holds, since further relabeling would no longer meaningfully move the training set.
        /// </summary>
        internal static bool HasConverged(int previousPositiveCount, int currentPositiveCount)
        {
            return Math.Abs(currentPositiveCount - previousPositiveCount)
                   < ConvergenceThreshold * Math.Max(1, previousPositiveCount);
        }

        /// <summary>
        /// Computes one PEP-derived q-value map per cross-fit fold. Each fold's map is built strictly
        /// from that fold's own held-out PEP predictions (a standard target/decoy q-value, monotonized),
        /// so the maps are pairwise disjoint and fold k is relabeled only from the model trained on
        /// folds != k. This is what preserves cross-fit isolation across training iterations.
        /// </summary>
        internal List<Dictionary<SpectralMatch, double>> ComputePerFoldPepQValues(List<SpectralMatchGroup> peptideGroups, List<int>[] peptideGroupIndices)
        {
            var perFold = new List<Dictionary<SpectralMatch, double>>(peptideGroupIndices.Length);
            foreach (List<int> foldIndices in peptideGroupIndices)
            {
                // Best match per full sequence in this fold - the same units CreatePsmData labels on.
                List<SpectralMatch> foldMatches = foldIndices
                    .SelectMany(idx => peptideGroups[idx].GetBestMatches())
                    .Where(psm => psm != null)
                    .OrderBy(psm => psm.GetFdrInfo(UsePeptideLevelQValueForTraining).PEP)
                    .ToList();

                // First pass, best PEP to worst: running raw target/decoy FDR.
                double cumulativeTarget = 0;
                double cumulativeDecoy = 0;
                var rawFdr = new double[foldMatches.Count];
                for (int i = 0; i < foldMatches.Count; i++)
                {
                    if (foldMatches[i].IsDecoy)
                    {
                        cumulativeDecoy++;
                    }
                    else
                    {
                        cumulativeTarget++;
                    }
                    rawFdr[i] = cumulativeDecoy / Math.Max(cumulativeTarget, 1);
                }

                // Second pass, worst to best: monotonize the raw FDR into q-values.
                var qValueByMatch = new Dictionary<SpectralMatch, double>(foldMatches.Count);
                double qValue = 1.0;
                for (int i = foldMatches.Count - 1; i >= 0; i--)
                {
                    qValue = Math.Min(qValue, rawFdr[i]);
                    qValueByMatch[foldMatches[i]] = qValue;
                }
                perFold.Add(qValueByMatch);
            }
            return perFold;
        }

        /// <summary>
        /// Sets the following static properties: ChargeStateMode, FileSpecificMedianFragmentMassErrors, FileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, FileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, and FileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE
        /// </summary>
        /// <param name="trainingData"> The PSMs that will be used for training </param>
        /// <param name="trainingVariables"> An array of training variables from PsmData.trainingInfos dictionary </param>
        public void BuildFileSpecificDictionaries(List<SpectralMatch> trainingData, string[] trainingVariables)
        {
            FileSpecificMedianFragmentMassErrors = GetFileSpecificMedianFragmentMassError(trainingData);
            ChargeStateMode = GetChargeStateMode(trainingData);

            if (trainingVariables.Contains("HydrophobicityZScore"))
            {
                FileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = ComputeRetentionTimeEquivalentValues(trainingData, false, RetentionTimePredictor);
                FileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = ComputeRetentionTimeEquivalentValues(trainingData, true, RetentionTimePredictor);
                FileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE = ComputeMobilityValues(trainingData);
            }
            if (trainingVariables.Contains("ChimeraCount"))
            {
                chimeraCountDictionary = trainingData.GroupBy(p => p.ChimeraIdString).ToDictionary(g => g.Key, g => g.Count());
            }
        }

        public static List<int>[] GetPeptideGroupIndices(List<SpectralMatchGroup> peptides, int numGroups)
        {
            List<int>[] groupsOfIndices = new List<int>[numGroups];

            List<int> targetIndices = new List<int>();
            List<int> decoyIndices = new List<int>();
            for (int i = 0; i < peptides.Count; i++)
            {
                if (peptides[i].BestMatch.IsDecoy)
                {
                    decoyIndices.Add(i);
                }
                else
                {
                    targetIndices.Add(i);
                }
            }

            var targetIndexGroups = DivideListIntoGroups(targetIndices, numGroups);
            var decoyIndexGroups = DivideListIntoGroups(decoyIndices, numGroups);

            for (int i = 0; i < numGroups; i++)
            {
                groupsOfIndices[i] = targetIndexGroups[i].Concat(decoyIndexGroups[i]).ToList();
            }

            return groupsOfIndices;
        }

        /// <summary>
        /// This takes in a list of ints, and partitions them into numGroups partitions,
        /// e.g., partition 1 = [0, 4, 8...], partition 2 = [1, 5, 9...], etc.
        /// </summary>
        /// <returns>A list containing numGroups partitions (lists of ints) </returns>
        static List<List<int>> DivideListIntoGroups(List<int> list, int numGroups)
        {
            var groups = new List<List<int>>();
            for (int i = 0; i < numGroups; i++)
            {
                groups.Add(new List<int>());
            }

            int mainIndex = 0;
            while (mainIndex < list.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && mainIndex < list.Count)
                {
                    groups[subIndex].Add(list[mainIndex]);

                    subIndex++;
                    mainIndex++;
                }
            }

            return groups;
        }


        /// <summary>
        /// Builds the labeled training rows for one cross-fit fold. A non-decoy match is labeled
        /// positive when its q-value is at or below <see cref="QValueCutoff"/>. On the first training
        /// iteration <paramref name="pepQValueForRelabeling"/> is null and the raw-score q-value is
        /// used (the historical rule); on later iterations it carries this fold's PEP-derived
        /// q-values so positives are re-selected from the model's own ranking.
        /// </summary>
        public IEnumerable<PsmData> CreatePsmData(string searchType,
            List<SpectralMatchGroup> peptideGroups, List<int> peptideGroupIndices,
            Dictionary<SpectralMatch, double> pepQValueForRelabeling = null)
        {
            List<PsmData> psmDataList = new List<PsmData>();

            for (int i = 0; i < peptideGroupIndices.Count; i++)
            {
                int modCount = 0;
                foreach (var psm in peptideGroups[peptideGroupIndices[i]].GetBestMatches().Where(psm => psm != null))
                {
                    // First iteration: label by raw-score q-value. Later iterations: label by this
                    // fold's PEP-derived q-value from the previous fit. A match missing from the
                    // relabel map (should not happen) is treated as mid-confidence and excluded.
                    double qValueForLabeling = pepQValueForRelabeling == null
                        ? psm.GetFdrInfo(UsePeptideLevelQValueForTraining).QValue
                        : (pepQValueForRelabeling.TryGetValue(psm, out double pepQ) ? pepQ : double.MaxValue);

                    PsmData newPsmData = new PsmData();
                    if (searchType == "crosslink" && ((CrosslinkSpectralMatch)psm)?.BetaPeptide != null)
                    {
                        CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)psm;

                        bool label;
                        if (csm.IsDecoy || csm.BetaPeptide.IsDecoy)
                        {
                            label = false;
                            newPsmData = CreateOnePsmDataEntry(searchType, csm, csm.BestMatchingBioPolymersWithSetMods.First(), label);
                        }
                        else if (!csm.IsDecoy && !csm.BetaPeptide.IsDecoy && qValueForLabeling <= QValueCutoff)
                        {
                            label = true;
                            newPsmData = CreateOnePsmDataEntry(searchType, csm, csm.BestMatchingBioPolymersWithSetMods.First(), label);
                        }
                        else
                        {
                            continue;
                        }
                        psmDataList.Add(newPsmData);
                    }
                    else
                    {
                        double bmp = 0;
                        foreach (SpectralMatchHypothesis bestMatch in psm.BestMatchingBioPolymersWithSetMods)
                        {
                            bool label;
                            double bmpc = psm.BestMatchingBioPolymersWithSetMods.Count();
                            if (bestMatch.SpecificBioPolymer.Parent.IsDecoy)
                            {
                                label = false;
                                newPsmData = CreateOnePsmDataEntry(searchType, psm, bestMatch, label);
                            }
                            else if (!bestMatch.SpecificBioPolymer.Parent.IsDecoy
                                && qValueForLabeling <= QValueCutoff)
                            {
                                label = true;
                                newPsmData = CreateOnePsmDataEntry(searchType, psm, bestMatch, label);
                            }
                            else
                            {
                                continue;
                            }
                            psmDataList.Add(newPsmData);

                            bmp += 1.0;
                        }
                    }
                    modCount++;
                }
            }

            return psmDataList;
        }

        public static string AggregateMetricsForOutput(List<CalibratedBinaryClassificationMetrics> allMetrics, int sumOfAllAmbiguousPeptidesResolved,
            int positiveTrainingCount, int negativeTrainingCount, double qValueCutoff, List<int> positivePoolPerIteration = null)

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
            s.AppendLine("*       Accuracy:  " + accuracy.Average());
            s.AppendLine("*       Area Under Curve:  " + areaUnderRocCurve.Average());
            s.AppendLine("*       Area under Precision recall Curve:  " + areaUnderPrecisionRecallCurve.Average());
            s.AppendLine("*       F1Score:  " + F1Score.Average());
            s.AppendLine("*       LogLoss:  " + logLossAverage);
            s.AppendLine("*       LogLossReduction:  " + logLossReductionAverage);
            s.AppendLine("*       PositivePrecision:  " + positivePrecision.Average());
            s.AppendLine("*       PositiveRecall:  " + positiveRecall.Average());
            s.AppendLine("*       NegativePrecision:  " + negativePrecision.Average());
            s.AppendLine("*       NegativeRecall:  " + negativeRecall.Average());
            s.AppendLine($"*       Count of Ambiguous {char.ToUpper(GlobalVariables.AnalyteType.GetUniqueFormLabel()[0]) + GlobalVariables.AnalyteType.GetUniqueFormLabel()[1..]}s Removed:  " + sumOfAllAmbiguousPeptidesResolved);
            s.AppendLine("*       Q-Value Cutoff for Training Targets:  " + qValueCutoff);
            s.AppendLine("*       Targets Used for Training:  " + positiveTrainingCount);
            s.AppendLine("*       Decoys Used for Training:  " + negativeTrainingCount);
            if (positivePoolPerIteration != null && positivePoolPerIteration.Count > 0)
            {
                s.AppendLine("*       Training Iterations Run:  " + positivePoolPerIteration.Count);
                s.AppendLine("*       Positive Training Pool Per Iteration:  " + string.Join(", ", positivePoolPerIteration));
            }
            s.AppendLine("************************************************************");
            return s.ToString();
        }

        /// <summary>
        /// Applies a trained cross-fit model to a held-out fold, writing each PSM's PEP. When
        /// <paramref name="removeAmbiguous"/> is true, ambiguous peptide hypotheses scoring well
        /// below the best are also pruned from the PSM. Pruning mutates the PSM, so iterative
        /// training defers it to the final pass; the PEP value itself is unaffected by the flag.
        /// </summary>
        public int Compute_PSM_PEP(List<SpectralMatchGroup> peptideGroups,
            List<int> peptideGroupIndices,
            MLContext mLContext, TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>> trainedModel, string searchType, string outputFolder, bool removeAmbiguous = true)
        {
            int maxThreads = FileSpecificParametersDictionary.Values.FirstOrDefault().MaxThreadsToUsePerFile;
            int ambiguousPeptidesResolved = 0;

            var predictionEnginePerThread =
                new ThreadLocal<PredictionEngine<PsmData, TruePositivePrediction>>(
                    () => mLContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(trainedModel),
                    trackAllValues: true);

            try
            {
                Parallel.ForEach(Partitioner.Create(0, peptideGroupIndices.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                    (range, loopState) =>
                    {
                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops) { return; }

                        // one prediction engine per thread, because the prediction engine is not thread-safe
                        var threadPredictionEngine = predictionEnginePerThread.Value;

                        int ambigousPeptidesRemovedinThread = 0;

                        List<int> indiciesOfPeptidesToRemove = new List<int>();
                        List<double> pepValuePredictions = new List<double>();
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            foreach (SpectralMatch psm in peptideGroups[peptideGroupIndices[i]])
                            {
                                // I'm not sure what's going one here vis-a-vis disambiguations, but I'm not going to touch it for now
                                if (psm != null)
                                {
                                    indiciesOfPeptidesToRemove.Clear();
                                    pepValuePredictions.Clear();

                                    //Here we compute the pepvalue predection for each ambiguous peptide in a PSM. Ambiguous peptides with lower pepvalue predictions are removed from the PSM.
                                    var bestMatchingBioPolymersWithSetMods = psm.BestMatchingBioPolymersWithSetMods.ToList();
                                    foreach (SpectralMatchHypothesis bestMatch in bestMatchingBioPolymersWithSetMods)
                                    {
                                        PsmData pd = CreateOnePsmDataEntry(searchType, psm, bestMatch, !bestMatch.IsDecoy);
                                        var pepValuePrediction = threadPredictionEngine.Predict(pd);
                                        pepValuePredictions.Add(pepValuePrediction.Probability);
                                        //A score is available using the variable pepvaluePrediction.Score
                                    }

                                    // PEP is 1 - the best hypothesis probability. Pruning only ever
                                    // drops hypotheses below that best, so it leaves the PEP unchanged
                                    // and can be safely skipped on non-final iterations.
                                    psm.PsmFdrInfo.PEP = 1 - pepValuePredictions.Max();
                                    psm.PeptideFdrInfo.PEP = 1 - pepValuePredictions.Max();

                                    if (removeAmbiguous)
                                    {
                                        GetIndiciesOfPeptidesToRemove(indiciesOfPeptidesToRemove, pepValuePredictions);
                                        RemoveBestMatchingPeptidesWithLowPEP(psm, indiciesOfPeptidesToRemove, bestMatchingBioPolymersWithSetMods, ref ambigousPeptidesRemovedinThread);
                                    }
                                }

                            }
                        }

                        Interlocked.Add(ref ambiguousPeptidesResolved, ambigousPeptidesRemovedinThread);
                    });
            }
            finally
            {
                foreach (var engine in predictionEnginePerThread.Values)
                {
                    engine?.Dispose();
                }

                predictionEnginePerThread.Dispose();
            }

            return ambiguousPeptidesResolved;
        }

        public PsmData CreateOnePsmDataEntry(string searchType, SpectralMatch psm, SpectralMatchHypothesis tentativeSpectralMatch, bool label)
        {
            double normalizationFactor = tentativeSpectralMatch.SpecificBioPolymer.BaseSequence.Length;
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
                totalMatchingFragmentCount = (float)(Math.Round(tentativeSpectralMatch.MatchedIons.Count(p => p.NeutralTheoreticalProduct.SecondaryProductType == null) / normalizationFactor * multiplier, 0));
                internalMatchingFragmentCount = (float)(Math.Round(tentativeSpectralMatch.MatchedIons.Count(p => p.NeutralTheoreticalProduct.SecondaryProductType != null) / normalizationFactor * multiplier, 0));
                intensity = (float)Math.Min(50, Math.Round((psm.Score - (int)psm.Score) / normalizationFactor * Math.Pow(multiplier, 2), 0));
                chargeDifference = -Math.Abs(ChargeStateMode - psm.ScanPrecursorCharge);
                deltaScore = (float)Math.Round(psm.DeltaScore / normalizationFactor * multiplier, 0);
                notch = tentativeSpectralMatch.Notch;
                modCount = Math.Min((float)tentativeSpectralMatch.SpecificBioPolymer.AllModsOneIsNterminus.Keys.Count(), 10);
                if (tentativeSpectralMatch.MatchedIons?.Count() > 0)
                {
                    absoluteFragmentMassError = (float)Math.Min(100.0, Math.Round(10.0 * Math.Abs(GetAverageFragmentMassError(tentativeSpectralMatch.MatchedIons) - FileSpecificMedianFragmentMassErrors[Path.GetFileName(psm.FullFilePath)])));
                }

                ambiguity = Math.Min((float)(psm.BestMatchingBioPolymersWithSetMods.Count() - 1), 10);
                //ambiguity = 10; // I'm pretty sure that you shouldn't train on ambiguity and its skewing the results
                longestSeq = (float)Math.Round(SpectralMatch.GetLongestIonSeriesBidirectional(tentativeSpectralMatch) / normalizationFactor * multiplier, 0);
                complementaryIonCount = (float)Math.Round(SpectralMatch.GetCountComplementaryIons(tentativeSpectralMatch) / normalizationFactor * multiplier, 0);
                isVariantPeptide = PeptideIsVariant(tentativeSpectralMatch.SpecificBioPolymer);
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

                if (psm.DigestionParams.DigestionAgent.Name != "top-down")
                {
                    missedCleavages = tentativeSpectralMatch.SpecificBioPolymer.MissedCleavages;
                    var fileName = Path.GetFileName(psm.FullFilePath);
                    bool fileIsCzeSeparationType = FileSpecificParametersDictionary.TryGetValue(fileName, out var fileParams) && fileParams.SeparationType == "CZE";

                    if (searchType != "RNA")
                    {
                        if (!fileIsCzeSeparationType)
                        {
                            var isUnmodified = tentativeSpectralMatch.SpecificBioPolymer.BaseSequence.Equals(tentativeSpectralMatch.SpecificBioPolymer.FullSequence);
                            var dict = isUnmodified
                                ? FileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified
                                : FileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified;
                            hydrophobicityZscore = (float)Math.Round(GetRetentionTimeEquivalentZscore(psm, tentativeSpectralMatch.SpecificBioPolymer, dict, RetentionTimePredictor) * 10.0, 0);
                        }
                        else
                        {
                            hydrophobicityZscore = (float)Math.Round(GetMobilityZScore(psm, tentativeSpectralMatch.SpecificBioPolymer) * 10.0, 0);
                        }
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
                var selectedAlphaPeptide = csm.BestMatchingBioPolymersWithSetMods.First();
                var selectedBetaPeptide = csm.BetaPeptide?.BestMatchingBioPolymersWithSetMods.First();

                float alphaNormalizationFactor = selectedAlphaPeptide.SpecificBioPolymer.BaseSequence.Length;
                float betaNormalizationFactor = selectedBetaPeptide == null ? (float)0 : selectedBetaPeptide.SpecificBioPolymer.BaseSequence.Length;
                float totalNormalizationFactor = alphaNormalizationFactor + betaNormalizationFactor;

                totalMatchingFragmentCount = (float)Math.Round(csm.XLTotalScore / totalNormalizationFactor * 10, 0);

                //Compute fragment mass error
                int alphaCount = 0;
                float alphaError = 0;
                if (selectedAlphaPeptide.MatchedIons?.Count > 0)
                {
                    alphaCount = selectedAlphaPeptide.MatchedIons.Count;
                    alphaError = Math.Abs(GetAverageFragmentMassError(selectedAlphaPeptide.MatchedIons));
                }
                int betaCount = 0;
                float betaError = 0;
                if (selectedBetaPeptide != null && selectedBetaPeptide.MatchedIons?.Count > 0)
                {
                    betaCount = selectedBetaPeptide.MatchedIons.Count;
                    betaError = Math.Abs(GetAverageFragmentMassError(selectedBetaPeptide.MatchedIons));
                }

                float averageError = 0;
                if ((alphaCount + betaCount) > 0)
                {
                    averageError = (alphaCount * alphaError + betaCount * betaError) / (alphaCount + betaCount);
                }

                absoluteFragmentMassError = (float)Math.Min(100, Math.Round(averageError - FileSpecificMedianFragmentMassErrors[Path.GetFileName(csm.FullFilePath)] * 10.0, 0));
                //End compute fragment mass error

                deltaScore = (float)Math.Round(csm.DeltaScore / totalNormalizationFactor * 10.0, 0);
                chargeDifference = -Math.Abs(ChargeStateMode - psm.ScanPrecursorCharge);
                alphaIntensity = (float)Math.Min(100, Math.Round((csm.Score - (int)csm.Score) / alphaNormalizationFactor * 100.0, 0));
                betaIntensity = csm.BetaPeptide == null ? (float)0 : (float)Math.Min(100.0, Math.Round((csm.BetaPeptide.Score - (int)csm.BetaPeptide.Score) / betaNormalizationFactor * 100.0, 0));
                longestFragmentIonSeries_Alpha = (float)Math.Round(SpectralMatch.GetLongestIonSeriesBidirectional(selectedAlphaPeptide) / alphaNormalizationFactor * 10.0, 0);
                longestFragmentIonSeries_Beta = selectedBetaPeptide == null ? (float)0 : SpectralMatch.GetLongestIonSeriesBidirectional(selectedBetaPeptide) / betaNormalizationFactor;
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

        public static void RemoveBestMatchingPeptidesWithLowPEP(SpectralMatch psm, List<int> indiciesOfPeptidesToRemove, List<SpectralMatchHypothesis> allPeptides, ref int ambiguousPeptidesRemovedCount)
        {
            int peptidesRemoved = 0;
            foreach (var toRemove in indiciesOfPeptidesToRemove)
            {
                psm.RemoveThisAmbiguousPeptide(allPeptides[toRemove - peptidesRemoved]);
                peptidesRemoved++;
            }
            ambiguousPeptidesRemovedCount += peptidesRemoved;
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

        #region Dictionary Builder Functions and Utilities

        /// <summary>
        /// Here we're getting the most common charge state for precursors that are Targets with q<=0.01.

        public int GetChargeStateMode(List<SpectralMatch> psms)
        {
            return psms.Where(p => p.IsDecoy != true && p.GetFdrInfo(UsePeptideLevelQValueForTraining).QValue <= 0.01).Select(p => p.ScanPrecursorCharge).GroupBy(n => n).OrderByDescending(g => g.Count()).Select(g => g.Key).FirstOrDefault();
        }

        public Dictionary<string, Dictionary<int, Tuple<double, double>>> ComputeRetentionTimeEquivalentValues(List<SpectralMatch> psms, bool computeHydrophobicitiesforModifiedPeptides, IRetentionTimePredictor predictor)
        {
            //TODO change the tuple so the values have names
            Dictionary<string, Dictionary<int, Tuple<double, double>>> rtHydrophobicityAvgDev = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            List<string> filenames = FileSpecificParametersDictionary.Select(kvp => Path.GetFileName(kvp.Key)).ToList();

            filenames = filenames.Distinct().ToList();

            foreach (string filename in filenames)
            {
                Dictionary<int, List<double>> hydrophobicities = new Dictionary<int, List<double>>();
                Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations = new Dictionary<int, Tuple<double, double>>();

                foreach (SpectralMatch psm in psms.Where(f => (f.FullFilePath == null || Path.GetFileName(f.FullFilePath) == filename) && f.FdrInfo.QValue <= 0.01 && !f.IsDecoy))
                {
                    List<string> fullSequences = new List<string>();
                    foreach (SpectralMatchHypothesis bestMatch in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        if (fullSequences.Contains(bestMatch.SpecificBioPolymer.FullSequence))
                        {
                            continue;
                        }
                        fullSequences.Add(bestMatch.SpecificBioPolymer.FullSequence);

                        double predictedHydrophobicity = bestMatch.SpecificBioPolymer is PeptideWithSetModifications pep ? predictor.PredictRetentionTimeEquivalent(pep, out _) ?? 0 : 0;

                        //here i'm grouping this in 2 minute increments becuase there are cases where you get too few data points to get a good standard deviation an average. This is for stability.
                        int possibleKey = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));

                        //First block of if statement is for modified peptides.
                        if (bestMatch.SpecificBioPolymer.AllModsOneIsNterminus.Any() && computeHydrophobicitiesforModifiedPeptides)
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
                        else if (!bestMatch.SpecificBioPolymer.AllModsOneIsNterminus.Any() && !computeHydrophobicitiesforModifiedPeptides)
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

        public Dictionary<string, Dictionary<int, Tuple<double, double>>> ComputeMobilityValues(List<SpectralMatch> psms)
        {
            Dictionary<string, Dictionary<int, Tuple<double, double>>> rtMobilityAvgDev = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            List<string> filenames = FileSpecificParametersDictionary.Select(kvp => Path.GetFileName(kvp.Key)).ToList();

            filenames = filenames.Distinct().ToList();

            foreach (string filename in filenames)
            {
                Dictionary<int, List<double>> mobilities = new Dictionary<int, List<double>>();
                Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations = new Dictionary<int, Tuple<double, double>>();

                foreach (SpectralMatch psm in psms.Where(f => (f.FullFilePath == null || Path.GetFileName(f.FullFilePath) == filename) && f.FdrInfo.QValue <= 0.01 && !f.IsDecoy))
                {
                    List<string> fullSequences = new List<string>();
                    foreach (SpectralMatchHypothesis bestMatch in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        if (fullSequences.Contains(bestMatch.SpecificBioPolymer.FullSequence))
                        {
                            continue;
                        }
                        fullSequences.Add(bestMatch.SpecificBioPolymer.FullSequence);

                        double predictedMobility = bestMatch.SpecificBioPolymer is PeptideWithSetModifications pep ? 100.0 * GetCifuentesMobility(pep) : 0;

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

        private static float GetRetentionTimeEquivalentZscore(SpectralMatch psm, IBioPolymerWithSetMods Peptide, Dictionary<string, Dictionary<int, Tuple<double, double>>> d, IRetentionTimePredictor predictor)
        {
            //Using SSRCalc3 but probably any number of different calculators could be used instead. One could also use the CE mobility.
            double hydrophobicityZscore = double.NaN;

            if (d.ContainsKey(Path.GetFileName(psm.FullFilePath)))
            {
                int time = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));
                if (d[Path.GetFileName(psm.FullFilePath)].Keys.Contains(time))
                {
                    double predictedHydrophobicity = Peptide is PeptideWithSetModifications pep ? predictor.PredictRetentionTimeEquivalent(pep, out _) ?? 0 : 0;

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

        private float GetMobilityZScore(SpectralMatch psm, IBioPolymerWithSetMods selectedPeptide)
        {
            double mobilityZScore = double.NaN;

            if (FileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE.ContainsKey(Path.GetFileName(psm.FullFilePath)))
            {
                int time = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));
                if (FileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[Path.GetFileName(psm.FullFilePath)].Keys.Contains(time))
                {
                    double predictedMobility = 100.0 * GetCifuentesMobility(selectedPeptide);

                    mobilityZScore = Math.Abs(FileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[Path.GetFileName(psm.FullFilePath)][time].Item1 - predictedMobility) / FileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[Path.GetFileName(psm.FullFilePath)][time].Item2;
                }
            }

            double maxMobilityZscore = 10; // each "Z" is one standard deviation. so, maxHydrophobicityZscore 10 is quite large
            if (double.IsNaN(mobilityZScore) || double.IsInfinity(mobilityZScore) || mobilityZScore > maxMobilityZscore)
            {
                mobilityZScore = maxMobilityZscore;
            }

            return (float)mobilityZScore;
        }

        private static bool PeptideIsVariant(IBioPolymerWithSetMods bpwsm)
        {
            if (bpwsm is not PeptideWithSetModifications pwsm)
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
                    foreach (var bestMatch in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        if (bestMatch.MatchedIons is { Count: > 0 })
                        {
                            averageMassErrors.Add(GetAverageFragmentMassError(bestMatch.MatchedIons));
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

        #endregion
    }
}
