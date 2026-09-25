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
using System.Runtime.CompilerServices;
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

        /// <summary>
        /// The cap on training rounds used when iterative training is switched on
        /// (<c>SearchParameters.IterativePepTraining</c>). A safety stop, not the number of rounds that
        /// will run: <see cref="TrainingImprovementTolerance"/> decides that.
        /// </summary>
        public const int IterativeTrainingRoundCap = 10;

        /// <summary>
        /// Upper bound on training rounds. The default of 1 trains once, on labels from the search-score
        /// q-value, which is the algorithm this engine had before iteration existed. Values above 1 turn on
        /// semi-supervised iteration, and the stopping criterion then decides how many rounds actually run.
        /// </summary>
        public int MaxTrainingRounds { get; set; } = 1;

        /// <summary>
        /// Stop when a round grows the accepted target count by less than this fraction. 0.001 = 0.1%.
        /// </summary>
        public double TrainingImprovementTolerance { get; set; } = 0.001;

        /// <summary>
        /// The PEP q-value below which a target peptide counts as "accepted" when judging whether a round
        /// helped. The count is computed the way <see cref="FdrAnalysisEngine"/> computes the reported
        /// peptide PEP q-value (see <see cref="CountAcceptedPeptides"/>), so the criterion tracks what users see.
        /// </summary>
        public double AcceptanceQValueCutoff { get; set; } = 0.01;

        /// <summary>
        /// Picks a fold's positives for the next round from its own model's PEPs on its own training
        /// matches. <see cref="SelectPositives"/> at <see cref="QValueCutoff"/>. Settable only so tests can
        /// reach the path where a later round has no positives.
        /// </summary>
        internal Func<List<SpectralMatch>, double[], HashSet<SpectralMatch>> SelectNextRoundPositives { get; set; }

        /// <summary>
        /// Feature vectors, computed once and reused for every training round.
        /// <remarks>
        /// EVERY field of <see cref="PsmData"/> except Label is round-invariant: each one derives from
        /// the match, the hypothesis, or a dictionary built once in the constructor
        /// (FileSpecificMedianFragmentMassErrors, ChargeStateMode, chimeraCountDictionary, the
        /// hydrophobicity tables). Retention-time prediction is deterministic from those fixed inputs.
        /// So recomputing them per round re-derives bit-identical values -- and RT prediction is the
        /// single most expensive thing this engine does (~67% of an entire search before #2841).
        ///
        /// This is only CORRECT because PEP no longer prunes. `Ambiguity` reads
        /// BestMatchingBioPolymersWithSetMods.Count(), which the old elimination step mutated
        /// mid-loop; caching across rounds would have frozen a stale value. Removing the pruning is
        /// what makes the cache sound.
        ///
        /// The cached instances are never handed out: callers get a copy via <see cref="PsmData.WithLabel"/>.
        /// Memory: one PsmData (~200 B) per hypothesis for the engine's lifetime. That is tens of MB for a
        /// few hundred thousand matches and some hundreds of MB for multi-million-match runs. The entry
        /// count is written to pep_training_rounds.txt.
        /// </remarks>
        /// </summary>
        private readonly ConcurrentDictionary<SpectralMatchHypothesis, PsmData> _featureCache
            = new ConcurrentDictionary<SpectralMatchHypothesis, PsmData>(new HypothesisReferenceComparer());

        /// <summary>
        /// Identity comparer for hypotheses used as cache keys.
        /// <remarks>
        /// SpectralMatchHypothesis's OWN equality cannot be used here: Equals compares
        /// MatchedIons.Count while GetHashCode hashes the MatchedIons reference, so two instances can
        /// be Equal yet hash differently. That violates the hash/equals contract and would make
        /// dictionary lookups miss unpredictably. Identity is what we actually want anyway -- each
        /// hypothesis instance belongs to exactly one match, so the reference identifies the pair.
        /// </remarks>
        /// </summary>
        private sealed class HypothesisReferenceComparer : IEqualityComparer<SpectralMatchHypothesis>
        {
            public bool Equals(SpectralMatchHypothesis x, SpectralMatchHypothesis y) => ReferenceEquals(x, y);

            public int GetHashCode(SpectralMatchHypothesis obj) => RuntimeHelpers.GetHashCode(obj);
        }

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
        /// When true, PEP also removes ambiguous match hypotheses predicted well below the best one. Only for
        /// callers with no DisambiguationEngine downstream (glyco, crosslink, and nonspecific or semi-specific searches); a classic or modern SearchTask leaves it false.
        /// </summary>
        public bool PruneAmbiguousHypotheses { get; }
        private int _ambiguousHypothesesRemoved;

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

        public PepAnalysisEngine(List<SpectralMatch> psms, string searchType, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, string outputFolder, IRetentionTimePredictor? rtPredictor = null, bool pruneAmbiguousHypotheses = false)
        {
            // This creates a new list of PSMs, but does not clone the Psms themselves.
            // This allows the PSMs to be modified and the order to be preserved
            AllPsms = psms.OrderByDescending(p => p).ToList();
            PruneAmbiguousHypotheses = pruneAmbiguousHypotheses;
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
            SelectNextRoundPositives = (matches, pep) => SelectPositives(matches, pep, QValueCutoff);
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
            // Pruning deletes hypotheses irreversibly, so a later, better round could not bring them back.
            // Callers that prune (glyco, crosslink) therefore train once, exactly as before iteration existed.
            int maxRounds = PruneAmbiguousHypotheses ? 1 : Math.Max(1, MaxTrainingRounds);

            // Timings go to pep_training_rounds.txt, never to the results block (see below).
            WriteRoundProgress($"=== PEP {DateTime.Now:yyyy-MM-dd HH:mm:ss}  search type {SearchType}, " +
                $"digestion {string.Join("+", AllPsms.Select(p => p.DigestionParams?.DigestionAgent?.Name).Distinct())}, " +
                $"{AllPsms.Count} matches, {peptideGroups.Count} training groups " +
                $"({(UsePeptideLevelQValueForTraining ? "peptide" : "PSM")} level), max rounds {maxRounds} ===");
            var roundClock = System.Diagnostics.Stopwatch.StartNew();

            // Round 0's labels come from the search-score q-value. They depend on no model and no fold, so the
            // four groups are built once, in parallel, before anything trains -- as before iteration existed.
            IEnumerable<PsmData>[] roundZeroData = new IEnumerable<PsmData>[numGroups];
            bool allGroupsHavePositiveAndNegativeTrainingExamples = true;
            Parallel.ForEach(
                Enumerable.Range(0, numGroups),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                group =>
                {
                    roundZeroData[group] = CreatePsmData(SearchType, peptideGroups, peptideGroupIndices[group]);
                    if (!HasPositiveAndNegativeExamples(roundZeroData[group]))
                    {
                        allGroupsHavePositiveAndNegativeTrainingExamples = false;
                    }
                });
            // Checked for every group, before any fold trains or scores, so a failure leaves no PEP half-assigned.
            // It also guarantees that each held-out fold, which is evaluated against these same labels, has both classes.
            if (!allGroupsHavePositiveAndNegativeTrainingExamples)
            {
                return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
            }
            WriteRoundProgress($"round 0 features built  ({roundClock.Elapsed.TotalSeconds:F1} s)");
            roundClock.Restart();

            MLContext mlContext = new MLContext(seed: _randomSeed);
            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(BGDTreeOptions);
            var pipeline = mlContext.Transforms.Concatenate("Features", TrainingVariables)
                .Append(trainer);

            // What each fold's model trains on: the three groups it does not score.
            var trainingGroupIndices = new List<int>[numGroups];
            var trainingData = new List<PsmData>[numGroups];
            for (int fold = 0; fold < numGroups; fold++)
            {
                trainingGroupIndices[fold] = Enumerable.Range(0, numGroups).Where(g => g != fold)
                    .SelectMany(g => peptideGroupIndices[g]).ToList();
                //concat doesn't work in a loop, therefore I had to hard code the concat to group 3 out of 4 lists. if the const int numGroups value is changed, then the concat has to be changed accordingly.
                var others = Enumerable.Range(0, numGroups).Where(g => g != fold).ToList();
                trainingData[fold] = roundZeroData[others[0]].Concat(roundZeroData[others[1]].Concat(roundZeroData[others[2]])).ToList();
            }

            List<CalibratedBinaryClassificationMetrics> allMetrics = new List<CalibratedBinaryClassificationMetrics>();
            int positiveTrainingCount = roundZeroData.SelectMany(p => p).Count(p => p.Label);
            int negativeTrainingcount = roundZeroData.SelectMany(p => p).Count(p => !p.Label);
            int roundsRun = 0;
            var roundLog = new StringBuilder();
            int previousAccepted = 0;
            // Snapshot so a round that makes things worse can be undone. One double per PSM.
            double[] bestPepSnapshot = null;

            for (int round = 0; round < maxRounds; round++)
            {
                if (round > 0)
                {
                    // FOLD-LOCAL, MODEL-LOCAL LABELS. Each fold's next positives come from ITS OWN previous
                    // model's scores on ITS OWN training folds, as in mokapot's brew(). No model that saw a
                    // fold's held-out data takes part in that fold's labels, and no held-out score is read.
                    // All four label sets are derived before any model is refitted or any PEP overwritten,
                    // so the result does not depend on the order the folds run in.
                    bool foldStarved = false;
                    var nextTrainingData = new List<PsmData>[numGroups];
                    for (int fold = 0; fold < numGroups && !foldStarved; fold++)
                    {
                        var trainingMatches = trainingGroupIndices[fold]
                            .SelectMany(i => peptideGroups[i].GetBestMatches())
                            .Where(m => m != null)
                            .ToList();
                        double[] ownPep = PredictPep(mlContext, trainedModels[fold], trainingMatches, maxThreads);
                        HashSet<SpectralMatch> positives = SelectNextRoundPositives(trainingMatches, ownPep);
                        nextTrainingData[fold] = CreatePsmData(SearchType, peptideGroups, trainingGroupIndices[fold], positives.Contains);
                        foldStarved = !HasPositiveAndNegativeExamples(nextTrainingData[fold]);
                    }

                    if (foldStarved)
                    {
                        // Nothing has been refitted or rescored this round, so the previous round's PEPs stand.
                        WriteRoundProgress($"round {round}: a fold has no positive or no negative examples; keeping round {round - 1}");
                        break;
                    }

                    trainingData = nextTrainingData;
                }

                var roundMetrics = new List<CalibratedBinaryClassificationMetrics>();
                for (int fold = 0; fold < numGroups; fold++)
                {
                    trainedModels[fold] = pipeline.Fit(mlContext.Data.LoadFromEnumerable(trainingData[fold]));

                    // Evaluated against the round-0 labels in every round: a fixed yardstick, the one master used,
                    // and one that the all-groups check above guarantees has both classes.
                    var myPredictions = trainedModels[fold].Transform(mlContext.Data.LoadFromEnumerable(roundZeroData[fold]));
                    roundMetrics.Add(mlContext.BinaryClassification.Evaluate(data: myPredictions, labelColumnName: "Label", scoreColumnName: "Score"));

                    //model is trained on peptides but here we can use that to compute PEP for all PSMs
                    Compute_PSM_PEP(peptideGroups, peptideGroupIndices[fold], mlContext, trainedModels[fold], SearchType, OutputFolder);
                }

                // Only after every fold of this round has scored. REPORTING and the stopping decision only --
                // never labels.
                int accepted = CountAcceptedPeptides(AllPsms, AcceptanceQValueCutoff);
                var verdict = JudgeRound(round, accepted, previousAccepted, TrainingImprovementTolerance, maxRounds);
                WriteRoundProgress($"round {round}: accepted {accepted}  {verdict}  ({roundClock.Elapsed.TotalSeconds:F1} s)");
                roundClock.Restart();

                if (verdict == RoundVerdict.RevertAndStop)
                {
                    // No better than the round before. Keep the better one -- an unguarded loop can make a
                    // search worse and nothing downstream would notice.
                    RestorePepValues(bestPepSnapshot);
                    break;
                }

                if (maxRounds > 1)
                {
                    bestPepSnapshot = SnapshotPepValues();
                }
                allMetrics = roundMetrics;
                if (round > 0)
                {
                    // Averaged over the folds: each peptide appears in the training set of 3 of the 4 models.
                    positiveTrainingCount = trainingData.Sum(d => d.Count(p => p.Label)) / (numGroups - 1);
                    negativeTrainingcount = trainingData.Sum(d => d.Count(p => !p.Label)) / (numGroups - 1);
                }
                roundsRun = round + 1;
                previousAccepted = accepted;
                // NO TIMINGS in the results block. It is compared for byte equality across two runs of
                // the same data by PepAnalysisEngineHasReproducibleOutput, and a wall clock is not
                // reproducible. Timings go to pep_training_rounds.txt, which nothing asserts on.
                roundLog.AppendLine($"*         round {round}:  accepted {accepted}");

                if (verdict == RoundVerdict.KeepAndStop)
                {
                    break;
                }
            }

            WriteRoundProgress($"done: {roundsRun} round(s) kept; feature cache holds {_featureCache.Count} vectors");

            return AggregateMetricsForOutput(allMetrics, positiveTrainingCount, negativeTrainingcount, QValueCutoff,
                PruneAmbiguousHypotheses ? _ambiguousHypothesesRemoved : null, roundsRun, previousAccepted, roundLog.ToString());
        }

        /// <summary>
        /// What to do after a training round has scored every fold.
        /// </summary>
        internal enum RoundVerdict
        {
            /// <summary>Keep this round and train another.</summary>
            Continue,
            /// <summary>Keep this round and stop: the cap is reached or the gain fell below the tolerance.</summary>
            KeepAndStop,
            /// <summary>This round accepted no more than the last one. Restore the last one and stop.</summary>
            RevertAndStop
        }

        /// <summary>
        /// The stopping rule, as a pure function so that every branch can be tested without training a model.
        /// A round is kept only if it accepts MORE targets than the round before, and the loop continues only
        /// while the relative gain is at least <paramref name="tolerance"/> and the cap allows another round.
        /// </summary>
        internal static RoundVerdict JudgeRound(int round, int accepted, int previousAccepted, double tolerance, int maxRounds)
        {
            if (round > 0 && accepted <= previousAccepted)
            {
                return RoundVerdict.RevertAndStop;
            }

            if (round + 1 >= maxRounds)
            {
                return RoundVerdict.KeepAndStop;
            }

            double improvement = round == 0 || previousAccepted == 0
                ? double.PositiveInfinity
                : (accepted - previousAccepted) / (double)previousAccepted;
            return improvement < tolerance ? RoundVerdict.KeepAndStop : RoundVerdict.Continue;
        }

        private static bool HasPositiveAndNegativeExamples(IEnumerable<PsmData> data)
        {
            return data.Any(p => p.Label) && data.Any(p => !p.Label);
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

        /// <summary>
        /// Appends one line of per-round progress to <c>pep_training_rounds.txt</c> in the output
        /// folder, flushed immediately so it can be watched while the engine is still running.
        /// Best-effort: progress reporting must never take a search down.
        /// </summary>
        private void WriteRoundProgress(string line)
        {
            if (string.IsNullOrEmpty(OutputFolder))
            {
                return;
            }

            try
            {
                File.AppendAllText(Path.Combine(OutputFolder, "pep_training_rounds.txt"),
                    $"{DateTime.Now:HH:mm:ss}  {line}{Environment.NewLine}");
            }
            catch (IOException)
            {
                // A locked or unwritable output folder is not a reason to fail the analysis.
            }
            catch (UnauthorizedAccessException)
            {
            }
        }

        /// <summary>
        /// Captures the PEP currently assigned to every match, so a round that turns out to be worse
        /// than its predecessor can be undone. One double per PSM.
        /// </summary>
        private double[] SnapshotPepValues()
        {
            var snapshot = new double[AllPsms.Count];
            for (int i = 0; i < AllPsms.Count; i++)
            {
                snapshot[i] = AllPsms[i].PsmFdrInfo.PEP;
            }

            return snapshot;
        }

        private void RestorePepValues(double[] snapshot)
        {
            if (snapshot == null)
            {
                return;
            }

            for (int i = 0; i < AllPsms.Count && i < snapshot.Length; i++)
            {
                AllPsms[i].PsmFdrInfo.PEP = snapshot[i];
                AllPsms[i].PeptideFdrInfo.PEP = snapshot[i];
            }
        }

        /// <summary>
        /// Round 0's rule for a positive training example: the SEARCH-SCORE q-value, which is all this
        /// engine used before iteration existed. Later rounds use <see cref="SelectPositives"/> instead.
        /// </summary>
        private bool IsPositiveBySearchScore(SpectralMatch psm)
        {
            return psm.GetFdrInfo(UsePeptideLevelQValueForTraining).QValue <= QValueCutoff;
        }

        /// <summary>
        /// PEP (1 - the best hypothesis's predicted probability) of each match under the given model, written
        /// to an array and NOT to the match. Used to relabel a model's own training folds, which must never
        /// disturb the PEPs other folds have assigned.
        /// </summary>
        private double[] PredictPep(MLContext mLContext,
            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>> trainedModel,
            List<SpectralMatch> matches, int maxThreads)
        {
            var pep = new double[matches.Count];
            var predictionEnginePerThread =
                new ThreadLocal<PredictionEngine<PsmData, TruePositivePrediction>>(
                    () => mLContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(trainedModel),
                    trackAllValues: true);
            try
            {
                Parallel.ForEach(Partitioner.Create(0, matches.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                    range =>
                    {
                        var threadPredictionEngine = predictionEnginePerThread.Value;
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            double best = double.MinValue;
                            foreach (SpectralMatchHypothesis hypothesis in matches[i].BestMatchingBioPolymersWithSetMods)
                            {
                                best = Math.Max(best, threadPredictionEngine.Predict(GetFeatures(SearchType, matches[i], hypothesis)).Probability);
                            }
                            pep[i] = 1 - best;
                        }
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

            return pep;
        }

        /// <summary>
        /// The next round's positive training examples: walk the matches in order of the model's own PEP,
        /// accumulate a target-decoy q-value, and take every target ranked at or above the last target
        /// whose q-value is within <paramref name="qValueCutoff"/>.
        /// <remarks>
        /// The cut is by RANK, not by PEP value. A threshold (PEP &lt;= t) would also admit every target
        /// tied at t, including tied targets the sweep had ordered past the point where q crosses the cutoff;
        /// ties are common where the calibrated output saturates near 0. Within a tie the order is the
        /// default SpectralMatch comparer's, as in <see cref="FdrAnalysisEngine"/>.
        ///
        /// The population is GetBestMatches(), one per full sequence: the unit the model is trained on.
        /// SpectralMatchGroup.BestMatch (one per BASE sequence) would silently drop every modified variant.
        /// </remarks>
        /// </summary>
        internal static HashSet<SpectralMatch> SelectPositives(List<SpectralMatch> matches, double[] pep, double qValueCutoff)
        {
            int[] order = Enumerable.Range(0, matches.Count)
                .OrderBy(i => pep[i])
                .ThenByDescending(i => matches[i])
                .ToArray();

            // q must be monotone in rank, so the running ratio is swept once forward and then minimised
            // from the bottom up -- the same shape as FdrAnalysisEngine.QValueInverted.
            var q = new double[order.Length];
            double cumulativeTarget = 0;
            double cumulativeDecoy = 0;
            for (int r = 0; r < order.Length; r++)
            {
                if (matches[order[r]].IsDecoy)
                {
                    cumulativeDecoy++;
                }
                else
                {
                    cumulativeTarget++;
                }

                q[r] = cumulativeDecoy / Math.Max(cumulativeTarget, 1);
            }

            double best = double.PositiveInfinity;
            for (int r = order.Length - 1; r >= 0; r--)
            {
                best = Math.Min(best, q[r]);
                q[r] = best;
            }

            int lastRank = -1;
            for (int r = 0; r < order.Length; r++)
            {
                if (!matches[order[r]].IsDecoy && q[r] <= qValueCutoff)
                {
                    lastRank = r;
                }
            }

            var positives = new HashSet<SpectralMatch>(ReferenceEqualityComparer.Instance);
            for (int r = 0; r <= lastRank; r++)
            {
                if (!matches[order[r]].IsDecoy)
                {
                    positives.Add(matches[order[r]]);
                }
            }

            return positives;
        }

        /// <summary>
        /// Target peptides whose peptide-level PEP q-value is below <paramref name="qValueCutoff"/>, computed
        /// the way <see cref="FdrAnalysisEngine"/> computes the reported one: the PEP-best match per full
        /// sequence (ties broken by the default comparer), fractional target and decoy counts per hypothesis
        /// as in CalculateQValue, and q = min over the tail of (decoys + 1) / targets as in PepQValueInverted.
        /// Reads the matches' current PEP. For reporting and the stopping decision only -- never labels.
        /// <remarks>
        /// Not written through CalculateQValue itself because that method stores its counts on the matches,
        /// and FdrAnalysisEngine leaves some of those fields untouched after PEP, so writing them here
        /// would change output.
        /// </remarks>
        /// </summary>
        internal static int CountAcceptedPeptides(IEnumerable<SpectralMatch> psms, double qValueCutoff)
        {
            var peptides = psms
                .Where(p => p != null)
                .OrderBy(p => p.FdrInfo.PEP)
                .ThenByDescending(p => p)
                .GroupBy(p => p.FullSequence)
                .Select(g => g.First())
                .ToList();

            var q = new double[peptides.Count];
            double cumulativeTarget = 0;
            double cumulativeDecoy = 0;
            for (int i = 0; i < peptides.Count; i++)
            {
                double totalHits = peptides[i].BestMatchingBioPolymersWithSetMods.Count();
                double targetHits = peptides[i].BestMatchingBioPolymersWithSetMods.Count(h => !h.IsDecoy);
                cumulativeTarget += targetHits / totalHits;
                cumulativeDecoy += (totalHits - targetHits) / totalHits;
                q[i] = (cumulativeDecoy + 1) / cumulativeTarget;
            }

            int accepted = 0;
            double best = double.PositiveInfinity;
            for (int i = peptides.Count - 1; i >= 0; i--)
            {
                best = Math.Min(best, q[i]);
                if (!peptides[i].IsDecoy && best < qValueCutoff)
                {
                    accepted++;
                }
            }

            return accepted;
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
        /// Training data for the given groups, labelled by the search-score q-value (round 0's rule).
        /// </summary>
        public IEnumerable<PsmData> CreatePsmData(string searchType,
            List<SpectralMatchGroup> peptideGroups, List<int> peptideGroupIndices)
        {
            return CreatePsmData(searchType, peptideGroups, peptideGroupIndices, IsPositiveBySearchScore);
        }

        /// <summary>
        /// Training data for the given groups. Decoys are negatives; a target is a positive when
        /// <paramref name="isPositive"/> says so and is left out otherwise. The rule is a parameter, not
        /// engine state, so the labels a caller gets never depend on what ran before.
        /// </summary>
        private List<PsmData> CreatePsmData(string searchType,
            List<SpectralMatchGroup> peptideGroups, List<int> peptideGroupIndices, Func<SpectralMatch, bool> isPositive)
        {
            List<PsmData> psmDataList = new List<PsmData>();

            for (int i = 0; i < peptideGroupIndices.Count; i++)
            {
                int modCount = 0;
                foreach (var psm in peptideGroups[peptideGroupIndices[i]].GetBestMatches().Where(psm => psm != null))
                {
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
                        else if (!csm.IsDecoy && !csm.BetaPeptide.IsDecoy && isPositive(csm))
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
                                && isPositive(psm))
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

        public static string AggregateMetricsForOutput(List<CalibratedBinaryClassificationMetrics> allMetrics,
            int positiveTrainingCount, int negativeTrainingCount, double qValueCutoff, int? ambiguousHypothesesRemoved = null,
            int trainingRounds = 1, int acceptedTargets = 0, string roundLog = null)

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
            if (ambiguousHypothesesRemoved.HasValue)
                s.AppendLine($"*       Count of Ambiguous {char.ToUpper(GlobalVariables.AnalyteType.GetUniqueFormLabel()[0]) + GlobalVariables.AnalyteType.GetUniqueFormLabel()[1..]}s Removed:  " + ambiguousHypothesesRemoved.Value);
            s.AppendLine("*       Q-Value Cutoff for Training Targets:  " + qValueCutoff);
            s.AppendLine("*       Targets Used for Training:  " + positiveTrainingCount);
            s.AppendLine("*       Decoys Used for Training:  " + negativeTrainingCount);
            s.AppendLine("*       Training Rounds Run:  " + trainingRounds);
            if (!string.IsNullOrEmpty(roundLog))
            {
                s.Append(roundLog);
            }
            if (acceptedTargets > 0)
            {
                s.AppendLine($"*       Target {GlobalVariables.AnalyteType.GetUniqueFormLabel()}s Accepted After Training:  " + acceptedTargets);
            }
            s.AppendLine("************************************************************");
            return s.ToString();
        }

        /// <summary>
        /// Assigns a PEP to every spectral match in the given groups. Unless <see cref="PruneAmbiguousHypotheses"/>
        /// is set, this method scores; it does not prune. Removing ambiguous match hypotheses is the
        /// DisambiguationEngine's job.
        /// </summary>
        public void Compute_PSM_PEP(List<SpectralMatchGroup> peptideGroups,
            List<int> peptideGroupIndices,
            MLContext mLContext, TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>> trainedModel, string searchType, string outputFolder)
        {
            int maxThreads = FileSpecificParametersDictionary.Values.FirstOrDefault().MaxThreadsToUsePerFile;

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

                        int ambiguousRemovedInThread = 0;
                        List<double> pepValuePredictions = new List<double>();
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            foreach (SpectralMatch psm in peptideGroups[peptideGroupIndices[i]])
                            {
                                if (psm != null)
                                {
                                    pepValuePredictions.Clear();

                                    // One prediction per ambiguous match hypothesis. The PSM keeps the best of them;
                                    // the others are left in place for the DisambiguationEngine to judge.
                                    var hypotheses = psm.BestMatchingBioPolymersWithSetMods.ToList();
                                    foreach (SpectralMatchHypothesis bestMatch in hypotheses)
                                    {
                                        PsmData pd = CreateOnePsmDataEntry(searchType, psm, bestMatch, !bestMatch.IsDecoy);
                                        var pepValuePrediction = threadPredictionEngine.Predict(pd);
                                        pepValuePredictions.Add(pepValuePrediction.Probability);
                                        //A score is available using the variable pepvaluePrediction.Score
                                    }

                                    ambiguousRemovedInThread += AssignPep(psm, hypotheses, pepValuePredictions, PruneAmbiguousHypotheses);
                                }

                            }
                        }
                        Interlocked.Add(ref _ambiguousHypothesesRemoved, ambiguousRemovedInThread);
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
        }

        /// <summary>
        /// Sets the PSM's PEP from the best of its hypotheses' predictions. When
        /// <paramref name="pruneAmbiguousHypotheses"/> is true, also removes every hypothesis predicted more than
        /// <see cref="AbsoluteProbabilityThatDistinguishesPeptides"/> below the best, and returns how many were
        /// removed. The PEP is the same either way, because the best prediction is never removed.
        /// </summary>
        public static int AssignPep(SpectralMatch psm, List<SpectralMatchHypothesis> hypotheses, List<double> pepValuePredictions, bool pruneAmbiguousHypotheses)
        {
            int removed = 0;
            if (pruneAmbiguousHypotheses)
            {
                List<int> indicesOfPeptidesToRemove = new List<int>();
                GetIndicesOfPeptidesToRemove(indicesOfPeptidesToRemove, pepValuePredictions);
                RemoveBestMatchingPeptidesWithLowPEP(psm, indicesOfPeptidesToRemove, hypotheses, ref removed);
            }

            psm.PsmFdrInfo.PEP = 1 - pepValuePredictions.Max();
            psm.PeptideFdrInfo.PEP = 1 - pepValuePredictions.Max();
            return removed;
        }

        public PsmData CreateOnePsmDataEntry(string searchType, SpectralMatch psm, SpectralMatchHypothesis tentativeSpectralMatch, bool label)
        {
            // A COPY, never the cached instance. Training and prediction both come through here,
            // ML.NET enumerates its training set lazily, and handing out one shared object would
            // let prediction mutate a Label inside a training list a later fold is about to fit on.
            psm.PsmData_forPEPandPercolator = GetFeatures(searchType, psm, tentativeSpectralMatch).WithLabel(label);
            return psm.PsmData_forPEPandPercolator;
        }

        /// <summary>
        /// The feature vector of one hypothesis, from <see cref="_featureCache"/>. Only Label changes between
        /// training rounds. Everything else is a pure function of inputs that are fixed for the lifetime of
        /// this engine, and computing it is dominated by retention-time prediction. The returned instance is
        /// the cached one: read it, never modify it. Its Label is meaningless.
        /// </summary>
        private PsmData GetFeatures(string searchType, SpectralMatch psm, SpectralMatchHypothesis tentativeSpectralMatch)
        {
            return _featureCache.GetOrAdd(tentativeSpectralMatch, _ => ComputeFeatures(searchType, psm, tentativeSpectralMatch));
        }

        private PsmData ComputeFeatures(string searchType, SpectralMatch psm, SpectralMatchHypothesis tentativeSpectralMatch)
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

            return new PsmData
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

                Label = false,

                SpectralAngle = spectralAngle,
                HasSpectralAngle = hasSpectralAngle,
                PeaksInPrecursorEnvelope = peaksInPrecursorEnvelope,
                ChimeraCount = chimeraCount,
                MostAbundantPrecursorPeakIntensity = mostAbundantPrecursorPeakIntensity,
                PrecursorFractionalIntensity = fractionalIntensity,
                InternalIonCount = internalMatchingFragmentCount,
                PrecursorDeconvolutionScore = (float)psm.PrecursorScanDeconvolutionScore,
            };
        }

        /// <summary>
        /// Removes the given ambiguous match hypotheses from the PSM.
        /// <remarks>
        /// Called only when <see cref="PruneAmbiguousHypotheses"/> is set (glyco, crosslink and nonspecific searches,
        /// which have no DisambiguationEngine downstream). Otherwise PEP scores and does not prune:
        /// disambiguation-by-PEP belongs in <see cref="SpectrumMatch.DisambiguationEngine"/>, whose
        /// own summary already names "PEPAnalysisEngine -> By PEP" as a site to consolidate there.
        /// </remarks>
        /// </summary>
        public static void RemoveBestMatchingPeptidesWithLowPEP(SpectralMatch psm, List<int> indicesOfPeptidesToRemove, List<SpectralMatchHypothesis> allPeptides, ref int ambiguousPeptidesRemovedCount)
        {
            int peptidesRemoved = 0;
            foreach (var toRemove in indicesOfPeptidesToRemove)
            {
                psm.RemoveThisAmbiguousPeptide(allPeptides[toRemove - peptidesRemoved]);
                peptidesRemoved++;
            }
            ambiguousPeptidesRemovedCount += peptidesRemoved;
        }

        /// <summary>
        /// Given a set of PEP values, this method will find the indices of BestMatchingBioPolymersWithSetMods that are not within the required tolerance
        /// This method will also remove the low scoring predictions from the set.
        /// <remarks>
        /// Called only when pruning -- see <see cref="RemoveBestMatchingPeptidesWithLowPEP"/>.
        /// Note that it never drops the maximum (max - max = 0 is not &gt; the threshold), which is why
        /// pruning or not leaves every assigned PEP unchanged: PEP is 1 - pepValuePredictions.Max().
        /// </remarks>
        /// </summary>
        public static void GetIndicesOfPeptidesToRemove(List<int> indicesOfPeptidesToRemove, List<double> pepValuePredictions)
        {
            double highestPredictedPEPValue = pepValuePredictions.Max();
            for (int i = 0; i < pepValuePredictions.Count; i++)
            {
                if ((highestPredictedPEPValue - pepValuePredictions[i]) > AbsoluteProbabilityThatDistinguishesPeptides)
                {
                    indicesOfPeptidesToRemove.Add(i);
                }
            }

            foreach (int i in indicesOfPeptidesToRemove.OrderByDescending(p => p))
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
