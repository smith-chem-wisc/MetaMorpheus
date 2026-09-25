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
        /// Upper bound on training rounds. This is a safety cap, NOT the number of rounds that will
        /// run -- <see cref="TrainingImprovementTolerance"/> decides that. A fixed count tuned on one
        /// dataset is wrong on the next, so the stopping rule is a criterion and this is only a stop.
        /// </summary>
        public int MaxTrainingRounds { get; set; } = 10;

        /// <summary>
        /// Stop when a round grows the accepted target count by less than this fraction. 0.001 = 0.1%.
        /// </summary>
        public double TrainingImprovementTolerance { get; set; } = 0.001;

        /// <summary>
        /// The q-value below which a target counts as "accepted" when judging whether a round helped.
        /// Matches the threshold results are reported at, so the criterion optimises the reported number.
        /// </summary>
        public double AcceptanceQValueCutoff { get; set; } = 0.01;

        /// <summary>
        /// The peptide groups whose best match qualifies as a POSITIVE training example this round.
        /// Null for the first round, which falls back to the search-score q-value exactly as before.
        /// Rounds after the first re-derive this from the model's own output -- the semi-supervised
        /// step this engine previously lacked.
        /// </summary>
        private double _positiveTrainingPepThreshold = double.NaN;

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
            IEnumerable<PsmData>[] PSMDataGroups = new IEnumerable<PsmData>[numGroups];
            int maxThreads = FileSpecificParametersDictionary.Values.FirstOrDefault().MaxThreadsToUsePerFile;

            MLContext mlContext = new MLContext(seed: _randomSeed);
            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(BGDTreeOptions);
            var pipeline = mlContext.Transforms.Concatenate("Features", TrainingVariables)
                .Append(trainer);

            List<CalibratedBinaryClassificationMetrics> allMetrics = new List<CalibratedBinaryClassificationMetrics>();
            int positiveTrainingCount = 0;
            int negativeTrainingcount = 0;
            int roundsRun = 0;
            // Per-round trace. Without it a long run is completely opaque: the metrics block is only
            // emitted at the end, which is the same complaint AggregateMetricsForOutput already earns.
            var roundLog = new StringBuilder();
            var roundClock = System.Diagnostics.Stopwatch.StartNew();
            int previousAccepted = 0;
            // Snapshot so a round that makes things worse can be undone. One double per PSM.
            double[] bestPepSnapshot = null;

            // Pruning deletes hypotheses irreversibly, so a later, better round could not bring them back.
            // Callers that prune (glyco, crosslink) therefore train once, exactly as before iteration existed.
            int maxRounds = PruneAmbiguousHypotheses ? 1 : Math.Max(1, MaxTrainingRounds);
            for (int round = 0; round < maxRounds; round++)
            {
                var roundMetrics = new List<CalibratedBinaryClassificationMetrics>();
                bool foldStarved = false;
                int roundPositives = 0;
                int roundNegatives = 0;

                // Round 0's labels come from the search score and do not depend on the fold, so its data is
                // built once, up front, before any fold's model has scored anything -- as before iteration
                // existed. Building it per fold instead would let a pruning caller's earlier folds delete
                // hypotheses from the data a later fold trains on.
                IEnumerable<PsmData>[] roundZeroData = null;
                if (round == 0)
                {
                    _positiveTrainingPepThreshold = double.NaN;
                    roundZeroData = new IEnumerable<PsmData>[numGroups];
                    for (int group = 0; group < numGroups; group++)
                    {
                        roundZeroData[group] = CreatePsmData(SearchType, peptideGroups, peptideGroupIndices[group]);
                    }
                }

                for (int groupIndexNumber = 0; groupIndexNumber < numGroups; groupIndexNumber++)
                {
                    List<int> allGroupIndexes = Enumerable.Range(0, numGroups).ToList();
                    allGroupIndexes.RemoveAt(groupIndexNumber);

                    // FOLD-LOCAL LABELS. The threshold that decides which matches are positive
                    // training examples is derived ONLY from the folds this model trains on. Deriving
                    // it globally would let a held-out peptide's own score influence the labels of the
                    // peptides that train the model which then scores it -- a leak that compounds every
                    // round. mokapot avoids it by running its whole iteration inside each fold
                    // (brew() deep-copies a model per fold); this is the same guarantee.
                    // Round 0 has nothing to leak: its labels come from the search score, not from us.
                    _positiveTrainingPepThreshold = round == 0
                        ? double.NaN
                        : ComputeTrainingPepThreshold(peptideGroups,
                            allGroupIndexes.SelectMany(i => peptideGroupIndices[i]));

                    // Built per model, not per group, because the label of a given peptide now depends
                    // on WHICH model is being trained. Cheap only because the feature vectors are
                    // cached -- this is what _featureCache was for.
                    var trainingParts = new List<IEnumerable<PsmData>>();
                    foreach (int trainingGroup in allGroupIndexes)
                    {
                        trainingParts.Add(round == 0
                            ? roundZeroData[trainingGroup]
                            : CreatePsmData(SearchType, peptideGroups, peptideGroupIndices[trainingGroup]));
                    }

                    if (!trainingParts.Any(part => part.Any(p => p.Label)) || !trainingParts.Any(part => part.Any(p => !p.Label)))
                    {
                        foldStarved = true;
                        break;
                    }

                    //concat doesn't work in a loop, therefore I had to hard code the concat to group 3 out of 4 lists. if the const int numGroups value is changed, then the concat has to be changed accordingly.
                    IDataView dataView = mlContext.Data.LoadFromEnumerable(trainingParts[0].Concat(trainingParts[1].Concat(trainingParts[2])));
                    trainedModels[groupIndexNumber] = pipeline.Fit(dataView);

                    // The held-out fold, labelled by the SAME threshold, so the evaluation measures the
                    // model against the rule it was trained under rather than a different one.
                    var heldOut = round == 0
                        ? roundZeroData[groupIndexNumber]
                        : CreatePsmData(SearchType, peptideGroups, peptideGroupIndices[groupIndexNumber]);
                    var myPredictions = trainedModels[groupIndexNumber].Transform(mlContext.Data.LoadFromEnumerable(heldOut));
                    CalibratedBinaryClassificationMetrics metrics = mlContext.BinaryClassification.Evaluate(data: myPredictions, labelColumnName: "Label", scoreColumnName: "Score");

                    //model is trained on peptides but here we can use that to compute PEP for all PSMs
                    Compute_PSM_PEP(peptideGroups, peptideGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], SearchType, OutputFolder);

                    roundMetrics.Add(metrics);
                    roundPositives += trainingParts.Sum(part => part.Count(p => p.Label));
                    roundNegatives += trainingParts.Sum(part => part.Count(p => !p.Label));
                }

                if (foldStarved)
                {
                    if (round == 0)
                    {
                        return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
                    }

                    // A later round starved a fold. Keep what the previous round produced.
                    RestorePepValues(bestPepSnapshot);
                    break;
                }

                // Global count, for REPORTING and for the stopping decision only -- never for labels.
                // The one-bit-per-round channel this opens (whether to run again) is the same one
                // mokapot's max_iter and BestFeatureIsBetterError use.
                int accepted = CountAcceptedTargets(peptideGroups);

                if (round > 0 && accepted <= previousAccepted)
                {
                    // No better than the round before. Keep the better one and stop -- an unguarded loop
                    // can make a search worse and nothing downstream would notice.
                    RestorePepValues(bestPepSnapshot);
                    break;
                }

                bestPepSnapshot = SnapshotPepValues();
                allMetrics = roundMetrics;
                // Averaged over the folds: each peptide appears in the training set of 3 of the 4 models.
                positiveTrainingCount = roundPositives / (numGroups - 1);
                negativeTrainingcount = roundNegatives / (numGroups - 1);
                roundsRun = round + 1;
                // NO TIMINGS in the results block. It is compared for byte equality across two runs of
                // the same data by PepAnalysisEngineHasReproducibleOutput, and a wall clock is not
                // reproducible. Timings go to the streamed file below, which nothing asserts on.
                roundLog.AppendLine($"*         round {round}:  accepted {accepted}");
                string roundLine = $"round {round}:  accepted {accepted}  ({roundClock.Elapsed.TotalSeconds:F1} s)";
                // Also stream it to disk as it happens. The metrics block is only emitted when the
                // whole engine finishes, so without this a multi-round run on a large dataset shows
                // nothing at all for as long as it takes -- which makes it impossible to tell a slow
                // run from a hung one, or to see the trajectory before it is over.
                WriteRoundProgress(roundLine);
                roundClock.Restart();

                double improvement = previousAccepted == 0
                    ? double.PositiveInfinity
                    : (accepted - previousAccepted) / (double)previousAccepted;
                previousAccepted = accepted;

                if (round > 0 && improvement < TrainingImprovementTolerance)
                {
                    break; // converged by criterion, not by a hard-coded count
                }
            }

            return AggregateMetricsForOutput(allMetrics, positiveTrainingCount, negativeTrainingcount, QValueCutoff,
                roundsRun, previousAccepted, roundLog.ToString(), PruneAmbiguousHypotheses ? _ambiguousHypothesesRemoved : null);
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
        /// Does this match count as a positive training example for the model currently being trained?
        /// <remarks>
        /// Round 0 uses the SEARCH-SCORE q-value, which is all this engine ever used. Later rounds use
        /// a PEP threshold derived by <see cref="ComputeTrainingPepThreshold"/> from the folds THIS
        /// model trains on -- never from the fold it will score.
        ///
        /// That is the semi-supervised step: without it the model only ever learns from matches the
        /// raw score already trusted, so the population that rescoring exists to rescue is absent from
        /// every training set it sees.
        /// </remarks>
        /// </summary>
        private bool IsPositiveTrainingExample(SpectralMatch psm)
        {
            return double.IsNaN(_positiveTrainingPepThreshold)
                ? psm.GetFdrInfo(UsePeptideLevelQValueForTraining).QValue <= QValueCutoff
                : psm.GetFdrInfo(UsePeptideLevelQValueForTraining).PEP <= _positiveTrainingPepThreshold;
        }

        /// <summary>
        /// Walks a population in PEP order accumulating a target-decoy q-value, and returns the
        /// monotone q for each position alongside the matches.
        /// <remarks>
        /// q must be monotone in score, so the running ratio is swept once forward and then minimised
        /// from the bottom up -- the same shape as FdrAnalysisEngine.QValueInverted.
        ///
        /// The population is GetBestMatches(), one per full sequence: the unit the model is trained on.
        /// SpectralMatchGroup.BestMatch (one per BASE sequence) would silently drop every modified
        /// variant.
        /// </remarks>
        /// </summary>
        private (List<SpectralMatch> Ordered, double[] MonotoneQ) SweepByPep(IEnumerable<SpectralMatch> matches)
        {
            var ordered = matches
                .Where(m => m != null)
                .OrderBy(m => m.GetFdrInfo(UsePeptideLevelQValueForTraining).PEP)
                .ThenByDescending(m => m)
                .ToList();

            var runningQ = new double[ordered.Count];
            double cumulativeTarget = 0;
            double cumulativeDecoy = 0;
            for (int i = 0; i < ordered.Count; i++)
            {
                if (ordered[i].IsDecoy)
                {
                    cumulativeDecoy++;
                }
                else
                {
                    cumulativeTarget++;
                }

                runningQ[i] = cumulativeDecoy / Math.Max(cumulativeTarget, 1);
            }

            double best = double.PositiveInfinity;
            for (int i = ordered.Count - 1; i >= 0; i--)
            {
                best = Math.Min(best, runningQ[i]);
                runningQ[i] = best;
            }

            return (ordered, runningQ);
        }

        /// <summary>
        /// The PEP at which the target-decoy q-value of the given TRAINING population crosses
        /// <see cref="QValueCutoff"/>. Matches at or below it are the next round's positives.
        /// </summary>
        private double ComputeTrainingPepThreshold(List<SpectralMatchGroup> peptideGroups, IEnumerable<int> populationIndices)
        {
            var population = populationIndices.SelectMany(i => peptideGroups[i].GetBestMatches());
            var (ordered, monotoneQ) = SweepByPep(population);

            double threshold = double.NegativeInfinity;
            for (int i = 0; i < ordered.Count; i++)
            {
                if (!ordered[i].IsDecoy && monotoneQ[i] <= QValueCutoff)
                {
                    threshold = ordered[i].GetFdrInfo(UsePeptideLevelQValueForTraining).PEP;
                }
            }

            return threshold;
        }

        /// <summary>
        /// Target groups accepted at <see cref="AcceptanceQValueCutoff"/> across ALL folds.
        /// Reporting and the stopping decision only -- never labels.
        /// </summary>
        private int CountAcceptedTargets(List<SpectralMatchGroup> peptideGroups)
        {
            var (ordered, monotoneQ) = SweepByPep(peptideGroups.SelectMany(g => g.GetBestMatches()));

            int accepted = 0;
            for (int i = 0; i < ordered.Count; i++)
            {
                if (!ordered[i].IsDecoy && monotoneQ[i] < AcceptanceQValueCutoff)
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


        public IEnumerable<PsmData> CreatePsmData(string searchType,
            List<SpectralMatchGroup> peptideGroups, List<int> peptideGroupIndices)
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
                        else if (!csm.IsDecoy && !csm.BetaPeptide.IsDecoy && IsPositiveTrainingExample(csm))
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
                                && IsPositiveTrainingExample(psm))
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
            int positiveTrainingCount, int negativeTrainingCount, double qValueCutoff,
            int trainingRounds = 1, int acceptedTargets = 0, string roundLog = null, int? ambiguousHypothesesRemoved = null)

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
            // Only Label changes between training rounds -- see _featureCache. Everything below is a
            // pure function of inputs that are fixed for the lifetime of this engine, and computing
            // it is dominated by retention-time prediction.
            if (_featureCache.TryGetValue(tentativeSpectralMatch, out PsmData cached))
            {
                // A COPY, never the cached instance. Training and prediction both come through here,
                // ML.NET enumerates its training set lazily, and handing out one shared object would
                // let prediction mutate a Label inside a training list a later fold is about to fit on.
                psm.PsmData_forPEPandPercolator = cached.WithLabel(label);
                return psm.PsmData_forPEPandPercolator;
            }

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
                PrecursorDeconvolutionScore = (float)psm.PrecursorScanDeconvolutionScore,
            };

            _featureCache[tentativeSpectralMatch] = psm.PsmData_forPEPandPercolator;

            return psm.PsmData_forPEPandPercolator;
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
