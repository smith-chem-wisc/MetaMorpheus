using Omics;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.Localization
{
    /// <summary>
    /// Where on its peptide (and protein) a site is, for the prior over where a modification occurs. The first and last
    /// residues stand for their peptide termini as well, which fragment identically.
    /// </summary>
    public enum LocalizationPositionClass
    {
        Anywhere,
        PeptideNTerminal,
        PeptideCTerminal,
        ProteinNTerminal,
        ProteinCTerminal
    }

    /// <summary>
    /// Computes localization probabilities from the <see cref="LocalizationEvidence"/> the <see cref="LocalizationEngine"/>
    /// collected, once every spectra file has been read (An et al., PTMiner, MCP 2019):
    ///
    /// P(placement | spectrum) is proportional to P(spectrum | placement) x P(placement).
    ///
    /// The likelihood takes, over the peaks some placements match and others do not, F(I) for each peak a placement
    /// explains and 1 - F(I) for each it does not, where F is <see cref="MatchedPeakIntensityModel"/> learned per spectra
    /// file from confident unmodified matches. Peaks every placement matches cancel. That is PTMiner's likelihood
    /// (<see cref="LocalizationLikelihood.PtminerCdf"/>); it counts explaining a peak fainter than the median matched
    /// peak against a placement. <see cref="LocalizationLikelihood.IntensityRatio"/> instead scores each predicted
    /// fragment as a likelihood ratio against chance (<see cref="PeakIntensityRatioModel"/>); <see cref="Likelihood"/>
    /// chooses.
    ///
    /// The prior is the product, over a placement's sites, of how often that modification occurs on that residue in
    /// that position (anywhere, peptide or protein terminal), learned by expectation maximization from the confident
    /// target matches, separately for each chemistry (and each cluster of open-search mass differences). A chemistry
    /// seen on too few matches keeps a flat prior.
    ///
    /// Each match gets <see cref="SpectralMatch.LocalizationProbabilities"/>; confident targets with more than one
    /// placement also get a localization q-value.
    /// </summary>
    public class LocalizationProbabilityEngine : MetaMorpheusEngine
    {
        public const int DefaultMinTrainingPeaksPerFile = 1000;
        public const int DefaultMinTrainingPeaks = 50;
        public const int DefaultMinMatchesForPrior = 5;
        public const double DefaultPseudocount = 1;
        public const double DefaultConvergence = 0.01;
        public const int DefaultMaxIterations = 100;
        public const double DefaultMassDifferenceClusterTolerance = 0.01;

        private readonly List<SpectralMatch> AllMatches;
        private readonly Func<SpectralMatch, bool> IsConfident;

        /// <param name="allMatches">The matches of every spectra file.</param>
        /// <param name="isConfident">Which matches are confident enough to learn from and to rank for the localization
        /// q-value; by default, those at PSM q-value 0.01 or less. Decoys are never learned from.</param>
        public LocalizationProbabilityEngine(IEnumerable<SpectralMatch> allMatches, CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds,
            Func<SpectralMatch, bool> isConfident = null)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllMatches = allMatches.ToList();
            IsConfident = isConfident ?? (psm => psm.PsmFdrInfo != null && psm.PsmFdrInfo.QValue <= 0.01);
        }

        public int MinTrainingPeaksPerFile { get; init; } = DefaultMinTrainingPeaksPerFile;
        public int MinTrainingPeaks { get; init; } = DefaultMinTrainingPeaks;
        public int MinMatchesForPrior { get; init; } = DefaultMinMatchesForPrior;
        public double Pseudocount { get; init; } = DefaultPseudocount;
        public double Convergence { get; init; } = DefaultConvergence;
        public int MaxIterations { get; init; } = DefaultMaxIterations;
        public double MassDifferenceClusterTolerance { get; init; } = DefaultMassDifferenceClusterTolerance;
        public double MinRelativeIntensity { get; init; } = MatchedPeakIntensityModel.DefaultMinRelativeIntensity;
        public MatchedPeakIntensityDistribution IntensityDistribution { get; init; } = MatchedPeakIntensityDistribution.Empirical;

        /// <summary>
        /// The likelihood a candidate's spectrum is scored with. <see cref="LocalizationLikelihood.IntensityRatio"/> falls
        /// back to <see cref="LocalizationLikelihood.PtminerCdf"/>, with a warning, when the evidence carries no
        /// unexplained peaks or predicted-fragment counts to learn it from.
        /// </summary>
        public LocalizationLikelihood Likelihood { get; init; } = DefaultLikelihood;

        public const LocalizationLikelihood DefaultLikelihood = LocalizationLikelihood.IntensityRatio;

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            var matches = AllMatches.Where(m => m.LocalizationEvidence != null).ToList();
            var results = new LocalizationProbabilityResults(this);
            if (matches.Count == 0)
            {
                return results;
            }

            Status("Computing localization probabilities...");

            // F (and the intensity-ratio model), per spectra file where the file has enough confident unmodified peaks,
            // else from every file
            var models = FitIntensityModels(matches, out FileModels pooled);
            results.TrainingPeaks = pooled.Cdf.TrainingPeakCount;
            if (pooled.Cdf.TrainingPeakCount < MinTrainingPeaks)
            {
                Warn($"Localization probabilities were not computed: only {pooled.Cdf.TrainingPeakCount} matched peaks of confident unmodified matches to learn peak intensities from, and at least {MinTrainingPeaks} are needed.");
                return results;
            }

            LocalizationLikelihood likelihood = Likelihood;
            if (likelihood == LocalizationLikelihood.IntensityRatio && !pooled.Ratio.IsTrained)
            {
                Warn("The intensity-ratio likelihood could not be learned (the confident unmodified matches recorded no unexplained peaks or predicted fragments); localization probabilities use PTMiner's likelihood instead.");
                likelihood = LocalizationLikelihood.PtminerCdf;
            }
            results.Likelihood = likelihood;
            if (likelihood == LocalizationLikelihood.IntensityRatio)
            {
                results.FragmentObservationProbability = pooled.Ratio.ObservationProbability;
                results.MeanChanceCoverage = pooled.Ratio.MeanCoverage;
            }

            // each candidate's log likelihood, and the prior class of each of its sites
            var siteClasses = new SiteClassRegistry(matches, MassDifferenceClusterTolerance);
            var scored = new ScoredMatch[matches.Count];
            Parallel.For(0, matches.Count, new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile }, i =>
            {
                SpectralMatch psm = matches[i];
                FileModels model = models.TryGetValue(psm.FullFilePath ?? "", out var fileModel) ? fileModel : pooled;
                double[] logLikelihoods = likelihood == LocalizationLikelihood.IntensityRatio
                    ? LogLikelihoodRatios(psm.LocalizationEvidence, model.Ratio)
                    : LogLikelihoods(psm.LocalizationEvidence, model.Cdf);
                scored[i] = new ScoredMatch(psm, logLikelihoods, IsConfident(psm) && !psm.IsDecoy);
            });
            foreach (ScoredMatch match in scored)
            {
                match.Classes = siteClasses.ClassesOf(match.Psm.LocalizationEvidence);
            }

            double[] logPrior = LearnPrior(scored, siteClasses, results);

            foreach (ScoredMatch match in scored)
            {
                match.Psm.LocalizationProbabilities = new LocalizationProbabilities(match.Psm.LocalizationEvidence.Placements, Posterior(match, logPrior));
            }

            AssignQValues(scored, results);
            return results;
        }

        /// <summary>The peak-intensity models of one spectra file, or of every file pooled.</summary>
        private sealed record FileModels(MatchedPeakIntensityModel Cdf, PeakIntensityRatioModel Ratio);

        /// <summary>What confident unmodified matches of one spectra file offer to learn from.</summary>
        private sealed class TrainingData
        {
            public List<double> MatchedPeaks { get; } = new();
            public List<IReadOnlyList<float>> NoiseCoverages { get; } = new();
            public long PredictedFragments { get; set; }
            public long ObservedFragments { get; set; }
        }

        /// <summary>
        /// Fits F, and the intensity-ratio model, to confident target matches that carry no modification but fixed ones
        /// and whose peptide explains the precursor: F and f_m to every peak all their candidates match, f_n to the
        /// peaks they left unexplained, and d to the fraction of their predicted fragments matched at or above the floor.
        /// </summary>
        private Dictionary<string, FileModels> FitIntensityModels(List<SpectralMatch> matches, out FileModels pooled)
        {
            var dataByFile = new Dictionary<string, TrainingData>();
            foreach (SpectralMatch psm in matches.Where(m => !m.IsDecoy && IsConfident(m) && IsUnmodified(m)))
            {
                LocalizationEvidence evidence = psm.LocalizationEvidence;
                string file = psm.FullFilePath ?? "";
                if (!dataByFile.TryGetValue(file, out var data))
                {
                    dataByFile[file] = data = new TrainingData();
                }
                int observed = 0;
                for (int p = 0; p < evidence.PeakCount; p++)
                {
                    if (evidence.IsSharedByAllCandidates(p))
                    {
                        data.MatchedPeaks.Add(evidence.PeakRelativeIntensity(p));
                        if (evidence.PeakRelativeIntensity(p) >= MinRelativeIntensity)
                        {
                            observed++;
                        }
                    }
                }
                if (evidence.NoiseCoverage != null)
                {
                    data.NoiseCoverages.Add(evidence.NoiseCoverage);
                }
                if (evidence.HasPredictedFragmentCounts && evidence.CandidateCount == 1)
                {
                    data.PredictedFragments += evidence.PredictedFragmentCount(0);
                    data.ObservedFragments += Math.Min(observed, evidence.PredictedFragmentCount(0));
                }
            }

            FileModels Fit(IEnumerable<TrainingData> data) => new(
                new MatchedPeakIntensityModel(data.SelectMany(d => d.MatchedPeaks), MinRelativeIntensity, IntensityDistribution),
                new PeakIntensityRatioModel(data.SelectMany(d => d.MatchedPeaks), data.SelectMany(d => d.NoiseCoverages),
                    data.Sum(d => d.PredictedFragments), data.Sum(d => d.ObservedFragments), MinRelativeIntensity));

            pooled = Fit(dataByFile.Values);
            var models = new Dictionary<string, FileModels>();
            foreach (var (file, data) in dataByFile)
            {
                FileModels model = Fit(new[] { data });
                if (model.Cdf.TrainingPeakCount >= MinTrainingPeaksPerFile)
                {
                    // a file whose ratio model cannot be learned on its own borrows the pooled one
                    models[file] = model.Ratio.IsTrained ? model : model with { Ratio = pooled.Ratio };
                }
            }
            return models;
        }

        private bool IsUnmodified(SpectralMatch psm)
        {
            IBioPolymerWithSetMods peptide = psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer;
            return peptide.NumVariableMods == 0
                && CommonParameters.PrecursorMassTolerance.Within(peptide.MonoisotopicMass + psm.LocalizationEvidence.MassDifference, peptide.MonoisotopicMass);
        }

        /// <summary>
        /// Each candidate's log likelihood: over the peaks some candidates match and others do not, log F(I) where the
        /// candidate matches and log(1 - F(I)) where it does not. Peaks every candidate matches, and peaks fainter than
        /// the model's minimum, are left out.
        /// </summary>
        internal static double[] LogLikelihoods(LocalizationEvidence evidence, MatchedPeakIntensityModel model)
        {
            var logLikelihoods = new double[evidence.CandidateCount];
            for (int p = 0; p < evidence.PeakCount; p++)
            {
                double intensity = evidence.PeakRelativeIntensity(p);
                if (intensity < model.MinRelativeIntensity || evidence.IsSharedByAllCandidates(p))
                {
                    continue;
                }

                double f = model.CumulativeProbability(intensity);
                double logMatched = Math.Log(f);
                double logUnmatched = Math.Log(1 - f);
                for (int c = 0; c < evidence.CandidateCount; c++)
                {
                    logLikelihoods[c] += evidence.Matches(c, p) ? logMatched : logUnmatched;
                }
            }
            return logLikelihoods;
        }

        /// <summary>
        /// Each candidate's log likelihood ratio under <see cref="LocalizationLikelihood.IntensityRatio"/>: over the peaks
        /// at or above the floor that some candidates match and others do not, log((1 - d) + d f_m(I) / (q f_n(I))) for
        /// each candidate that explains the peak; and log(1 - d) for each fragment a candidate predicts within the scan's
        /// range but does not match at or above the floor. Peaks every candidate matches count the same for each and are
        /// left out of the first sum. q is the match's own chance-peak coverage (see <see cref="PeakIntensityRatioModel"/>).
        /// </summary>
        internal static double[] LogLikelihoodRatios(LocalizationEvidence evidence, PeakIntensityRatioModel model)
        {
            int candidates = evidence.CandidateCount;
            var logLikelihoods = new double[candidates];
            var observed = new int[candidates];
            double coverage = model.CoverageOf(evidence);
            for (int p = 0; p < evidence.PeakCount; p++)
            {
                double intensity = evidence.PeakRelativeIntensity(p);
                if (intensity < model.MinRelativeIntensity)
                {
                    continue;
                }

                bool shared = evidence.IsSharedByAllCandidates(p);
                double explained = shared ? 0 : model.LogExplainedPeak(intensity, coverage);
                for (int c = 0; c < candidates; c++)
                {
                    if (evidence.Matches(c, p))
                    {
                        observed[c]++;
                        logLikelihoods[c] += explained;
                    }
                }
            }

            if (evidence.HasPredictedFragmentCounts)
            {
                double logMissing = model.LogMissingFragment;
                for (int c = 0; c < candidates; c++)
                {
                    logLikelihoods[c] += Math.Max(0, evidence.PredictedFragmentCount(c) - observed[c]) * logMissing;
                }
            }
            return logLikelihoods;
        }

        /// <summary>
        /// Learns the log prior of each site class by expectation maximization over the confident target matches with more
        /// than one candidate: each class's new weight is its expected count of sites under the current posteriors, plus
        /// a pseudocount, normalized within its chemistry. A chemistry on fewer than <see cref="MinMatchesForPrior"/> such
        /// matches keeps a flat prior.
        /// </summary>
        private double[] LearnPrior(ScoredMatch[] scored, SiteClassRegistry siteClasses, LocalizationProbabilityResults results)
        {
            int classCount = siteClasses.ClassCount;
            var training = scored.Where(m => m.IsTraining && m.LogLikelihoods.Length > 1).ToList();

            // which chemistries have enough matches to learn from
            var matchesPerGroup = new int[siteClasses.GroupCount];
            foreach (ScoredMatch match in training)
            {
                foreach (int group in match.Classes.SelectMany(c => c).Select(k => siteClasses.GroupOf(k)).Distinct())
                {
                    matchesPerGroup[group]++;
                }
            }
            bool[] learned = matchesPerGroup.Select(n => n >= MinMatchesForPrior).ToArray();
            results.ChemistriesLearned = learned.Count(l => l);
            results.Chemistries = siteClasses.GroupCount;

            // start flat within each chemistry
            var prior = new double[classCount];
            for (int k = 0; k < classCount; k++)
            {
                prior[k] = 1.0 / siteClasses.ClassesInGroup(siteClasses.GroupOf(k));
            }

            for (int iteration = 1; iteration <= MaxIterations; iteration++)
            {
                double[] logPrior = prior.Select(p => Math.Log(p)).ToArray();
                var expected = new double[classCount];
                foreach (ScoredMatch match in training)
                {
                    double[] posterior = Posterior(match, logPrior);
                    for (int c = 0; c < posterior.Length; c++)
                    {
                        foreach (int k in match.Classes[c])
                        {
                            expected[k] += posterior[c];
                        }
                    }
                }

                var groupTotals = new double[siteClasses.GroupCount];
                for (int k = 0; k < classCount; k++)
                {
                    groupTotals[siteClasses.GroupOf(k)] += expected[k] + Pseudocount;
                }

                var change = new double[siteClasses.GroupCount];
                for (int k = 0; k < classCount; k++)
                {
                    int group = siteClasses.GroupOf(k);
                    if (!learned[group])
                    {
                        continue;
                    }
                    double updated = (expected[k] + Pseudocount) / groupTotals[group];
                    change[group] += Math.Abs(updated - prior[k]);
                    prior[k] = updated;
                }

                results.Iterations = iteration;
                if (change.All(d => d < Convergence))
                {
                    break;
                }
            }

            results.Prior = siteClasses.Describe(prior, learned);
            return prior.Select(p => Math.Log(p)).ToArray();
        }

        /// <summary>
        /// The posterior probability of each candidate: its likelihood times the prior of its sites, normalized.
        /// </summary>
        private static double[] Posterior(ScoredMatch match, double[] logPrior)
        {
            int n = match.LogLikelihoods.Length;
            var logPosterior = new double[n];
            double max = double.NegativeInfinity;
            for (int c = 0; c < n; c++)
            {
                double value = match.LogLikelihoods[c];
                foreach (int k in match.Classes[c])
                {
                    value += logPrior[k];
                }
                logPosterior[c] = value;
                max = Math.Max(max, value);
            }

            // Without a pseudocount a learned prior can be zero for every candidate a match has (PTMiner sets the prior
            // of classes it never sees to zero); the spectrum alone then decides.
            if (double.IsNegativeInfinity(max))
            {
                Array.Copy(match.LogLikelihoods, logPosterior, n);
                max = logPosterior.Max();
            }

            double sum = 0;
            for (int c = 0; c < n; c++)
            {
                logPosterior[c] = Math.Exp(logPosterior[c] - max);
                sum += logPosterior[c];
            }
            for (int c = 0; c < n; c++)
            {
                logPosterior[c] /= sum;
            }
            return logPosterior;
        }

        /// <summary>
        /// Ranks the confident targets with more than one placement by best probability. The false localization rate of
        /// the top k is the mean of one minus their best probabilities; a match's q-value is the lowest such rate of any
        /// set that includes it, and tied probabilities share a rate.
        /// </summary>
        private static void AssignQValues(ScoredMatch[] scored, LocalizationProbabilityResults results)
        {
            var ranked = scored
                .Where(m => m.IsTraining && m.Psm.LocalizationProbabilities.IsLocalizationInQuestion)
                .Select(m => m.Psm.LocalizationProbabilities)
                .OrderByDescending(p => p.BestProbability)
                .ToList();

            var rates = new double[ranked.Count];
            double sumOfErrors = 0;
            int i = 0;
            while (i < ranked.Count)
            {
                int end = i;
                while (end < ranked.Count && ranked[end].BestProbability == ranked[i].BestProbability)
                {
                    sumOfErrors += 1 - ranked[end].BestProbability;
                    end++;
                }
                for (int t = i; t < end; t++)
                {
                    rates[t] = sumOfErrors / end;
                }
                i = end;
            }

            double lowest = double.PositiveInfinity;
            for (int r = ranked.Count - 1; r >= 0; r--)
            {
                lowest = Math.Min(lowest, rates[r]);
                ranked[r].QValue = lowest;
            }

            results.MatchesInQuestion = ranked.Count;
            results.MatchesAtOnePercentFlr = ranked.Count(p => p.QValue <= 0.01);
        }

        /// <summary>
        /// One match's candidate log likelihoods and the prior classes of each candidate's sites.
        /// </summary>
        private sealed class ScoredMatch
        {
            public ScoredMatch(SpectralMatch psm, double[] logLikelihoods, bool isTraining)
            {
                Psm = psm;
                LogLikelihoods = logLikelihoods;
                IsTraining = isTraining;
            }

            public SpectralMatch Psm { get; }
            public double[] LogLikelihoods { get; }
            public bool IsTraining { get; }
            public int[][] Classes { get; set; }
        }

        /// <summary>
        /// Numbers the prior classes: (chemistry, residue, position class). A chemistry is an identified modification's
        /// formula (or mass, when it has none), or a cluster of open-search mass differences.
        /// </summary>
        private sealed class SiteClassRegistry
        {
            private readonly Dictionary<(int Group, char Residue, LocalizationPositionClass Position), int> _classIndex = new();
            private readonly List<(int Group, char Residue, LocalizationPositionClass Position)> _classes = new();
            private readonly Dictionary<string, int> _groupIndex = new();
            private readonly List<string> _groupNames = new();
            private readonly List<int> _classesPerGroup = new();
            private readonly double[] _clusterStarts;

            public SiteClassRegistry(List<SpectralMatch> matches, double clusterTolerance)
            {
                // single-linkage clusters of the mass differences of open-search placements
                var massDifferences = matches
                    .Select(m => m.LocalizationEvidence)
                    .Where(e => e.Placements.Any(p => p.Sites.Any(s => s.Modification == null)))
                    .Select(e => e.MassDifference)
                    .OrderBy(d => d)
                    .ToList();
                var starts = new List<double>();
                for (int i = 0; i < massDifferences.Count; i++)
                {
                    if (i == 0 || massDifferences[i] - massDifferences[i - 1] > clusterTolerance)
                    {
                        starts.Add(massDifferences[i]);
                    }
                }
                _clusterStarts = starts.ToArray();

                foreach (SpectralMatch match in matches)
                {
                    ClassesOf(match.LocalizationEvidence);
                }
            }

            public int ClassCount => _classes.Count;
            public int GroupCount => _groupNames.Count;
            public int GroupOf(int classIndex) => _classes[classIndex].Group;
            public int ClassesInGroup(int group) => _classesPerGroup[group];

            public int[][] ClassesOf(LocalizationEvidence evidence)
            {
                var classes = new int[evidence.CandidateCount][];
                for (int c = 0; c < evidence.CandidateCount; c++)
                {
                    var sites = evidence.Placements[c].Sites;
                    classes[c] = new int[sites.Count];
                    for (int s = 0; s < sites.Count; s++)
                    {
                        int group = GroupIndex(ChemistryOf(sites[s].Modification, evidence.MassDifference));
                        var key = (group, sites[s].Residue, PositionOf(sites[s].ResidueIndex, evidence));
                        if (!_classIndex.TryGetValue(key, out int index))
                        {
                            index = _classes.Count;
                            _classIndex.Add(key, index);
                            _classes.Add(key);
                            _classesPerGroup[group]++;
                        }
                        classes[c][s] = index;
                    }
                }
                return classes;
            }

            /// <summary>The prior of each learned chemistry's classes, most probable first, for the log.</summary>
            public string Describe(double[] prior, bool[] learned)
            {
                var sb = new StringBuilder();
                for (int group = 0; group < GroupCount; group++)
                {
                    if (!learned[group])
                    {
                        continue;
                    }
                    var top = Enumerable.Range(0, ClassCount)
                        .Where(k => _classes[k].Group == group)
                        .OrderByDescending(k => prior[k])
                        .Take(5)
                        .Select(k => $"{_classes[k].Residue}{PositionLabel(_classes[k].Position)} {prior[k].ToString("F3", CultureInfo.InvariantCulture)}");
                    sb.AppendLine($"  {_groupNames[group]}: {string.Join(", ", top)}");
                }
                return sb.ToString();
            }

            private static string PositionLabel(LocalizationPositionClass position) => position switch
            {
                LocalizationPositionClass.Anywhere => "",
                _ => $" ({position})"
            };

            private string ChemistryOf(Modification modification, double massDifference)
            {
                if (modification != null)
                {
                    return modification.ChemicalFormula?.Formula
                        ?? "mass " + (modification.MonoisotopicMass ?? 0).ToString("F3", CultureInfo.InvariantCulture);
                }

                int cluster = Array.BinarySearch(_clusterStarts, massDifference);
                if (cluster < 0)
                {
                    cluster = ~cluster - 1;
                }
                return "mass difference " + _clusterStarts[Math.Max(cluster, 0)].ToString("F3", CultureInfo.InvariantCulture);
            }

            private int GroupIndex(string chemistry)
            {
                if (!_groupIndex.TryGetValue(chemistry, out int group))
                {
                    group = _groupNames.Count;
                    _groupIndex.Add(chemistry, group);
                    _groupNames.Add(chemistry);
                    _classesPerGroup.Add(0);
                }
                return group;
            }

            private static LocalizationPositionClass PositionOf(int residueIndex, LocalizationEvidence evidence)
            {
                if (residueIndex == 0)
                {
                    return evidence.IsProteinNTerminal ? LocalizationPositionClass.ProteinNTerminal : LocalizationPositionClass.PeptideNTerminal;
                }
                if (residueIndex == evidence.PeptideLength - 1)
                {
                    return evidence.IsProteinCTerminal ? LocalizationPositionClass.ProteinCTerminal : LocalizationPositionClass.PeptideCTerminal;
                }
                return LocalizationPositionClass.Anywhere;
            }
        }
    }

    internal class LocalizationProbabilityResults : MetaMorpheusEngineResults
    {
        public LocalizationProbabilityResults(MetaMorpheusEngine s) : base(s)
        {
        }

        public int TrainingPeaks { get; set; }
        public LocalizationLikelihood Likelihood { get; set; }
        public double FragmentObservationProbability { get; set; } = double.NaN;
        public double MeanChanceCoverage { get; set; } = double.NaN;
        public int Chemistries { get; set; }
        public int ChemistriesLearned { get; set; }
        public int Iterations { get; set; }
        public int MatchesInQuestion { get; set; }
        public int MatchesAtOnePercentFlr { get; set; }
        public string Prior { get; set; } = "";

        public override string ToString()
        {
            var sb = new StringBuilder(base.ToString());
            sb.AppendLine();
            sb.AppendLine($"Matched peaks learned from: {TrainingPeaks}");
            sb.AppendLine($"Likelihood: {Likelihood}");
            if (Likelihood == LocalizationLikelihood.IntensityRatio)
            {
                sb.AppendLine($"Fragment observation probability d: {FragmentObservationProbability.ToString("F3", CultureInfo.InvariantCulture)}; mean chance-peak coverage q: {MeanChanceCoverage.ToString("G3", CultureInfo.InvariantCulture)}");
            }
            sb.AppendLine($"Chemistries with a learned prior: {ChemistriesLearned} of {Chemistries}, after {Iterations} iterations");
            sb.AppendLine($"Confident matches with more than one placement: {MatchesInQuestion}; localized at 1% FLR: {MatchesAtOnePercentFlr}");
            sb.Append(Prior);
            return sb.ToString();
        }
    }
}
