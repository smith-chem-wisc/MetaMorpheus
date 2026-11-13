using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Localization
{
    /// <summary>
    /// Which likelihood <see cref="LocalizationProbabilityEngine"/> scores a candidate placement's spectrum with.
    /// </summary>
    public enum LocalizationLikelihood
    {
        /// <summary>
        /// PTMiner's (An et al., MCP 2019, Eq. 5): over the peaks some placements match and others do not, F(I) for a
        /// placement that explains the peak and 1 - F(I) for one that does not, F being the cumulative distribution of
        /// matched-peak intensities. Explaining a peak fainter than the median matched peak counts against a placement.
        /// </summary>
        PtminerCdf,

        /// <summary>
        /// A likelihood ratio against chance (<see cref="PeakIntensityRatioModel"/>): explaining a peak is worth how much
        /// likelier a real fragment is than a chance peak to sit in the predicted window at that intensity, and each
        /// predicted fragment left unobserved costs the probability that a real fragment goes unobserved.
        /// </summary>
        IntensityRatio
    }

    /// <summary>
    /// The model behind <see cref="LocalizationLikelihood.IntensityRatio"/>, learned from confident unmodified matches.
    ///
    /// Each fragment a placement predicts within the scan's range either is observed or is not. If it is real it is
    /// observed with probability d, at a relative intensity with density f_m; a window that holds no real fragment holds
    /// a chance (unexplained) peak with probability q, at a relative intensity with density f_n. Against the hypothesis
    /// that the window holds no fragment of the placement, a predicted fragment
    /// - matched to a peak of relative intensity I counts (1 - d) + d f_m(I) / (q f_n(I)), never less than 1 - d;
    /// - left unobserved counts 1 - d.
    /// A placement's log likelihood ratio is the sum over its predicted fragments; peaks every placement explains count
    /// the same for each and are left out.
    ///
    /// d is the fraction of predicted fragments matched at or above the intensity floor; f_m is the density, over log
    /// intensity, of the matched peaks; f_n the density of the peaks the identified peptide left unexplained, each
    /// weighted by the width of the tolerance window around it (the chance that a predicted fragment lands on it);
    /// q, per spectrum, the fraction of the fragment mass range that unexplained peaks at or above the floor cover.
    /// Densities are histograms over <see cref="BinCount"/> bins of <see cref="BinWidth"/> decades from 10^<see cref="LogMin"/>,
    /// smoothed and mixed with a little of the uniform density so that no ratio is 0 or infinite.
    /// </summary>
    public sealed class PeakIntensityRatioModel
    {
        /// <summary>log10 of the lower edge of the first intensity bin.</summary>
        public const double LogMin = -4;

        /// <summary>The width of an intensity bin, in decades.</summary>
        public const double BinWidth = 0.25;

        /// <summary>Bins cover relative intensities 10^-4 to 1.</summary>
        public const int BinCount = 16;

        /// <summary>The weight of the uniform density each learned density is mixed with.</summary>
        private const double UniformWeight = 0.01;

        /// <summary>A spectrum's chance-peak coverage is kept within these, so that one with no unexplained peaks does
        /// not make every match certain.</summary>
        internal const double MinCoverage = 1e-4;
        internal const double MaxCoverage = 0.5;

        private readonly int _firstBin;
        private readonly double[] _logRatio = new double[BinCount];

        /// <param name="matchedRelativeIntensities">Relative intensities of the peaks confident unmodified matches explain.</param>
        /// <param name="noiseCoverages">For confident unmodified matches, the chance-peak coverage per intensity bin
        /// (<see cref="LocalizationEvidence.NoiseCoverage"/>).</param>
        /// <param name="predictedFragments">The fragments those matches predicted within their scans' ranges.</param>
        /// <param name="observedFragments">How many of them matched a peak at or above the floor.</param>
        public PeakIntensityRatioModel(IEnumerable<double> matchedRelativeIntensities, IEnumerable<IReadOnlyList<float>> noiseCoverages,
            long predictedFragments, long observedFragments, double minRelativeIntensity = MatchedPeakIntensityModel.DefaultMinRelativeIntensity)
        {
            MinRelativeIntensity = minRelativeIntensity;
            _firstBin = FirstBinAtOrAbove(minRelativeIntensity);

            var matched = new double[BinCount];
            foreach (double intensity in matchedRelativeIntensities.Where(i => i >= minRelativeIntensity))
            {
                matched[BinOf(intensity)]++;
                TrainingPeakCount++;
            }

            var noise = new double[BinCount];
            double coverageSum = 0;
            foreach (IReadOnlyList<float> coverage in noiseCoverages)
            {
                for (int b = _firstBin; b < BinCount; b++)
                {
                    noise[b] += coverage[b];
                }
                coverageSum += CoverageAtOrAbove(coverage, _firstBin);
                NoiseTrainingMatches++;
            }
            MeanCoverage = NoiseTrainingMatches > 0 ? Math.Clamp(coverageSum / NoiseTrainingMatches, MinCoverage, MaxCoverage) : double.NaN;

            ObservationProbability = predictedFragments > 0
                ? Math.Clamp((double)observedFragments / predictedFragments, 0.01, 0.99)
                : double.NaN;

            IsTrained = TrainingPeakCount > 0 && noise.Sum() > 0 && predictedFragments > 0;
            if (!IsTrained)
            {
                return;
            }

            double[] fm = Density(matched, _firstBin);
            double[] fn = Density(noise, _firstBin);
            for (int b = _firstBin; b < BinCount; b++)
            {
                _logRatio[b] = Math.Log(fm[b]) - Math.Log(fn[b]);
            }
        }

        public double MinRelativeIntensity { get; }

        /// <summary>Whether there were matched peaks, unexplained peaks and predicted fragments to learn from.</summary>
        public bool IsTrained { get; }

        /// <summary>The number of matched peaks f_m was learned from.</summary>
        public int TrainingPeakCount { get; }

        /// <summary>The number of matches f_n was learned from.</summary>
        public int NoiseTrainingMatches { get; }

        /// <summary>d: the probability that a predicted fragment of the right peptide is matched at or above the floor.</summary>
        public double ObservationProbability { get; }

        /// <summary>The mean chance-peak coverage q of the training spectra, used for a match that recorded none.</summary>
        public double MeanCoverage { get; }

        /// <summary>log(1 - d), what each predicted but unobserved fragment counts.</summary>
        public double LogMissingFragment => Math.Log(1 - ObservationProbability);

        /// <summary>
        /// log(f_m(I) / f_n(I)), interpolated linearly in log intensity between bin centres and constant beyond the
        /// outermost centres.
        /// </summary>
        public double LogDensityRatio(double relativeIntensity)
        {
            double position = (Math.Log10(Math.Max(relativeIntensity, 1e-300)) - LogMin) / BinWidth - 0.5; // in bin-centre units
            if (position <= _firstBin)
            {
                return _logRatio[_firstBin];
            }
            if (position >= BinCount - 1)
            {
                return _logRatio[BinCount - 1];
            }
            int lower = (int)Math.Floor(position);
            double t = position - lower;
            return (1 - t) * _logRatio[lower] + t * _logRatio[lower + 1];
        }

        /// <summary>
        /// What a predicted fragment matched to a peak of relative intensity I counts, given the spectrum's chance-peak
        /// coverage q: log((1 - d) + d f_m(I) / (q f_n(I))).
        /// </summary>
        public double LogExplainedPeak(double relativeIntensity, double coverage)
        {
            double d = ObservationProbability;
            return Math.Log((1 - d) + d * Math.Exp(LogDensityRatio(relativeIntensity)) / coverage);
        }

        /// <summary>
        /// q for one match: the fraction of its fragment mass range covered by unexplained peaks at or above the floor,
        /// or the training mean when the match recorded none.
        /// </summary>
        public double CoverageOf(LocalizationEvidence evidence)
        {
            if (evidence.NoiseCoverage == null)
            {
                return MeanCoverage;
            }
            return Math.Clamp(CoverageAtOrAbove(evidence.NoiseCoverage, _firstBin), MinCoverage, MaxCoverage);
        }

        /// <summary>The bin of a relative intensity, clamped to the bins there are.</summary>
        public static int BinOf(double relativeIntensity)
        {
            if (!(relativeIntensity > 0))
            {
                return 0;
            }
            int bin = (int)Math.Floor((Math.Log10(relativeIntensity) - LogMin) / BinWidth);
            return Math.Clamp(bin, 0, BinCount - 1);
        }

        /// <summary>The first bin whose lower edge is at or above the floor (to within rounding).</summary>
        private static int FirstBinAtOrAbove(double minRelativeIntensity)
        {
            if (!(minRelativeIntensity > 0))
            {
                return 0;
            }
            int bin = (int)Math.Ceiling((Math.Log10(minRelativeIntensity) - LogMin) / BinWidth - 1e-9);
            return Math.Clamp(bin, 0, BinCount - 1);
        }

        private static double CoverageAtOrAbove(IReadOnlyList<float> coverage, int firstBin)
        {
            double sum = 0;
            for (int b = firstBin; b < BinCount; b++)
            {
                sum += coverage[b];
            }
            return sum;
        }

        /// <summary>
        /// Normalizes a histogram over the bins from <paramref name="firstBin"/>, smooths it with a (1/4, 1/2, 1/4)
        /// kernel, and mixes it with the uniform density.
        /// </summary>
        private static double[] Density(double[] counts, int firstBin)
        {
            int used = BinCount - firstBin;
            double total = 0;
            for (int b = firstBin; b < BinCount; b++)
            {
                total += counts[b];
            }

            var smoothed = new double[BinCount];
            for (int b = firstBin; b < BinCount; b++)
            {
                double sum = 0, weight = 0;
                for (int k = -1; k <= 1; k++)
                {
                    int j = b + k;
                    if (j < firstBin || j >= BinCount)
                    {
                        continue;
                    }
                    double w = k == 0 ? 0.5 : 0.25;
                    sum += w * counts[j] / total;
                    weight += w;
                }
                smoothed[b] = sum / weight;
            }

            double smoothedTotal = 0;
            for (int b = firstBin; b < BinCount; b++)
            {
                smoothedTotal += smoothed[b];
            }
            var density = new double[BinCount];
            for (int b = firstBin; b < BinCount; b++)
            {
                density[b] = (1 - UniformWeight) * smoothed[b] / smoothedTotal + UniformWeight / used;
            }
            return density;
        }
    }
}
