using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer
{
    /// <summary>
    /// The Andromeda score (Cox et al. 2011, J Proteome Res 10:1794), for putting a MetaMorpheus
    /// PSM next to a MaxQuant one. It is an output column only: nothing in FDR or PEP reads it.
    ///
    /// For a peak depth q, the spectrum keeps the q most intense peaks in each 100 Th window, so a
    /// theoretical fragment matches by chance with p = q/100. With k of n theoretical fragments
    /// matched, the score is -10 log10 of the binomial tail
    ///     P(X >= k) = sum_{j=k..n} C(n,j) p^j (1-p)^(n-j),
    /// computed for every q in 1..maxPeaksPer100Th, and the highest score is reported.
    ///
    /// Approximation: fragments are matched as singly charged m/z against the centroided peaks, with
    /// the search's product tolerance; Andromeda's own peak preprocessing (deisotoping and charge
    /// reduction of high-resolution spectra) is not reproduced.
    /// </summary>
    public static class AndromedaScoring
    {
        public const int DefaultMaxPeaksPer100Th = 10;
        private const double WindowWidthTh = 100.0;

        /// <summary>
        /// -10 log10 of the probability of matching at least <paramref name="matched"/> of
        /// <paramref name="theoretical"/> fragments when each matches by chance with
        /// p = <paramref name="peaksPer100Th"/>/100. Never negative.
        /// </summary>
        public static double BinomialScore(int theoretical, int matched, int peaksPer100Th)
        {
            if (theoretical < 0 || matched < 0 || matched > theoretical)
                throw new ArgumentOutOfRangeException(nameof(matched), $"need 0 <= matched <= theoretical, got matched={matched}, theoretical={theoretical}");
            if (peaksPer100Th < 1 || peaksPer100Th >= WindowWidthTh)
                throw new ArgumentOutOfRangeException(nameof(peaksPer100Th), $"need 1 <= q < 100, got {peaksPer100Th}");
            if (matched == 0)
                return 0; // P(X >= 0) = 1

            double p = peaksPer100Th / WindowWidthTh;
            double logP = Math.Log(p);
            double logQ = Math.Log(1 - p);

            // log C(n,j) built incrementally from log C(n,0) = 0, so large n does not overflow
            var logTerms = new double[theoretical - matched + 1];
            double logChoose = 0;
            for (int j = 1; j <= theoretical; j++)
            {
                logChoose += Math.Log(theoretical - j + 1) - Math.Log(j);
                if (j >= matched)
                    logTerms[j - matched] = logChoose + j * logP + (theoretical - j) * logQ;
            }

            // log-sum-exp
            double max = logTerms.Max();
            double sum = 0;
            foreach (double t in logTerms)
                sum += Math.Exp(t - max);
            double logTail = max + Math.Log(sum);

            // the tail is a probability; rounding can push it a hair above 0 in log space
            return Math.Max(0, -10.0 * logTail / Math.Log(10));
        }

        /// <summary>
        /// The Andromeda score of one hypothesis against one spectrum: the best binomial score over
        /// peak depths q = 1..<paramref name="maxPeaksPer100Th"/>.
        /// </summary>
        public static double Score(MzSpectrum spectrum, IReadOnlyCollection<Product> theoreticalProducts, Tolerance productTolerance, int maxPeaksPer100Th = DefaultMaxPeaksPer100Th)
        {
            int n = theoreticalProducts.Count;
            if (n == 0 || spectrum == null || spectrum.Size == 0)
                return 0;

            // Rank each peak by intensity within its 100 Th window (0 = most intense). A peak
            // survives depth q when its rank is < q, so one pass covers every q.
            int[] rankInWindow = RankPeaksWithinWindows(spectrum.XArray, spectrum.YArray);

            // For each theoretical fragment, the best (lowest) rank of any peak inside tolerance:
            // it is matched at depth q exactly when that rank is < q.
            var matchedAtDepth = new int[maxPeaksPer100Th + 1];
            foreach (Product product in theoreticalProducts)
            {
                double mz = product.NeutralMass.ToMz(1);
                int bestRank = int.MaxValue;
                int i = Array.BinarySearch(spectrum.XArray, productTolerance.GetMinimumValue(mz));
                if (i < 0) i = ~i;
                double maxMz = productTolerance.GetMaximumValue(mz);
                for (; i < spectrum.XArray.Length && spectrum.XArray[i] <= maxMz; i++)
                    bestRank = Math.Min(bestRank, rankInWindow[i]);

                if (bestRank < maxPeaksPer100Th)
                    matchedAtDepth[bestRank + 1]++; // counted at depth bestRank+1 and every deeper q
            }

            double best = 0;
            int k = 0;
            for (int q = 1; q <= maxPeaksPer100Th; q++)
            {
                k += matchedAtDepth[q];
                best = Math.Max(best, BinomialScore(n, k, q));
            }
            return best;
        }

        /// <summary>
        /// Fills <see cref="SpectralMatch.AndromedaScore"/> for every PSM, targets and decoys alike,
        /// whichever engine produced it. Ambiguous PSMs report the best score among their tied
        /// hypotheses.
        /// </summary>
        public static void ScorePsms(IEnumerable<SpectralMatch> psms, Ms2ScanWithSpecificMass[] scansSortedByMass, CommonParameters commonParameters, int maxPeaksPer100Th = DefaultMaxPeaksPer100Th)
        {
            SpectralMatch[] toScore = psms.Where(p => p != null).ToArray();
            Parallel.For(0, toScore.Length, new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, i =>
            {
                SpectralMatch psm = toScore[i];
                Ms2ScanWithSpecificMass scan = scansSortedByMass[psm.ScanIndex];
                DissociationType dissociationType = commonParameters.DissociationType == DissociationType.Autodetect
                    ? scan.TheScan.DissociationType ?? DissociationType.Unknown
                    : commonParameters.DissociationType;

                double best = 0;
                var products = new List<Product>();
                foreach (IBioPolymerWithSetMods hypothesis in psm.BestMatchingBioPolymersWithSetMods.Select(h => h.SpecificBioPolymer))
                {
                    products.Clear();
                    hypothesis.Fragment(dissociationType, commonParameters.DigestionParams.FragmentationTerminus, products, commonParameters.FragmentationParameters);
                    best = Math.Max(best, Score(scan.TheScan.MassSpectrum, products, commonParameters.ProductMassTolerance, maxPeaksPer100Th));
                }
                psm.AndromedaScore = best;
            });
        }

        private static int[] RankPeaksWithinWindows(double[] mz, double[] intensity)
        {
            var rank = new int[mz.Length];
            int start = 0;
            while (start < mz.Length)
            {
                double window = Math.Floor(mz[start] / WindowWidthTh);
                int end = start;
                while (end < mz.Length && Math.Floor(mz[end] / WindowWidthTh) == window)
                    end++;

                int[] order = Enumerable.Range(start, end - start)
                    .OrderByDescending(j => intensity[j]).ThenBy(j => j).ToArray();
                for (int r = 0; r < order.Length; r++)
                    rank[order[r]] = r;
                start = end;
            }
            return rank;
        }
    }
}
