using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;


namespace EngineLayer.DIA
{
    /// <summary>
    /// PrecursorFragmentPair class represents a pair of precursor XIC and fragment XIC, along with their relationship metrics such as correlation and overlap.
    /// </summary>
    public class PrecursorFragmentPair
    {
        public ExtractedIonChromatogram PrecursorXic { get; set; }
        public ExtractedIonChromatogram FragmentXic { get; set; }
        public double? Correlation { get; set; }
        public double? Overlap { get; set; }
        public int? PrecursorRank { get; set; }
        public int? FragmentRank { get; set; }

        public PrecursorFragmentPair(ExtractedIonChromatogram precursorXic, ExtractedIonChromatogram fragmentXic, double? correlation = null, double? overlap = null)
        {
            PrecursorXic = precursorXic;
            FragmentXic = fragmentXic;
            Correlation = correlation;
            Overlap = overlap;
        }

        public static void SetPrecursorRankForPfPairs(IEnumerable<PrecursorFragmentPair> allPfPairs)
        {
            if (allPfPairs.IsNullOrEmpty()) return;
            var fragmentPairMap = new Dictionary<ExtractedIonChromatogram, List<PrecursorFragmentPair>>();
            foreach (var pair in allPfPairs)
            {
                if (!fragmentPairMap.ContainsKey(pair.FragmentXic))
                {
                    fragmentPairMap[pair.FragmentXic] = new List<PrecursorFragmentPair>();
                }
                fragmentPairMap[pair.FragmentXic].Add(pair);
            }
            foreach (var fragmentXic in fragmentPairMap.Keys)
            {
                fragmentPairMap[fragmentXic].Sort((a, b) => b.Correlation.Value.CompareTo(a.Correlation.Value));
                for (int i = 0; i < fragmentPairMap[fragmentXic].Count; i++)
                {
                    fragmentPairMap[fragmentXic][i].PrecursorRank = i + 1;
                }
            }
        }

        public static double CalculateSharedXIC(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            if (xic1.EndScanIndex <= xic2.StartScanIndex || xic2.EndScanIndex <= xic1.StartScanIndex)
            {
                return 0;
            }

            var overlapStart = Math.Max(xic1.StartScanIndex, xic2.StartScanIndex);
            var overlapEnd = Math.Min(xic1.EndScanIndex, xic2.EndScanIndex);
            var maxLength = overlapEnd - overlapStart + 1;

            var overlapArea = new List<(double, double)>();
            var scanCycles1 = xic1.Peaks.Select(p => p.ZeroBasedScanIndex).ToArray();
            var scanCycles2 = xic2.Peaks.Select(p => p.ZeroBasedScanIndex).ToArray();
            var index1 = Array.BinarySearch(scanCycles1, overlapStart);
            var index2 = Array.BinarySearch(scanCycles2, overlapStart);

            for (int i = 0; i < maxLength - 1; i++)
            {
                double diff0 = xic1.NormalizedPeakIntensities[index1 + i] - xic2.NormalizedPeakIntensities[index2 + i];
                double diff1 = xic1.NormalizedPeakIntensities[index1 + i + 1] - xic2.NormalizedPeakIntensities[index2 + i + 1];

                overlapArea.Add((overlapStart + i, Math.Min(xic1.NormalizedPeakIntensities[index1 + i], xic2.NormalizedPeakIntensities[index2 + i])));
                if (diff0 * diff1 < 0)
                {
                    double slope = xic1.NormalizedPeakIntensities[index1 + i + 1] - xic1.NormalizedPeakIntensities[index1 + i];
                    double y = xic1.NormalizedPeakIntensities[index1 + i] + slope * Math.Abs(diff0 / (diff1 - diff0));
                    overlapArea.Add((overlapStart + i + Math.Abs(diff0 / (diff1 - diff0)), y));
                }
            }
            double overlapAUC = CalculateNormalizedArea(overlapArea);
            return overlapAUC;
        }

        public static double CalculateNormalizedArea(List<(double, double)> data)
        {
            double area = data[0].Item2 / 2;
            for (int i = 1; i < data.Count; i++)
            {
                double x1 = data[i - 1].Item1;
                double x2 = data[i].Item1;
                double y1 = data[i - 1].Item2;
                double y2 = data[i].Item2;
                area += (x2 - x1) * (y1 + y2) / 2;
            }
            area += data[data.Count - 1].Item2 / 2;
            return area;
        }
    }
}

