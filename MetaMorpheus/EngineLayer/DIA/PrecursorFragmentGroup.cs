using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PrecursorFragmentGroup
    {
        public ExtractedIonChromatogram PrecursorXic { get; set; }

        public List<PrecursorFragmentPair> PFpairs { get; set; }

        public PrecursorFragmentGroup(ExtractedIonChromatogram precursorXic, List<PrecursorFragmentPair> pfPairs = null)
        {
            PrecursorXic = precursorXic;
            PFpairs = pfPairs ?? new List<PrecursorFragmentPair>();
        }

        public static double CalculateXicCorrelationXYData(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            if (xic1.XYData == null || xic2.XYData == null)
            {
                var xy1 = xic1.Peaks.Select(p => ((double)p.ZeroBasedScanIndex, p.Intensity)).ToArray();
                var xy2 = xic2.Peaks.Select(p => ((double)p.ZeroBasedScanIndex, p.Intensity)).ToArray();
                return CalculateCorrelation(xy1, xy2);
            }
            double corr = CalculateCorrelation(xic1.XYData, xic2.XYData);
            return corr;
        }

        public static double CalculateCorrelation((double, double)[] xy1, (double, double)[] xy2)
        {
            if (xy1 == null || xy2 == null)
            {
                return double.NaN;
            }

            double start = Math.Max(xy1[0].Item1, xy2[0].Item1);
            double end = Math.Min(xy1[xy1.Length - 1].Item1, xy2[xy2.Length - 1].Item1);

            var validxy1 = xy1.Where(p => p.Item1 >= start && p.Item1 <= end).ToArray();
            var validxy2 = xy2.Where(p => p.Item1 >= start && p.Item1 <= end).ToArray();
            int numPoints = Math.Min(validxy1.Length, validxy2.Length);
            if (numPoints < 3)
            {
                return double.NaN;
            }
            var xy = validxy1.Take(numPoints).Zip(validxy2.Take(numPoints), (a, b) => (a.Item2, b.Item2)).ToArray();
            var y1 = xy.Select(p => p.Item1).ToArray();
            var y2 = xy.Select(p => p.Item2).ToArray();
            double corr = MathNet.Numerics.Statistics.Correlation.Pearson(y1, y2);

            return corr;
        }
    }
}
