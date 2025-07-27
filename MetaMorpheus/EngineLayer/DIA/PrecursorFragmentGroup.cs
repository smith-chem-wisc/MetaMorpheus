using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class PrecursorFragmentsGroup
    {
        public ExtractedIonChromatogram PrecursorXic { get; set; }

        public List<PrecursorFragmentPair> PFpairs { get; set; }
        public int PFgroupIndex { get; set; } 

        public PrecursorFragmentsGroup(ExtractedIonChromatogram precursorXic, List<PrecursorFragmentPair> pfPairs)
        {
            PrecursorXic = precursorXic;
            PFpairs = pfPairs;
        }

        

        public static double CalculateXicCorrelationXYData(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            if (xic1.XYData == null || xic2.XYData == null)
            {
                var xy1 = xic1.Peaks.Select(p => ((float)p.ZeroBasedScanIndex, p.Intensity)).ToArray();
                var xy2 = xic2.Peaks.Select(p => ((float)p.ZeroBasedScanIndex, p.Intensity)).ToArray();
                return CalculateCorrelation(xy1, xy2);
            }
            double corr = CalculateCorrelation(xic1.XYData, xic2.XYData);
            return corr;
        }

        public static double CalculateCorrelation((float, float)[] xy1, (float, float)[] xy2)
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
            double corr = PearsonCorrelation(y1, y2);

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
            double corr = Correlation.Pearson(y1, y2);

            return corr;
        }

        public static float PearsonCorrelation(float[] x, float[] y)
        {
            if (x.Length != y.Length || x.Length == 0)
                throw new ArgumentException("Arrays must have the same non-zero length.");

            float meanX = x.Average();
            float meanY = y.Average();

            float numerator = 0f;
            float sumSqX = 0f;
            float sumSqY = 0f;

            for (int i = 0; i < x.Length; i++)
            {
                float diffX = x[i] - meanX;
                float diffY = y[i] - meanY;
                numerator += diffX * diffY;
                sumSqX += diffX * diffX;
                sumSqY += diffY * diffY;
            }

            float denominator = (float)Math.Sqrt(sumSqX * sumSqY);
            if (denominator == 0f) return 0f;

            return numerator / denominator;
        }

        public static double CalculateXicOverlapRatio(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            if (xic1.EndScanIndex < xic2.StartScanIndex || xic2.EndScanIndex < xic1.StartScanIndex)
            {
                return 0;
            }
            var start = Math.Min(xic1.StartScanIndex, xic2.StartScanIndex);
            var end = Math.Max(xic1.EndScanIndex, xic2.EndScanIndex);
            var overlapStart = Math.Max(xic1.StartScanIndex, xic2.StartScanIndex);
            var overlapEnd = Math.Min(xic1.EndScanIndex, xic2.EndScanIndex);

            double overlap = (overlapEnd - overlapStart) / (double)(end - start);
            return overlap;
        }

        public static Ms2ScanWithSpecificMass GetPseudoMs2ScanFromPfGroup(PrecursorFragmentsGroup pfGroup, PseudoMs2ConstructionType pseudoMs2ConstructionType, CommonParameters commonParameters, string dataFilePath)
        {
            pfGroup.PFpairs.Sort((a, b) => b.FragmentXic.AveragedM.CompareTo(a.FragmentXic.AveragedM));
            var mzs = pfGroup.PFpairs.Select(pf => pf.FragmentXic.AveragedM).ToArray();
            var intensities = pfGroup.PFpairs.Select(pf => pf.FragmentXic.Peaks.Max(p => (double)p.Intensity)).ToArray();
            var newMs2Scan = new MsDataScan(new MzSpectrum(mzs, intensities, false), pfGroup.PFgroupIndex, 2, true, Polarity.Positive, pfGroup.PrecursorXic.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null, MZAnalyzerType.Unknown, intensities.Sum(), null, null, null, oneBasedPrecursorScanNumber: pfGroup.PFgroupIndex);

            IsotopicEnvelope[] neutralExperimentalFragments = null;
            switch (pseudoMs2ConstructionType)
            {
                case PseudoMs2ConstructionType.MzPeak:
                    neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(newMs2Scan, commonParameters);
                    break;
                case PseudoMs2ConstructionType.Mass:
                    neutralExperimentalFragments = pfGroup.PFpairs.Select(pf => new IsotopicEnvelope(1,
                            new List<(double mz, double intensity)> { (1, 1) }, pf.FragmentXic.AveragedM, pf.FragmentXic.Peaks.Cast<IndexedMass>().First().Charge, 1, 0)).ToArray();
                    break;
                default:
                    throw new ArgumentException("Invalid pseudo MS2 construction type specified.");
            }
            var charge = pfGroup.PrecursorXic.Peaks.Cast<IndexedMass>().First().Charge;
            var monoMz = pfGroup.PrecursorXic.Peaks.Cast<IndexedMass>().First().M.ToMz(charge);
            Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, monoMz, charge, dataFilePath,
                commonParameters, neutralExperimentalFragments);

            return scanWithprecursor;
        }

        
    }
}
