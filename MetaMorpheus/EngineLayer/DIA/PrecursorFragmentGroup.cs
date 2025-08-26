using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    /// <summary>
    /// PrecursorFragmentsGroup represents a group of precursor-fragment pairs belonging to the same precursor XIC.
    /// <summary>
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

        public static double CalculateXicCorrXYData_Umpire(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2, int NoPointPerInterval)
        {
            double start = Math.Max(xic1.XYData[0].Item1, xic2.XYData[0].Item1);

            //int num = Math.Min(curve1.SmoothedData.Count, curve2.SmoothedData.Count) / 2;
            int num = Math.Max(xic1.XYData.Count(), xic2.XYData.Count()) / 2;
            double timeInterval = (double)2 / NoPointPerInterval;

            if (num < 6)
            {
                return 0f;
            }
            double[] arrayA = new double[num];
            double[] arrayB = new double[num];

            int i = 0;
            double low = start;
            double up = start + timeInterval;

            for (int j = 0; j < xic1.XYData.Length; j++)
            {
                while (xic1.XYData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (xic1.XYData[j].Item1 >= low && xic1.XYData[j].Item1 < up)
                {
                    if (xic1.XYData[j].Item2 > arrayA[i])
                    {
                        arrayA[i] = xic1.XYData[j].Item2;
                    }
                }
            }

            i = 0;
            low = start;
            up = start + timeInterval;

            for (int j = 0; j < xic2.XYData.Length; j++)
            {
                while (xic2.XYData[j].Item1 > up)
                {
                    i++;
                    low = up;
                    up = low + timeInterval;
                }
                if (i >= num)
                {
                    break;
                }
                if (xic2.XYData[j].Item1 >= low && xic2.XYData[j].Item1 < up)
                {
                    if (xic2.XYData[j].Item2 > arrayB[i])
                    {
                        arrayB[i] = xic2.XYData[j].Item2;
                    }
                }
            }

            for (int idx = 1; idx < num - 1; idx++)
            {
                if (arrayA[idx] == 0f)
                {
                    arrayA[idx] = (arrayA[idx - 1] + arrayA[idx + 1]) / 2;
                }
                if (arrayB[idx] == 0f)
                {
                    arrayB[idx] = (arrayB[idx - 1] + arrayB[idx + 1]) / 2;
                }
            }

            double corr = Correlation.Pearson(arrayA, arrayB);
            return corr;
        }

        public static double CalculateCorrelation<T> ((T, T)[] xy1, (T, T)[] xy2) where T : INumber<T>
        {
            if (xy1 == null || xy2 == null)
            {
                return double.NaN;
            }

            T start = T.Max(xy1[0].Item1, xy2[0].Item1);
            T end = T.Min(xy1[xy1.Length - 1].Item1, xy2[xy2.Length - 1].Item1);

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

        public static double PearsonCorrelation<T>(T[] x, T[] y) where T : INumber<T>
        {
            if (x == null || y == null) throw new ArgumentNullException();
            if (x.Length != y.Length) throw new ArgumentException("Arrays must have the same length.");

            int n = x.Length;
            T sumX = T.Zero, sumY = T.Zero, sumXY = T.Zero, sumX2 = T.Zero, sumY2 = T.Zero;
            for (int i = 0; i < n; i++)
            {
                sumX += x[i];
                sumY += y[i];
                sumXY += x[i] * y[i];
                sumX2 += x[i] * x[i];
                sumY2 += y[i] * y[i];
            }

            double meanX = double.CreateChecked(sumX) / n;
            double meanY = double.CreateChecked(sumY) / n;

            double numerator = double.CreateChecked(sumXY) - n * meanX * meanY;
            double denominatorX = double.CreateChecked(sumX2) - n * meanX * meanX;
            double denominatorY = double.CreateChecked(sumY2) - n * meanY * meanY;
            double denominator = Math.Sqrt(denominatorX * denominatorY);

            if (denominator == 0) return 0; // Avoid division by zero

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

        public static double CalculateXicOverlapRatio_Umpire(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            double overlap = 0;
            var ms1rtrange = xic1.EndRT - xic1.StartRT;
            var ms2rtrange = xic2.EndRT - xic2.StartRT;
            if (xic1.StartRT >= xic2.StartRT && xic1.StartRT <= xic2.EndRT && xic1.EndRT >= xic2.EndRT)
            {
                overlap = (xic2.EndRT - xic1.StartRT) / ms1rtrange;
            }
            else if (xic1.EndRT >= xic2.StartRT && xic1.EndRT <= xic2.EndRT && xic1.StartRT <= xic2.StartRT)
            {
                overlap = (xic1.EndRT - xic2.StartRT) / ms1rtrange;
            }
            else if (xic1.StartRT <= xic2.StartRT && xic1.EndRT >= xic2.EndRT)
            {
                overlap = ms2rtrange / ms1rtrange;
            }
            else if (xic1.StartRT >= xic2.StartRT && xic1.EndRT <= xic2.EndRT)
            {
                overlap = 1;
            }
            return overlap;
        }

        public static Ms2ScanWithSpecificMass GetPseudoMs2ScanFromPfGroup(PrecursorFragmentsGroup pfGroup, PseudoMs2ConstructionType pseudoMs2ConstructionType, CommonParameters commonParameters, string dataFilePath)
        {
            pfGroup.PFpairs.Sort((a, b) => a.FragmentXic.AveragedM.CompareTo(b.FragmentXic.AveragedM));
            var mzs = pfGroup.PFpairs.Select(pf => (double)pf.FragmentXic.Peaks.First().M).ToArray();
            var intensities = pfGroup.PFpairs.Select(pf => (double)pf.FragmentXic.Peaks.First().Intensity).ToArray();
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

        public void SetFragmentRankForPfPairs()
        {
            if (PFpairs.IsNullOrEmpty()) return;
            PFpairs.Sort((a, b) => b.Correlation.Value.CompareTo(a.Correlation.Value));
            for (int i = 0; i < PFpairs.Count; i++)
            {
                PFpairs[i].FragmentRank = i + 1;
            }
        }
    }
}
