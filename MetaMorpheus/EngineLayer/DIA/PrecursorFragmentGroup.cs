using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace EngineLayer.DIA
{
    /// <summary>
    /// PrecursorFragmentsGroup represents a group of precursor-fragment pairs belonging to the same precursor XIC.
    /// </summary>
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

        public static double CalculateXicCorrelation(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            if (xic1.XYData == null || xic2.XYData == null)
            {
                return CalculateXicCorrelationRawPeaks(xic1, xic2);
            }
            double corr = CalculateCorrelationXYData(xic1.XYData, xic2.XYData);
            return corr;
        }

        public static double CalculateXicCorrelationRawPeaks(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            var start = Math.Max(xic1.StartScanIndex, xic2.StartScanIndex);
            var end = Math.Min(xic1.EndScanIndex, xic2.EndScanIndex);

            if (end - start + 1 < 3) // Ideally we need at least 3 points to calculate correlation
            {
                return double.NaN;
            }
            var y1 = new float[end - start + 1];
            var y2 = new float[end - start + 1];
            var peakScanIndices1 = xic1.Peaks.Select(p => p.ZeroBasedScanIndex).ToArray();
            var peakScanIndices2 = xic2.Peaks.Select(p => p.ZeroBasedScanIndex).ToArray();
            for (int i = 0; i < y1.Length; i++)
            {
                int index1 = Array.BinarySearch(peakScanIndices1, i + start);
                int index2 = Array.BinarySearch(peakScanIndices2, i + start);
                if (index1 >= 0)
                {
                    y1[i] = xic1.Peaks[index1].Intensity;
                }
                else
                {
                    y1[i] = 0;
                }
                if (index2 >= 0)
                {
                    y2[i] = xic2.Peaks[index2].Intensity;
                }
                else
                {
                    y2[i] = 0;
                }
            }
            return PearsonCorrelation(y1, y2);
        }

        public static double CalculateCorrelationXYData<T> ((T, T)[] xy1, (T, T)[] xy2) where T : INumber<T>
        {
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
            if (xic1.EndScanIndex <= xic2.StartScanIndex || xic2.EndScanIndex <= xic1.StartScanIndex)
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
            //sort all fragment XICs by their m/z value
            pfGroup.PFpairs.Sort((a, b) => a.FragmentXic.AveragedMassOrMz.CompareTo(b.FragmentXic.AveragedMassOrMz));
            // This is currently taking the mz value of the first peak as the representative mz for the fragment XIC; it is only for testing purposes
            // It will be changed to the mz value of the highest peak or the averaged mz value of the XIC when there is a release for the updated mzLib.
            var mzs = pfGroup.PFpairs.Select(pf => (double)pf.FragmentXic.Peaks.First().M).ToArray();
            var intensities = pfGroup.PFpairs.Select(pf => (double)pf.FragmentXic.Peaks.First().Intensity).ToArray();
            var newMs2Scan = new MsDataScan(new MzSpectrum(mzs, intensities, false), pfGroup.PFgroupIndex, 2, true, Polarity.Positive, pfGroup.PrecursorXic.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null, MZAnalyzerType.Unknown, intensities.Sum(), null, null, null, oneBasedPrecursorScanNumber: pfGroup.PFgroupIndex);

            // assign neutral experimental fragments for the new MS2 scan
            IsotopicEnvelope[] neutralExperimentalFragments = null;
            switch (pseudoMs2ConstructionType)
            {
                // If the fragment XIC is mzPeak-based, we deconvolute the scan as normal
                case PseudoMs2ConstructionType.MzPeak:
                    neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(newMs2Scan, commonParameters);
                    break;
                // If the fragment XIC is mass-based, we create an isotopic envelope as the neutral mass fragment for each fragment XIC
                case PseudoMs2ConstructionType.Mass:
                    neutralExperimentalFragments = pfGroup.PFpairs.Select(pf => new IsotopicEnvelope(1,
                            new List<(double mz, double intensity)> { (1, 1) }, pf.FragmentXic.ApexPeak.M, pf.FragmentXic.Peaks.Cast<IndexedMass>().First().Charge, 1, 0)).ToArray();
                    break;
                default:
                    throw new MetaMorpheusException("Invalid pseudo MS2 construction type specified.");
            }
            // add precursor information
            var charge = pfGroup.PrecursorXic.Peaks.First() is IndexedMass im ? im.Charge : 1;
            var monoMz = pfGroup.PrecursorXic.Peaks.First() is IndexedMass im2 ? im2.M.ToMz(charge) : pfGroup.PrecursorXic.ApexPeak.M;
            Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, monoMz, charge, dataFilePath,
                commonParameters, neutralExperimentalFragments);

            return scanWithprecursor;
        }
    }
}
