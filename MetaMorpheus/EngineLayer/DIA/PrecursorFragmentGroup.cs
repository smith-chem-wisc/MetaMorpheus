using EngineLayer.DIA.Enums;
using MassSpectrometry;
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

        public PrecursorFragmentsGroup(ExtractedIonChromatogram precursorXic, List<PrecursorFragmentPair> pfPairs = null)
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

        public static double CalculateXicOverlapRatio(ExtractedIonChromatogram xic1, ExtractedIonChromatogram xic2)
        {
            if ((xic1.StartScanIndex < xic2.StartScanIndex && xic1.EndScanIndex < xic2.StartScanIndex) ||
                (xic2.StartScanIndex < xic1.StartScanIndex && xic2.EndScanIndex < xic1.StartScanIndex))
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

        public static Ms2ScanWithSpecificMass ConstructNewMs2Scans(PrecursorFragmentsGroup pfGroup, CommonParameters commonParameters, PseudoMs2ConstructionType pseudoMs2Type, string dataFilePath)
        {
            switch (pseudoMs2Type)
            {

                default: return null;
            }
        }

        //public static Ms2ScanWithSpecificMass GetPseudoMs2Scan_mzPeak(PrecursorFragmentsGroup pfGroup, CommonParameters commonParameters, string dataFilePath)
        //{
        //    //var mzs = pfGroup.PFpairs.Select(pf => pf.FragmentXic.AveragedM).ToArray();
        //    //var intensities = pfGroup.PFpairs.Select(pf => pf.FragmentXic.AveragedM).ToArray();
        //    //var spectrum = new MzSpectrum(mzs, intensities, false);
        //    //var newMs2Scan = new MsDataScan(spectrum, pfGroup.PFgroupIndex, 2, true, Polarity.Positive, pfGroup.PrecursorXic.ApexRT, new MzRange(mzs.Min(), mzs.Max()), null, MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, null, oneBasedPrecursorScanNumber: pfGroup.PrecursorXic.Index);
        //    //var neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(newMs2Scan, commonParameters);
        //    //var charge = pfGroup.PrecursorXic.Charge;
        //    //var monoMz = pfGroup.PrecursorPeakCurve.MonoisotopicMass.ToMz(charge);
        //    ////should highestPeakMz used in ms2withmass???
        //    //Ms2ScanWithSpecificMass scanWithprecursor = new Ms2ScanWithSpecificMass(newMs2Scan, monoMz, charge, dataFilePath,
        //    //    commonParameters, neutralExperimentalFragments);

        //    return scanWithprecursor;
        //}
    }
}
