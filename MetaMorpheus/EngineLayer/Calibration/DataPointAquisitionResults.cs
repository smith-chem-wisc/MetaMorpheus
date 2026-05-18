using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.Statistics;
using System.IO;

namespace EngineLayer.Calibration
{
    /// <summary>
    /// Returns PSMs that can be used for calibration based on tolerance limits passed from the calibration task
    /// </summary>
    public class DataPointAquisitionResults : MetaMorpheusEngineResults
    {
        public DataPointAquisitionResults(
            MetaMorpheusEngine dataPointAcquisitionEngine,
            List<SpectralMatch> psms,
            List<LabeledDataPoint> ms1List,
            List<LabeledDataPoint> ms2List,
            int numMs1MassChargeCombinationsConsidered,
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
            int numMs2MassChargeCombinationsConsidered,
            int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks)
            : base(dataPointAcquisitionEngine)
        {
            Psms = psms;

            Ms1List = ms1List;
            Ms2List = ms2List;

            var ms1Range = Ms1List.Select(b => b.ExperimentalMz - b.TheoreticalMz).ToArray();
            var ms2Range = Ms2List.Select(b => b.ExperimentalMz - b.TheoreticalMz).ToArray();

            var ms1PpmRange = Ms1List.Select(b => (b.ExperimentalMz - b.TheoreticalMz) / b.TheoreticalMz).ToArray();
            var ms2PpmRange = Ms2List.Select(b => (b.ExperimentalMz - b.TheoreticalMz) / b.TheoreticalMz).ToArray();

            var precursorErrors = psms.Select(p => (p.ScanPrecursorMass - p.BioPolymerWithSetModsMonoisotopicMass.Value) / p.BioPolymerWithSetModsMonoisotopicMass.Value * 1e6).ToList();
            PsmPrecursorIqrPpmError = precursorErrors.InterquartileRange();
            PsmPrecursorMedianPpmError = precursorErrors.Median();

            var productErrors = psms.Where(p => p.MatchedFragmentIons != null).SelectMany(p => p.MatchedFragmentIons).Select(p => (p.Mz.ToMass(p.Charge) - p.NeutralTheoreticalProduct.NeutralMass) / p.NeutralTheoreticalProduct.NeutralMass * 1e6).ToList();
            PsmProductIqrPpmError = productErrors.InterquartileRange();
            PsmProductMedianPpmError = productErrors.Median();
        }

        public List<LabeledDataPoint> Ms1List { get; }
        public List<LabeledDataPoint> Ms2List { get; }

        public readonly double PsmPrecursorMedianPpmError;
        public readonly double PsmProductMedianPpmError;
        public readonly double PsmPrecursorIqrPpmError;
        public readonly double PsmProductIqrPpmError;
        public readonly List<SpectralMatch> Psms;

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("MS1 calibration datapoint count: " + Ms1List.Count);
            sb.AppendLine("MS1 ppm error median: " + Math.Round(PsmPrecursorMedianPpmError, 3));
            sb.AppendLine("MS1 ppm error interquartile range: " + Math.Round(PsmPrecursorIqrPpmError, 3));

            sb.AppendLine("MS2 calibration datapoint count: " + Ms2List.Count);
            sb.AppendLine("MS2 ppm error median: " + Math.Round(PsmProductMedianPpmError, 3));
            sb.AppendLine("MS2 ppm error interquartile range: " + Math.Round(PsmProductIqrPpmError, 3));
            return sb.ToString();
        }
    }
}