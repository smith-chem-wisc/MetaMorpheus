using Chemistry;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Calibration
{
    public class DataPointAquisitionResults : MetaMorpheusEngineResults
    {
        public DataPointAquisitionResults(
            MetaMorpheusEngine dataPointAcquisitionEngine,
            List<PeptideSpectralMatch> psms,
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

            Ms1InfoTh = Ms1List.Select(b => b.ExperimentalMz - b.TheoreticalMz).MeanStandardDeviation();
            Ms2InfoTh = Ms2List.Select(b => b.ExperimentalMz - b.TheoreticalMz).MeanStandardDeviation();

            Ms1InfoPpm = Ms1List.Select(b => (b.ExperimentalMz - b.TheoreticalMz) / b.TheoreticalMz).MeanStandardDeviation();
            Ms2InfoPpm = Ms2List.Select(b => (b.ExperimentalMz - b.TheoreticalMz) / b.TheoreticalMz).MeanStandardDeviation();

            NumMs1MassChargeCombinationsConsidered = numMs1MassChargeCombinationsConsidered;
            NumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
            NumMs2MassChargeCombinationsConsidered = numMs2MassChargeCombinationsConsidered;
            NumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

            var precursorErrors = psms.Select(p => (p.ScanPrecursorMass - p.PeptideMonisotopicMass.Value) / p.PeptideMonisotopicMass.Value * 1e6).ToList();
            PsmPrecursorIqrPpmError = Statistics.InterquartileRange(precursorErrors);
            PsmPrecursorMedianPpmError = Statistics.Median(precursorErrors);

            var productErrors = psms.Where(p => p.MatchedFragmentIons != null).SelectMany(p => p.MatchedFragmentIons).Select(p => (p.Mz.ToMass(p.Charge) - p.NeutralTheoreticalProduct.NeutralMass) / p.NeutralTheoreticalProduct.NeutralMass * 1e6).ToList();
            PsmProductIqrPpmError = Statistics.InterquartileRange(productErrors);
            PsmProductMedianPpmError = Statistics.Median(productErrors);
        }

        public Tuple<double, double> Ms1InfoTh { get; }
        public Tuple<double, double> Ms2InfoTh { get; }
        public Tuple<double, double> Ms1InfoPpm { get; }
        public Tuple<double, double> Ms2InfoPpm { get; }

        public int NumMs1MassChargeCombinationsConsidered { get; }
        public int NumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks { get; }
        public int NumMs2MassChargeCombinationsConsidered { get; }
        public int NumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks { get; }

        public List<LabeledDataPoint> Ms1List { get; }
        public List<LabeledDataPoint> Ms2List { get; }

        public readonly double PsmPrecursorMedianPpmError;
        public readonly double PsmProductMedianPpmError;
        public readonly double PsmPrecursorIqrPpmError;
        public readonly double PsmProductIqrPpmError;
        public readonly List<PeptideSpectralMatch> Psms;

        public int Count { get { return Ms1List.Count + Ms2List.Count; } }

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