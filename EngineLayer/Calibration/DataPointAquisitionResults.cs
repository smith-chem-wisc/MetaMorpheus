using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Calibration
{
    public class DataPointAquisitionResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public DataPointAquisitionResults(
            MetaMorpheusEngine dataPointAcquisitionEngine,
            List<LabeledDataPoint> ms1List,
            List<LabeledDataPoint> ms2List,
            int numMs1MassChargeCombinationsConsidered,
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
            int numMs2MassChargeCombinationsConsidered,
            int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks)
            : base(dataPointAcquisitionEngine)
        {
            Ms1List = ms1List;
            Ms2List = ms2List;

            Ms1InfoTh = Ms1List.Select(b => b.experimentalMz - b.theoreticalMz).MeanStandardDeviation();
            Ms2InfoTh = Ms2List.Select(b => b.experimentalMz - b.theoreticalMz).MeanStandardDeviation();

            Ms1InfoPpm = Ms1List.Select(b => (b.experimentalMz - b.theoreticalMz) / b.theoreticalMz).MeanStandardDeviation();
            Ms2InfoPpm = Ms2List.Select(b => (b.experimentalMz - b.theoreticalMz) / b.theoreticalMz).MeanStandardDeviation();

            NumMs1MassChargeCombinationsConsidered = numMs1MassChargeCombinationsConsidered;
            NumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
            NumMs2MassChargeCombinationsConsidered = numMs2MassChargeCombinationsConsidered;
            NumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
        }

        #endregion Public Constructors

        #region Public Properties

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

        public int Count { get { return Ms1List.Count + Ms2List.Count; } }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("Ms1List.Count: " + Ms1List.Count + " Ms2List.Count: " + Ms2List.Count);
            sb.AppendLine("Ms1ppm mean: " + Ms1InfoPpm.Item1 + " Ms1ppm sd: " + Ms1InfoPpm.Item2);
            sb.AppendLine("Ms1th mean: " + Ms1InfoTh.Item1 + " Ms1th sd: " + Ms1InfoTh.Item2);
            sb.AppendLine("Ms2ppm mean: " + Ms2InfoPpm.Item1 + " Ms2ppm sd: " + Ms2InfoPpm.Item2);
            sb.AppendLine("Ms2th mean: " + Ms2InfoTh.Item1 + " Ms2th sd: " + Ms2InfoTh.Item2);
            sb.AppendLine("NumMs1MassChargeCombinationsConsidered: " + NumMs1MassChargeCombinationsConsidered);
            sb.AppendLine("NumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks: " + NumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            sb.AppendLine("NumMs2MassChargeCombinationsConsidered: " + NumMs2MassChargeCombinationsConsidered);
            sb.Append("NumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks: " + NumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}