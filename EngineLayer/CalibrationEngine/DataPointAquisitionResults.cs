using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Calibration
{
    public class DataPointAquisitionResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public DataPointAquisitionResults(
            MetaMorpheusEngine dataPointAcquisitionEngine,
            List<LabeledMs1DataPoint> ms1List,
            List<LabeledMs2DataPoint> ms2List,
            int numMs1MassChargeCombinationsConsidered,
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
            int numMs2MassChargeCombinationsConsidered,
            int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks) : base(dataPointAcquisitionEngine)
        {
            Ms1List = ms1List;
            Ms2List = ms2List;

            NumMs1MassChargeCombinationsConsidered = numMs1MassChargeCombinationsConsidered;
            NumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
            NumMs2MassChargeCombinationsConsidered = numMs2MassChargeCombinationsConsidered;
            NumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

            Ms1InfoTh = Ms1List.Select(b => b.LabelTh).MeanStandardDeviation();
            Ms2InfoTh = Ms2List.Select(b => b.LabelTh).MeanStandardDeviation();

            Ms1InfoPPM = Ms1List.Select(b => b.LabelPPM).MeanStandardDeviation();
            Ms2InfoPPM = Ms2List.Select(b => b.LabelPPM).MeanStandardDeviation();
        }

        #endregion Public Constructors

        #region Public Properties

        public int NumMs1MassChargeCombinationsConsidered { get; }
        public int NumMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks { get; }
        public int NumMs2MassChargeCombinationsConsidered { get; }
        public int NumMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks { get; }
        public Tuple<double, double> Ms1InfoTh { get; }
        public Tuple<double, double> Ms2InfoTh { get; }
        public Tuple<double, double> Ms1InfoPPM { get; }
        public Tuple<double, double> Ms2InfoPPM { get; }
        public List<LabeledMs1DataPoint> Ms1List { get; }
        public List<LabeledMs2DataPoint> Ms2List { get; }

        public int Count { get { return Ms1List.Count + Ms2List.Count; } }

        #endregion Public Properties
    }
}