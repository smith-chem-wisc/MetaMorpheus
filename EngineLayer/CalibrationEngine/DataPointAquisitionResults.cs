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
            List<LabeledMs2DataPoint> ms2List) : base(dataPointAcquisitionEngine)
        {
            Ms1List = ms1List;
            Ms2List = ms2List;

            Ms1InfoTh = Ms1List.Select(b => b.LabelTh).MeanStandardDeviation();
            Ms2InfoTh = Ms2List.Select(b => b.LabelTh).MeanStandardDeviation();

            Ms1InfoPpm = Ms1List.Select(b => b.LabelPpm).MeanStandardDeviation();
            Ms2InfoPpm = Ms2List.Select(b => b.LabelPpm).MeanStandardDeviation();
        }

        #endregion Public Constructors

        #region Public Properties

        public Tuple<double, double> Ms1InfoTh { get; }
        public Tuple<double, double> Ms2InfoTh { get; }
        public Tuple<double, double> Ms1InfoPpm { get; }
        public Tuple<double, double> Ms2InfoPpm { get; }
        public List<LabeledMs1DataPoint> Ms1List { get; }
        public List<LabeledMs2DataPoint> Ms2List { get; }

        public int Count { get { return Ms1List.Count + Ms2List.Count; } }

        #endregion Public Properties
    }
}