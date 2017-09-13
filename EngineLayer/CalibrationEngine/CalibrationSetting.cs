using SharpLearning.Common.Interfaces;
using SharpLearning.Metrics.Regression;
using System;
using System.Collections.Generic;

namespace EngineLayer.Calibration
{
    public class CalibrationSetting
    {
        #region Public Properties

        public List<ILearner<double>> Learners { get; set; }

        public (Func<DataPointAquisitionResults, DataPointAquisitionResults, bool>, string) ContinueLoop { get; set; }

        public List<ILearner<double>> DoFirst { get; set; }
        public IRegressionMetric Metric { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return " continueLoopString= " + ContinueLoop.Item2 + " DoFirst " + DoFirst == null ? "null" : string.Join(",", DoFirst) + " Learners " + string.Join(",", Learners) + " Metric " + Metric;
        }

        #endregion Public Methods
    }
}