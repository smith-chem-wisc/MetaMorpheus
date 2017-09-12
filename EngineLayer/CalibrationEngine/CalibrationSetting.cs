using SharpLearning.Common.Interfaces;
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

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return " continueLoopString= " + ContinueLoop.Item2 + " DoFirst " + DoFirst == null ? "null" : string.Join(",", DoFirst) + " Learners " + string.Join(",", Learners);
        }

        #endregion Public Methods
    }
}