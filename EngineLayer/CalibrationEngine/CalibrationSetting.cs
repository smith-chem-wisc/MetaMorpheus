using SharpLearning.Common.Interfaces;
using System;

namespace EngineLayer.Calibration
{
    public class CalibrationSetting
    {
        #region Public Fields

        public ILearner<double> learner = new IdentityCalibrationFunction();

        #endregion Public Fields

        #region Public Properties

        public (Func<DataPointAquisitionResults, DataPointAquisitionResults, bool>, string) ContinueLoop { get; set; } = ((DataPointAquisitionResults a, DataPointAquisitionResults b) => false, "");

        public bool DoCalibration { get; set; } = true;

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return " learner= " + learner.ToString() + " continueLoopString= " + ContinueLoop.Item2;
        }

        #endregion Public Methods
    }
}