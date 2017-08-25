﻿using MzLibUtil;

namespace EngineLayer
{
    public class CalibrationParameters
    {
        #region Public Constructors

        public CalibrationParameters()
        {
            PrecursorMassTolerance = new PpmTolerance(10);
            NonLinearCalibration = true;
            WriteIntermediateFiles = false;
        }

        #endregion Public Constructors

        #region Public Properties

        public Tolerance PrecursorMassTolerance { get; set; }
        public bool NonLinearCalibration { get; set; }
        public bool WriteIntermediateFiles { get; set; }

        #endregion Public Properties
    }
}