using System.Collections.Generic;

namespace EngineLayer.Calibration
{
    internal class DataPointAquisitionResults
    {

        #region Public Fields

        public int numMs1MassChargeCombinationsConsidered;
        public int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
        public int numMs2MassChargeCombinationsConsidered;
        public int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

        #endregion Public Fields

        #region Public Properties

        public List<LabeledMs1DataPoint> Ms1List { get; set; }
        public List<LabeledMs2DataPoint> Ms2List { get; set; }

        public int Count { get { return Ms1List.Count + Ms2List.Count; } }

        #endregion Public Properties

    }
}