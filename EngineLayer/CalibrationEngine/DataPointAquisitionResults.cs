using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Calibration
{
    public class DataPointAquisitionResults
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

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(" Ms1List.Count " + Ms1List.Count);
            sb.AppendLine(" Ms2List.Count " + Ms2List.Count);
            sb.AppendLine(" numMs1MassChargeCombinationsConsidered " + numMs1MassChargeCombinationsConsidered);
            sb.AppendLine(" numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks " + numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            sb.AppendLine(" numMs2MassChargeCombinationsConsidered " + numMs2MassChargeCombinationsConsidered);
            sb.Append(" numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks " + numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}