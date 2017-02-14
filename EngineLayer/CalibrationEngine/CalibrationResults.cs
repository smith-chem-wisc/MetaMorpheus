using MassSpectrometry;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Calibration
{
    public class CalibrationResults : MyResults
    {

        #region Private Fields

        private readonly List<int> numMs1MassChargeCombinationsConsideredList;
        private readonly List<int> numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList;
        private readonly List<int> numMs2MassChargeCombinationsConsideredList;
        private readonly List<int> numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList;
        private readonly List<int> countList;
        private readonly List<SeparateCalibrationFunction> calibrationFunctions;

        #endregion Private Fields

        #region Public Constructors

        public CalibrationResults(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile, CalibrationEngine s) : base(s)
        {
            this.MyMSDataFile = myMSDataFile;
            calibrationFunctions = new List<SeparateCalibrationFunction>();
            numMs1MassChargeCombinationsConsideredList = new List<int>();
            numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList = new List<int>();
            numMs2MassChargeCombinationsConsideredList = new List<int>();
            numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList = new List<int>();
            countList = new List<int>();
        }

        #endregion Public Constructors

        #region Public Properties

        public IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> MyMSDataFile { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            for (int i = 0; i < countList.Count; i++)
            {
                sb.AppendLine("\t\tRound " + (i + 1));
                sb.AppendLine("\t\t\tTraining points: " + countList[i]);
                sb.AppendLine("\t\t\tMs1MassChargeSeen: " + numMs1MassChargeCombinationsConsideredList[i]);
                sb.AppendLine("\t\t\tMs1MassChargeSeenAndIgnoredBecause too many: " + numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList[i]);
                sb.AppendLine("\t\t\tMs2MassChargeSeen: " + numMs2MassChargeCombinationsConsideredList[i]);
                sb.AppendLine("\t\t\tMs2MassChargeSeenAndIgnoredBecause too many: " + numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList[i]);

                if (i < calibrationFunctions.Count)
                {
                    sb.AppendLine("\t\t\tMs1Calibration function: " + calibrationFunctions[i].CalibrationFunction1.ToString());
                    sb.AppendLine("\t\t\tMs2Calibration function: " + calibrationFunctions[i].CalibrationFunction2.ToString());
                }
            }
            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal void Add(int numMs1MassChargeCombinationsConsidered, int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks, int count, int numMs2MassChargeCombinationsConsidered, int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks)
        {
            numMs1MassChargeCombinationsConsideredList.Add(numMs1MassChargeCombinationsConsidered);
            numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList.Add(numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            numMs2MassChargeCombinationsConsideredList.Add(numMs2MassChargeCombinationsConsidered);
            numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList.Add(numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            countList.Add(count);
        }

        internal void Add(CalibrationFunction combinedCalibration)
        {
            calibrationFunctions.Add((SeparateCalibrationFunction)combinedCalibration);
        }

        #endregion Internal Methods

    }
}