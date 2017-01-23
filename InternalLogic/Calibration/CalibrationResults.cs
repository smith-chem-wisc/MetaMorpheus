using InternalLogicEngineLayer;
using MassSpectrometry;
using Spectra;
using System.Collections.Generic;
using System.Text;

namespace InternalLogicCalibration
{
    public class CalibrationResults : MyResults
    {

        #region Private Fields

        private List<int> calibrationRoundList;
        private List<int> numMs1MassChargeCombinationsConsideredList;
        private List<int> numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList;
        private List<int> numMs2MassChargeCombinationsConsideredList;
        private List<int> numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList;
        private List<int> countList;

        #endregion Private Fields

        #region Public Constructors

        public CalibrationResults(IMsDataFile<IMzSpectrum<MzPeak>> myMSDataFile, CalibrationEngine s) : base(s)
        {
            this.MyMSDataFile = myMSDataFile;
            calibrationRoundList = new List<int>();
            numMs1MassChargeCombinationsConsideredList = new List<int>();
            numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList = new List<int>();
            numMs2MassChargeCombinationsConsideredList = new List<int>();
            numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList = new List<int>();
            countList = new List<int>();
        }

        #endregion Public Constructors

        #region Public Properties

        public IMsDataFile<IMzSpectrum<MzPeak>> MyMSDataFile { get; private set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string StringForOutput
        {
            get
            {
                var sb = new StringBuilder();
                for (int i = 0; i < calibrationRoundList.Count; i++)
                {
                    sb.AppendLine("\t\tRound " + calibrationRoundList[i]);
                    sb.AppendLine("\t\t\tTraining points: " + countList[i]);
                    sb.AppendLine("\t\t\tMs2MassChargeSeen: " + numMs2MassChargeCombinationsConsideredList[i]);
                    sb.AppendLine("\t\t\tMs2MassChargeSeenAndIgnoredBecause too many: " + numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList[i]);

                    sb.AppendLine("\t\t\tMs1MassChargeSeen: " + numMs1MassChargeCombinationsConsideredList[i]);
                    sb.AppendLine("\t\t\tMs1MassChargeSeenAndIgnoredBecause too many: " + numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList[i]);
                }
                return sb.ToString();
            }
        }

        #endregion Protected Properties

        #region Internal Methods

        internal void Add(int calibrationRound, int numMs1MassChargeCombinationsConsidered, int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks, int count, int numMs2MassChargeCombinationsConsidered, int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks)
        {
            calibrationRoundList.Add(calibrationRound);
            numMs1MassChargeCombinationsConsideredList.Add(numMs1MassChargeCombinationsConsidered);
            numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList.Add(numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            numMs2MassChargeCombinationsConsideredList.Add(numMs2MassChargeCombinationsConsidered);
            numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList.Add(numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
            countList.Add(count);
        }

        #endregion Internal Methods

    }
}