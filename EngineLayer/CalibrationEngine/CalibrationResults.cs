using MassSpectrometry;
using Spectra;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Calibration
{
    public class CalibrationResults : MyResults
    {

        #region Private Fields

        private List<int> numMs1MassChargeCombinationsConsideredList;
        private List<int> numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList;
        private List<int> numMs2MassChargeCombinationsConsideredList;
        private List<int> numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList;
        private List<int> countList;
        private List<SeparateCalibrationFunction> calibrationFunctions;

        #endregion Private Fields

        #region Public Constructors

        public CalibrationResults(IMsDataFile<IMzSpectrum<MzPeak>> myMSDataFile, CalibrationEngine s) : base(s)
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

        public IMsDataFile<IMzSpectrum<MzPeak>> MyMSDataFile { get; private set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string StringForOutput
        {
            get
            {
                var sb = new StringBuilder();
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
        }

        #endregion Protected Properties

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