using InternalLogicEngineLayer;
using MassSpectrometry;
using Spectra;
using System.Text;

namespace InternalLogicCalibration
{
    public class CalibrationResults : MyResults
    {
        #region Public Constructors

        public CalibrationResults(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, CalibrationEngine s) : base(s)
        {
            this.myMsDataFile = myMsDataFile;
        }

        #endregion Public Constructors

        #region Public Properties

        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile { get; private set; }

        #endregion Public Properties

        #region Protected Methods

        protected override string GetStringForOutput()
        {
            StringBuilder sb = new StringBuilder();
            return sb.ToString();
        }

        #endregion Protected Methods
    }
}