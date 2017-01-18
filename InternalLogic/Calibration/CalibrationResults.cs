using InternalLogicEngineLayer;
using MassSpectrometry;
using Spectra;
using System.Text;

namespace InternalLogicCalibration
{
    public class CalibrationResults : MyResults
    {

        #region Public Constructors

        public CalibrationResults(IMsDataFile<IMzSpectrum<MzPeak>> myMSDataFile, CalibrationEngine s) : base(s)
        {
            this.MyMSDataFile = myMSDataFile;
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
                return sb.ToString();
            }
        }

        #endregion Protected Properties

    }
}