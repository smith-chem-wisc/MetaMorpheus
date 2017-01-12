using System;
using InternalLogicEngineLayer;
using MassSpectrometry;
using Spectra;

namespace InternalLogicCalibration
{
    public class CalibrationResults : MyResults
    {
        public CalibrationResults(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, CalibrationEngine s) : base(s)
        {
            this.myMsDataFile = myMsDataFile;
        }

        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile { get; private set; }

        protected override string GetStringForOutput()
        {
            throw new NotImplementedException();
        }
    }
}