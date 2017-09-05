using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using Spectra;
using System;

namespace Test
{
    internal class TestScan : IMsDataScan<IMzSpectrum<IMzPeak>>
    {
        #region Public Constructors

        public TestScan(int OneBasedScanNumber, double RetentionTime, MzmlMzSpectrum MassSpectrum, double InjectionTime)
        {
            MsnOrder = 1;
            this.OneBasedScanNumber = OneBasedScanNumber;
            this.RetentionTime = RetentionTime;
            this.MassSpectrum = MassSpectrum;
            this.InjectionTime = InjectionTime;
        }

        #endregion Public Constructors

        #region Public Properties

        public string Id
        {
            get
            {
                return OneBasedScanNumber.ToString();
            }
        }

        public double? InjectionTime { get; protected set; }

        public bool IsCentroid
        {
            get
            {
                return true;
            }
        }

        public MzmlMzSpectrum MassSpectrum { get; protected set; }

        public int MsnOrder { get; protected set; }

        public MZAnalyzerType MzAnalyzer
        {
            get
            {
                throw new NotImplementedException();
            }
        }

        public double[,] NoiseData
        {
            get
            {
                throw new NotImplementedException();
            }
        }

        public int OneBasedScanNumber { get; protected set; }

        public Polarity Polarity
        {
            get
            {
                return Polarity.Positive;
            }
        }

        public double RetentionTime { get; protected set; }

        public string ScanFilter
        {
            get
            {
                return "FTMS" + Id;
            }
        }

        public MzRange ScanWindowRange
        {
            get
            {
                return new MzRange(0, 10000);
            }
        }

        public double TotalIonCurrent
        {
            get
            {
                return MassSpectrum.SumOfAllY;
            }
        }

        IMzSpectrum<IMzPeak> IMsDataScan<IMzSpectrum<IMzPeak>>.MassSpectrum
        {
            get
            {
                return MassSpectrum;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public byte[] Get64BitNoiseDataBaseline()
        {
            throw new NotImplementedException();
        }

        public byte[] Get64BitNoiseDataMass()
        {
            throw new NotImplementedException();
        }

        public byte[] Get64BitNoiseDataNoise()
        {
            throw new NotImplementedException();
        }

        public void TransformByApplyingFunctionToSpectra(Func<IPeak, double> convertorForSpectrum)
        {
            MassSpectrum.ReplaceXbyApplyingFunction(convertorForSpectrum);
        }

        #endregion Public Methods
    }
}