using MassSpectrometry;
using Spectra;
using System;

namespace Test
{
    internal class TestScan : IMsDataScan<DefaultMzSpectrum>
    {

        #region Private Fields

        private readonly double isolationMZ;
        private readonly int precursorOneBasedScanNumber;
        private double selectedIonGuessMonoisotopicIntensity;

        private int selectedIonGuessChargeStateGuess;

        private double selectedIonGuessMonoisotopicMZ;

        #endregion Private Fields

        #region Public Constructors

        public TestScan(int OneBasedScanNumber, double RetentionTime, DefaultMzSpectrum MassSpectrum, double selectedIonGuessMonoisotopicMZ, int selectedIonGuessChargeStateGuess, double selectedIonGuessMonoisotopicIntensity, double isolationMZ, double InjectionTime, int precursorOneBasedScanNumber)
        {
            MsnOrder = 2;
            this.OneBasedScanNumber = OneBasedScanNumber;
            this.RetentionTime = RetentionTime;
            this.MassSpectrum = MassSpectrum;

            this.selectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
            this.selectedIonGuessChargeStateGuess = selectedIonGuessChargeStateGuess;
            this.selectedIonGuessMonoisotopicIntensity = selectedIonGuessMonoisotopicIntensity;
            this.isolationMZ = isolationMZ;
            this.InjectionTime = InjectionTime;
            this.precursorOneBasedScanNumber = precursorOneBasedScanNumber;
        }

        public TestScan(int OneBasedScanNumber, double RetentionTime, DefaultMzSpectrum MassSpectrum, double InjectionTime)
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

        public double InjectionTime { get; private set; }

        public bool IsCentroid
        {
            get
            {
                return true;
            }
        }

        public DefaultMzSpectrum MassSpectrum { get; private set; }

        public int MsnOrder { get; private set; }

        public MZAnalyzerType MzAnalyzer
        {
            get
            {
                throw new NotImplementedException();
            }
        }

        public int OneBasedScanNumber { get; private set; }

        public Polarity Polarity
        {
            get
            {
                return Polarity.Positive;
            }
        }

        public double RetentionTime { get; private set; }

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

        #endregion Public Properties

        #region Public Methods

        public void TranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(Func<MzPeak, double> convertorForSpectrum, double selectedIonGuessMZ, double newSelectedIonGuessMonoisotopicMZ)
        {
            MassSpectrum.replaceXbyApplyingFunction(convertorForSpectrum);
            selectedIonGuessMonoisotopicMZ = newSelectedIonGuessMonoisotopicMZ;
        }

        public bool TryGetDissociationType(out DissociationType DissociationType)
        {
            if (MsnOrder == 2)
            {
                DissociationType = DissociationType.HCD;
                return true;
            }
            DissociationType = DissociationType.Unknown;
            return false;
        }

        public bool TryGetIsolationMZ(out double IsolationMZ)
        {
            if (MsnOrder == 2)
            {
                IsolationMZ = isolationMZ;
                return true;
            }
            IsolationMZ = double.NaN;
            return false;
        }

        public bool TryGetIsolationRange(out MzRange IsolationRange)
        {
            if (MsnOrder == 2)
            {
                IsolationRange = new MzRange(isolationMZ - 3, isolationMZ + 3);
                return true;
            }
            IsolationRange = null;
            return false;
        }

        public bool TryGetIsolationWidth(out double IsolationWidth)
        {
            if (MsnOrder == 2)
            {
                IsolationWidth = 6;
                return true;
            }
            IsolationWidth = double.NaN;
            return false;
        }

        public bool TryGetPrecursorID(out string PrecursorID)
        {
            if (MsnOrder == 2)
            {
                PrecursorID = 1.ToString();
                return true;
            }
            PrecursorID = null;
            return false;
        }

        public bool TryGetPrecursorOneBasedScanNumber(out int PrecursorOneBasedScanNumber)
        {
            if (MsnOrder == 2)
            {
                PrecursorOneBasedScanNumber = precursorOneBasedScanNumber;
                return true;
            }
            PrecursorOneBasedScanNumber = 0;
            return false;
        }

        public bool TryGetSelectedIonGuessChargeStateGuess(out int SelectedIonGuessChargeStateGuess)
        {
            if (MsnOrder == 2)
            {
                SelectedIonGuessChargeStateGuess = selectedIonGuessChargeStateGuess;
                return true;
            }
            SelectedIonGuessChargeStateGuess = 0;
            return false;
        }

        public bool TryGetSelectedIonGuessIntensity(out double SelectedIonGuessIntensity)
        {
            if (MsnOrder == 2)
            {
                SelectedIonGuessIntensity = selectedIonGuessMonoisotopicIntensity;
                return true;
            }
            SelectedIonGuessIntensity = double.NaN;
            return false;
        }

        public bool TryGetSelectedIonGuessMonoisotopicIntensity(out double SelectedIonGuessMonoisotopicIntensity)
        {
            if (MsnOrder == 2)
            {
                SelectedIonGuessMonoisotopicIntensity = selectedIonGuessMonoisotopicIntensity;
                return true;
            }
            SelectedIonGuessMonoisotopicIntensity = double.NaN;
            return false;
        }

        public bool TryGetSelectedIonGuessMonoisotopicMZ(out double SelectedIonGuessMonoisotopicMZ)
        {
            if (MsnOrder == 2)
            {
                SelectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
                return true;
            }
            SelectedIonGuessMonoisotopicMZ = double.NaN;
            return false;
        }

        public bool TryGetSelectedIonGuessMZ(out double SelectedIonGuessMZ)
        {
            if (MsnOrder == 2)
            {
                SelectedIonGuessMZ = selectedIonGuessMonoisotopicMZ;
                return true;
            }
            SelectedIonGuessMZ = double.NaN;
            return false;
        }

        #endregion Public Methods

    }
}