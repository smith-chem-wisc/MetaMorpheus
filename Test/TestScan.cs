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

        public string id
        {
            get
            {
                throw new NotImplementedException();
            }
        }

        public double InjectionTime { get; private set; }

        public bool isCentroid
        {
            get
            {
                throw new NotImplementedException();
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
                throw new NotImplementedException();
            }
        }

        public double RetentionTime { get; private set; }

        public string ScanFilter
        {
            get
            {
                throw new NotImplementedException();
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

        public void tranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(Func<MzPeak, double> convertorForSpectrum, double selectedIonGuessMZ, double newSelectedIonGuessMonoisotopicMZ)
        {
            MassSpectrum.replaceXbyApplyingFunction(convertorForSpectrum);
            selectedIonGuessMonoisotopicMZ = newSelectedIonGuessMonoisotopicMZ;
        }

        public bool TryGetDissociationType(out DissociationType DissociationType)
        {
            throw new NotImplementedException();
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
            throw new NotImplementedException();
        }

        public bool TryGetIsolationWidth(out double IsolationWidth)
        {
            throw new NotImplementedException();
        }

        public bool TryGetPrecursorID(out string PrecursorID)
        {
            throw new NotImplementedException();
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