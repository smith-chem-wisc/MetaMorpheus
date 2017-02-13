using System;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;

namespace Test
{
    internal class TestScanWithPrecursor : TestScan, IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>
    {

        #region Private Fields

        private readonly double isolationMZ;
        private readonly int precursorOneBasedScanNumber;
        private double selectedIonGuessMonoisotopicIntensity;

        private int? selectedIonGuessChargeStateGuess;

        private double selectedIonGuessMonoisotopicMZ;

        #endregion Private Fields
        public TestScanWithPrecursor(int OneBasedScanNumber, double RetentionTime, MzmlMzSpectrum MassSpectrum, double selectedIonGuessMonoisotopicMZ, int selectedIonGuessChargeStateGuess, double selectedIonGuessMonoisotopicIntensity, double isolationMZ, double InjectionTime, int precursorOneBasedScanNumber) : base(OneBasedScanNumber, RetentionTime, MassSpectrum, InjectionTime)
        {
            MsnOrder = 2;
            this.selectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
            this.selectedIonGuessChargeStateGuess = selectedIonGuessChargeStateGuess;
            this.selectedIonGuessMonoisotopicIntensity = selectedIonGuessMonoisotopicIntensity;
            this.isolationMZ = isolationMZ;
            this.precursorOneBasedScanNumber = precursorOneBasedScanNumber;

        }

        public void TranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(Func<IMzPeak, double> convertorForSpectrum, double selectedIonGuessMZ, double newSelectedIonGuessMonoisotopicMZ)
        {
            MassSpectrum.ReplaceXbyApplyingFunction(convertorForSpectrum);
            selectedIonGuessMonoisotopicMZ = newSelectedIonGuessMonoisotopicMZ;
        }

        public DissociationType DissociationType
        {
            get
            {
                return DissociationType.HCD;
            }
        }

        public double IsolationMz
        {
            get
            {
                return isolationMZ;
            }
        }

        public MzRange IsolationRange
        {
            get
            {
                return new MzRange(isolationMZ - IsolationWidth / 2, isolationMZ + IsolationWidth / 2);
            }
        }

        public double IsolationWidth
        {
            get
            {
                return 2;
            }
        }

        public int OneBasedPrecursorScanNumber
        {
            get
            {
                return precursorOneBasedScanNumber;
            }
        }

        public string PrecursorID
        {
            get
            {
                return precursorOneBasedScanNumber.ToString();
            }
        }

        public int? SelectedIonGuessChargeStateGuess
        {
            get
            {
                return selectedIonGuessChargeStateGuess;
            }
        }

        public double SelectedIonGuessIntensity
        {
            get
            {
                return selectedIonGuessMonoisotopicIntensity;
            }
        }

        public double SelectedIonGuessMonoisotopicIntensity
        {
            get
            {
                return selectedIonGuessMonoisotopicIntensity;
            }
        }

        public double SelectedIonGuessMonoisotopicMZ
        {
            get
            {
                return selectedIonGuessMonoisotopicMZ;
            }
        }

        public double SelectedIonGuessMZ
        {
            get
            {
                return selectedIonGuessMonoisotopicMZ;
            }
        }
    }
}