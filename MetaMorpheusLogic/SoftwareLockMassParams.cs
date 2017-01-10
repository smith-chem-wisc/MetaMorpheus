using MassSpectrometry;
using Spectra;
using System.Collections.Generic;

namespace MetaMorpheusLogic
{
    public class SoftwareLockMassParams
    {
        #region isotopologue parameters

        // THIS PARAMETER IS FRAGILE!!!
        // TUNED TO CORRESPOND TO SPECTROMETER OUTPUT
        // BETTER SPECTROMETERS WOULD HAVE BETTER (LOWER) RESOLUIONS
        // Parameter for isotopolouge distribution searching
        public double fineResolution = 0.1;

        #endregion isotopologue parameters

        public HashSet<int> MS2spectraToWatch;
        public HashSet<int> MS1spectraToWatch;
        public DoubleRange mzRange;

        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
        public List<NewPsmWithFDR> identifications;

        public delegate void PostProcessing(SoftwareLockMassParams p);

        public PostProcessing postProcessing;

        public delegate string GetFormulaFromDictionary(string dictionary, string acession);

        public bool calibrateSpectra = true;
        internal int randomSeed;
        public string paramString = "";
        public int minMS2 = 2;
        public int minMS1 = 3;
        public double toleranceInMZforMS2Search;
        internal double toleranceInMZforMS1Search = 0.01;
        public HashSet<int> matchesToExclude;
        internal AllTasksParams po;
        internal MyTaskResults myTaskResults;

        #region Constructors

        public SoftwareLockMassParams(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, int randomSeed, double toleranceInMZforMS2Search, MyTaskResults myTaskResults)
        {
            this.myMsDataFile = myMsDataFile;
            MS1spectraToWatch = new HashSet<int>();
            MS2spectraToWatch = new HashSet<int>();
            this.randomSeed = randomSeed;
            this.toleranceInMZforMS2Search = toleranceInMZforMS2Search;
            this.myTaskResults = myTaskResults;
        }

        #endregion Constructors
    }
}