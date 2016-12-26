using MassSpectrometry;
using Spectra;
using System;
using System.Collections.Generic;

namespace mzCal
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

        public event EventHandler<string> outputHandler;

        public event EventHandler<int> progressHandler;

        public event EventHandler<string> watchHandler;

        public event EventHandler<string> finishedFileHandler;

        public HashSet<int> MS2spectraToWatch;
        public HashSet<int> MS1spectraToWatch;
        public DoubleRange mzRange;

        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
        public Identifications identifications;

        public delegate void PostProcessing(SoftwareLockMassParams p);

        public PostProcessing postProcessing;

        public delegate string GetFormulaFromDictionary(string dictionary, string acession);

        public GetFormulaFromDictionary getFormulaFromDictionary;
        public bool calibrateSpectra = true;
        internal int randomSeed;
        public string paramString = "";
        public int minMS2 = 2;
        public int minMS1 = 3;
        public double toleranceInMZforMS2Search;
        internal double toleranceInMZforMS1Search = 0.01;
        public HashSet<int> matchesToExclude;

        #region Constructors

        public SoftwareLockMassParams(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, int randomSeed, bool deconvolute, double toleranceInMZforMS2Search)
        {
            this.myMsDataFile = myMsDataFile;
            MS1spectraToWatch = new HashSet<int>();
            MS2spectraToWatch = new HashSet<int>();
            this.randomSeed = randomSeed;
            this.toleranceInMZforMS2Search = toleranceInMZforMS2Search;
        }

        #endregion Constructors

        public virtual void OnFileFinished(string e)
        {
            finishedFileHandler?.Invoke(this, e);
        }

        public virtual void OnOutput(string e)
        {
            outputHandler?.Invoke(this, e);
        }

        public virtual void OnProgress(int e)
        {
            progressHandler?.Invoke(this, e);
        }

        public virtual void OnWatch(string e)
        {
            watchHandler?.Invoke(this, e);
        }
    }
}