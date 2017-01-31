using MassSpectrometry;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net;
using System.Reflection;

namespace EngineLayer
{
    public abstract class MyEngine
    {

        #region Internal Fields

        internal readonly int Level;

        #endregion Internal Fields

        #region Private Fields

        private static readonly string elementsLocation = Path.Combine("Data", @"elements.dat");
        private static readonly string unimodLocation = Path.Combine("Data", @"unimod_tables.xml");
        private static readonly string uniprotLocation = Path.Combine("Data", @"ptmlist.txt");

        private static readonly double[] mms = new double[] { 1.0025, 2.005, 3.0075, 4.010 };

        #endregion Private Fields

        #region Public Constructors

        static MyEngine()
        {
            try
            {
                UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
                unimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation);
                uniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation);
            }
            catch (WebException)
            {
            }

            MetaMorpheusVersion = Assembly.GetExecutingAssembly().GetName().Version.ToString();
        }

        #endregion Public Constructors

        #region Protected Constructors

        protected MyEngine(int Level)
        {
            this.Level = Level;
        }

        #endregion Protected Constructors

        #region Public Events

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        #endregion Public Events

        #region Public Properties

        public static string MetaMorpheusVersion { get; private set; }
        public static UsefulProteomicsDatabases.Generated.unimod unimodDeserialized { get; private set; }
        public static Dictionary<int, ChemicalFormulaModification> uniprotDeseralized { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public static IEnumerable<LocalMS2Scan> GetMs2Scans(IMsDataFile<IMzSpectrum<MzPeak>> myMSDataFile)
        {
            foreach (var heh in myMSDataFile)
            {
                if (heh.MsnOrder == 2)
                {
                    int kk;
                    heh.TryGetPrecursorOneBasedScanNumber(out kk);
                    var uu = myMSDataFile.GetOneBasedScan(kk);
                    double isolationMz;
                    heh.TryGetIsolationMZ(out isolationMz);
                    int? monoisotopicPrecursorChargehere;
                    int mc;
                    if (heh.TryGetSelectedIonGuessChargeStateGuess(out monoisotopicPrecursorChargehere) && monoisotopicPrecursorChargehere.HasValue)
                        mc = monoisotopicPrecursorChargehere.Value;
                    else
                        mc = GuessCharge(uu.MassSpectrum.NewSpectrumExtract(isolationMz - 2.1, isolationMz + 2.1));
                    yield return new LocalMS2Scan(heh, mc);
                }
            }
        }

        public MyResults Run()
        {
            startingSingleEngine();
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            var myResults = RunSpecific();
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            finishedSingleEngine(myResults);
            return myResults;
        }

        #endregion Public Methods

        #region Protected Methods

        protected void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v));
        }

        protected void Status(string v)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MyResults RunSpecific();

        #endregion Protected Methods

        #region Private Methods

        private static int GuessCharge(IMzSpectrum<MzPeak> mzSpectrum)
        {
            double tolHere = 0.01;
            int[] chargeCount = new int[4]; // charges 1,2,3,4
            for (int i = 0; i < mzSpectrum.Count; i++)
                for (int j = i + 1; j < mzSpectrum.Count; j++)
                {
                    for (int charge = 1; charge <= 4; charge++)
                    {
                        for (int isotope = 0; isotope < 4; isotope++)
                        {
                            if (Math.Abs(mzSpectrum.XArray[j] - mzSpectrum.XArray[i] - mms[isotope] / charge) < tolHere)
                            {
                                chargeCount[charge - 1]++;
                            }
                        }
                    }
                }
            return Array.IndexOf(chargeCount, chargeCount.Max()) + 1;
        }

        private void startingSingleEngine()
        {
            StartingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void finishedSingleEngine(MyResults myResults)
        {
            FinishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }

        #endregion Private Methods

    }
}