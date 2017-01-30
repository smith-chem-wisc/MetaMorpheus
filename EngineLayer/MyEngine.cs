using MassSpectrometry;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
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
                int? chargeState;
                if (heh.TryGetSelectedIonGuessChargeStateGuess(out chargeState) && chargeState.HasValue)
                {
                    yield return new LocalMS2Scan(heh);
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