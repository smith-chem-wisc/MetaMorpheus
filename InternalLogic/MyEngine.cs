using MassSpectrometry;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Net;
using System.Reflection;

namespace InternalLogicEngineLayer
{
    public abstract class MyEngine
    {

        #region Internal Fields

        internal readonly int Level;

        #endregion Internal Fields

        #region Private Fields

        private const string elementsLocation = @"elements.dat";
        private const string unimodLocation = @"unimod_tables.xml";
        private const string uniprotLocation = @"ptmlist.txt";

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
                int chargeState = 0;
                if (heh.TryGetSelectedIonGuessChargeStateGuess(out chargeState) && chargeState > 0)
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