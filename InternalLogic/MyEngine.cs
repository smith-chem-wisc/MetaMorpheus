using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Net;
using System.Reflection;

namespace InternalLogicEngineLayer
{
    public abstract class MyEngine
    {
        #region Public Fields

        public static readonly string MetaMorpheusVersion;
        public static UsefulProteomicsDatabases.Generated.unimod unimodDeserialized;
        public static Dictionary<int, ChemicalFormulaModification> uniprotDeseralized;

        #endregion Public Fields

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

        public static event EventHandler<SingleEngineEventArgs> startingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> finishedSingleEngineHandler;

        public static event EventHandler<string> outLabelStatusHandler;

        public static event EventHandler<ProgressEventArgs> outProgressHandler;

        #endregion Public Events

        #region Public Methods

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

        protected void status(string v)
        {
            outLabelStatusHandler?.Invoke(this, v);
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            outProgressHandler?.Invoke(this, v);
        }

        protected abstract MyResults RunSpecific();

        #endregion Protected Methods

        #region Private Methods

        private void startingSingleEngine()
        {
            startingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void finishedSingleEngine(MyResults myResults)
        {
            finishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }

        #endregion Private Methods
    }
}