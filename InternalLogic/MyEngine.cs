using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace InternalLogicEngineLayer
{
    public abstract class MyEngine
    {
        protected static UsefulProteomicsDatabases.Generated.unimod unimodDeserialized;

        protected static Dictionary<int, ChemicalFormulaModification> uniprotDeseralized;

        public static event EventHandler<SingleEngineEventArgs> startingSingleEngineHander;

        public static event EventHandler<SingleEngineEventArgs> finishedSingleEngineHandler;

        public static event EventHandler<string> outLabelStatusHandler;

        public static event EventHandler<ProgressEventArgs> outProgressHandler;

        public static event EventHandler<string> outRichTextBoxHandler;

        static MyEngine()
        {
            unimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(@"unimod_tables.xml");
            uniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(@"ptmlist.txt");
        }

        protected void status(string v)
        {
            outLabelStatusHandler?.Invoke(this, v);
        }

        private void startingSingleEngine()
        {
            startingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void finishedSingleEngine()
        {
            finishedSingleEngineHandler?.Invoke(this, new SingleEngineEventArgs(this));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            outProgressHandler?.Invoke(this, v);
        }

        protected void output(string v)
        {
            outRichTextBoxHandler?.Invoke(this, v);
        }

        public MyResults Run()
        {
            startingSingleEngine();
            ValidateParams();
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            var myResults = RunSpecific();
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            finishedSingleEngine();
            return myResults;
        }

        public abstract void ValidateParams();
        protected abstract MyResults RunSpecific();
    }
}