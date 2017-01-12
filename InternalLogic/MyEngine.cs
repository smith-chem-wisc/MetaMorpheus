using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace InternalLogicEngineLayer
{
    public abstract class MyEngine
    {
        public static UsefulProteomicsDatabases.Generated.unimod unimodDeserialized;

        public static Dictionary<int, ChemicalFormulaModification> uniprotDeseralized;

        internal readonly int Level;

        public static event EventHandler<SingleEngineEventArgs> startingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> finishedSingleEngineHandler;

        public static event EventHandler<string> outLabelStatusHandler;

        public static event EventHandler<ProgressEventArgs> outProgressHandler;

        protected MyEngine(int Level)
        {
            this.Level = Level;
        }

        protected void status(string v)
        {
            outLabelStatusHandler?.Invoke(this, v);
        }

        private void startingSingleEngine()
        {
            startingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            outProgressHandler?.Invoke(this, v);
        }

        private void finishedSingleEngine(MyResults myResults)
        {
            finishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
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
            finishedSingleEngine(myResults);
            return myResults;
        }

        public abstract void ValidateParams();
        protected abstract MyResults RunSpecific();
    }
}