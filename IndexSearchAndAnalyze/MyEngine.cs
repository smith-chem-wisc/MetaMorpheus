using System;
using System.Diagnostics;

namespace IndexSearchAndAnalyze
{
    public abstract class MyEngine
    {
        public MyResults Run()
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            var myResults = RunSpecific();
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            return myResults;
        }

        protected abstract MyResults RunSpecific();
    }
}