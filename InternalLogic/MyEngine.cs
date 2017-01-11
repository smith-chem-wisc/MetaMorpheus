using System.Diagnostics;

namespace InternalLogic
{
    public abstract class MyEngine
    {
        protected MyParams myParams;

        public MyResults Run()
        {
            myParams.Validate();
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