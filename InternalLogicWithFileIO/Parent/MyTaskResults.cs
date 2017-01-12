using InternalLogicEngineLayer;
using System.Collections.Generic;

namespace InternalLogicTaskLayer
{
    public abstract class MyTaskResults : MyResults
    {
        public List<string> newSpectra;
        public List<string> newDatabases;

        public MyTaskResults(MyTaskEngine s) : base(s)
        {
        }
    }
}