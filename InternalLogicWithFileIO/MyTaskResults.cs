using InternalLogic;
using System.Collections.Generic;

namespace InternalLogicWithFileIO
{
    public class MyTaskResults : MyResults
    {
        public List<string> newSpectra;
        public List<string> newDatabases;

        public MyTaskResults(MyTaskEngine s) : base(s)
        {
        }
    }
}