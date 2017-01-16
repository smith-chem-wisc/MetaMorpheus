using InternalLogicEngineLayer;
using System.Collections.Generic;

namespace InternalLogicTaskLayer
{
    public abstract class MyTaskResults : MyResults
    {
        #region Public Fields

        public List<string> newSpectra;
        public List<string> newDatabases;

        #endregion Public Fields

        #region Public Constructors

        protected MyTaskResults(MyTaskEngine s) : base(s)
        {
        }

        #endregion Public Constructors
    }
}