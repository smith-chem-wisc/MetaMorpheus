using InternalLogicEngineLayer;
using System.Collections.Generic;

namespace InternalLogicTaskLayer
{
    public abstract class MyTaskResults : MyResults
    {

        #region Public Fields

        public List<string> newSpectra;
        public List<DbForTask> newDatabases;

        #endregion Public Fields

        #region Protected Constructors

        protected MyTaskResults(MyTaskEngine s) : base(s)
        {
        }

        #endregion Protected Constructors

    }
}