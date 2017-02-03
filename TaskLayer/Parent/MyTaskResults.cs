using EngineLayer;
using System.Collections.Generic;

namespace TaskLayer
{
    public abstract class MyTaskResults : MyResults
    {

        #region Public Fields

        public List<string> newSpectra;
        public List<DbForTask> newDatabases;

        #endregion Public Fields

        #region Protected Constructors

        protected MyTaskResults(MetaMorpheusTask s) : base(s)
        {
        }

        #endregion Protected Constructors

    }
}