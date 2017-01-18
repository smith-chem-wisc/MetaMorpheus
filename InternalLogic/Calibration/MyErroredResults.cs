using InternalLogicEngineLayer;

namespace InternalLogicCalibration
{
    internal class MyErroredResults : MyResults
    {

        #region Private Fields

        private string v;

        #endregion Private Fields

        #region Public Constructors

        public MyErroredResults(MyEngine s, string v) : base(s)
        {
            this.v = v;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override string StringForOutput
        {
            get
            {
                return v;
            }
        }

        #endregion Protected Methods

    }
}