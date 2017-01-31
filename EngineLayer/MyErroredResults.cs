namespace EngineLayer
{
    public class MyErroredResults : MyResults
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

        #region Protected Properties

        protected override string StringForOutput
        {
            get
            {
                return "\t\t" + v;
            }
        }

        #endregion Protected Properties

    }
}