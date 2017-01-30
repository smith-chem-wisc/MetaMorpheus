namespace MetaMorpheusGUI
{
    public class FinishedFile
    {
        #region Public Constructors

        public FinishedFile(string filepath)
        {
            this.filepath = filepath;
        }

        #endregion Public Constructors

        #region Public Properties

        public string filepath { get; private set; }

        #endregion Public Properties
    }
}