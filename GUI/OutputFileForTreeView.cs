namespace MetaMorpheusGUI
{
    public class OutputFileForTreeView : ForTreeView
    {
        #region Public Constructors

        public OutputFileForTreeView(string fullPath, string displayName) : base(displayName)
        {
            FullPath = fullPath;
        }

        #endregion Public Constructors

        #region Public Properties

        public string FullPath { get; }

        #endregion Public Properties
    }
}