namespace MetaMorpheusGUI
{
    public class OutputFileForTreeView : ForTreeView
    {
        #region Public Fields

        public readonly string fullPath;

        #endregion Public Fields

        #region Public Constructors

        public OutputFileForTreeView(string fullPath, string displayName) : base(displayName)
        {
            this.fullPath = fullPath;
        }

        #endregion Public Constructors
    }
}