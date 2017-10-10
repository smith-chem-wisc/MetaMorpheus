using System.IO;

namespace MetaMorpheusGUI
{
    public class OutputFileForTreeView : ForTreeView
    {
        public readonly string fullPath;
        #region Public Constructors

        public OutputFileForTreeView(string fullPath, string displayName) : base(displayName)
        {
            this.fullPath = fullPath;
        }

        #endregion Public Constructors
    }
}