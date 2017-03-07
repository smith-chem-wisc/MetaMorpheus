using System.IO;

namespace MetaMorpheusGUI
{
    public class OutputFileForTreeView : ForTreeView
    {

        #region Public Constructors

        public OutputFileForTreeView(string fullPath) : base(fullPath, Path.GetFileName(fullPath))
        {
        }

        #endregion Public Constructors

    }
}