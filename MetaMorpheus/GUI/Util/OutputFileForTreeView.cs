namespace MetaMorpheusGUI
{
    public class OutputFileForTreeView : ForTreeView
    {
        public OutputFileForTreeView(string fullPath, string displayName) : base(displayName, fullPath)
        {
            FullPath = fullPath;
        }

        public string FullPath { get; }
    }
}