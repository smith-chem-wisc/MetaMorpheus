namespace MetaMorpheusGUI
{
    public class FinishedFileForDataGrid
    {
        public FinishedFileForDataGrid(string filePath)
        {
            FilePath = filePath;
        }

        public string FilePath { get; private set; }
    }
}