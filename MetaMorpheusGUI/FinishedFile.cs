namespace MetaMorpheusGUI
{
    internal class FinishedFile
    {
        public string filepath { get; private set; }

        public FinishedFile(string filepath)
        {
            this.filepath = filepath;
        }
    }
}