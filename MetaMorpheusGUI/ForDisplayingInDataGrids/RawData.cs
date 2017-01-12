namespace MetaMorpheusGUI
{
    public class RawData
    {
        public bool Use { get; set; }
        public string FileName { get; private set; }

        public RawData(string FileName)
        {
            this.FileName = FileName;
            if (FileName != null)
                Use = true;
        }

        public void AddFilePath(string FileName)
        {
            this.FileName = FileName;
            Use = true;
        }
    }
}