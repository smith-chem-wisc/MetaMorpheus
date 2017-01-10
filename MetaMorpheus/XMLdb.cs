namespace MetaMorpheus
{
    public class XMLdb
    {
        public bool Use { get; set; }
        public string FileName { get; private set; }

        public XMLdb(string FileName)
        {
            Use = true;
            this.FileName = FileName;
        }
    }
}