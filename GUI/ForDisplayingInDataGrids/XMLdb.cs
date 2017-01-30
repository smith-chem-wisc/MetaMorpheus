namespace MetaMorpheusGUI
{
    public class XMLdb
    {

        #region Public Constructors

        public XMLdb(string FileName)
        {
            Use = true;
            this.FileName = FileName;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Use { get; set; }
        public bool Contaminant { get; set; }
        public string FileName { get; private set; }

        #endregion Public Properties

    }
}