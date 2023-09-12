namespace MetaMorpheusGUI
{
    internal class ModListForSearchTask
    {
        #region Public Constructors

        public ModListForSearchTask(string filePath)
        {
            this.FileName = filePath;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Fixed { get; set; }
        public bool Variable { get; set; }
        public bool Localize { get; set; }
        public bool AlwaysKeep { get; set; }

        public string FileName { get; }

        #endregion Public Properties
    }
}