namespace MetaMorpheusGUI
{
    internal class ModListForGptmdTask
    {
        #region Public Constructors

        public ModListForGptmdTask(string filePath)
        {
            this.FileName = filePath;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Fixed { get; set; }
        public bool Variable { get; set; }
        public bool Localize { get; set; }
        public bool Gptmd { get; set; }

        public string FileName { get; }

        #endregion Public Properties
    }
}