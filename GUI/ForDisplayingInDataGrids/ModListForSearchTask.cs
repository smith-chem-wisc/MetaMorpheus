namespace MetaMorpheusGUI
{
    internal class ModListForSearchTask
    {
        public ModListForSearchTask(string filePath)
        {
            FileName = filePath;
        }

        public bool Fixed { get; set; }
        public bool Variable { get; set; }
        public bool Localize { get; set; }
        public bool AlwaysKeep { get; set; }

        public string FileName { get; }
    }
}