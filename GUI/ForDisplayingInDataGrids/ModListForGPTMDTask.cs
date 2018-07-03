namespace MetaMorpheusGUI
{
    internal class ModListForGptmdTask
    {
        public ModListForGptmdTask(string filePath)
        {
            FileName = filePath;
        }

        public bool Fixed { get; set; }
        public bool Variable { get; set; }
        public bool Localize { get; set; }
        public bool Gptmd { get; set; }

        public string FileName { get; }
    }
}