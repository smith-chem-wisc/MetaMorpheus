namespace MetaMorpheusGUI
{
    internal class ModListForCalibrationTask
    {
        public ModListForCalibrationTask(string filePath)
        {
            this.FileName = filePath;
        }

        public bool Fixed { get; set; }
        public bool Variable { get; set; }
        public bool Localize { get; set; }

        public string FileName { get; }
    }
}