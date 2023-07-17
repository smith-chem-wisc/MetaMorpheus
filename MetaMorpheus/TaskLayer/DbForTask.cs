using EngineLayer;
using System.IO;

namespace TaskLayer
{
    public class DbForTask
    {
        public DbForTask(string filePath, bool isContaminant)
        {
            FilePath = filePath;
            IsContaminant = isContaminant;
            FileName = System.IO.Path.GetFileName(filePath);
            IsSpectralLibrary = GlobalVariables.GetFileExtension(filePath).ToLowerInvariant() == ".msp";
        }

        public bool IsSpectralLibrary { get; }
        public string FilePath { get; }
        public bool IsContaminant { get; }
        public string FileName { get; }

        public void CloseConnections()
        {
            if(this.IsSpectralLibrary)
            {
                //this.FileName.Close();
            }
        }
    }
}