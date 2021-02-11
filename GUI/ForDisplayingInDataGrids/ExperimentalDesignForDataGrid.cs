using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MetaMorpheusGUI
{
    public class ExperimentalDesignForDataGrid
    {
        public ExperimentalDesignForDataGrid(string fullFilePathWithExtension)
        {
            FullFilePathWithExtension = fullFilePathWithExtension;
            FileNameWithExtension = Path.GetFileName(fullFilePathWithExtension);
        }

        public string FullFilePathWithExtension { get; private set; }
        public string FileNameWithExtension { get; private set; }
        public string Condition { get; set; }
        public string Biorep { get; set; }
        public string Fraction { get; set; }
        public string Techrep { get; set; }
    }
}
