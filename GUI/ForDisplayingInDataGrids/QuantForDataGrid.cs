using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MetaMorpheusGUI
{
    public class QuantForDataGrid
    {
        public QuantForDataGrid(string path)
        {
            FileName = Path.GetFileName(path);
            FilePath = path;
        }

        public string FileName { get; private set; }      
        public string Qcondition { get; set; }
        public string QbioRep { get; set; }
        public string Qfraction { get; set; }
        public string QtechRep { get; set; }
        public string FilePath { get; private set; }

        public void SetQconditionText(string text)
        {
            Qcondition = text;
        }
        public void SetQbioRepText(string text)
        {
            QbioRep = text;
        }
        public void SetQfractionText(string text)
        {
            Qfraction = text;
        }
        public void SetQtechRepText(string text)
        {
            QtechRep = text;
        }

    }
}
