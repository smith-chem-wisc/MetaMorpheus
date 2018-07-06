using System.IO;

namespace MetaMorpheusGUI
{
    public class ExperimentalDesignForDataGrid
    {
        public ExperimentalDesignForDataGrid(string filename)
        {
            FileName = Path.GetFileNameWithoutExtension(filename);
        }

        public string FileName { get; private set; }
        public string Condition { get; set; }
        public string Biorep { get; set; }
        public string Fraction { get; set; }
        public string Techrep { get; set; }

        public void SetQconditionText(string text)
        {
            Condition = text;
        }

        public void SetQbioRepText(string text)
        {
            Biorep = text;
        }

        public void SetQfractionText(string text)
        {
            Fraction = text;
        }

        public void SetQtechRepText(string text)
        {
            Techrep = text;
        }
    }
}