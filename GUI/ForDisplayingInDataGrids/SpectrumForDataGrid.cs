using System.IO;

namespace MetaMorpheusGUI
{
    class SpectrumForDataGrid
    {
        public SpectrumForDataGrid(int scanNum, string fullSequence)
        {
            ScanNum = scanNum;
            FullSequence = fullSequence;
        }

        public int ScanNum { get; set; }

        public string FullSequence { get; set; }
    }
}
