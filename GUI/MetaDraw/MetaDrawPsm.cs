using MassSpectrometry;
using System.Collections.Generic;

namespace MetaMorpheusGUI
{
    public class MetaDrawPsm
    {
        public int ScanNum { get; private set; }
        public string FullSequence { get; private set; }
        public string FileName { get; private set; }
        public List<TheoreticalFragmentIon> FragmentIons { get; private set; }

        public MetaDrawPsm(int oneBasedScanNumber, string fileName, string fullSequence, List<TheoreticalFragmentIon> fragmentIons)
        {
            this.ScanNum = oneBasedScanNumber;
            this.FileName = fileName;
            this.FullSequence = fullSequence;
            this.FragmentIons = fragmentIons;
        }
    }
}
