using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    public class MetaDrawPsm
    {
        public int ScanNum { get; private set; }
        public string FullSequence { get; private set; }
        public PeptideWithSetModifications Peptide { get; private set; }

        public MetaDrawPsm(int oneBasedScanNumber, PeptideWithSetModifications peptide)
        {
            this.ScanNum = oneBasedScanNumber;
            this.FullSequence = peptide.Sequence;
            this.Peptide = peptide;
        }
    }
}
