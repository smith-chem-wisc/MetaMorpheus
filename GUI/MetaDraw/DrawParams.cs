using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    public class DrawParams
    {
        public DrawParams()
        {
            ProductMassTolerance = new PpmTolerance(20);
            BIons = true;
            YIons = true;
            ZdotIons = false;
            CIons = false;
            TopNpeaks = 200;
            MinRatio = 0.01;
            TrimMs1Peaks = false;
            TrimMsMsPeaks = true;
        }

        public bool BIons { get; set; }
        public bool YIons { get; set; }
        public bool ZdotIons { get; set; }
        public bool CIons { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public int TopNpeaks { get; set; }
        public double MinRatio { get; set; }
        public bool TrimMs1Peaks { get; set; }
        public bool TrimMsMsPeaks { get; set; }
    }
}
