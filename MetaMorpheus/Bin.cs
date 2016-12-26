using System;
using System.Collections.Generic;
using System.Linq;

namespace MetaMorpheus
{
    public class Bin
    {
        public string UnimodId = "-";
        public string psimodID = "-";
        public string uniprotID = "-";
        public string UnimodFormulas = "-";
        public string AA = "-";
        public double fracOfTotal;
        public string combos = "-";

        public double MassShift { get; private set; }

        public Bin(double MassShift)
        {
            this.MassShift = MassShift;
            uniquePSMs = new Dictionary<string, PeptideSpectrumMatch>();
        }

        public double Count
        {
            get
            {
                return uniquePSMs.Count;
            }
        }

        public double CountDecoy
        {
            get
            {
                return uniquePSMs.Values.Where(b => b.isDecoy).Count();
            }
        }

        public double CountTarget
        {
            get
            {
                return uniquePSMs.Values.Where(b => !b.isDecoy).Count();
            }
        }

        public double LocalizeableTarget
        {
            get
            {
                return uniquePSMs.Values.Where(b => !b.isDecoy && b.LocalizedScores.Max() >= b.MetaMorpheusScore + 1).Count();
            }
        }

        public Dictionary<string, PeptideSpectrumMatch> uniquePSMs;
        public Dictionary<char, int> residueCount;
        public double NlocCount;
        public double ClocCount;
        public Dictionary<string, int> mods;

        internal void Add(PeptideSpectrumMatch ok)
        {
            var seq = ok.Peptide.ExtendedSequence;
            if (uniquePSMs.ContainsKey(seq))
            {
                var current = uniquePSMs[seq];
                if (current.MetaMorpheusScore < ok.MetaMorpheusScore)
                    uniquePSMs[seq] = ok;
            }
            else
                uniquePSMs.Add(seq, ok);
        }

        public double ComputeZ(double v)
        {
            return Math.Sqrt(Count) * (CountDecoy / Count - v) / (v * (1 - v));
        }
    }
}