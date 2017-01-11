using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public class Bin
    {
        public string UnimodId = "-";
        public string psimodID = "-";
        public string uniprotID = "-";
        public string UnimodFormulas = "-";
        public string AA = "-";
        public string combos = "-";
        public double MassShift { get; private set; }
        public Dictionary<char, int> residueCount;
        public int NlocCount;
        public int ClocCount;

        public Bin(double MassShift)
        {
            this.MassShift = MassShift;
            uniquePSMs = new Dictionary<string, Tuple<string, string, NewPsmWithFDR>>();
        }

        public int Count
        {
            get
            {
                return uniquePSMs.Count;
            }
        }

        public int CountDecoy
        {
            get
            {
                return uniquePSMs.Values.Where(b => b.Item3.isDecoy).Count();
            }
        }

        public int CountTarget
        {
            get
            {
                return Count - CountDecoy;
            }
        }

        public int LocalizeableTarget
        {
            get
            {
                return uniquePSMs.Values.Where(b => !b.Item3.isDecoy && b.Item3.thisPSM.LocalizedScores.Max() >= b.Item3.thisPSM.Score + 1).Count();
            }
        }

        public string mine { get; internal set; }
        public Dictionary<char, int> AAsInCommon { get; internal set; }

        public Dictionary<string, Tuple<string, string, NewPsmWithFDR>> uniquePSMs;

        internal void Add(NewPsmWithFDR ok)
        {
            if (uniquePSMs.ContainsKey(ok.thisPSM.FullSequence))
            {
                var current = uniquePSMs[ok.thisPSM.FullSequence];
                if (current.Item3.thisPSM.Score < ok.thisPSM.Score)
                    uniquePSMs[ok.thisPSM.FullSequence] = new Tuple<string, string, NewPsmWithFDR>(ok.thisPSM.BaseSequence, ok.thisPSM.FullSequence, ok);
            }
            else
                uniquePSMs.Add(ok.thisPSM.FullSequence, new Tuple<string, string, NewPsmWithFDR>(ok.thisPSM.BaseSequence, ok.thisPSM.FullSequence, ok));
        }

        public Dictionary<string, int> modsInCommon;

        public double ComputeZ(double v)
        {
            return Math.Sqrt(Count) * ((double)CountDecoy / Count - v) / (v * (1 - v));
        }
    }
}