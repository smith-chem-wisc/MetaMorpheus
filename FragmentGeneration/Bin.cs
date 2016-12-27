using System;
using System.Collections.Generic;
using System.Linq;

namespace FragmentGeneration
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
        public double NlocCount;
        public double ClocCount;

        public Bin(double MassShift)
        {
            this.MassShift = MassShift;
            uniquePSMs = new Dictionary<string, Tuple<string, string, NewPsmWithFDR>>();
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
                return uniquePSMs.Values.Where(b => b.Item3.isDecoy).Count();
            }
        }

        public double CountTarget
        {
            get
            {
                return Count - CountDecoy;
            }
        }

        public double LocalizeableTarget
        {
            get
            {
                return uniquePSMs.Values.Where(b => !b.Item3.isDecoy && b.Item3.thisPSM.LocalizedScores.Max() >= b.Item3.thisPSM.ScoreFromSearch + 1).Count();
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
                if (current.Item3.thisPSM.ScoreFromSearch < ok.thisPSM.ScoreFromSearch)
                    uniquePSMs[ok.thisPSM.FullSequence] = new Tuple<string, string, NewPsmWithFDR>(ok.thisPSM.BaseSequence, ok.thisPSM.FullSequence, ok);
            }
            else
                uniquePSMs.Add(ok.thisPSM.FullSequence, new Tuple<string, string, NewPsmWithFDR>(ok.thisPSM.BaseSequence, ok.thisPSM.FullSequence, ok));
        }

        public Dictionary<string, int> modsInCommon;

        public double ComputeZ(double v)
        {
            return Math.Sqrt(Count) * (CountDecoy / Count - v) / (v * (1 - v));
        }
    }
}