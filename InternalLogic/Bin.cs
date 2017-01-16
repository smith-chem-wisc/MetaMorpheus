using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public class Bin
    {
        #region Public Fields

        public string UnimodId = "-";
        public string psimodID = "-";
        public string uniprotID = "-";
        public string UnimodFormulas = "-";
        public string AA = "-";
        public string combos = "-";
        public Dictionary<char, int> residueCount;
        public int NlocCount;
        public int ClocCount;
        public Dictionary<string, Tuple<string, string, NewPsmWithFDR>> uniquePSMs;
        public Dictionary<string, int> modsInCommon;

        #endregion Public Fields

        #region Public Constructors

        public Bin(double MassShift)
        {
            this.MassShift = MassShift;
            uniquePSMs = new Dictionary<string, Tuple<string, string, NewPsmWithFDR>>();
        }

        #endregion Public Constructors

        #region Public Properties

        public double MassShift { get; private set; }

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
                return uniquePSMs.Values.Count(b => b.Item3.isDecoy);
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
                return uniquePSMs.Values.Count(b => !b.Item3.isDecoy && b.Item3.thisPSM.LocalizedScores.Max() >= b.Item3.thisPSM.Score + 1);
            }
        }

        public string mine { get; internal set; }
        public Dictionary<char, int> AAsInCommon { get; internal set; }

        #endregion Public Properties

        #region Public Methods

        public double ComputeZ(double v)
        {
            return Math.Sqrt(Count) * ((double)CountDecoy / Count - v) / (v * (1 - v));
        }

        #endregion Public Methods

        #region Internal Methods

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

        #endregion Internal Methods
    }
}