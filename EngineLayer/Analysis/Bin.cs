using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Analysis
{
    public class Bin
    {

        #region Public Fields

        public string UnimodId = "-";
        public string psimodID = "-";
        public string uniprotID = "-";
        public string UnimodFormulas = "-";
        public string UnimodDiffs = "-";
        public string AA = "-";
        public string combos = "-";
        public Dictionary<char, int> residueCount;
        public int pepNlocCount;
        public int pepClocCount;
        public int protNlocCount;
        public int protClocCount;
        public Dictionary<string, Tuple<string, string, Psm>> uniquePSMs;
        public Dictionary<string, int> modsInCommon;

        #endregion Public Fields

        #region Public Constructors

        public Bin(double massShift)
        {
            this.MassShift = massShift;
            uniquePSMs = new Dictionary<string, Tuple<string, string, Psm>>();
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
                return uniquePSMs.Values.Count(b => b.Item3.MostProbableProteinInfo.IsDecoy);
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
                return uniquePSMs.Values.Count(b => !b.Item3.MostProbableProteinInfo.IsDecoy && b.Item3.LocalizationResults.LocalizedScores.Max() >= b.Item3.Score + 1);
            }
        }

        public string Mine { get; internal set; }
        public Dictionary<char, int> AAsInCommon { get; internal set; }
        public int Overlapping { get; internal set; }
        public double FracWithSingle { get; set; }
        public double MedianLength { get; internal set; }

        #endregion Public Properties

        #region Public Methods

        public double ComputeZ(double v)
        {
            return Math.Sqrt(Count) * ((double)CountDecoy / Count - v) / (v * (1 - v));
        }

        #endregion Public Methods

        #region Internal Methods

        internal void Add(Psm ok)
        {
            if (uniquePSMs.ContainsKey(ok.MostProbableProteinInfo.FullSequence))
            {
                var current = uniquePSMs[ok.MostProbableProteinInfo.FullSequence];
                if (current.Item3.Score < ok.Score)
                    uniquePSMs[ok.MostProbableProteinInfo.FullSequence] = new Tuple<string, string, Psm>(ok.MostProbableProteinInfo.BaseSequence, ok.MostProbableProteinInfo.FullSequence, ok);
            }
            else
                uniquePSMs.Add(ok.MostProbableProteinInfo.FullSequence, new Tuple<string, string, Psm>(ok.MostProbableProteinInfo.BaseSequence, ok.MostProbableProteinInfo.FullSequence, ok));
        }

        #endregion Internal Methods

    }
}