using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Analysis
{
    public class Bin
    {
        #region Public Fields

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

        public string UnimodDiffs { get; private set; } = "-";
        public string UniprotID { get; private set; } = "-";
        public string UnimodFormulas { get; private set; } = "-";
        public string UnimodId { get; private set; } = "-";
        public double MassShift { get; }

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
                return uniquePSMs.Values.Count(b => b.Item3.IsDecoy);
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
                return uniquePSMs.Values.Where(b => b.Item3.LocalizedScores != null).Count(b => !b.Item3.IsDecoy && b.Item3.LocalizedScores.Max() >= b.Item3.Score + 1);
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

        public void IdentifyUniprotBins(double v)
        {
            var ok = new HashSet<string>();
            foreach (var hm in GlobalEngineLevelSettings.UniprotDeseralized)
            {
                if (hm is ModificationWithMass theMod && Math.Abs(theMod.monoisotopicMass - MassShift) <= v)
                    ok.Add(hm.id);
            }
            UniprotID = string.Join(" or ", ok);
        }

        public void IdentifyAA(double v)
        {
            var ok = new HashSet<string>();
            for (char c = 'A'; c <= 'Z'; c++)
            {
                if (Residue.TryGetResidue(c, out Residue residue))
                {
                    if (Math.Abs(residue.MonoisotopicMass - MassShift) <= v)
                        ok.Add("Add " + residue.Name);
                    if (Math.Abs(residue.MonoisotopicMass + MassShift) <= v)
                        ok.Add("Remove " + residue.Name);
                    for (char cc = 'A'; cc <= 'Z'; cc++)
                    {
                        if (Residue.TryGetResidue(cc, out Residue residueCC))
                        {
                            if (Math.Abs(residueCC.MonoisotopicMass + residue.MonoisotopicMass - MassShift) <= v)
                                ok.Add("Add (" + residue.Name + "+" + residueCC.Name + ")");
                            if (Math.Abs(residueCC.MonoisotopicMass + residue.MonoisotopicMass + MassShift) <= v)
                                ok.Add("Remove (" + residue.Name + "+" + residueCC.Name + ")");
                        }
                    }
                }
            }
            AA = string.Join(" or ", ok);
        }

        public void IdentifyUnimodBins(double v)
        {
            var ok = new HashSet<string>();
            var okformula = new HashSet<string>();
            var okDiff = new HashSet<double>();
            foreach (var hm in GlobalEngineLevelSettings.UnimodDeserialized)
            {
                var theMod = hm as ModificationWithMassAndCf;
                if (Math.Abs(theMod.monoisotopicMass - MassShift) <= v)
                {
                    ok.Add(hm.id);
                    okformula.Add(theMod.chemicalFormula.Formula);
                    okDiff.Add(theMod.monoisotopicMass - MassShift);
                }
            }
            UnimodId = string.Join(" or ", ok);
            UnimodFormulas = string.Join(" or ", okformula);
            UnimodDiffs = string.Join(" or ", okDiff);
        }

        #endregion Public Methods

        #region Internal Methods

        internal void Add(Psm ok)
        {
            if (ok.FullSequence != null)
            {
                if (uniquePSMs.ContainsKey(ok.FullSequence))
                {
                    var current = uniquePSMs[ok.FullSequence];
                    if (current.Item3.Score < ok.Score)
                        uniquePSMs[ok.FullSequence] = new Tuple<string, string, Psm>(ok.BaseSequence, ok.FullSequence, ok);
                }
                else
                    uniquePSMs.Add(ok.FullSequence, new Tuple<string, string, Psm>(ok.BaseSequence, ok.FullSequence, ok));
            }
        }

        #endregion Internal Methods
    }
}