using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace EngineLayer.Analysis
{
    public class Bin
    {
        #region Public Fields

        public string AA = "-";
        public Dictionary<char, int> residueCount;
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

        public int PepNlocCount { get; private set; }
        public int PepClocCount { get; private set; }
        public int ProtNlocCount { get; private set; }
        public int ProtClocCount { get; private set; }
        public string Combos { get; private set; } = "-";

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

        public void IdentifyResidues()
        {
            residueCount = new Dictionary<char, int>();
            foreach (var hehe in uniquePSMs.Values.Where(b => b.Item3.LocalizedScores != null))
            {
                double bestScore = hehe.Item3.LocalizedScores.Max();
                if (bestScore >= hehe.Item3.Score + 1 && !hehe.Item3.IsDecoy)
                {
                    for (int i = 0; i < hehe.Item1.Count(); i++)
                        if (bestScore - hehe.Item3.LocalizedScores[i] < 0.5)
                            if (residueCount.ContainsKey(hehe.Item1[i]))
                                residueCount[hehe.Item1[i]]++;
                            else
                                residueCount.Add(hehe.Item1[i], 1);
                    if (hehe.Item3.LocalizedScores.Max() - hehe.Item3.LocalizedScores[0] < 0.5)
                    {
                        PepNlocCount++;
                        if (hehe.Item3.OneBasedStartResidueInProtein.HasValue && hehe.Item3.OneBasedStartResidueInProtein.Value <= 2)
                            ProtNlocCount++;
                    }
                    if (hehe.Item3.LocalizedScores.Max() - hehe.Item3.LocalizedScores.Last() < 0.5)
                    {
                        PepClocCount++;
                        if (hehe.Item3.OneBasedEndResidueInProtein.HasValue && hehe.Item3.ProteinLength.HasValue && hehe.Item3.OneBasedEndResidueInProtein.Value == hehe.Item3.ProteinLength.Value)
                            ProtClocCount++;
                    }
                }
            }
        }

        public void IdentifyCombos(double v, HashSet<Tuple<double, double, double>> ok)
        {
            var okk = new HashSet<string>();
            foreach (var hm in ok)
                if (Math.Abs(hm.Item1 + hm.Item2 - MassShift) <= v && CountTarget < hm.Item3)
                    okk.Add("Combo " + Math.Min(hm.Item1, hm.Item2).ToString("F3", CultureInfo.InvariantCulture) + " and " + Math.Max(hm.Item1, hm.Item2).ToString("F3", CultureInfo.InvariantCulture));
            Combos = string.Join(" or ", okk);
        }

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