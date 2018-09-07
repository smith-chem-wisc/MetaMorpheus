using Proteomics;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace EngineLayer.HistogramAnalysis
{
    public class Bin
    {
        public string AA = "-";
        public Dictionary<char, int> ResidueCount;
        public Dictionary<string, Tuple<string, string, PeptideSpectralMatch>> UniquePSMs;
        public Dictionary<string, int> ModsInCommon;

        public Bin(double massShift)
        {
            this.MassShift = massShift;
            UniquePSMs = new Dictionary<string, Tuple<string, string, PeptideSpectralMatch>>();
        }

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
                return UniquePSMs.Count;
            }
        }

        public int CountDecoy
        {
            get
            {
                return UniquePSMs.Values.Count(b => b.Item3.IsDecoy);
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
                return UniquePSMs.Values.Where(b => b.Item3.LocalizedScores != null).Count(b => !b.Item3.IsDecoy && b.Item3.LocalizedScores.Max() >= b.Item3.Score + 1);
            }
        }

        public string Mine { get; internal set; }
        public Dictionary<char, int> AAsInCommon { get; internal set; }
        public int Overlapping { get; internal set; }
        public double FracWithSingle { get; set; }
        public double MedianLength { get; internal set; }

        public void IdentifyResidues()
        {
            ResidueCount = new Dictionary<char, int>();
            foreach (var hehe in UniquePSMs.Values.Where(b => b.Item3.LocalizedScores != null))
            {
                double bestScore = hehe.Item3.LocalizedScores.Max();
                if (bestScore >= hehe.Item3.Score + 1 && !hehe.Item3.IsDecoy)
                {
                    for (int i = 0; i < hehe.Item1.Count(); i++)
                        if (bestScore - hehe.Item3.LocalizedScores[i] < 0.5)
                            if (ResidueCount.ContainsKey(hehe.Item1[i]))
                                ResidueCount[hehe.Item1[i]]++;
                            else
                                ResidueCount.Add(hehe.Item1[i], 1);
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
            Combos = string.Join("|", okk);
        }

        public double ComputeZ(double v)
        {
            return Math.Sqrt(Count) * ((double)CountDecoy / Count - v) / (v * (1 - v));
        }

        public void IdentifyUniprotBins(double v)
        {
            var modIdOptions = new HashSet<string>();
            foreach (Modification mod in GlobalVariables.UniprotDeseralized)
            {
                if (mod.MonoisotopicMass.HasValue && Math.Abs(mod.MonoisotopicMass.Value - MassShift) <= v)
                {
                    modIdOptions.Add(mod.IdWithMotif);
                }
            }
            UniprotID = string.Join("|", modIdOptions);
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
            AA = string.Join("|", ok);
        }

        public void IdentifyUnimodBins(double v)
        {
            var ok = new HashSet<string>();
            var okformula = new HashSet<string>();
            var okDiff = new HashSet<double>();
            foreach (var hm in GlobalVariables.UnimodDeserialized)
            {
                var theMod = hm as Modification;
                if (Math.Abs(theMod.MonoisotopicMass.Value - MassShift) <= v)
                {
                    ok.Add(hm.IdWithMotif);
                    okformula.Add(theMod.ChemicalFormula.Formula);
                    okDiff.Add(theMod.MonoisotopicMass.Value - MassShift);
                }
            }
            UnimodId = string.Join("|", ok);
            UnimodFormulas = string.Join("|", okformula);
            UnimodDiffs = string.Join("|", okDiff);
        }

        internal void Add(PeptideSpectralMatch ok)
        {
            if (ok.FullSequence != null)
            {
                if (UniquePSMs.ContainsKey(ok.FullSequence))
                {
                    var current = UniquePSMs[ok.FullSequence];
                    if (current.Item3.Score < ok.Score)
                        UniquePSMs[ok.FullSequence] = new Tuple<string, string, PeptideSpectralMatch>(ok.BaseSequence, ok.FullSequence, ok);
                }
                else
                    UniquePSMs.Add(ok.FullSequence, new Tuple<string, string, PeptideSpectralMatch>(ok.BaseSequence, ok.FullSequence, ok));
            }
        }
    }
}