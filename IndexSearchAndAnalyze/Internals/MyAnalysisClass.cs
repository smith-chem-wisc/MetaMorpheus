using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using UsefulProteomicsDatabases.Generated;

namespace IndexSearchAndAnalyze
{
    internal class MyAnalysisClass
    {
        public static void IdentifyAAsInCommon(BinTreeStructure myTreeStructure)
        {
            foreach (Bin bin in myTreeStructure.finalBins)
            {
                bin.AAsInCommon = new Dictionary<char, int>();
                foreach (var hehe in bin.uniquePSMs.Values.Where(b => !b.Item3.isDecoy))
                {
                    HashSet<char> chars = new HashSet<char>();
                    for (int i = 0; i < hehe.Item1.Count(); i++)
                    {
                        chars.Add(hehe.Item1[i]);
                    }
                    foreach (var ch in chars)
                    {
                        if (bin.AAsInCommon.ContainsKey(ch))
                        {
                            bin.AAsInCommon[ch]++;
                        }
                        else
                            bin.AAsInCommon.Add(ch, 1);
                    }
                }
            }
        }

        public static void IdentifyMods(BinTreeStructure myTreeStructure)
        {
            foreach (Bin bin in myTreeStructure.finalBins)
            {
                bin.modsInCommon = new Dictionary<string, int>();
                foreach (var hehe in bin.uniquePSMs.Values.Where(b => !b.Item3.isDecoy))
                {
                    int inModLevel = 0;
                    string currentMod = "";
                    for (int i = 0; i < hehe.Item2.Count(); i++)
                    {
                        char ye = hehe.Item2[i];
                        if (ye.Equals('('))
                        {
                            inModLevel++;
                            if (inModLevel == 1)
                            {
                                continue;
                            }
                        }
                        else if (ye.Equals(')'))
                        {
                            inModLevel--;
                            if (inModLevel == 0)
                            {
                                if (bin.modsInCommon.ContainsKey(currentMod))
                                    bin.modsInCommon[currentMod]++;
                                else
                                    bin.modsInCommon.Add(currentMod, 1);
                                currentMod = "";
                            }
                            continue;
                        }
                        if (inModLevel > 0)
                        {
                            currentMod += ye;
                        }
                    }
                }
            }
        }

        internal static void IdentifyResidues(BinTreeStructure myTreeStructure)
        {
            foreach (Bin bin in myTreeStructure.finalBins)
            {
                bin.residueCount = new Dictionary<char, int>();
                foreach (var hehe in bin.uniquePSMs.Values)
                {
                    double bestScore = hehe.Item3.thisPSM.LocalizedScores.Max();
                    if (bestScore >= hehe.Item3.thisPSM.Score + 1 && !hehe.Item3.isDecoy)
                    {
                        for (int i = 0; i < hehe.Item1.Count(); i++)
                        {
                            if (bestScore - hehe.Item3.thisPSM.LocalizedScores[i] < 0.5)
                            {
                                if (bin.residueCount.ContainsKey(hehe.Item1[i]))
                                    bin.residueCount[hehe.Item1[i]]++;
                                else
                                    bin.residueCount.Add(hehe.Item1[i], 1);
                            }
                        }
                        if (hehe.Item3.thisPSM.LocalizedScores.Max() - hehe.Item3.thisPSM.LocalizedScores[0] < 0.5)
                            bin.NlocCount++;
                        if (hehe.Item3.thisPSM.LocalizedScores.Max() - hehe.Item3.thisPSM.LocalizedScores.Last() < 0.5)
                            bin.ClocCount++;
                    }
                }
            }
        }

        internal static void IdentifyUnimodBins(BinTreeStructure myTreeStructure, double v, unimod unimodDeserialized)
        {
            foreach (var bin in myTreeStructure.finalBins)
            {
                HashSet<string> ok = new HashSet<string>();
                HashSet<string> okformula = new HashSet<string>();
                foreach (var hm in unimodDeserialized.modifications)
                {
                    if (Math.Abs(hm.mono_mass - bin.MassShift) <= v)
                    {
                        ok.Add(hm.full_name);
                        okformula.Add(hm.composition);
                    }
                }
                bin.UnimodId = string.Join(" or ", ok);
                bin.UnimodFormulas = string.Join(" or ", okformula);
            }
        }

        internal static void IdentifyUniprotBins(BinTreeStructure myTreeStructure, double v, Dictionary<int, ChemicalFormulaModification> uniprotDeseralized)
        {
            foreach (var bin in myTreeStructure.finalBins)
            {
                HashSet<string> ok = new HashSet<string>();
                foreach (var hm in uniprotDeseralized)
                {
                    if (Math.Abs(hm.Value.MonoisotopicMass - bin.MassShift) <= v)
                    {
                        ok.Add(hm.Value.NameAndSites);
                    }
                }
                bin.uniprotID = string.Join(" or ", ok);
            }
        }

        internal static void IdentifyCombos(BinTreeStructure myTreeStructure, double v)
        {
            double totalTargetCount = myTreeStructure.finalBins.Select(b => b.CountTarget).Sum();
            HashSet<Tuple<double, double, double>> ok = new HashSet<Tuple<double, double, double>>();
            foreach (var bin in myTreeStructure.finalBins.Where(b => Math.Abs(b.MassShift) > v))
                foreach (var bin2 in myTreeStructure.finalBins.Where(b => Math.Abs(b.MassShift) > v))
                    if (bin.CountTarget * bin2.CountTarget >= totalTargetCount * 3)
                        ok.Add(new Tuple<double, double, double>(bin.MassShift, bin2.MassShift, Math.Min(bin.CountTarget, bin2.CountTarget)));

            foreach (var bin in myTreeStructure.finalBins)
            {
                HashSet<string> okk = new HashSet<string>();
                foreach (var hm in ok)
                {
                    if (Math.Abs(hm.Item1 + hm.Item2 - bin.MassShift) <= v && bin.CountTarget < hm.Item3)
                    {
                        okk.Add("Combo " + Math.Min(hm.Item1, hm.Item2).ToString("F3", CultureInfo.InvariantCulture) + " and " + Math.Max(hm.Item1, hm.Item2).ToString("F3", CultureInfo.InvariantCulture));
                    }
                }
                bin.combos = string.Join(" or ", okk);
            }
        }

        internal static void IdentifyAA(BinTreeStructure myTreeStructure, double v)
        {
            foreach (var bin in myTreeStructure.finalBins)
            {
                HashSet<string> ok = new HashSet<string>();
                for (char c = 'A'; c <= 'Z'; c++)
                {
                    AminoAcid residue;
                    if (AminoAcid.TryGetResidue(c, out residue))
                    {
                        if (Math.Abs(residue.MonoisotopicMass - bin.MassShift) <= v)
                        {
                            ok.Add("Add " + residue.Name);
                        }
                        if (Math.Abs(residue.MonoisotopicMass + bin.MassShift) <= v)
                        {
                            ok.Add("Remove " + residue.Name);
                        }
                        for (char cc = 'A'; cc <= 'Z'; cc++)
                        {
                            AminoAcid residueCC;
                            if (AminoAcid.TryGetResidue(cc, out residueCC))
                            {
                                if (Math.Abs(residueCC.MonoisotopicMass + residue.MonoisotopicMass - bin.MassShift) <= v)
                                {
                                    ok.Add("Add (" + residue.Name + "+" + residueCC.Name + ")");
                                }
                                if (Math.Abs(residueCC.MonoisotopicMass + residue.MonoisotopicMass + bin.MassShift) <= v)
                                {
                                    ok.Add("Remove (" + residue.Name + "+" + residueCC.Name + ")");
                                }
                            }
                        }
                    }
                }
                bin.AA = string.Join(" or ", ok);
            }
        }

        internal static void IdentifyMine(BinTreeStructure myTreeStructure, double v)
        {
            List<MyInfo> myInfos = new List<MyInfo>();
            myInfos.Add(new MyInfo(0, "Exact match!"));
            myInfos.Add(new MyInfo(-48.128629, "Phosphorylation-Lysine: Probably reverse is the correct match"));
            myInfos.Add(new MyInfo(-76.134779, "Phosphorylation-Arginine: Probably reverse is the correct match"));
            myInfos.Add(new MyInfo(1.003, "1 MM"));
            myInfos.Add(new MyInfo(2.005, "2 MM"));
            myInfos.Add(new MyInfo(3.008, "3 MM"));
            myInfos.Add(new MyInfo(173.051055, "Acetylation + Methionine: Usually on protein N terminus"));
            myInfos.Add(new MyInfo(-91.009185, "neg Carbamidomethylation - H2S: Usually on cysteine."));
            myInfos.Add(new MyInfo(-32.008456, "oxidation and then loss of oxidized M side chain"));
            myInfos.Add(new MyInfo(-79.966331, "neg Phosphorylation. Probably real thing does not have it, but somehow matched! Might want to exclude."));
            myInfos.Add(new MyInfo(189.045969, "Carboxymethylated + Methionine. Usually on protein N terminus"));
            myInfos.Add(new MyInfo(356.20596, "Lysine+V+E or Lysine+L+D"));
            myInfos.Add(new MyInfo(239.126988, "Lysine+H(5) C(5) N O(2), possibly Nmethylmaleimide"));
            foreach (Bin bin in myTreeStructure.finalBins)
            {
                bin.mine = "";
                foreach (MyInfo myInfo in myInfos)
                {
                    if (Math.Abs(myInfo.MassShift - bin.MassShift) <= v)
                    {
                        bin.mine = myInfo.infostring;
                    }
                }
            }
        }
    }
}