using MathNet.Numerics.Distributions;
using MetaMorpheus;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace NotchAnalysisCmd
{
    public class Program
    {
        public static UsefulProteomicsDatabases.Generated.unimod unimodDeserialized;
        public static UsefulProteomicsDatabases.Generated.obo psimodDeserialized;
        public static Dictionary<int, ChemicalFormulaModification> uniprotDeseralized;

        public static string unimodLocation = @"unimod_tables.xml";
        public static string psimodLocation = @"PSI-MOD.obo.xml";
        public static string elementsLocation = @"elements.dat";
        public static string uniprotLocation = @"ptmlist.txt";

        private static void Main(string[] args)
        {
            UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
            unimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation);
            psimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(psimodLocation);
            uniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation);
            Console.WriteLine(args[0]);

            IEnumerable<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>> psms_with_fdr = ReadPsmsWithFdrs(args[0]);

            MyTreeStructure myTreeStructure = new MyTreeStructure();
            Console.WriteLine("Generating bins...");
            myTreeStructure.GenerateBins(psms_with_fdr.Where(b => (b.QValue <= 0.01 && !b.Identification.isDecoy)).Select(b => b.Identification), 0.003);

            Console.WriteLine("Adding to bins...");
            myTreeStructure.AddToBins(psms_with_fdr.Where(b => b.QValue <= 0.01).Select(b => b.Identification));

            Console.WriteLine("Identifying bins...");
            IdentifyUnimodBins(myTreeStructure, 0.003);
            IdentifyUniprotBins(myTreeStructure, 0.003);
            IdentifyAA(myTreeStructure, 0.003);

            Console.WriteLine("Identifying combos...");
            IdentifyCombos(myTreeStructure, 0.003);

            Console.WriteLine("Extracting residues from localizeable...");
            IdentifyResidues(myTreeStructure);

            Console.WriteLine("Identifying mods...");
            IdentifyMods(myTreeStructure);

            Console.WriteLine("Writing bins...");
            using (StreamWriter output = new StreamWriter("okk.tsv"))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tmods\tUnimodID\tUnimodFormulas\tAA\tCombos\tres\tNlocCount\tClocCount\tUniprot");
                foreach (Bin bin in myTreeStructure.finalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.Count == 0 ? double.NaN : bin.CountDecoy / bin.Count).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : bin.LocalizeableTarget / bin.CountTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + string.Join(",", bin.mods.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + b.Value / bin.CountTarget))
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.AA
                        + "\t" + bin.combos
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : bin.NlocCount / bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : bin.ClocCount / bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.uniprotID);
                }
            }
        }

        public static void IdentifyMods(MyTreeStructure myTreeStructure)
        {
            foreach (Bin bin in myTreeStructure.finalBins)
            {
                bin.mods = new Dictionary<string, int>();
                foreach (PeptideSpectrumMatch hehe in bin.uniquePSMs.Values.Where(b => !b.isDecoy))
                {
                    int inModLevel = 0;
                    string currentMod = "";
                    for (int i = 0; i < hehe.Peptide.ExtendedSequence.Count(); i++)
                    {
                        char ye = hehe.Peptide.ExtendedSequence[i];
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
                                if (bin.mods.ContainsKey(currentMod))
                                    bin.mods[currentMod]++;
                                else
                                    bin.mods.Add(currentMod, 1);
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

        public static void IdentifyResidues(MyTreeStructure myTreeStructure)
        {
            foreach (Bin bin in myTreeStructure.finalBins)
            {
                bin.residueCount = new Dictionary<char, int>();
                foreach (PeptideSpectrumMatch hehe in bin.uniquePSMs.Values)
                {
                    double bestScore = hehe.LocalizedScores.Max();
                    if (bestScore >= hehe.MetaMorpheusScore + 1 && !hehe.isDecoy)
                    {
                        for (int i = 0; i < hehe.Peptide.BaseSequence.Count(); i++)
                        {
                            if (bestScore - hehe.LocalizedScores[i] < 0.5)
                            {
                                if (bin.residueCount.ContainsKey(hehe.Peptide.BaseSequence[i]))
                                    bin.residueCount[hehe.Peptide.BaseSequence[i]]++;
                                else
                                    bin.residueCount.Add(hehe.Peptide.BaseSequence[i], 1);
                            }
                        }
                        if (hehe.LocalizedScores.Max() - hehe.LocalizedScores[0] < 0.5)
                            bin.NlocCount++;
                        if (hehe.LocalizedScores.Max() - hehe.LocalizedScores.Last() < 0.5)
                            bin.ClocCount++;
                    }
                }
            }
        }

        public static void IdentifyCombos(MyTreeStructure myTreeStructure, double v)
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
                        okk.Add("Combo " + Math.Min(hm.Item1, hm.Item2) + " and " + Math.Max(hm.Item1, hm.Item2));
                    }
                }
                bin.combos = string.Join(" or ", okk);
            }
        }

        public static void IdentifyAA(MyTreeStructure myTreeStructure, double v)
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

        public static void IdentifyUniprotBins(MyTreeStructure myTreeStructure, double v)
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

        public static void IdentifyUnimodBins(MyTreeStructure myTreeStructure, double v)
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

        public static IEnumerable<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>> ReadPsmsWithFdrs(string v)
        {
            foreach (string line in File.ReadLines(v, Encoding.UTF8).Skip(1))
            {
                string[] ok = line.Split('\t');
                var hm = ok[29].Substring(1, ok[29].Length - 2).Split(',').Select(double.Parse);
                //PeptideWithSetModifications peptide = new PeptideWithSetModifications(ok[11], ok[12]);
                //PeptideSpectrumMatch identification = new PeptideSpectrumMatch(Convert.ToBoolean(ok[33]), Convert.ToDouble(ok[19]), hm.ToList(), Convert.ToDouble(ok[26]), peptide);
                //yield return new IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>(identification, Convert.ToInt32(ok[34]), Convert.ToInt32(ok[35]), Convert.ToDouble(ok[36]) / 100);

                yield return null;
            }
        }
    }
}