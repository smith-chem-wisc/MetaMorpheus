using MathNet.Numerics.Statistics;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace EngineLayer.Analysis
{
    public class BinTreeStructure
    {

        #region Private Fields

        private const int minAdditionalPsmsInBin = 1;

        #endregion Private Fields

        #region Public Properties

        public List<Bin> FinalBins { get; private set; }

        #endregion Public Properties

        #region Internal Methods

        internal void GenerateBins(List<NewPsmWithFdr> targetAndDecoyMatches, double dc)
        {
            List<double> listOfMassShifts = targetAndDecoyMatches.Select(b => b.thisPSM.ScanPrecursorMass - b.thisPSM.PeptideMonoisotopicMass).OrderBy(b => b).ToList();
            double minMassShift = listOfMassShifts.Min();
            double maxMassShift = listOfMassShifts.Max();

            int[] p = new int[listOfMassShifts.Count];

            int firstIndex = 0;
            int lastIndex = 0;
            for (int i = 0; i < listOfMassShifts.Count; i++)
            {
                var thisMassShift = listOfMassShifts[i];

                while (thisMassShift - listOfMassShifts[firstIndex] > dc)
                    firstIndex++;
                while (lastIndex + 1 < listOfMassShifts.Count && listOfMassShifts[lastIndex + 1] - thisMassShift <= dc)
                    lastIndex++;

                p[i] = lastIndex - firstIndex;
            }

            int maxP = p.Max();
            double[] sigma = new double[listOfMassShifts.Count];

            for (int i = 0; i < listOfMassShifts.Count; i++)
            {
                var thisMassShift = listOfMassShifts[i];
                var thisP = p[i];
                if (thisP == maxP)
                    sigma[i] = Math.Max(maxMassShift - thisMassShift, thisMassShift - minMassShift);
                else
                {
                    // SIGMA IS THE DISTANCE TO THE CLOSEST MASS SHIFT THAT HAS A HIGHER P VALUE THAN ITSELF

                    sigma[i] = getSigma(thisMassShift, thisP, i, listOfMassShifts, p);
                }
            }

            var listokbin = new List<OkBin>();
            for (int i = 0; i < sigma.Count(); i++)
                listokbin.Add(new OkBin(listOfMassShifts[i], sigma[i], p[i]));

            var prelimBins = new HashSet<double>();
            foreach (OkBin okbin in listokbin.OrderByDescending(b => b.p))
            {
                if (okbin.sigma < dc || okbin.p < minAdditionalPsmsInBin)
                    continue;
                bool add = true;
                foreach (double a in prelimBins)
                {
                    if (Math.Abs(okbin.massShift - a) <= dc)
                    {
                        add = false;
                        break;
                    }
                }
                if (add)
                    prelimBins.Add(okbin.massShift);
            }

            var forFinalBins = new Dictionary<double, List<double>>();
            foreach (double ok in prelimBins)
                forFinalBins.Add(ok, new List<double>());
            foreach (double a in listOfMassShifts)
                foreach (double b in prelimBins)
                    if (Math.Abs(a - b) <= dc)
                        forFinalBins[b].Add(a);

            FinalBins = forFinalBins.Select(b => new Bin(b.Value.Average())).ToList();

            for (int i = 0; i < targetAndDecoyMatches.Count; i++)
            {
                foreach (Bin bin in FinalBins)
                    if (Math.Abs(targetAndDecoyMatches[i].thisPSM.ScanPrecursorMass - targetAndDecoyMatches[i].thisPSM.PeptideMonoisotopicMass - bin.MassShift) <= dc)
                        bin.Add(targetAndDecoyMatches[i]);
            }

            FinalBins = FinalBins.Where(b => b.Count > 1).ToList();

            IdentifyUnimodBins(dc);
            IdentifyUniprotBins(dc);
            IdentifyAA(dc);

            IdentifyCombos(dc);

            IdentifyResidues();

            IdentifyMods();

            IdentifyAAsInCommon();

            IdentifyMine(dc);

            OverlappingIonSequences();

            IdentifyFracWithSingle();
            IdentifyMedianLength();
        }

        #endregion Internal Methods

        #region Private Methods

        private static double getSigma(double thisMassShift, int thisP, int i, List<double> listOfMassShifts, int[] p)
        {
            int currentDown = i - 1;
            int currentUp = i + 1;
            int lookingAtP;
            double distDown, distUp;
            while (true)
            {
                distDown = currentDown == -1 ? double.MaxValue : thisMassShift - listOfMassShifts[currentDown];
                distUp = currentUp == listOfMassShifts.Count ? double.MaxValue : listOfMassShifts[currentUp] - thisMassShift;
                if (distDown < distUp)
                {
                    lookingAtP = p[currentDown];
                    if (lookingAtP > thisP)
                        return distDown;
                    currentDown--;
                }
                else
                {
                    lookingAtP = p[currentUp];
                    if (lookingAtP > thisP)
                        return distUp;
                    currentUp++;
                }
            }
        }

        private void OverlappingIonSequences()
        {
            foreach (Bin bin in FinalBins)
                foreach (var hm in bin.uniquePSMs.Where(b => !b.Value.Item3.IsDecoy))
                {
                    var ya = hm.Value.Item3.thisPSM.newPsm.matchedIonsList;
                    if (ya.ContainsKey(ProductType.B)
                        && ya.ContainsKey(ProductType.Y)
                        && ya[ProductType.B].Any(b => b > 0)
                        && ya[ProductType.Y].Any(b => b > 0)
                        && ya[ProductType.B].Last(b => b > 0) + ya[ProductType.Y].Last(b => b > 0) > hm.Value.Item3.thisPSM.PeptideMonoisotopicMass)
                        bin.Overlapping++;
                }
        }

        private void IdentifyFracWithSingle()
        {
            foreach (Bin bin in FinalBins)
            {
                var numTarget = bin.uniquePSMs.Values.Count(b => !b.Item3.IsDecoy);
                if (numTarget > 0)
                    bin.FracWithSingle = (double)bin.uniquePSMs.Values.Count(b => !b.Item3.IsDecoy && b.Item3.thisPSM.peptidesWithSetModifications.Count == 1) / numTarget;
            }
        }

        private void IdentifyMedianLength()
        {
            foreach (Bin bin in FinalBins)
            {
                var numTarget = bin.uniquePSMs.Values.Count(b => !b.Item3.IsDecoy);
                if (numTarget > 0)
                    bin.MedianLength = Statistics.Median(bin.uniquePSMs.Values.Where(b => !b.Item3.IsDecoy).Select(b => (double)b.Item3.thisPSM.BaseSequence.Length));
            }
        }

        private void IdentifyAAsInCommon()
        {
            foreach (Bin bin in FinalBins)
            {
                bin.AAsInCommon = new Dictionary<char, int>();
                foreach (var hehe in bin.uniquePSMs.Values.Where(b => !b.Item3.IsDecoy))
                {
                    var chars = new HashSet<char>();
                    for (int i = 0; i < hehe.Item1.Count(); i++)
                        chars.Add(hehe.Item1[i]);
                    foreach (var ch in chars)
                        if (bin.AAsInCommon.ContainsKey(ch))
                            bin.AAsInCommon[ch]++;
                        else
                            bin.AAsInCommon.Add(ch, 1);
                }
            }
        }

        private void IdentifyMods()
        {
            foreach (Bin bin in FinalBins)
            {
                bin.modsInCommon = new Dictionary<string, int>();
                foreach (var hehe in bin.uniquePSMs.Values.Where(b => !b.Item3.IsDecoy))
                {
                    int inModLevel = 0;
                    string currentMod = "";
                    for (int i = 0; i < hehe.Item2.Count(); i++)
                    {
                        char ye = hehe.Item2[i];
                        if (ye.Equals('['))
                        {
                            inModLevel++;
                            if (inModLevel == 1)
                                continue;
                        }
                        else if (ye.Equals(']'))
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
                            currentMod += ye;
                    }
                }
            }
        }

        private void IdentifyResidues()
        {
            foreach (Bin bin in FinalBins)
            {
                bin.residueCount = new Dictionary<char, int>();
                foreach (var hehe in bin.uniquePSMs.Values)
                {
                    double bestScore = hehe.Item3.thisPSM.LocalizedScores.Max();
                    if (bestScore >= hehe.Item3.thisPSM.Score + 1 && !hehe.Item3.IsDecoy)
                    {
                        for (int i = 0; i < hehe.Item1.Count(); i++)
                            if (bestScore - hehe.Item3.thisPSM.LocalizedScores[i] < 0.5)
                                if (bin.residueCount.ContainsKey(hehe.Item1[i]))
                                    bin.residueCount[hehe.Item1[i]]++;
                                else
                                    bin.residueCount.Add(hehe.Item1[i], 1);
                        if (hehe.Item3.thisPSM.LocalizedScores.Max() - hehe.Item3.thisPSM.LocalizedScores[0] < 0.5)
                            bin.NlocCount++;
                        if (hehe.Item3.thisPSM.LocalizedScores.Max() - hehe.Item3.thisPSM.LocalizedScores.Last() < 0.5)
                            bin.ClocCount++;
                    }
                }
            }
        }

        private void IdentifyUnimodBins(double v)
        {
            foreach (var bin in FinalBins)
            {
                var ok = new HashSet<string>();
                var okformula = new HashSet<string>();
                foreach (var hm in MyEngine.UnimodDeserialized)
                {
                    var theMod = hm as ModificationWithMassAndCf;
                    if (Math.Abs(theMod.monoisotopicMass - bin.MassShift) <= v)
                    {
                        ok.Add(hm.id);
                        okformula.Add(theMod.chemicalFormula.Formula);
                    }
                }
                bin.UnimodId = string.Join(" or ", ok);
                bin.UnimodFormulas = string.Join(" or ", okformula);
            }
        }

        private void IdentifyUniprotBins(double v)
        {
            foreach (var bin in FinalBins)
            {
                var ok = new HashSet<string>();
                foreach (var hm in MyEngine.UniprotDeseralized)
                {
                    var theMod = hm as ModificationWithMass;
                    if (theMod != null && Math.Abs(theMod.monoisotopicMass - bin.MassShift) <= v)
                        ok.Add(hm.id);
                }
                bin.uniprotID = string.Join(" or ", ok);
            }
        }

        private void IdentifyCombos(double v)
        {
            double totalTargetCount = FinalBins.Select(b => b.CountTarget).Sum();
            var ok = new HashSet<Tuple<double, double, double>>();

            // For every non-zero bin
            foreach (var bin in FinalBins.Where(b => Math.Abs(b.MassShift) > v))
                foreach (var bin2 in FinalBins.Where(b => Math.Abs(b.MassShift) > v))
                    if (bin.CountTarget * bin2.CountTarget >= totalTargetCount)
                        ok.Add(new Tuple<double, double, double>(bin.MassShift, bin2.MassShift, Math.Min(bin.CountTarget, bin2.CountTarget)));

            foreach (var bin in FinalBins)
            {
                var okk = new HashSet<string>();
                foreach (var hm in ok)
                    if (Math.Abs(hm.Item1 + hm.Item2 - bin.MassShift) <= v && bin.CountTarget < hm.Item3)
                        okk.Add("Combo " + Math.Min(hm.Item1, hm.Item2).ToString("F3", CultureInfo.InvariantCulture) + " and " + Math.Max(hm.Item1, hm.Item2).ToString("F3", CultureInfo.InvariantCulture));
                bin.combos = string.Join(" or ", okk);
            }
        }

        private void IdentifyAA(double v)
        {
            foreach (var bin in FinalBins)
            {
                var ok = new HashSet<string>();
                for (char c = 'A'; c <= 'Z'; c++)
                {
                    Residue residue;
                    if (Residue.TryGetResidue(c, out residue))
                    {
                        if (Math.Abs(residue.MonoisotopicMass - bin.MassShift) <= v)
                            ok.Add("Add " + residue.Name);
                        if (Math.Abs(residue.MonoisotopicMass + bin.MassShift) <= v)
                            ok.Add("Remove " + residue.Name);
                        for (char cc = 'A'; cc <= 'Z'; cc++)
                        {
                            Residue residueCC;
                            if (Residue.TryGetResidue(cc, out residueCC))
                            {
                                if (Math.Abs(residueCC.MonoisotopicMass + residue.MonoisotopicMass - bin.MassShift) <= v)
                                    ok.Add("Add (" + residue.Name + "+" + residueCC.Name + ")");
                                if (Math.Abs(residueCC.MonoisotopicMass + residue.MonoisotopicMass + bin.MassShift) <= v)
                                    ok.Add("Remove (" + residue.Name + "+" + residueCC.Name + ")");
                            }
                        }
                    }
                }
                bin.AA = string.Join(" or ", ok);
            }
        }

        private void IdentifyMine(double v)
        {
            var myInfos = new List<MyInfo>();
            myInfos.Add(new MyInfo(0, "Exact match!"));
            myInfos.Add(new MyInfo(-48.128629, "Phosphorylation-Lysine: Probably reverse is the correct match"));
            myInfos.Add(new MyInfo(-76.134779, "Phosphorylation-Arginine: Probably reverse is the correct match"));
            myInfos.Add(new MyInfo(1.003, "1 MM"));
            myInfos.Add(new MyInfo(2.005, "2 MM"));
            myInfos.Add(new MyInfo(3.008, "3 MM"));
            myInfos.Add(new MyInfo(173.051055, "Acetylation + Methionine: Usually on protein N terminus"));
            myInfos.Add(new MyInfo(-91.009185, "neg Carbamidomethylation - H2S: Usually on cysteine."));
            myInfos.Add(new MyInfo(-32.008456, "oxidation and then loss of oxidized M side chain"));
            myInfos.Add(new MyInfo(-79.966331, "neg Phosphorylation."));
            myInfos.Add(new MyInfo(189.045969, "Carboxymethylated + Methionine. Usually on protein N terminus"));
            myInfos.Add(new MyInfo(356.20596, "Lysine+V+E or Lysine+L+D"));
            myInfos.Add(new MyInfo(239.126988, "Lysine+H(5) C(5) N O(2), possibly Nmethylmaleimide"));
            foreach (Bin bin in FinalBins)
            {
                bin.Mine = "";
                foreach (MyInfo myInfo in myInfos)
                    if (Math.Abs(myInfo.MassShift - bin.MassShift) <= v)
                        bin.Mine = myInfo.infostring;
            }
        }

        #endregion Private Methods

    }
}