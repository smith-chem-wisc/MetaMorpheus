using MassSpectrometry;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Proteomics.Fragmentation;

namespace EngineLayer.HistogramAnalysis
{
    public class BinTreeStructure
    {
        private const int MinAdditionalPsmsInBin = 1;

        public List<Bin> FinalBins { get; private set; }

        public void GenerateBins(List<PeptideSpectralMatch> targetAndDecoyMatches, double dc)
        {
            List<double> listOfMassShifts = targetAndDecoyMatches.Where(b => b.PeptideMonisotopicMass.HasValue).Select(b => b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value).OrderBy(b => b).ToList();
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

                    sigma[i] = GetSigma(thisMassShift, thisP, i, listOfMassShifts, p);
                }
            }

            var listokbin = new List<OkBin>();
            for (int i = 0; i < sigma.Count(); i++)
                listokbin.Add(new OkBin(listOfMassShifts[i], sigma[i], p[i]));

            var prelimBins = new HashSet<double>();
            foreach (OkBin okbin in listokbin.OrderByDescending(b => b.P))
            {
                if (okbin.Sigma < dc || okbin.P < MinAdditionalPsmsInBin)
                    continue;
                bool add = true;
                foreach (double a in prelimBins)
                {
                    if (Math.Abs(okbin.MassShift - a) <= dc)
                    {
                        add = false;
                        break;
                    }
                }
                if (add)
                    prelimBins.Add(okbin.MassShift);
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
                    if (targetAndDecoyMatches[i].PeptideMonisotopicMass.HasValue && Math.Abs(targetAndDecoyMatches[i].ScanPrecursorMass - targetAndDecoyMatches[i].PeptideMonisotopicMass.Value - bin.MassShift) <= dc)
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
            
            IdentifyFracWithSingle();
            IdentifyMedianLength();
        }

        private static double GetSigma(double thisMassShift, int thisP, int i, List<double> listOfMassShifts, int[] p)
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
        
        private void IdentifyFracWithSingle()
        {
            foreach (Bin bin in FinalBins)
            {
                var numTarget = bin.UniquePSMs.Values.Count(b => !b.Item3.IsDecoy);
                if (numTarget > 0)
                    bin.FracWithSingle = (double)bin.UniquePSMs.Values.Count(b => !b.Item3.IsDecoy && b.Item3.NumDifferentMatchingPeptides == 1) / numTarget;
            }
        }

        private void IdentifyMedianLength()
        {
            foreach (Bin bin in FinalBins)
            {
                var numTarget = bin.UniquePSMs.Values.Count(b => !b.Item3.IsDecoy);
                if (numTarget > 0)
                    bin.MedianLength = Statistics.Median(bin.UniquePSMs.Values.Where(b => !b.Item3.IsDecoy).Where(b => b.Item3.PeptideLength.HasValue).Select(b => (double)b.Item3.PeptideLength.Value));
            }
        }

        private void IdentifyAAsInCommon()
        {
            foreach (Bin bin in FinalBins)
            {
                bin.AAsInCommon = new Dictionary<char, int>();
                foreach (var hehe in bin.UniquePSMs.Values.Where(b => !b.Item3.IsDecoy))
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
                bin.ModsInCommon = new Dictionary<string, int>();
                foreach (var hehe in bin.UniquePSMs.Values.Where(b => !b.Item3.IsDecoy))
                {
                    int inModLevel = 0;
                    StringBuilder currentMod = new StringBuilder();
                    HashSet<string> modsHere = new HashSet<string>();
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
                                if (!currentMod.ToString().StartsWith("Common Fixed:"))
                                    modsHere.Add(currentMod.ToString());
                                currentMod.Clear();
                            }
                            continue;
                        }
                        if (inModLevel > 0)
                            currentMod.Append(ye);
                    }
                    foreach (var modInHS in modsHere)
                    {
                        if (bin.ModsInCommon.ContainsKey(modInHS))
                            bin.ModsInCommon[modInHS]++;
                        else
                            bin.ModsInCommon.Add(modInHS, 1);
                    }
                }
            }
        }

        private void IdentifyResidues()
        {
            foreach (Bin bin in FinalBins)
            {
                bin.IdentifyResidues();
            }
        }

        private void IdentifyUnimodBins(double v)
        {
            foreach (var bin in FinalBins)
            {
                bin.IdentifyUnimodBins(v);
            }
        }

        private void IdentifyUniprotBins(double v)
        {
            foreach (var bin in FinalBins)
            {
                bin.IdentifyUniprotBins(v);
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
                bin.IdentifyCombos(v, ok);
            }
        }

        private void IdentifyAA(double v)
        {
            foreach (var bin in FinalBins)
            {
                bin.IdentifyAA(v);
            }
        }

        private void IdentifyMine(double v)
        {
            var myInfos = new List<MyInfo>
                    {
                        new MyInfo(0, "Exact match!"),
                        new MyInfo(-48.128629, "Phosphorylation-Lysine: Probably reverse is the correct match"),
                        new MyInfo(-76.134779, "Phosphorylation-Arginine: Probably reverse is the correct match"),
                        new MyInfo(1.0029, "1 MM"),
                        new MyInfo(2.0052, "2 MM"),
                        new MyInfo(3.0077, "3 MM"),
                        new MyInfo(173.051055, "Acetylation + Methionine: Usually on protein N terminus"),
                        new MyInfo(-91.009185, "neg Carbamidomethylation - H2S: Usually on cysteine."),
                        new MyInfo(-32.008456, "oxidation and then loss of oxidized M side chain"),
                        new MyInfo(-79.966331, "neg Phosphorylation."),
                        new MyInfo(189.045969, "Carboxymethylated + Methionine. Usually on protein N terminus"),
                        new MyInfo(356.20596, "Lysine+V+E or Lysine+L+D"),
                        new MyInfo(239.126988, "Lysine+H(5) C(5) N O(2), possibly Nmethylmaleimide"),
                        new MyInfo(-105.02484, "Methionine loss then acetaldehyde"),
                        new MyInfo(52.911464, "Fe[III]"),
                        new MyInfo(71.000729, "H C2 N O2"),
                        new MyInfo(50.000394, "H2 O3")
                    };
            foreach (Bin bin in FinalBins)
            {
                bin.Mine = "";
                foreach (MyInfo myInfo in myInfos)
                    if (Math.Abs(myInfo.MassShift - bin.MassShift) <= v)
                        bin.Mine = myInfo.infostring;
            }
        }
    }
}