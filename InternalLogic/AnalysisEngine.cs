using MassSpectrometry;
using OldInternalLogic;
using Proteomics;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

namespace InternalLogicEngineLayer
{
    public class AnalysisEngine : MyEngine
    {

        #region Private Fields

        private readonly ParentSpectrumMatch[][] newPsms;

        private readonly Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching;

        private readonly List<Protein> proteinList;

        private readonly List<MorpheusModification> variableModifications;

        private readonly List<MorpheusModification> fixedModifications;

        private readonly List<MorpheusModification> localizeableModifications;

        private readonly Protease protease;

        private readonly List<SearchMode> searchModes;

        private readonly IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;

        private readonly Tolerance fragmentTolerance;

        private readonly Action<BinTreeStructure, string> action1;

        private readonly Action<List<NewPsmWithFDR>, string> action2;

        private readonly Action<List<ProteinGroup>, string> action3;

        private readonly bool doParsimony;

        #endregion Private Fields

        #region Public Constructors

        public AnalysisEngine(ParentSpectrumMatch[][] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease, List<SearchMode> searchModes, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, Tolerance fragmentTolerance, Action<BinTreeStructure, string> action1, Action<List<NewPsmWithFDR>, string> action2, Action<List<ProteinGroup>, string> action3, bool doParsimony) : base(2)
        {
            this.doParsimony = doParsimony;
            this.newPsms = newPsms;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.localizeableModifications = localizeableModifications;
            this.protease = protease;
            this.searchModes = searchModes;
            this.myMsDataFile = myMsDataFile;
            this.fragmentTolerance = fragmentTolerance;
            this.action1 = action1;
            this.action2 = action2;
            this.action3 = action3;
        }

        #endregion Public Constructors

        #region Public Methods

        public static Dictionary<CompactPeptide, PeptideWithSetModifications> GetSingleMatchDictionary(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching)
        {
            // Right now very stupid, add the first decoy one, and if no decoy, add the first one
            var outDict = new Dictionary<CompactPeptide, PeptideWithSetModifications>();
            foreach (var kvp in fullSequenceToProteinPeptideMatching)
            {
                var k = kvp.Key;
                var val = kvp.Value;
                bool sawDecoy = false;
                foreach (var entry in val)
                {
                    if (entry.protein.isDecoy)
                    {
                        outDict[kvp.Key] = entry;
                        sawDecoy = true;
                        break;
                    }
                }
                if (sawDecoy == false)
                    outDict[kvp.Key] = kvp.Value.First();
            }

            return outDict;
        }

        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> ApplyProteinParsimony(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching)
        {
            // TODO**: protein group may have a shared peptide that is only shared within the group (peptide is unique to the group but shared within the group)
            // TODO**: how to handle indistinguishable proteins?
            // TODO**: rewrite parsimony function using protein group objects (more efficient)
            // TODO**: add proteins with unique peptides first, makes algo faster, use parsimony on remainder

            // finds unique peptides first (not really used for parsimony but used later for protein groups)
            foreach (var kvp in fullSequenceToProteinPeptideMatching)
            {
                List<Protein> proteinList = new List<Protein>();

                foreach(var peptide in kvp.Value)
                {
                    if(!proteinList.Contains(peptide.protein))
                    {
                        proteinList.Add(peptide.protein);
                    }
                }

                if(proteinList.Count == 1)
                {
                    kvp.Key.isUnique = true;
                }
            }

            // makes dictionary with proteins as keys and list of associated peptides as the value
            Dictionary<Protein, HashSet<CompactPeptide>> newDict = new Dictionary<Protein, HashSet<CompactPeptide>>();

            foreach (var kvp in fullSequenceToProteinPeptideMatching)
            {
                foreach (var virtualPeptide in kvp.Value)
                {
                    // find peptides associated with new protein
                    if (!newDict.ContainsKey(virtualPeptide.protein))
                    {
                        HashSet<CompactPeptide> peptides = new HashSet<CompactPeptide>();

                        foreach (var kvp1 in fullSequenceToProteinPeptideMatching)
                        {
                            foreach (var virtualPeptide2 in kvp1.Value)
                            {
                                // if the virtual peptides come from the same protein, add the peptide to the protein's list
                                if (virtualPeptide.protein == virtualPeptide2.protein)
                                {
                                    peptides.Add(kvp1.Key);
                                }
                            }
                        }

                        newDict.Add(virtualPeptide.protein, peptides);
                    }
                }
            }

            // greedy algorithm adds the next protein that will account for the most unaccounted-for peptides
            Dictionary<Protein, HashSet<CompactPeptide>> parsimonyDict = new Dictionary<Protein, HashSet<CompactPeptide>>();
            HashSet<CompactPeptide> usedPeptides = new HashSet<CompactPeptide>();
            HashSet<CompactPeptide> bestProteinPeptideList = new HashSet<CompactPeptide>();
            Protein bestProtein = null;
            int currentBestNumNewPeptides = -1;
            bool bestProteinHasOnePeptide = false;

            // as long as there are peptides that have not been accounted for, keep going
            while (currentBestNumNewPeptides != 0)
            {
                currentBestNumNewPeptides = 0;

                if (bestProteinHasOnePeptide == false)
                {
                    // attempt to find protein that best accounts for unaccounted-for peptides
                    foreach (var kvp in newDict)
                    {
                        int comparisonProteinNewPeptides = 0;

                        // determines number of unaccounted-for peptides for the current protein
                        foreach (CompactPeptide peptide in kvp.Value)
                        {
                            if (!usedPeptides.Contains(peptide))
                            {
                                comparisonProteinNewPeptides++;
                            }
                        }

                        // if the current protein is better than the best so far, current protein is the new best protein
                        if (comparisonProteinNewPeptides > currentBestNumNewPeptides)
                        {
                            bestProtein = kvp.Key;
                            bestProteinPeptideList = kvp.Value;
                            currentBestNumNewPeptides = comparisonProteinNewPeptides;
                        }
                    }

                    if (currentBestNumNewPeptides == 1)
                    {
                        bestProteinHasOnePeptide = true;
                    }

                    // adds the best protein if algo found unaccounted-for peptides
                    if (currentBestNumNewPeptides > 1)
                    {
                        parsimonyDict.Add(bestProtein, bestProteinPeptideList);

                        // the new best protein's peptides have been accounted for
                        foreach (CompactPeptide peptide in bestProteinPeptideList)
                        {
                            if (!usedPeptides.Contains(peptide))
                                usedPeptides.Add(peptide);
                        }
                    }
                }

                // best protein has one unaccounted-for peptide (stop searching for more than that peptide)
                else
                {
                    foreach (var kvp in newDict)
                    {
                        if (!parsimonyDict.ContainsKey(kvp.Key))
                        {
                            foreach (CompactPeptide peptide in kvp.Value)
                            {
                                if (!usedPeptides.Contains(peptide))
                                {
                                    bestProtein = kvp.Key;
                                    bestProteinPeptideList = kvp.Value;

                                    parsimonyDict.Add(bestProtein, bestProteinPeptideList);

                                    foreach (CompactPeptide peptide2 in bestProteinPeptideList)
                                    {
                                        if (!usedPeptides.Contains(peptide))
                                            usedPeptides.Add(peptide2);
                                    }

                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // add indistinguishable proteins to the protein group

            // build protein group after parsimony and match it to its peptide
            Dictionary<CompactPeptide, HashSet<Protein>> peptideProteinListMatch = new Dictionary<CompactPeptide, HashSet<Protein>>();

            foreach (var kvp in parsimonyDict)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!peptideProteinListMatch.ContainsKey(peptide))
                    {
                        HashSet<Protein> proteinList = new HashSet<Protein>();

                        foreach (var kvp1 in parsimonyDict)
                        {
                            if (kvp1.Value.Contains(peptide))
                            {
                                proteinList.Add(kvp1.Key);
                            }
                        }

                        peptideProteinListMatch.Add(peptide, proteinList);
                    }
                }
            }

            // constructs return dictionary (only use parsimony proteins for the new virtual peptide list)
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> answer = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            foreach (var kvp in parsimonyDict)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!answer.ContainsKey(peptide))
                    {
                        HashSet<PeptideWithSetModifications> oldVirtualPeptides = new HashSet<PeptideWithSetModifications>();
                        HashSet<PeptideWithSetModifications> newVirtualPeptides = new HashSet<PeptideWithSetModifications>();
                        HashSet<Protein> proteinList = new HashSet<Protein>();

                        // get the peptide's protein group after parsimony
                        peptideProteinListMatch.TryGetValue(peptide, out proteinList);

                        // find peptide's original (unparsimonious) virtual peptide matches
                        fullSequenceToProteinPeptideMatching.TryGetValue(peptide, out oldVirtualPeptides);

                        // get the virtual peptides that belong to the post-parsimony protein(s) only
                        foreach (var virtualPeptide in oldVirtualPeptides)
                        {
                            if (proteinList.Contains(virtualPeptide.protein))
                            {
                                newVirtualPeptides.Add(virtualPeptide);
                            }
                        }

                        // make new dictionary using only virtual peptides from parsimonious protein list
                        answer.Add(peptide, newVirtualPeptides);
                    }
                }
            }

            return answer;
        }

        public List<ProteinGroup> BuildProteinGroupsAndDoProteinFDR(List<NewPsmWithFDR> psmList, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching)
        {
            List<ProteinGroup> proteinGroups = new List<ProteinGroup>();

            Dictionary<Protein, List<NewPsmWithFDR>> dictionary1 = new Dictionary<Protein, List<NewPsmWithFDR>>();
            Dictionary<NewPsmWithFDR, HashSet<Protein>> dictionary2 = new Dictionary<NewPsmWithFDR, HashSet<Protein>>();

            HashSet<Protein> proteinGroupProteinList = new HashSet<Protein>();
            List<NewPsmWithFDR> proteinGroupPsmList = new List<NewPsmWithFDR>();

            // build protein list (get proteins associated with psm using the parsimony dictionary)
            foreach (var psm in psmList)
            {
                proteinGroupProteinList = new HashSet<Protein>();
                CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications);

                HashSet<PeptideWithSetModifications> virtualPeptides = new HashSet<PeptideWithSetModifications>();
                compactPeptideToProteinPeptideMatching.TryGetValue(peptide, out virtualPeptides);

                foreach (var virtualPeptide in virtualPeptides)
                {
                    if (!proteinGroupProteinList.Contains(virtualPeptide.protein))
                    {
                        proteinGroupProteinList.Add(virtualPeptide.protein);
                    }
                }

                dictionary2.Add(psm, proteinGroupProteinList);
            }

            // build psm list
            foreach (var kvp in dictionary2)
            {
                foreach (var protein in kvp.Value)
                {
                    if (!dictionary1.ContainsKey(protein))
                    {
                        proteinGroupPsmList = new List<NewPsmWithFDR>();

                        foreach (var kvp1 in dictionary2)
                        {
                            // if the kvp's have the same protein groups, add the newly found psm
                            if (kvp1.Value.SetEquals(kvp.Value))
                            {
                                proteinGroupPsmList.Add(kvp1.Key);
                                dictionary1.Add(protein, proteinGroupPsmList);
                            }
                        }

                        // build protein group from the dictionaries
                        proteinGroups.Add(new ProteinGroup(kvp.Value, proteinGroupPsmList, variableModifications, localizeableModifications));
                    }
                }
            }

            // protein groups are now properly constructed and scored

            // do fdr

            return proteinGroups;
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            status("Running analysis engine!");
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy
            status("Adding observed peptides to dictionary...");
            AddObservedPeptidesToDictionary();

            if (doParsimony)
            {
                status("Getting protein parsimony dictionary...");
                ApplyProteinParsimony(compactPeptideToProteinPeptideMatching);
            }

            //status("Getting single match just for FDR purposes...");
            //var fullSequenceToProteinSingleMatch = GetSingleMatchDictionary(compactPeptideToProteinPeptideMatching);

            List<NewPsmWithFDR>[] yeah = new List<NewPsmWithFDR>[searchModes.Count];
            List<ProteinGroup> proteinGroups = null;
            for (int j = 0; j < searchModes.Count; j++)
            {
                if (newPsms[j] != null)
                {
                    PSMwithTargetDecoyKnown[] psmsWithTargetDecoyKnown = new PSMwithTargetDecoyKnown[newPsms[0].Length];
                    Parallel.ForEach(Partitioner.Create(0, newPsms[0].Length), fff =>
                    {
                        for (int i = fff.Item1; i < fff.Item2; i++)
                        {
                            var huh = newPsms[j][i];
                            if (huh != null)
                                if (huh.Score >= 1)
                                    psmsWithTargetDecoyKnown[i] = new PSMwithTargetDecoyKnown(huh, compactPeptideToProteinPeptideMatching[huh.GetCompactPeptide(variableModifications, localizeableModifications)], fragmentTolerance, myMsDataFile);
                        }
                    });

                    var orderedPsmsWithPeptides = psmsWithTargetDecoyKnown.Where(b => b != null).OrderByDescending(b => b.Score);

                    status("Running FDR analysis...");
                    var orderedPsmsWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides);
                    var limitedpsms_with_fdr = orderedPsmsWithFDR.Where(b => (b.QValue <= 0.01)).ToList();
                    if (limitedpsms_with_fdr.Count(b => !b.isDecoy) > 0)
                    {
                        status("Running histogram analysis...");
                        var hm = MyAnalysis(limitedpsms_with_fdr);
                        action1(hm, searchModes[j].FileNameAddition);
                    }

                    if (doParsimony)
                    {
                        status("Building protein groups and doing FDR...");
                        proteinGroups = BuildProteinGroupsAndDoProteinFDR(orderedPsmsWithFDR, compactPeptideToProteinPeptideMatching);
                    }

                    action2(orderedPsmsWithFDR, searchModes[j].FileNameAddition);
                    action3(proteinGroups, searchModes[j].FileNameAddition);

                    yeah[j] = orderedPsmsWithFDR;
                }
            }

            return new AnalysisResults(this, yeah, proteinGroups);
        }

        #endregion Protected Methods

        #region Private Methods

        private static void IdentifyAAsInCommon(BinTreeStructure myTreeStructure)
        {
            foreach (Bin bin in myTreeStructure.finalBins)
            {
                bin.AAsInCommon = new Dictionary<char, int>();
                foreach (var hehe in bin.uniquePSMs.Values.Where(b => !b.Item3.isDecoy))
                {
                    var chars = new HashSet<char>();
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

        private static void IdentifyMods(BinTreeStructure myTreeStructure)
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

        private static void IdentifyResidues(BinTreeStructure myTreeStructure)
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

        private static void IdentifyUnimodBins(BinTreeStructure myTreeStructure, double v)
        {
            foreach (var bin in myTreeStructure.finalBins)
            {
                var ok = new HashSet<string>();
                var okformula = new HashSet<string>();
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

        private static void IdentifyUniprotBins(BinTreeStructure myTreeStructure, double v)
        {
            foreach (var bin in myTreeStructure.finalBins)
            {
                var ok = new HashSet<string>();
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

        private static void IdentifyCombos(BinTreeStructure myTreeStructure, double v)
        {
            double totalTargetCount = myTreeStructure.finalBins.Select(b => b.CountTarget).Sum();
            var ok = new HashSet<Tuple<double, double, double>>();
            foreach (var bin in myTreeStructure.finalBins.Where(b => Math.Abs(b.MassShift) > v))
                foreach (var bin2 in myTreeStructure.finalBins.Where(b => Math.Abs(b.MassShift) > v))
                    if (bin.CountTarget * bin2.CountTarget >= totalTargetCount * 3)
                        ok.Add(new Tuple<double, double, double>(bin.MassShift, bin2.MassShift, Math.Min(bin.CountTarget, bin2.CountTarget)));

            foreach (var bin in myTreeStructure.finalBins)
            {
                var okk = new HashSet<string>();
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

        private static void IdentifyAA(BinTreeStructure myTreeStructure, double v)
        {
            foreach (var bin in myTreeStructure.finalBins)
            {
                var ok = new HashSet<string>();
                for (char c = 'A'; c <= 'Z'; c++)
                {
                    Residue residue;
                    if (Residue.TryGetResidue(c, out residue))
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
                            Residue residueCC;
                            if (Residue.TryGetResidue(cc, out residueCC))
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

        private static void IdentifyMine(BinTreeStructure myTreeStructure, double v)
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

        private static BinTreeStructure MyAnalysis(List<NewPsmWithFDR> limitedpsms_with_fdr)
        {
            var myTreeStructure = new BinTreeStructure();
            myTreeStructure.GenerateBins(limitedpsms_with_fdr, 0.003);
            myTreeStructure.AddToBins(limitedpsms_with_fdr);

            IdentifyUnimodBins(myTreeStructure, 0.003);
            IdentifyUniprotBins(myTreeStructure, 0.003);
            IdentifyAA(myTreeStructure, 0.003);

            IdentifyCombos(myTreeStructure, 0.003);

            IdentifyResidues(myTreeStructure);

            IdentifyMods(myTreeStructure);

            IdentifyAAsInCommon(myTreeStructure);

            IdentifyMine(myTreeStructure, 0.003);

            return myTreeStructure;
        }

        private static List<NewPsmWithFDR> DoFalseDiscoveryRateAnalysis(IEnumerable<PSMwithTargetDecoyKnown> items)
        {
            var ids = new List<NewPsmWithFDR>();

            int cumulative_target = 0;
            int cumulative_decoy = 0;
            foreach (PSMwithTargetDecoyKnown item in items)
            {
                var isDecoy = item.isDecoy;
                if (isDecoy)
                    cumulative_decoy++;
                else
                    cumulative_target++;
                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                ids.Add(new NewPsmWithFDR(item, cumulative_target, cumulative_decoy, temp_q_value));
            }

            double min_q_value = double.PositiveInfinity;
            for (int i = ids.Count - 1; i >= 0; i--)
            {
                NewPsmWithFDR id = ids[i];
                if (id.QValue > min_q_value)
                    id.QValue = min_q_value;
                else if (id.QValue < min_q_value)
                    min_q_value = id.QValue;
            }

            return ids;
        }

        private void AddObservedPeptidesToDictionary()
        {
            status("Adding new keys to peptide dictionary...");
            foreach (var psmListForAspecificSerchMode in newPsms)
            {
                if (psmListForAspecificSerchMode != null)
                    foreach (var psm in psmListForAspecificSerchMode)
                    {
                        if (psm != null)
                        {
                            var cp = psm.GetCompactPeptide(variableModifications, localizeableModifications);
                            if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp))
                                compactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());
                        }
                    }
            }

            int totalProteins = proteinList.Count;

            status("Adding possible sources to peptide dictionary...");

            int proteinsSeen = 0;
            int old_progress = 0;

            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> local = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());

                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, 2, InitiatorMethionineBehavior.Variable).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length == 1 || peptide.Length > 252)
                            continue;

                        peptide.SetFixedModifications(fixedModifications);

                        var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(variableModifications, 4096, 3).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            HashSet<PeptideWithSetModifications> v;
                            if (local.TryGetValue(new CompactPeptide(yyy, variableModifications, localizeableModifications), out v))
                            {
                                v.Add(yyy);
                            }
                        }
                    }
                }
                lock (compactPeptideToProteinPeptideMatching)
                {
                    foreach (var ye in local)
                    {
                        HashSet<PeptideWithSetModifications> v;
                        if (compactPeptideToProteinPeptideMatching.TryGetValue(ye.Key, out v))
                        {
                            foreach (var huh in ye.Value)
                                v.Add(huh);
                        }
                    }
                    proteinsSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In adding possible sources to peptide dictionary loop"));
                        old_progress = new_progress;
                    }
                }
            });
        }

        #endregion Private Methods

    }
}