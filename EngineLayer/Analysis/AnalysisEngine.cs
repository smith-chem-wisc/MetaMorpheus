using MassSpectrometry;
using Proteomics;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Analysis
{
    public class AnalysisEngine : MyEngine
    {

        #region Private Fields

        private const int max_mods_for_peptide = 3;
        private readonly double binTol;
        private readonly int maximumMissedCleavages;
        private readonly int maxModIsoforms;
        private readonly PsmParent[][] newPsms;
        private readonly List<Protein> proteinList;
        private readonly List<MetaMorpheusModification> variableModifications;
        private readonly List<MetaMorpheusModification> fixedModifications;
        private readonly List<MetaMorpheusModification> localizeableModifications;
        private readonly Protease protease;
        private readonly List<SearchMode> searchModes;
        private readonly IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
        private readonly Tolerance fragmentTolerance;
        private readonly Action<BinTreeStructure, string> writeHistogramPeaksAction;
        private readonly Action<List<NewPsmWithFdr>, string> writePsmsAction;
        private readonly Action<List<ProteinGroup>, string> action3;
        private readonly bool doParsimony;
        private readonly bool doHistogramAnalysis;
        private readonly List<ProductType> lp;
        private Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching;

        #endregion Private Fields

        #region Public Constructors

        public AnalysisEngine(PsmParent[][] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<MetaMorpheusModification> variableModifications, List<MetaMorpheusModification> fixedModifications, List<MetaMorpheusModification> localizeableModifications, Protease protease, List<SearchMode> searchModes, IMsDataFile<IMzSpectrum<MzPeak>> myMSDataFile, Tolerance fragmentTolerance, Action<BinTreeStructure, string> action1, Action<List<NewPsmWithFdr>, string> action2, Action<List<ProteinGroup>, string> action3, bool doParsimony, int maximumMissedCleavages, int maxModIsoforms, bool doHistogramAnalysis, List<ProductType> lp, double binTol) : base(2)
        {
            this.doParsimony = doParsimony;
            this.doHistogramAnalysis = doHistogramAnalysis;
            this.newPsms = newPsms;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.localizeableModifications = localizeableModifications;
            this.protease = protease;
            this.searchModes = searchModes;
            this.myMsDataFile = myMSDataFile;
            this.fragmentTolerance = fragmentTolerance;
            this.writeHistogramPeaksAction = action1;
            this.writePsmsAction = action2;
            this.action3 = action3;
            this.maximumMissedCleavages = maximumMissedCleavages;
            this.maxModIsoforms = maxModIsoforms;
            this.lp = lp;
            this.binTol = binTol;
        }

        #endregion Public Constructors

        #region Public Methods

        public Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> ApplyProteinParsimony(out List<ProteinGroup> proteinGroups)
        {
            Status("Applying protein parsimony...");

            var uniquePeptides = new HashSet<CompactPeptide>();
            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<CompactPeptide>>();
            var parsimonyDict = new Dictionary<Protein, HashSet<CompactPeptide>>();
            var proteinsWithUniquePeptides = new HashSet<Protein>();
            var usedBaseSequences = new HashSet<string>();

            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                // if a peptide is associated with a decoy protein, remove all target protein associations with the peptide
                var decoyPeptidesHere = new HashSet<PeptideWithSetModifications>(kvp.Value.Where(p => p.Protein.IsDecoy));

                if (decoyPeptidesHere.Any())
                {
                    var pepWithSetModsToRemove = new HashSet<PeptideWithSetModifications>(kvp.Value.Except(decoyPeptidesHere));

                    foreach (var pepWithSetMods in pepWithSetModsToRemove)
                        kvp.Value.Remove(pepWithSetMods);
                }

                // finds unique peptides (peptides that can belong to only one protein)
                var proteinListHere = new HashSet<Protein>(kvp.Value.Select(p => p.Protein));

                if (proteinListHere.Count == 1)
                    uniquePeptides.Add(kvp.Key);
            }

            // makes dictionary with proteins as keys and list of associated peptides as the value (makes parsimony algo easier)
            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!proteinToPeptidesMatching.ContainsKey(peptide.Protein))
                    {
                        var peptides = new HashSet<CompactPeptide>();
                        peptides.Add(kvp.Key);
                        proteinToPeptidesMatching.Add(peptide.Protein, peptides);
                    }
                    else
                    {
                        HashSet<CompactPeptide> peptides;
                        proteinToPeptidesMatching.TryGetValue(peptide.Protein, out peptides);
                        peptides.Add(kvp.Key);
                    }
                }
            }

            // add proteins with unique peptides to the parsimony dictionary before applying parsimony algorithm (more efficient)
            foreach (var kvp in proteinToPeptidesMatching)
            {
                var uniquePeptidesHere = kvp.Value.Intersect(uniquePeptides);

                if (uniquePeptidesHere.Any())
                {
                    parsimonyDict.Add(kvp.Key, kvp.Value);
                    proteinsWithUniquePeptides.Add(kvp.Key);

                    foreach (var peptide in kvp.Value)
                    {
                        string peptideBaseSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));
                        usedBaseSequences.Add(peptideBaseSequence);
                    }
                }
            }

            // build protein list for each peptide before parsimony has been applied
            var peptideBaseSeqProteinListMatch = new Dictionary<string, HashSet<Protein>>();
            foreach (var kvp in proteinToPeptidesMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    var proteinListHere = new HashSet<Protein>();
                    string peptideBaseSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    if (!peptideBaseSeqProteinListMatch.ContainsKey(peptideBaseSequence))
                    {
                        proteinListHere.Add(kvp.Key);
                        peptideBaseSeqProteinListMatch.Add(peptideBaseSequence, proteinListHere);
                    }
                    else
                    {
                        peptideBaseSeqProteinListMatch.TryGetValue(peptideBaseSequence, out proteinListHere);
                        proteinListHere.Add(kvp.Key);
                    }
                }
            }

            // greedy algorithm adds the next protein that will account for the most unaccounted-for peptides
            bool currentBestPeptidesIsOne = false;
            int currentBestNumNewPeptides;

            do
            {
                currentBestNumNewPeptides = 0;
                Protein bestProtein = null;

                if (!currentBestPeptidesIsOne)
                {
                    // attempt to find protein that best accounts for unaccounted-for peptides
                    foreach (var kvp in proteinToPeptidesMatching)
                    {
                        if (!parsimonyDict.ContainsKey(kvp.Key))
                        {
                            var baseSeqs = new HashSet<string>(kvp.Value.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)));
                            var comparisonProteinNewPeptideBaseSeqs = new HashSet<string>(baseSeqs.Except(usedBaseSequences));
                            int comparisonProteinNewPeptides = comparisonProteinNewPeptideBaseSeqs.Count;

                            if (comparisonProteinNewPeptides >= currentBestNumNewPeptides)
                            {
                                // if the current protein is better than the best so far, current protein is the new best protein
                                bestProtein = kvp.Key;
                                currentBestNumNewPeptides = comparisonProteinNewPeptides;
                            }
                        }
                    }

                    if (currentBestNumNewPeptides < 2)
                        currentBestPeptidesIsOne = true;

                    // adds the best protein for this iteration
                    if (!currentBestPeptidesIsOne)
                    {
                        HashSet<CompactPeptide> bestProteinPeptideList;
                        proteinToPeptidesMatching.TryGetValue(bestProtein, out bestProteinPeptideList);
                        var baseSeqs = new HashSet<string>(bestProteinPeptideList.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)));
                        var bestProteinNewPeptideBaseSeqs = new HashSet<string>(baseSeqs.Except(usedBaseSequences));
                        var proteinsWithNewSeqs = new HashSet<Protein>();

                        // find all proteins that have the new base sequences
                        foreach (var newBaseSeq in bestProteinNewPeptideBaseSeqs)
                        {
                            HashSet<Protein> proteinsWithThisBaseSeq;
                            peptideBaseSeqProteinListMatch.TryGetValue(newBaseSeq, out proteinsWithThisBaseSeq);
                            
                            foreach(var protein in proteinsWithThisBaseSeq)
                            {
                                HashSet<CompactPeptide> t;
                                proteinToPeptidesMatching.TryGetValue(protein, out t);
                                var temp = new HashSet<string>(t.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)));
                                if (baseSeqs.IsSubsetOf(temp))
                                    proteinsWithNewSeqs.Add(protein);
                            }
                        }

                        // multiple proteins have the same new base seqs - pick the one with the most total peptides
                        if (proteinsWithNewSeqs.Count > 1)
                        {
                            var proteinsWithNumBaseSeqs = new Dictionary<Protein, int>();

                            // count how many peptides each protein has
                            foreach (var protein in proteinsWithNewSeqs)
                            {
                                HashSet<CompactPeptide> peps;
                                proteinToPeptidesMatching.TryGetValue(protein, out peps);
                                baseSeqs = new HashSet<string>(peps.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)));
                                proteinsWithNumBaseSeqs.Add(protein, baseSeqs.Count);
                            }

                            // pick the protein with the most peptides
                            var temp = proteinsWithNumBaseSeqs.OrderByDescending(b => b.Value);
                            bestProtein = temp.First().Key;
                        }

                        proteinToPeptidesMatching.TryGetValue(bestProtein, out bestProteinPeptideList);
                        parsimonyDict.Add(bestProtein, bestProteinPeptideList);

                        foreach (var peptide in bestProteinPeptideList)
                        {
                            string peptideBaseSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                            usedBaseSequences.Add(peptideBaseSequence);
                        }
                    }
                }

                // best protein has one unaccounted-for peptide (stop searching for more than that peptide)
                else
                {
                    foreach (var kvp in proteinToPeptidesMatching)
                    {
                        if (!parsimonyDict.ContainsKey(kvp.Key))
                        {
                            foreach (var peptide in kvp.Value)
                            {
                                string peptideBaseSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                                if (!usedBaseSequences.Contains(peptideBaseSequence))
                                {
                                    bestProtein = kvp.Key;

                                    HashSet<Protein> proteins;
                                    peptideBaseSeqProteinListMatch.TryGetValue(peptideBaseSequence, out proteins);

                                    // multiple proteins have the same new base seq - pick the one with the most total peptides
                                    if (proteins.Count > 1)
                                    {
                                        var proteinsWithNumBaseSeqs = new Dictionary<Protein, int>();

                                        // count how many peptides each protein has
                                        foreach (var protein in proteins)
                                        {
                                            HashSet<CompactPeptide> peps;
                                            proteinToPeptidesMatching.TryGetValue(protein, out peps);
                                            var baseSeqs = new HashSet<string>(peps.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)));
                                            proteinsWithNumBaseSeqs.Add(protein, baseSeqs.Count);
                                        }

                                        // pick the protein with the most peptides
                                        var temp = proteinsWithNumBaseSeqs.OrderByDescending(b => b.Value);
                                        bestProtein = temp.First().Key;
                                    }

                                    parsimonyDict.Add(bestProtein, kvp.Value);
                                    usedBaseSequences.Add(peptideBaseSequence);
                                    break;
                                }
                            }
                        }
                    }
                }
            } while (currentBestNumNewPeptides != 0);

            // build protein group after parsimony (each group only has 1 protein at this point)
            proteinGroups = new List<ProteinGroup>();
            foreach (var kvp in parsimonyDict)
            {
                var uniquePeptidesHere = new HashSet<CompactPeptide>(kvp.Value.Intersect(uniquePeptides));
                proteinGroups.Add(new ProteinGroup(new HashSet<Protein>() { kvp.Key }, kvp.Value, uniquePeptidesHere));
            }

            // grab indistinguishable proteins ("if" conditions are to narrow search space)
            foreach (var proteinGroup in proteinGroups)
            {
                if (!proteinGroup.TotalUniquePeptideList.Any())
                {
                    foreach (var kvp in proteinToPeptidesMatching)
                    {
                        if (!proteinsWithUniquePeptides.Contains(kvp.Key))
                        {
                            // prevents looking at itself
                            if (!parsimonyDict.ContainsKey(kvp.Key))
                            {
                                if (kvp.Value.Count == proteinGroup.TotalPeptideList.Count)
                                {
                                    if (kvp.Value.SetEquals(proteinGroup.TotalPeptideList))
                                    {
                                        proteinGroup.Proteins.Add(kvp.Key);
                                        parsimonyDict.Add(kvp.Key, kvp.Value);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // build protein list for each peptide after parsimony has been applied (helps make return dictionary)
            var peptideProteinListMatch = new Dictionary<CompactPeptide, HashSet<Protein>>();
            foreach (var kvp in parsimonyDict)
            {
                foreach (var peptide in kvp.Value)
                {
                    var proteinListHere = new HashSet<Protein>();

                    if (!peptideProteinListMatch.ContainsKey(peptide))
                    {
                        proteinListHere.Add(kvp.Key);
                        peptideProteinListMatch.Add(peptide, proteinListHere);
                    }
                    else
                    {
                        peptideProteinListMatch.TryGetValue(peptide, out proteinListHere);
                        proteinListHere.Add(kvp.Key);
                    }
                }
            }

            // constructs return dictionary (only use parsimony proteins for the new PeptideWithSetModifications list)
            var answer = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            foreach (var kvp in parsimonyDict)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!answer.ContainsKey(peptide))
                    {
                        // find CompactPeptide's original (unparsimonious) peptide matches
                        HashSet<PeptideWithSetModifications> oldPepsWithSetMods;
                        compactPeptideToProteinPeptideMatching.TryGetValue(peptide, out oldPepsWithSetMods);

                        // get the CompactPeptide's protein list after parsimony
                        HashSet<Protein> proteinListHere;
                        peptideProteinListMatch.TryGetValue(peptide, out proteinListHere);

                        // get the peptides that belong to the post-parsimony protein(s) only
                        var newPeptides = new HashSet<PeptideWithSetModifications>(oldPepsWithSetMods.Where(p => proteinListHere.Contains(p.Protein)));

                        answer.Add(peptide, newPeptides);
                    }
                }
            }

            Status("Finished Parsimony");

            compactPeptideToProteinPeptideMatching = answer;

            // returns for test class (TODO**: remove)
            return answer;
        }

        public void ScoreProteinGroups(List<ProteinGroup> proteinGroups, List<NewPsmWithFdr> psmList)
        {
            Status("Scoring protein groups...");

            Dictionary<CompactPeptide, HashSet<ProteinGroup>> peptideToProteinGroupMatching = new Dictionary<CompactPeptide, HashSet<ProteinGroup>>();
            HashSet<CompactPeptide> allRazorPeptides = new HashSet<CompactPeptide>();
            HashSet<ProteinGroup> proteinGroupsToRemove = new HashSet<ProteinGroup>();

            // add each protein group's psm's
            Dictionary<CompactPeptide, NewPsmWithFdr> psmToCompactPeptideMatching = new Dictionary<CompactPeptide, NewPsmWithFdr>();
            foreach (var psm in psmList)
            {
                CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications, fixedModifications);
                if (!psmToCompactPeptideMatching.ContainsKey(peptide))
                    psmToCompactPeptideMatching.Add(peptide, psm);
            }

            foreach (var proteinGroup in proteinGroups)
            {
                foreach (var peptide in proteinGroup.TotalPeptideList)
                {
                    // build PSM list for scoring
                    NewPsmWithFdr psm;
                    psmToCompactPeptideMatching.TryGetValue(peptide, out psm);
                    if (psm != null)
                        proteinGroup.TotalPsmList.Add(psm);

                    // build PeptideWithSetMod list to calc sequence coverage
                    HashSet<PeptideWithSetModifications> peptidesWithSetMods = null;
                    compactPeptideToProteinPeptideMatching.TryGetValue(peptide, out peptidesWithSetMods);
                    foreach (var pep in peptidesWithSetMods)
                    {
                        proteinGroup.TotalPeptideWithSetModsList.Add(pep);
                    }
                }
            }

            // score the group and get number of protein groups per peptide
            foreach (var proteinGroup in proteinGroups)
            {
                // score the group (scoring algorithm defined in the ProteinGroup class)
                proteinGroup.ScoreThisProteinGroup(variableModifications, localizeableModifications, fixedModifications);

                // for finding razor peptides later
                foreach (var peptide in proteinGroup.StrictPeptideList)
                {
                    HashSet<ProteinGroup> proteinGroupsHere = new HashSet<ProteinGroup>();
                    if (peptideToProteinGroupMatching.ContainsKey(peptide))
                    {
                        peptideToProteinGroupMatching.TryGetValue(peptide, out proteinGroupsHere);
                        proteinGroupsHere.Add(proteinGroup);
                    }
                    else
                    {
                        proteinGroupsHere.Add(proteinGroup);
                        peptideToProteinGroupMatching.Add(peptide, proteinGroupsHere);
                    }
                }
            }

            /*
            // TODO**
            // merge indistinguishable protein groups after scoring
            var pg = proteinGroups.OrderByDescending(p => p.proteinGroupScore).ToArray();
            for(int i = 0; i < (pg.Length - 1); i++)
            {
                if (pg[i].proteinGroupScore == pg[i + 1].proteinGroupScore)
                {
                    if (pg[i].StrictPeptideList.SetEquals(pg[i + 1].StrictPeptideList))
                    {
                        pg[i].MergeProteinGroupWith(pg[i + 1]);
                    }
                }
            }
            */

            foreach (var proteinGroup in proteinGroups)
            {
                if (proteinGroup.proteinGroupScore == 0)
                    proteinGroupsToRemove.Add(proteinGroup);
            }

            // remove empty protein groups (peptides were too poor quality and group doesn't exist anymore)
            foreach (var proteinGroup in proteinGroupsToRemove)
            {
                proteinGroups.Remove(proteinGroup);
            }

            // build razor peptide list (peptides that have >1 protein groups in the final, scored protein group list)
            foreach (var kvp in peptideToProteinGroupMatching)
            {
                if (kvp.Value.Count > 1)
                {
                    allRazorPeptides.Add(kvp.Key);
                }
            }

            foreach (var proteinGroup in proteinGroups)
            {
                foreach (var peptide in proteinGroup.TotalPeptideList)
                {
                    // build razor peptide list for each protein group
                    if (allRazorPeptides.Contains(peptide))
                    {
                        // TODO**
                        // if the razor pep is associated with >1 protein group, it's a razor only for the group with the most ID'd peptides
                        // var sortedProteinGroups =
                        proteinGroup.StrictRazorPeptideList.Add(peptide);
                    }


                }

                // calculate sequence coverage for each protein in the group
                proteinGroup.CalculateSequenceCoverage();
            }
        }

        public List<ProteinGroup> DoProteinFdr(List<ProteinGroup> proteinGroups)
        {
            Status("Calculating protein FDR...");

            // order protein groups by score
            var sortedProteinGroups = proteinGroups.OrderByDescending(b => b.proteinGroupScore).ToList();

            // do fdr
            int cumulativeTarget = 0;
            int cumulativeDecoy = 0;
            foreach (var proteinGroup in sortedProteinGroups)
            {
                if (proteinGroup.isDecoy)
                {
                    cumulativeDecoy++;
                }
                else
                {
                    cumulativeTarget++;
                }

                proteinGroup.cumulativeTarget = cumulativeTarget;
                proteinGroup.cumulativeDecoy = cumulativeDecoy;
                proteinGroup.QValue = ((double)cumulativeDecoy / (cumulativeTarget + cumulativeDecoy));
            }

            return sortedProteinGroups;
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            Status("Running analysis engine!");
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy
            Status("Adding observed peptides to dictionary...");
            AddObservedPeptidesToDictionary();
            List<ProteinGroup>[] proteinGroups = new List<ProteinGroup>[searchModes.Count];
            if (doParsimony)
            {
                // TODO**: make this faster (only apply parsimony once but make multiple instances of the same ProteinGroups
                for (int i = 0; i < searchModes.Count; i++)
                    ApplyProteinParsimony(out proteinGroups[i]);
            }

            //status("Getting single match just for FDR purposes...");
            //var fullSequenceToProteinSingleMatch = GetSingleMatchDictionary(compactPeptideToProteinPeptideMatching);
            List<NewPsmWithFdr>[] allResultingIdentifications = new List<NewPsmWithFdr>[searchModes.Count];

            for (int j = 0; j < searchModes.Count; j++)
            {
                if (newPsms[j] != null)
                {
                    PsmWithMultiplePossiblePeptides[] psmsWithProteinHashSet = new PsmWithMultiplePossiblePeptides[newPsms[0].Length];
                    for (int i = 0; i < newPsms[0].Length; i++)
                    {
                        var huh = newPsms[j][i];
                        if (huh != null && huh.score >= 1)
                            psmsWithProteinHashSet[i] = new PsmWithMultiplePossiblePeptides(huh, compactPeptideToProteinPeptideMatching[huh.GetCompactPeptide(variableModifications, localizeableModifications, fixedModifications)], fragmentTolerance, myMsDataFile, lp);
                    }

                    var orderedPsmsWithPeptides = psmsWithProteinHashSet.Where(b => b != null).OrderByDescending(b => b.Score);

                    Status("Running FDR analysis...");
                    var orderedPsmsWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides, searchModes[j]);
                    writePsmsAction(orderedPsmsWithFDR, searchModes[j].FileNameAddition);

                    // This must come before the unique peptide identifications..
                    // Or instead!!!
                    if (doHistogramAnalysis)
                    {
                        var limitedpsms_with_fdr = orderedPsmsWithFDR.Where(b => (b.qValue <= 0.01)).ToList();
                        if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                        {
                            Status("Running histogram analysis...");
                            var hm = MyAnalysis(limitedpsms_with_fdr, binTol);
                            writeHistogramPeaksAction(hm, searchModes[j].FileNameAddition);
                        }
                    }
                    else
                    {
                        Status("Running FDR analysis on unique peptides...");
                        var uniquePsmsWithFdr = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides.Distinct(new SequenceComparer()), searchModes[j]);
                        writePsmsAction(uniquePsmsWithFdr, "uniquePeptides" + searchModes[j].FileNameAddition);
                    }

                    if (doParsimony)
                    {
                        ScoreProteinGroups(proteinGroups[j], orderedPsmsWithFDR);
                        proteinGroups[j] = DoProteinFdr(proteinGroups[j]);

                        action3(proteinGroups[j], searchModes[j].FileNameAddition);
                    }

                    allResultingIdentifications[j] = orderedPsmsWithFDR;
                }
            }
            return new AnalysisResults(this, allResultingIdentifications, proteinGroups);
        }

        #endregion Protected Methods

        #region Private Methods

        private static void OverlappingIonSequences(BinTreeStructure myTreeStructure)
        {
            foreach (Bin bin in myTreeStructure.FinalBins)
            {
                foreach (var hm in bin.uniquePSMs.Where(b => !b.Value.Item3.IsDecoy))
                {
                    var ya = hm.Value.Item3.thisPSM.newPsm.matchedIonsList;
                    if (ya.ContainsKey(ProductType.B) && ya.ContainsKey(ProductType.Y) && ya[ProductType.B].Any(b => b > 0) && ya[ProductType.Y].Any(b => b > 0))
                        if (ya[ProductType.B].Last(b => b > 0) + ya[ProductType.Y].Last(b => b > 0) > hm.Value.Item3.thisPSM.PeptideMonoisotopicMass)
                            bin.Overlapping++;
                }
            }
        }

        private static void IdentifyFracWithSingle(BinTreeStructure myTreeStructure)
        {
            foreach (Bin bin in myTreeStructure.FinalBins)
            {
                bin.FracWithSingle = (double)bin.uniquePSMs.Values.Count(b => !b.Item3.IsDecoy && b.Item3.thisPSM.peptidesWithSetModifications.Count == 1) / bin.uniquePSMs.Values.Count(b => !b.Item3.IsDecoy);
            }
        }

        private static void IdentifyAAsInCommon(BinTreeStructure myTreeStructure)
        {
            foreach (Bin bin in myTreeStructure.FinalBins)
            {
                bin.AAsInCommon = new Dictionary<char, int>();
                foreach (var hehe in bin.uniquePSMs.Values.Where(b => !b.Item3.IsDecoy))
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
            foreach (Bin bin in myTreeStructure.FinalBins)
            {
                bin.modsInCommon = new Dictionary<string, int>();
                foreach (var hehe in bin.uniquePSMs.Values.Where(b => !b.Item3.IsDecoy))
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
            foreach (Bin bin in myTreeStructure.FinalBins)
            {
                bin.residueCount = new Dictionary<char, int>();
                foreach (var hehe in bin.uniquePSMs.Values)
                {
                    double bestScore = hehe.Item3.thisPSM.LocalizedScores.Max();
                    if (bestScore >= hehe.Item3.thisPSM.Score + 1 && !hehe.Item3.IsDecoy)
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
            foreach (var bin in myTreeStructure.FinalBins)
            {
                var ok = new HashSet<string>();
                var okformula = new HashSet<string>();
                foreach (var hm in UnimodDeserialized.modifications)
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
            foreach (var bin in myTreeStructure.FinalBins)
            {
                var ok = new HashSet<string>();
                foreach (var hm in UniprotDeseralized)
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
            double totalTargetCount = myTreeStructure.FinalBins.Select(b => b.CountTarget).Sum();
            var ok = new HashSet<Tuple<double, double, double>>();
            foreach (var bin in myTreeStructure.FinalBins.Where(b => Math.Abs(b.MassShift) > v))
                foreach (var bin2 in myTreeStructure.FinalBins.Where(b => Math.Abs(b.MassShift) > v))
                    if (bin.CountTarget * bin2.CountTarget >= totalTargetCount)
                        ok.Add(new Tuple<double, double, double>(bin.MassShift, bin2.MassShift, Math.Min(bin.CountTarget, bin2.CountTarget)));

            foreach (var bin in myTreeStructure.FinalBins)
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
            foreach (var bin in myTreeStructure.FinalBins)
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
            foreach (Bin bin in myTreeStructure.FinalBins)
            {
                bin.Mine = "";
                foreach (MyInfo myInfo in myInfos)
                {
                    if (Math.Abs(myInfo.MassShift - bin.MassShift) <= v)
                    {
                        bin.Mine = myInfo.infostring;
                    }
                }
            }
        }

        private static BinTreeStructure MyAnalysis(List<NewPsmWithFdr> limitedpsms_with_fdr, double binTol)
        {
            var myTreeStructure = new BinTreeStructure();
            myTreeStructure.GenerateBins(limitedpsms_with_fdr, binTol);

            IdentifyUnimodBins(myTreeStructure, binTol);
            IdentifyUniprotBins(myTreeStructure, binTol);
            IdentifyAA(myTreeStructure, binTol);

            IdentifyCombos(myTreeStructure, binTol);

            IdentifyResidues(myTreeStructure);

            IdentifyMods(myTreeStructure);

            IdentifyAAsInCommon(myTreeStructure);

            IdentifyMine(myTreeStructure, binTol);

            OverlappingIonSequences(myTreeStructure);

            IdentifyFracWithSingle(myTreeStructure);

            return myTreeStructure;
        }

        private static List<NewPsmWithFdr> DoFalseDiscoveryRateAnalysis(IEnumerable<PsmWithMultiplePossiblePeptides> items, SearchMode sm)
        {
            var ids = new List<NewPsmWithFdr>();
            foreach (PsmWithMultiplePossiblePeptides item in items)
            {
                ids.Add(new NewPsmWithFdr(item));
            }

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            int[] cumulative_target_per_notch = new int[sm.NumNotches];
            int[] cumulative_decoy_per_notch = new int[sm.NumNotches];

            for (int i = 0; i < ids.Count; i++)
            {
                var item = ids[i];
                var isDecoy = item.IsDecoy;
                int notch = item.thisPSM.newPsm.notch;
                if (isDecoy)
                    cumulative_decoy++;
                else
                    cumulative_target++;

                if (isDecoy)
                    cumulative_decoy_per_notch[notch]++;
                else
                    cumulative_target_per_notch[notch]++;

                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                double temp_q_value_for_notch = (double)cumulative_decoy_per_notch[notch] / (cumulative_target_per_notch[notch] + cumulative_decoy_per_notch[notch]);
                item.SetValues(cumulative_target, cumulative_decoy, temp_q_value, cumulative_target_per_notch[notch], cumulative_decoy_per_notch[notch], temp_q_value_for_notch);
            }

            double min_q_value = double.PositiveInfinity;
            double[] min_q_value_notch = new double[sm.NumNotches];
            for (int i = 0; i < sm.NumNotches; i++)
                min_q_value_notch[i] = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                NewPsmWithFdr id = ids[i];
                if (id.qValue > min_q_value)
                    id.qValue = min_q_value;
                else if (id.qValue < min_q_value)
                    min_q_value = id.qValue;

                int notch = id.thisPSM.newPsm.notch;
                if (id.qValueNotch > min_q_value_notch[notch])
                    id.qValueNotch = min_q_value_notch[notch];
                else if (id.qValueNotch < min_q_value_notch[notch])
                    min_q_value_notch[notch] = id.qValueNotch;
            }

            return ids;
        }

        private void AddObservedPeptidesToDictionary()
        {
            Status("Adding new keys to peptide dictionary...");
            foreach (var psmListForAspecificSerchMode in newPsms)
            {
                if (psmListForAspecificSerchMode != null)
                    foreach (var psm in psmListForAspecificSerchMode)
                    {
                        if (psm != null)
                        {
                            var cp = psm.GetCompactPeptide(variableModifications, localizeableModifications, fixedModifications);
                            if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp))
                                compactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());
                        }
                    }
            }

            int totalProteins = proteinList.Count;

            Status("Adding possible sources to peptide dictionary...");

            int proteinsSeen = 0;
            int old_progress = 0;

            var obj = new object();

            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> local = compactPeptideToProteinPeptideMatching.ToDictionary(b => b.Key, b => new HashSet<PeptideWithSetModifications>());

                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, InitiatorMethionineBehavior.Variable).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length == 1 || peptide.Length > byte.MaxValue - 2) // 2 is for indexing terminal modifications
                            continue;

                        var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(variableModifications, maxModIsoforms, max_mods_for_peptide, fixedModifications).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            HashSet<PeptideWithSetModifications> v;
                            if (local.TryGetValue(new CompactPeptide(yyy, variableModifications, localizeableModifications, fixedModifications), out v))
                            {
                                v.Add(yyy);
                            }
                        }
                    }
                }
                lock (obj)
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

        #region Private Classes

        private class SequenceComparer : IEqualityComparer<PsmWithMultiplePossiblePeptides>
        {

            #region Public Methods

            public bool Equals(PsmWithMultiplePossiblePeptides x, PsmWithMultiplePossiblePeptides y)
            {
                return x.FullSequence.Equals(y.FullSequence);
            }

            public int GetHashCode(PsmWithMultiplePossiblePeptides obj)
            {
                return obj.FullSequence.GetHashCode();
            }

            #endregion Public Methods

        }

        #endregion Private Classes

    }
}