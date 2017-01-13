using MassSpectrometry;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace InternalLogicEngineLayer
{
    public class AnalysisEngine : MyEngine
    {

        #region Private Fields

        private ParentSpectrumMatch[][] newPsms;

        private Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching;

        private List<Protein> proteinList;

        private List<MorpheusModification> variableModifications;

        private List<MorpheusModification> fixedModifications;

        private List<MorpheusModification> localizeableModifications;

        private Protease protease;

        private List<SearchMode> searchModes;

        private IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;

        private Tolerance fragmentTolerance;

        private Action<BinTreeStructure, string> action1;

        private Action<List<NewPsmWithFDR>, string> action2;

        private bool doParsimony;

        #endregion Private Fields

        #region Public Constructors

        public AnalysisEngine(ParentSpectrumMatch[][] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<Protein> proteinList, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, Protease protease, List<SearchMode> searchModes, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, Tolerance fragmentTolerance, Action<BinTreeStructure, string> action1, Action<List<NewPsmWithFDR>, string> action2, bool doParsimony) : base(2)
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
        }

        #endregion Public Constructors

        #region Public Methods

        public static Dictionary<CompactPeptide, PeptideWithSetModifications> GetSingleMatchDictionary(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> fullSequenceToProteinPeptideMatching)
        {
            // Right now very stupid, add the first decoy one, and if no decoy, add the first one
            Dictionary<CompactPeptide, PeptideWithSetModifications> outDict = new Dictionary<CompactPeptide, PeptideWithSetModifications>();
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
            /* TODO: KeyedCollection instead of dictionary? Faster lookup */

            // makes dictionary with proteins as keys and list of associated peptides as the value (swaps input parameter dictionary keys/values)
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

            // build each protein group after parsimony and match it to its peptide
            Dictionary<CompactPeptide, HashSet<Protein>> peptideProteinGroupMatch = new Dictionary<CompactPeptide, HashSet<Protein>>();

            foreach (var kvp in parsimonyDict)
            {
                foreach (var peptide in kvp.Value)
                {
                    if (!peptideProteinGroupMatch.ContainsKey(peptide))
                    {
                        HashSet<Protein> proteinGroup = new HashSet<Protein>();

                        foreach (var kvp1 in parsimonyDict)
                        {
                            if (kvp1.Value.Contains(peptide))
                            {
                                proteinGroup.Add(kvp1.Key);
                            }
                        }

                        peptideProteinGroupMatch.Add(peptide, proteinGroup);
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
                        HashSet<Protein> proteinGroup = new HashSet<Protein>();

                        // get the peptide's protein group after parsimony
                        peptideProteinGroupMatch.TryGetValue(peptide, out proteinGroup);

                        // find peptide's original (unparsimonious) virtual peptide matches
                        fullSequenceToProteinPeptideMatching.TryGetValue(peptide, out oldVirtualPeptides);

                        // get the virtual peptides that belong to the post-parsimony protein(s) only
                        foreach (var virtualPeptide in oldVirtualPeptides)
                        {
                            if (proteinGroup.Contains(virtualPeptide.protein))
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

        #endregion Public Methods

        #region Protected Methods

        protected override void ValidateParams()
        {
            if (newPsms == null)
                throw new EngineValidationException("newPsms cannot be null");
        }

        protected override MyResults RunSpecific()
        {
            status("Running analysis engine!");
            //At this point have Spectrum-Sequence matching, without knowing which protein, and without know if target/decoy
            status("Adding observed peptides to dictionary...");
            AddObservedPeptidesToDictionary();

            //status("Getting single match just for FDR purposes...");
            //var fullSequenceToProteinSingleMatch = GetSingleMatchDictionary(compactPeptideToProteinPeptideMatching);

            List<NewPsmWithFDR>[] yeah = new List<NewPsmWithFDR>[searchModes.Count];
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

                    var orderedPsmsWithFDR = DoFalseDiscoveryRateAnalysis(orderedPsmsWithPeptides);
                    var limitedpsms_with_fdr = orderedPsmsWithFDR.Where(b => (b.QValue <= 0.01)).ToList();
                    if (limitedpsms_with_fdr.Where(b => !b.isDecoy).Count() > 0)
                    {
                        status("Running histogram analysis...");
                        var hm = MyAnalysis(limitedpsms_with_fdr);
                        action1(hm, searchModes[j].FileNameAddition);
                    }

                    action2(orderedPsmsWithFDR, searchModes[j].FileNameAddition);

                    yeah[j] = orderedPsmsWithFDR;
                }
            }
            if (doParsimony)
            {
                status("Getting protein parsimony dictionary...");
                compactPeptideToProteinPeptideMatching = ApplyProteinParsimony(compactPeptideToProteinPeptideMatching);

                // to do protein FDR, need peptides and their virtual peptide matches (for protein group analysis) with scores
                //Dictionary<Protein, List<NewPsmWithFDR>> aggregateProteinList = new Dictionary<Protein, List<NewPsmWithFDR>();
                //aggregateProteinList = BuildProteinGroupsAndDoProteinFDR(orderedPsmsWithFDR, compactPeptideToProteinPeptideMatching);
            }

            return new AnalysisResults(this, yeah, compactPeptideToProteinPeptideMatching);
        }

        #endregion Protected Methods

        #region Private Methods

        private static BinTreeStructure MyAnalysis(List<NewPsmWithFDR> limitedpsms_with_fdr)
        {
            BinTreeStructure myTreeStructure = new BinTreeStructure();
            myTreeStructure.GenerateBins(limitedpsms_with_fdr, 0.003);
            myTreeStructure.AddToBins(limitedpsms_with_fdr);

            MyAnalysisClass.IdentifyUnimodBins(myTreeStructure, 0.003, unimodDeserialized);
            MyAnalysisClass.IdentifyUniprotBins(myTreeStructure, 0.003, uniprotDeseralized);
            MyAnalysisClass.IdentifyAA(myTreeStructure, 0.003);

            MyAnalysisClass.IdentifyCombos(myTreeStructure, 0.003);

            MyAnalysisClass.IdentifyResidues(myTreeStructure);

            MyAnalysisClass.IdentifyMods(myTreeStructure);

            MyAnalysisClass.IdentifyAAsInCommon(myTreeStructure);

            MyAnalysisClass.IdentifyMine(myTreeStructure, 0.003);

            return myTreeStructure;
        }

        private static List<NewPsmWithFDR> DoFalseDiscoveryRateAnalysis(IEnumerable<PSMwithTargetDecoyKnown> items)
        {
            List<NewPsmWithFDR> ids = new List<NewPsmWithFDR>();

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

        private List<ProteinGroup> BuildProteinGroupsAndDoProteinFDR(List<NewPsmWithFDR> PSMs, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching)
        {
            return new List<ProteinGroup>();
        }

        private void AddObservedPeptidesToDictionary()
        {
            status("Adding new keys to peptide dictionary...");
            foreach (var ah in newPsms)
            {
                if (ah != null)
                    foreach (var fhh in ah)
                    {
                        if (fhh != null)
                        {
                            var cp = fhh.GetCompactPeptide(variableModifications, localizeableModifications);
                            if (!compactPeptideToProteinPeptideMatching.ContainsKey(cp))
                                compactPeptideToProteinPeptideMatching.Add(cp, new HashSet<PeptideWithSetModifications>());
                        }
                    }
            }

            int totalProteins = proteinList.Count;

            status("Adding possible sources to peptide dictionary...");

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

                        var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(variableModifications, 4098, 3, localizeableModifications).ToList();
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
                }
            });
        }

        #endregion Private Methods

    }
}