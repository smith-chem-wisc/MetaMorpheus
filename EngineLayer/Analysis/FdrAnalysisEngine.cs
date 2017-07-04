using Proteomics;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Analysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {

        #region Private Fields

        private readonly IEnumerable<PsmParent>[] newPsms;
        private readonly List<MassDiffAcceptor> searchModes;
        private readonly bool doParsimony;
        private readonly bool noOneHitWonders;
        private readonly bool treatModPeptidesAsDifferentPeptides;
        private Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching;

        #endregion Private Fields

        #region Public Constructors

        public FdrAnalysisEngine(List<PsmParent>[] newPsms, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching, List<MassDiffAcceptor> searchModes, bool doParsimony, bool noOneHitWonders, bool modPeptidesAreUnique, List<string> nestedIds) : base(nestedIds)
        {
            this.doParsimony = doParsimony;
            this.noOneHitWonders = noOneHitWonders;
            this.newPsms = newPsms;
            this.compactPeptideToProteinPeptideMatching = compactPeptideToProteinPeptideMatching;
            this.searchModes = searchModes;
            this.treatModPeptidesAsDifferentPeptides = modPeptidesAreUnique;
        }

        #endregion Public Constructors

        #region Public Methods

        public void ApplyProteinParsimony(out List<ProteinGroup> proteinGroups)
        {
            // digesting an XML database results in a non-mod-agnostic digestion; need to fix this if mod-agnostic parsimony enabled
            if (!treatModPeptidesAsDifferentPeptides)
            {
                var groups = compactPeptideToProteinPeptideMatching.GroupBy(p => p.Value.First().BaseSequence);

                foreach (var baseseq in groups)
                {
                    var proteinsForAllCompactPeptidesWithThisBaseSeq = baseseq.SelectMany(p => p.Value.Select(v => v.Protein)).Distinct().ToList();

                    foreach (var compactPeptide in baseseq)
                    {
                        var proteinsHere = compactPeptide.Value.Select(p => p.Protein).Distinct();

                        if (proteinsHere.Count() != proteinsForAllCompactPeptidesWithThisBaseSeq.Count)
                        {
                            var missingProteins = proteinsForAllCompactPeptidesWithThisBaseSeq.Except(proteinsHere);

                            foreach (var missingProtein in missingProteins)
                            {
                                var templatePepWithSetMod = baseseq.SelectMany(p => p.Value).Where(v => v.Protein.Equals(missingProtein)).First();
                                compactPeptideToProteinPeptideMatching[compactPeptide.Key].Add(new PeptideWithSetModifications(compactPeptide.Value.First(), templatePepWithSetMod));
                            }
                        }
                    }
                }
            }

            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<CompactPeptide>>();
            var parsimonyProteinList = new Dictionary<Protein, HashSet<CompactPeptide>>();
            var proteinsWithUniquePeptides = new Dictionary<Protein, HashSet<PeptideWithSetModifications>>();

            // peptide matched to fullseq (used depending on user preference)
            var compactPeptideToFullSeqMatch = compactPeptideToProteinPeptideMatching.ToDictionary(x => x.Key, x => x.Value.First().Sequence);

            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                // finds unique peptides (peptides that can belong to only one protein)
                HashSet<Protein> proteinsAssociatedWithThisPeptide = new HashSet<Protein>(kvp.Value.Select(p => p.Protein));
                if (proteinsAssociatedWithThisPeptide.Count() == 1)
                {
                    var peptides = new HashSet<PeptideWithSetModifications>();
                    if (!proteinsWithUniquePeptides.TryGetValue(kvp.Value.First().Protein, out peptides))
                        proteinsWithUniquePeptides.Add(kvp.Value.First().Protein, new HashSet<PeptideWithSetModifications>(kvp.Value));
                    else
                        peptides.UnionWith(kvp.Value);
                }

                // if a peptide is associated with a decoy protein, remove all target protein associations with the peptide
                if (kvp.Value.Where(p => p.Protein.IsDecoy).Any())
                    kvp.Value.RemoveWhere(p => !p.Protein.IsDecoy);

                // if a peptide is associated with a contaminant protein, remove all target protein associations with the peptide
                if (kvp.Value.Where(p => p.Protein.IsContaminant).Any())
                    kvp.Value.RemoveWhere(p => !p.Protein.IsContaminant);
            }

            // makes dictionary with proteins as keys and list of associated peptides as the value (makes parsimony algo easier)
            foreach (var kvp in compactPeptideToProteinPeptideMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    var peptides = new HashSet<CompactPeptide>();
                    if (!proteinToPeptidesMatching.TryGetValue(peptide.Protein, out peptides))
                        proteinToPeptidesMatching.Add(peptide.Protein, new HashSet<CompactPeptide>() { kvp.Key });
                    else
                        peptides.Add(kvp.Key);
                }
            }

            // build protein list for each peptide before parsimony has been applied
            var peptideSeqProteinListMatch = new Dictionary<string, HashSet<Protein>>();
            foreach (var kvp in proteinToPeptidesMatching)
            {
                foreach (var peptide in kvp.Value)
                {
                    var proteinListHere = new HashSet<Protein>();
                    string pepSequence;
                    if (!treatModPeptidesAsDifferentPeptides)
                        pepSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));
                    else
                        pepSequence = compactPeptideToFullSeqMatch[peptide];

                    if (!peptideSeqProteinListMatch.TryGetValue(pepSequence, out proteinListHere))
                        peptideSeqProteinListMatch.Add(pepSequence, new HashSet<Protein>() { kvp.Key });
                    else
                        proteinListHere.Add(kvp.Key);
                }
            }

            // dictionary associates proteins w/ unused base seqs (list will shrink over time)
            var algDictionary = new Dictionary<Protein, HashSet<string>>();
            foreach (var kvp in peptideSeqProteinListMatch)
            {
                foreach (var protein in kvp.Value)
                {
                    HashSet<string> newPeptideBaseSeqs;
                    if (algDictionary.TryGetValue(protein, out newPeptideBaseSeqs))
                        newPeptideBaseSeqs.Add(kvp.Key);
                    else
                        algDictionary.Add(protein, new HashSet<string> { kvp.Key });
                }
            }

            // dictionary associates proteins w/ unused base seqs (list will NOT shrink over time)
            var proteinToPepSeqMatch = algDictionary.ToDictionary(x => x.Key, x => x.Value);

            // *** main parsimony loop
            bool uniquePeptidesLeft = false;
            if (proteinsWithUniquePeptides.Any())
                uniquePeptidesLeft = true;
            int numNewSeqs = algDictionary.Max(p => p.Value.Count);

            while (numNewSeqs != 0)
            {
                var possibleBestProteinList = new List<KeyValuePair<Protein, HashSet<string>>>();

                if (uniquePeptidesLeft)
                {
                    var proteinsWithUniquePeptidesLeft = algDictionary.Where(p => proteinsWithUniquePeptides.ContainsKey(p.Key));
                    if (proteinsWithUniquePeptidesLeft.Any())
                        possibleBestProteinList.Add(proteinsWithUniquePeptidesLeft.First());
                    else
                        uniquePeptidesLeft = false;
                }

                // gets list of proteins with the most unaccounted-for peptide base sequences
                if (!uniquePeptidesLeft)
                    possibleBestProteinList = algDictionary.Where(p => p.Value.Count == numNewSeqs).ToList();

                Protein bestProtein = possibleBestProteinList.First().Key;
                HashSet<string> newSeqs = new HashSet<string>(algDictionary[bestProtein]);

                // may need to select different protein
                if (possibleBestProteinList.Count > 1)
                {
                    var proteinsWithTheseBaseSeqs = new HashSet<Protein>();

                    foreach (var kvp in possibleBestProteinList)
                    {
                        if (newSeqs.IsSubsetOf(kvp.Value))
                            proteinsWithTheseBaseSeqs.Add(kvp.Key);
                    }

                    if (proteinsWithTheseBaseSeqs.Count > 1)
                    {
                        var proteinsOrderedByTotalPeptideCount = new Dictionary<Protein, HashSet<string>>();
                        foreach (var protein in proteinsWithTheseBaseSeqs)
                            proteinsOrderedByTotalPeptideCount.Add(protein, proteinToPepSeqMatch[protein]);

                        bestProtein = proteinsOrderedByTotalPeptideCount.OrderByDescending(kvp => kvp.Value.Count).First().Key;
                    }
                }

                parsimonyProteinList.Add(bestProtein, proteinToPeptidesMatching[bestProtein]);

                // remove used peptides from their proteins
                foreach (var newBaseSeq in newSeqs)
                {
                    HashSet<Protein> proteinsWithThisPeptide = peptideSeqProteinListMatch[newBaseSeq];

                    foreach (var protein in proteinsWithThisPeptide)
                        algDictionary[protein].Remove(newBaseSeq);
                }

                algDictionary.Remove(bestProtein);
                if (algDictionary.Any())
                    numNewSeqs = algDictionary.Max(p => p.Value.Count);
                else
                    numNewSeqs = 0;
            }

            // *** done with parsimony

            // add indistinguishable proteins
            var proteinsGroupedByNumPeptides = proteinToPeptidesMatching.GroupBy(p => p.Value.Count);
            var parsimonyProteinsGroupedByNumPeptides = parsimonyProteinList.GroupBy(p => p.Value.Count);
            var indistinguishableProteins = new ConcurrentDictionary<Protein, HashSet<CompactPeptide>>();

            foreach (var group in proteinsGroupedByNumPeptides)
            {
                var parsimonyProteinsWithSameNumPeptides = parsimonyProteinsGroupedByNumPeptides.Where(p => p.Key == group.Key).FirstOrDefault();
                var list = group.ToList();
                if (parsimonyProteinsWithSameNumPeptides != null)
                {
                    Parallel.ForEach(Partitioner.Create(0, list.Count),
                        new ParallelOptions { MaxDegreeOfParallelism = -1 },
                        (range, loopState) =>
                        {
                            for (int i = range.Item1; i < range.Item2; i++)
                            {
                                foreach (var parsimonyProteinWithThisNumPeptides in parsimonyProteinsWithSameNumPeptides)
                                {
                                    if (parsimonyProteinWithThisNumPeptides.Key != list[i].Key)
                                        if (proteinToPeptidesMatching[parsimonyProteinWithThisNumPeptides.Key].SetEquals(proteinToPeptidesMatching[list[i].Key]))
                                            indistinguishableProteins.GetOrAdd(list[i].Key, proteinToPeptidesMatching[list[i].Key]);
                                }
                            }
                        }
                    );
                }
            }
            foreach (var protein in indistinguishableProteins)
                parsimonyProteinList.Add(protein.Key, protein.Value);

            foreach (var kvp in compactPeptideToProteinPeptideMatching)
                kvp.Value.RemoveWhere(p => !parsimonyProteinList.ContainsKey(p.Protein));

            proteinGroups = ConstructProteinGroups(new HashSet<PeptideWithSetModifications>(proteinsWithUniquePeptides.Values.SelectMany(p => p)), new HashSet<PeptideWithSetModifications>(compactPeptideToProteinPeptideMatching.Values.SelectMany(p => p)));

            Status("Finished Parsimony", nestedIds);
        }

        public List<ProteinGroup> ConstructProteinGroups(HashSet<PeptideWithSetModifications> uniquePeptides, HashSet<PeptideWithSetModifications> allPeptides)
        {
            List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<PeptideWithSetModifications>>();

            foreach (var peptide in allPeptides)
            {
                HashSet<PeptideWithSetModifications> peptidesHere;
                if (proteinToPeptidesMatching.TryGetValue(peptide.Protein, out peptidesHere))
                    peptidesHere.Add(peptide);
                else
                    proteinToPeptidesMatching.Add(peptide.Protein, new HashSet<PeptideWithSetModifications> { peptide });
            }

            // build protein group after parsimony (each group only has 1 protein at this point, indistinguishables will be merged during scoring)
            foreach (var kvp in proteinToPeptidesMatching)
            {
                var allPeptidesHere = proteinToPeptidesMatching[kvp.Key];
                var uniquePeptidesHere = new HashSet<PeptideWithSetModifications>(allPeptidesHere.Where(p => uniquePeptides.Contains(p)));

                proteinGroups.Add(new ProteinGroup(new HashSet<Protein>() { kvp.Key }, allPeptidesHere, uniquePeptidesHere));
            }

            foreach (var proteinGroup in proteinGroups)
            {
                proteinGroup.AllPeptides.RemoveWhere(p => !proteinGroup.Proteins.Contains(p.Protein));
                proteinGroup.DisplayModsOnPeptides = treatModPeptidesAsDifferentPeptides;
            }

            return proteinGroups;
        }

        public void ScoreProteinGroups(List<ProteinGroup> proteinGroups, IEnumerable<PsmParent> psmList)
        {
            Status("Scoring protein groups...", nestedIds);

            // add each protein groups PSMs
            var peptideToPsmMatching = new Dictionary<PeptideWithSetModifications, HashSet<PsmParent>>();
            foreach (var psm in psmList)
            {
                if (psm.FdrInfo.QValue <= 0.01)
                {
                    foreach (var pepWithSetMods in psm.Pli.PeptidesWithSetModifications)
                    {
                        var psmsForThisPeptide = new HashSet<PsmParent>();

                        if (!peptideToPsmMatching.TryGetValue(pepWithSetMods, out psmsForThisPeptide))
                            peptideToPsmMatching.Add(pepWithSetMods, new HashSet<PsmParent> { psm });
                        else
                            psmsForThisPeptide.Add(psm);
                    }
                }
            }

            foreach (var proteinGroup in proteinGroups)
            {
                List<PeptideWithSetModifications> pepsToRemove = new List<PeptideWithSetModifications>();
                foreach (var peptide in proteinGroup.AllPeptides)
                {
                    // build PSM list for scoring
                    HashSet<PsmParent> psms;
                    if (peptideToPsmMatching.TryGetValue(peptide, out psms))
                        proteinGroup.AllPsmsBelowOnePercentFDR.UnionWith(psms);
                    else
                        pepsToRemove.Add(peptide);
                }

                proteinGroup.AllPeptides.ExceptWith(pepsToRemove);
                proteinGroup.UniquePeptides.ExceptWith(pepsToRemove);
            }

            // score the group
            foreach (var proteinGroup in proteinGroups)
                proteinGroup.Score();

            // merge protein groups that are indistinguishable after scoring
            var pg = proteinGroups.OrderByDescending(p => p.ProteinGroupScore).ToList();
            for (int i = 0; i < (pg.Count - 1); i++)
            {
                if (pg[i].ProteinGroupScore == pg[i + 1].ProteinGroupScore && pg[i].ProteinGroupScore != 0)
                {
                    var pgsWithThisScore = pg.Where(p => p.ProteinGroupScore == pg[i].ProteinGroupScore).ToList();

                    // check to make sure they have the same peptides, then merge them
                    foreach (var p in pgsWithThisScore)
                    {
                        var seqs1 = new HashSet<string>(p.AllPeptides.Select(x => x.Sequence));
                        var seqs2 = new HashSet<string>(pg[i].AllPeptides.Select(x => x.Sequence));

                        if (p != pg[i] && seqs1.SetEquals(seqs2))
                        {
                            pg[i].MergeProteinGroupWith(p);
                        }
                    }
                }
            }

            // remove empty protein groups (peptides were too poor quality or group was merged)
            proteinGroups.RemoveAll(p => p.ProteinGroupScore == 0);

            // calculate sequence coverage
            foreach (var proteinGroup in proteinGroups)
                proteinGroup.CalculateSequenceCoverage();

            // distribute razor peptides
            var sharedPepWithProteinGroups = new Dictionary<string, HashSet<ProteinGroup>>();
            var seqToPeptide = new Dictionary<string, HashSet<PeptideWithSetModifications>>();
            foreach (var proteinGroup in proteinGroups)
            {
                var sharedPeps = proteinGroup.AllPeptides.Except(proteinGroup.UniquePeptides);
                foreach (var sharedPep in sharedPeps)
                {
                    HashSet<ProteinGroup> v;
                    if (sharedPepWithProteinGroups.TryGetValue(sharedPep.Sequence, out v))
                        v.Add(proteinGroup);
                    else
                        sharedPepWithProteinGroups.Add(sharedPep.Sequence, new HashSet<ProteinGroup> { proteinGroup });

                    HashSet<PeptideWithSetModifications> v1;
                    if (seqToPeptide.TryGetValue(sharedPep.Sequence, out v1))
                        v1.Add(sharedPep);
                    else
                        seqToPeptide.Add(sharedPep.Sequence, new HashSet<PeptideWithSetModifications> { sharedPep });
                }
            }

            foreach (var kvp in sharedPepWithProteinGroups)
            {
                int i = kvp.Value.Select(p => p.AllPeptides.Select(x => x.BaseSequence).Count()).Max();
                HashSet<ProteinGroup> t = new HashSet<ProteinGroup>(kvp.Value.Where(p => p.AllPeptides.Select(x => x.BaseSequence).Count() == i));
                foreach (var proteinGroup in t)
                {
                    var peptides = seqToPeptide[kvp.Key];
                    proteinGroup.RazorPeptides.UnionWith(peptides.Where(p => proteinGroup.Proteins.Contains(p.Protein)));
                }
            }
        }

        public List<ProteinGroup> DoProteinFdr(List<ProteinGroup> proteinGroups)
        {
            Status("Calculating protein FDR...", nestedIds);

            if (noOneHitWonders)
            {
                if (treatModPeptidesAsDifferentPeptides)
                    proteinGroups = proteinGroups.Where(p => p.isDecoy || new HashSet<string>(p.AllPeptides.Select(x => x.Sequence)).Count() > 1).ToList();
                else
                    proteinGroups = proteinGroups.Where(p => p.isDecoy || new HashSet<string>(p.AllPeptides.Select(x => x.BaseSequence)).Count() > 1).ToList();
            }

            // order protein groups by score
            var sortedProteinGroups = proteinGroups.OrderByDescending(b => b.ProteinGroupScore).ToList();

            // do fdr
            int cumulativeTarget = 0;
            int cumulativeDecoy = 0;
            foreach (var proteinGroup in sortedProteinGroups)
            {
                if (proteinGroup.isDecoy)
                    cumulativeDecoy++;
                else
                    cumulativeTarget++;

                proteinGroup.CumulativeTarget = cumulativeTarget;
                proteinGroup.CumulativeDecoy = cumulativeDecoy;
                proteinGroup.QValue = ((double)cumulativeDecoy / (cumulativeTarget + cumulativeDecoy));
            }

            return sortedProteinGroups;
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            FdrAnalysisResults myAnalysisResults = new FdrAnalysisResults(this);
            Status("Running analysis engine!", nestedIds);

            List<ProteinGroup>[] proteinGroups = null;
            if (doParsimony)
            {
                Status("Applying protein parsimony...", nestedIds);
                proteinGroups = new List<ProteinGroup>[searchModes.Count];
                // TODO**: make this faster (only apply parsimony once but make multiple instances of the same ProteinGroups
                for (int i = 0; i < searchModes.Count; i++)
                    ApplyProteinParsimony(out proteinGroups[i]);
            }

            for (int j = 0; j < searchModes.Count; j++)
            {
                Status("Running FDR analysis...", nestedIds);
                DoFalseDiscoveryRateAnalysis(newPsms[j], searchModes[j]);

                //if (!doParsimony && writeMzIdentmlAction != null)
                //    writeMzIdentmlAction.Invoke(newPsms[j], null, searchModes[j], searchModes[j].FileNameAddition, nestedIds);

                // write aggregate results
                if (doParsimony)
                {
                    ScoreProteinGroups(proteinGroups[j], newPsms[j]);
                    proteinGroups[j] = DoProteinFdr(proteinGroups[j]);

                    //writeProteinGroupsAction(proteinGroups[j], "allProteinGroups_" + searchModes[j].FileNameAddition, nestedIds);
                    //writeMzIdentmlAction(newPsms[j], proteinGroups[j], searchModes[j], searchModes[j].FileNameAddition, nestedIds);
                }
            }

            myAnalysisResults.ProteinGroups = proteinGroups;
            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private static List<PsmParent> DoFalseDiscoveryRateAnalysis(IEnumerable<PsmParent> items, MassDiffAcceptor sm)
        {
            var ids = new List<PsmParent>();
            foreach (PsmParent item in items)
                ids.Add(item);

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            int[] cumulative_target_per_notch = new int[sm.NumNotches];
            int[] cumulative_decoy_per_notch = new int[sm.NumNotches];

            for (int i = 0; i < ids.Count; i++)
            {
                var item = ids[i];
                var isDecoy = item.Pli.IsDecoy;
                int notch = item.Notch;
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
                PsmParent id = ids[i];
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;

                int notch = id.Notch;
                if (id.FdrInfo.QValueNotch > min_q_value_notch[notch])
                    id.FdrInfo.QValueNotch = min_q_value_notch[notch];
                else if (id.FdrInfo.QValueNotch < min_q_value_notch[notch])
                    min_q_value_notch[notch] = id.FdrInfo.QValueNotch;
            }

            return ids;
        }

        #endregion Private Methods

    }
}