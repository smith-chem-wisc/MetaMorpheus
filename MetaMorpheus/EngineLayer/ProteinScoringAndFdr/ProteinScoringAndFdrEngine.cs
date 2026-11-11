using EngineLayer.SpectrumMatch;
using System.Collections.Generic;
using System.Linq;
using Omics;

namespace EngineLayer
{
    public class ProteinScoringAndFdrEngine : MetaMorpheusEngine
    {
        private readonly IList<SpectralMatch> _FilteredPsms;
        private readonly bool NoOneHitWonders;
        private readonly bool TreatModPeptidesAsDifferentPeptides;
        private readonly bool MergeIndistinguishableProteinGroups;
        private readonly bool LabelPeptidesUniqueToOneGroupAfterFdr;
        private readonly List<ProteinGroup> ProteinGroups;
        private readonly HashSet<string> _decoyIdentifiers;
        private readonly FilterType _filterType;
        private readonly double _filterThreshold;
        private readonly bool _psmsAlreadyFiltered;

        /// <summary>
        /// Whether protein FDR ranks groups by their best peptide's q-value or by its PEP.
        /// </summary>
        public FilterType FilterType => _filterType;

        /// <summary>
        /// Takes every PSM from the search; the engine keeps those passing the q-value or PEP q-value threshold
        /// from <paramref name="commonParameters"/>, whichever is stricter.
        /// </summary>
        /// <param name="labelPeptidesUniqueToOneGroup">
        /// True for the search's own call: after FDR, label each group's peptides unique or shared among the groups it
        /// returns. False for a per-file call, whose groups are subsets of the search's and already carry its labels
        /// (see <see cref="ProteinGroup.ConstructSubsetProteinGroup"/>), so the per-file tables agree with the combined one.
        /// </param>
        public ProteinScoringAndFdrEngine(List<ProteinGroup> proteinGroups, IList<SpectralMatch> filteredPsms, bool noOneHitWonders, bool treatModPeptidesAsDifferentPeptides, bool mergeIndistinguishableProteinGroups,
            CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds,
            bool labelPeptidesUniqueToOneGroup = true)
            : this(proteinGroups, filteredPsms, psmsAlreadyFiltered: false, noOneHitWonders, treatModPeptidesAsDifferentPeptides, mergeIndistinguishableProteinGroups,
                commonParameters, fileSpecificParameters, nestedIds, labelPeptidesUniqueToOneGroup) { }

        /// <summary>
        /// When <paramref name="psmsAlreadyFiltered"/> is true, <paramref name="psms"/> are scored as-is with no q-value
        /// filtering. Parallel search relies on this to score every baseline PSM, since it aligns transient FDR to the
        /// baseline rather than recalculating it.
        /// </summary>
        public ProteinScoringAndFdrEngine(List<ProteinGroup> proteinGroups, IList<SpectralMatch> psms, bool psmsAlreadyFiltered, bool noOneHitWonders, bool treatModPeptidesAsDifferentPeptides, bool mergeIndistinguishableProteinGroups,
            CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds,
            bool labelPeptidesUniqueToOneGroup = true) : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            _FilteredPsms = psms;
            _psmsAlreadyFiltered = psmsAlreadyFiltered;
            ProteinGroups = proteinGroups;
            NoOneHitWonders = noOneHitWonders;
            TreatModPeptidesAsDifferentPeptides = treatModPeptidesAsDifferentPeptides;
            MergeIndistinguishableProteinGroups = mergeIndistinguishableProteinGroups;
            LabelPeptidesUniqueToOneGroupAfterFdr = labelPeptidesUniqueToOneGroup;
            _decoyIdentifiers = proteinGroups.SelectMany(p => p.Proteins.Where(b => b.IsDecoy).Select(b => b.Accession.Split('_')[0])).ToHashSet();

            // Pre-filtered callers keep their original tie-break: equal thresholds sort by PEP q-value.
            bool usePepQValue = psmsAlreadyFiltered
                ? !(CommonParameters.QValueThreshold < CommonParameters.PepQValueThreshold)
                : CommonParameters.PepQValueThreshold < CommonParameters.QValueThreshold;
            if (usePepQValue)
            {
                _filterType = FilterType.PepQValue;
                _filterThreshold = CommonParameters.PepQValueThreshold;
            }
            else
            {
                _filterType = FilterType.QValue;
                _filterThreshold = CommonParameters.QValueThreshold;
            }
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            ProteinScoringAndFdrResults myAnalysisResults = new ProteinScoringAndFdrResults(this);
            ScoreProteinGroups(ProteinGroups, _FilteredPsms);
            myAnalysisResults.SortedAndScoredProteinGroups = DoProteinFdr(ProteinGroups);

            if (LabelPeptidesUniqueToOneGroupAfterFdr)
            {
                LabelPeptidesUniqueToOneGroup(myAnalysisResults.SortedAndScoredProteinGroups, TreatModPeptidesAsDifferentPeptides);
            }

            return myAnalysisResults;
        }

        /// <summary>
        /// Sets each group's <see cref="ProteinGroup.UniquePeptides"/> to its peptides whose sequence no other group
        /// lists, so a peptide is unique to a protein group when that group is the only one containing it, however many
        /// proteins in the database carry it. Calmodulin (three genes, one sequence) forms one group, and its peptides
        /// are unique to that group.
        /// </summary>
        /// <remarks>
        /// Parsimony chose the proteins from database-level uniqueness, and that is untouched: this runs after the
        /// groups are formed and changes only what they report (the Unique and Shared Peptides columns, their counts,
        /// the per-file tables and the mzIdentML count). It is the meaning
        /// <see cref="Omics.BioPolymerGroup.IBioPolymerGroup.UniqueBioPolymersWithSetMods"/> documents ("not found in any
        /// other biopolymer group"), and the one quantification uses: FlashLFQ counts a peptide toward a group when its
        /// spectra belong to that group alone, and mzLib's isobaric quantification (from smith-chem-wisc/mzLib#1463) when that group alone lists it.
        /// Judged among every group scoring keeps, whatever its q-value, decoy or contaminant status, so a peptide is not
        /// called unique because the group it shares with is above the threshold. Those are the groups FlashLFQ and the
        /// isobaric engine quantify against: a group the two-peptide rule or the picked FDR drops does not count. Each
        /// group counts once, though the picked FDR can list a rescued group once per accession. The sequence is the one
        /// parsimony compares: the full sequence when modified peptides are treated as different, the base sequence
        /// otherwise, and, as in parsimony, per protease: a sequence another group lists only from a different protease
        /// does not make it shared. So with several proteases one sequence can be unique to two groups, one per protease,
        /// and the tables, which do not name the protease, show it unique in both. Quantification does not split by
        /// protease: FlashLFQ and the isobaric engine key a peptide on its sequence alone, so they treat that sequence as
        /// shared between the two groups.
        /// </remarks>
        internal static void LabelPeptidesUniqueToOneGroup(List<ProteinGroup> proteinGroups, bool treatModPeptidesAsDifferentPeptides)
        {
            (string, object) Sequence(IBioPolymerWithSetMods peptide) =>
                (treatModPeptidesAsDifferentPeptides ? peptide.FullSequence : peptide.BaseSequence, peptide.DigestionParams?.DigestionAgent);

            var groupsListingSequence = new Dictionary<(string, object), int>();
            foreach (var proteinGroup in proteinGroups.Distinct<ProteinGroup>(ReferenceEqualityComparer.Instance))
            {
                foreach (var sequence in proteinGroup.AllPeptides.Select(Sequence).Distinct())
                {
                    groupsListingSequence[sequence] = groupsListingSequence.GetValueOrDefault(sequence) + 1;
                }
            }

            foreach (var proteinGroup in proteinGroups)
            {
                proteinGroup.UniquePeptides = new HashSet<IBioPolymerWithSetMods>(
                    proteinGroup.AllPeptides.Where(peptide => groupsListingSequence[Sequence(peptide)] == 1));
            }
        }

        private static string StripDecoyIdentifier(string proteinGroupName, HashSet<string> decoyIdentifiers) //we're keeping only the better scoring protein group for each target/decoy pair. to do that we need to strip decoy from the name temporarily. this is the "top-picked" method
        {
            foreach (var ident in decoyIdentifiers.Where(proteinGroupName.Contains))
                return proteinGroupName.Replace($"{ident}_", "");

            return proteinGroupName;
        }

        protected void ScoreProteinGroups(List<ProteinGroup> proteinGroups, IEnumerable<SpectralMatch> psmList)
        {
            // add each protein groups PSMs
            var peptideToPsmMatching = new Dictionary<IBioPolymerWithSetMods, HashSet<SpectralMatch>>();
            var psmsToScore = _psmsAlreadyFiltered
                ? psmList
                : psmList.FilterByQValue(
                    includeHighQValuePsms: false,
                    qValueThreshold: _filterThreshold,
                    filterAtPeptideLevel: false,
                    filterType: _filterType);

            foreach (var psm in psmsToScore)
            {
                if ((TreatModPeptidesAsDifferentPeptides && psm.FullSequence != null) || (!TreatModPeptidesAsDifferentPeptides && psm.BaseSequence != null))
                {
                    foreach (var pepWithSetMods in psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer))
                    {
                        if (!peptideToPsmMatching.TryGetValue(pepWithSetMods, out HashSet<SpectralMatch> psmsForThisPeptide))
                            peptideToPsmMatching.Add(pepWithSetMods, new HashSet<SpectralMatch> { psm });
                        else
                            psmsForThisPeptide.Add(psm);
                    }
                }
            }

            foreach (var proteinGroup in proteinGroups)
            {
                List<IBioPolymerWithSetMods> pepsToRemove = new();
                foreach (var peptide in proteinGroup.AllPeptides)
                {
                    // build PSM list for scoring
                    if (peptideToPsmMatching.TryGetValue(peptide, out HashSet<SpectralMatch> psms))
                        proteinGroup.AllPsmsBelowOnePercentFDR.UnionWith(psms);
                    else
                        pepsToRemove.Add(peptide);
                }

                proteinGroup.AllPeptides.ExceptWith(pepsToRemove);
                proteinGroup.UniquePeptides.ExceptWith(pepsToRemove);
            }

            // score the group
            foreach (var proteinGroup in proteinGroups)
            {
                proteinGroup.Score();
            }

            if (MergeIndistinguishableProteinGroups)
            {
                // merge protein groups that are indistinguishable after scoring
                var pg = proteinGroups.OrderByDescending(p => p.ProteinGroupScore).ToList();

                // Within each run of groups that tie on score, every set of groups with the same peptides is merged into one.
                int idx = 0;
                while (idx < pg.Count - 1)
                {
                    double score = pg[idx].ProteinGroupScore;
                    if (score != 0 && score == pg[idx + 1].ProteinGroupScore)
                    {
                        // find the end of this score-tied run
                        int clusterEnd = idx + 1;
                        while (clusterEnd < pg.Count && pg[clusterEnd].ProteinGroupScore == score)
                            clusterEnd++;

                        var cluster = pg.GetRange(idx, clusterEnd - idx);

                        // a merge only joins groups with equal sequence sets, so each group's set is fixed for the whole run
                        var sequences = cluster.Select(g => new HashSet<string>(g.AllPeptides.Select(x => x.FullSequence + x.DigestionParams.DigestionAgent))).ToList();

                        // check to make sure they have the same peptides, then merge them
                        for (int i = 0; i < cluster.Count; i++)
                        {
                            if (cluster[i].ProteinGroupScore == 0) continue; // already merged into another group
                            for (int j = 0; j < cluster.Count; j++)
                            {
                                if (i == j) continue;
                                if (cluster[j].ProteinGroupScore == 0) continue; // already merged, skip to avoid reverse-merge zeroing the keeper
                                if (sequences[j].SetEquals(sequences[i]))
                                    cluster[i].MergeProteinGroupWith(cluster[j]);
                            }
                        }

                        idx = clusterEnd;
                    }
                    else
                    {
                        idx++;
                    }
                }
            }

            // remove empty protein groups (peptides were too poor quality or group was merged)
            proteinGroups.RemoveAll(p => p.ProteinGroupScore == 0);

            // calculate sequence coverage
            foreach (var proteinGroup in proteinGroups)
            {
                proteinGroup.CalculateSequenceCoverage();
            }
        }

        protected List<ProteinGroup> ApplyNoOneHitWondersFilter(List<ProteinGroup> proteinGroups)
        {
            if (!NoOneHitWonders)
            {
                return proteinGroups;
            }

            // GroupBy(...).Count() > 1 always builds a full dictionary over all peptides.
            // Distinct().Skip(1).Any() short-circuits as soon as a second distinct value is
            // found, avoiding the allocation for groups with many redundant sequences.
            // Semantics are identical: both return true iff there are >= 2 distinct values.
            if (TreatModPeptidesAsDifferentPeptides)
            {
                return proteinGroups.Where(p => p.AllPeptides.Select(x => x.FullSequence).Distinct().Skip(1).Any()).ToList();
            }

            return proteinGroups.Where(p => p.AllPeptides.Select(x => x.BaseSequence).Distinct().Skip(1).Any()).ToList();
        }

        protected static void PopulateBestPeptideMetrics(IEnumerable<ProteinGroup> proteinGroups)
        {
            foreach (var pg in proteinGroups)
            {
                if (pg.AllPsmsBelowOnePercentFDR.Count == 0)
                {
                    pg.BestPeptideScore = 0;
                    pg.BestPeptideQValue = 1;
                    pg.BestPeptidePEP = 1;
                    continue;
                }

                // Cast rather than OfType: an unexpected ISpectralMatch implementation should throw, not
                // vanish from the QValue and PEP minimums.
                var bestPeptidePsms = pg.AllPsmsBelowOnePercentFDR.Cast<SpectralMatch>().ToList();
                pg.BestPeptideScore = bestPeptidePsms.Max(psm => psm.Score);
                pg.BestPeptideQValue = bestPeptidePsms.Min(psm => psm.FdrInfo.QValueNotch);
                pg.BestPeptidePEP = bestPeptidePsms.Min(psm => psm.FdrInfo.PEP);
            }
        }

        protected virtual List<ProteinGroup> DoProteinFdr(List<ProteinGroup> proteinGroups)
        {
            proteinGroups = ApplyNoOneHitWondersFilter(proteinGroups);

            // Do Classic protein FDR (all targets, all decoys)
            // order protein groups based on filter type
            var sortedProteinGroups = SortProteinGroupsByFilterType(proteinGroups);
            AssignQValuesToProteins(sortedProteinGroups);

            // Do "Picked" protein FDR
            // adapted from "A Scalable Approach for Protein False Discovery Rate Estimation in Large Proteomic Data Sets" ~ MCP, 2015, Savitski
            // pair decoys and targets by accession
            // then use the best peptide metric (QValue or PEP) as the score for the protein group
            Dictionary<string, List<ProteinGroup>> accessionToProteinGroup = new Dictionary<string, List<ProteinGroup>>();
            foreach (var pg in proteinGroups)
            {
                foreach (var protein in pg.Proteins)
                {
                    string stippedAccession = StripDecoyIdentifier(protein.Accession, _decoyIdentifiers);

                    if (accessionToProteinGroup.TryGetValue(stippedAccession, out List<ProteinGroup> groups))
                    {
                        groups.Add(pg);
                    }
                    else
                    {
                        accessionToProteinGroup.Add(stippedAccession, new List<ProteinGroup> { pg });
                    }
                }

            }

            PopulateBestPeptideMetrics(proteinGroups);

            // pick the best for each paired accession based on filter type
            // this compares target-decoy pairs for each protein and saves the best scoring group
            List<ProteinGroup> rescuedProteins = new List<ProteinGroup>();
            // Collect all groups to remove across all accessions first, then do a single O(N) pass.
            // The original per-iteration Except(pgList).ToList() rebuilt the full list on every
            // loop iteration — O(M*N) total. Batching reduces this to O(M+N).
            // HashSet<ProteinGroup> uses reference equality by default (ProteinGroup does not
            // override GetHashCode/Equals(object)), which matches the original Except semantics.
            var toRemove = new HashSet<ProteinGroup>();
            foreach (var accession in accessionToProteinGroup)
            {
                if (accession.Value.Count > 1)
                {
                    var pgList = SortProteinGroupsByFilterType(accession.Value);
                    var pgToUse = pgList.First(); // pick best (lowest QValue or lowest PEP) and remove the rest
                    pgList.Remove(pgToUse);
                    rescuedProteins.AddRange(pgList); // save the remaining protein groups
                    foreach (var pg in pgList)
                        toRemove.Add(pg);
                }
            }
            proteinGroups = proteinGroups.Where(pg => !toRemove.Contains(pg)).ToList();

            sortedProteinGroups = SortProteinGroupsByFilterType(proteinGroups);
            AssignQValuesToProteins(sortedProteinGroups);

            // Rescue the removed TARGET proteins that have the classic protein fdr.
            // This isn't super transparent, but the "Picked" TDS (target-decoy strategy) does a good job of removing a lot of decoys from accumulating in large datasets.
            // It sounds biased, but the Picked TDS is actually necessary to keep the chance of a random assignment being assigned as a target or a decoy at 50:50.
            // The targets that we're re-adding have higher q-values (for their score) from the Classic TDS than the Picked TDS (the classic is conservative).
            // If we add the decoys, it will raise questions on if the FDR is being calculated correctly, 
            // because lots of decoys (which are out-competed in the Picked TDS) will be written with high(ish) scores
            // so really, we're only outputting targets for a cleanliness of output (but the decoys are still there for the classic TDS)
            // TL;DR 99% of the protein output is from the Picked TDS, but a small fraction is from the Classic TDS.
            // A group that lost under several of its accessions was saved once per loss; add each one back once.
            sortedProteinGroups.AddRange(rescuedProteins.Where(x => !x.IsDecoy).Distinct<ProteinGroup>(ReferenceEqualityComparer.Instance));

            return sortedProteinGroups.OrderBy(b => b.QValue).ToList();
        }
        /// <summary>
        /// Sorts protein groups based on the filter type.
        /// QValue: Sort by best peptide Q-value (ascending), then by best peptide score (descending) - higher scores are better
        /// PepQValue: Sort by best peptide PEP (ascending) - lower PEP is better
        /// </summary>
        private List<ProteinGroup> SortProteinGroupsByFilterType(IEnumerable<ProteinGroup> proteinGroups)
        {
            return _filterType switch
            {
                FilterType.PepQValue => proteinGroups
                    .OrderBy(p => p.BestPeptidePEP)
                    .ThenByDescending(p => p.BestPeptideScore)
                    .ToList(),
                _ => proteinGroups
                    .OrderBy(b => b.BestPeptideQValue)
                    .ThenByDescending(p => p.BestPeptideScore)
                    .ToList()
            };
        }

        private void AssignQValuesToProteins(List<ProteinGroup> sortedProteinGroups)
        {
            // sum targets and decoys
            int cumulativeTarget = 0;
            int cumulativeDecoy = 0;

            foreach (var proteinGroup in sortedProteinGroups)
            {
                if (proteinGroup.IsDecoy)
                {
                    cumulativeDecoy++;
                }
                else
                {
                    cumulativeTarget++;
                }
                proteinGroup.CumulativeTarget = cumulativeTarget;
                proteinGroup.CumulativeDecoy = cumulativeDecoy;
            }

            //calculate q-values, assuming that q-values can never decrease with decreasing score
            double maxQValue = double.PositiveInfinity;
            for (int i = sortedProteinGroups.Count - 1; i >= 0; i--)
            {
                ProteinGroup proteinGroup = sortedProteinGroups[i];
                double currentQValue = 1d * proteinGroup.CumulativeDecoy / proteinGroup.CumulativeTarget;
                if (currentQValue < maxQValue)
                {
                    maxQValue = currentQValue;
                }
                proteinGroup.QValue = maxQValue;
            }
        }
    }
}
