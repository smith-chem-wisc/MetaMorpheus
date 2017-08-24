﻿using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class ProteinScoringAndFdrEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly IEnumerable<Psm> newPsms;
        private readonly List<MassDiffAcceptor> searchModes;
        private readonly bool noOneHitWonders;
        private List<ProteinGroup> proteinGroups;

        #endregion Private Fields

        #region Public Constructors

        public ProteinScoringAndFdrEngine(List<ProteinGroup> proteinGroups, List<Psm> newPsms, List<MassDiffAcceptor> searchModes, bool noOneHitWonders, List<string> nestedIds) : base(nestedIds)
        {
            this.newPsms = newPsms;
            this.searchModes = searchModes;
            this.proteinGroups = proteinGroups;
            this.noOneHitWonders = noOneHitWonders;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            ProteinScoringAndFdrResults myAnalysisResults = new ProteinScoringAndFdrResults(this);
            Status("Running protein scoring and FDR engine!", nestedIds);

            ScoreProteinGroups(proteinGroups, newPsms);
            myAnalysisResults.sortedAndScoredProteinGroups = DoProteinFdr(proteinGroups);

            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private void ScoreProteinGroups(List<ProteinGroup> proteinGroups, IEnumerable<Psm> psmList)
        {
            Status("Scoring protein groups...", nestedIds);

            // add each protein groups PSMs
            var peptideToPsmMatching = new Dictionary<PeptideWithSetModifications, HashSet<Psm>>();
            foreach (var psm in psmList)
            {
                if (psm.FdrInfo.QValue <= 0.01)
                {
                    foreach (var pepWithSetMods in psm.MostProbableProteinInfo.PeptidesWithSetModifications)
                    {
                        if (!peptideToPsmMatching.TryGetValue(pepWithSetMods, out HashSet<Psm> psmsForThisPeptide))
                            peptideToPsmMatching.Add(pepWithSetMods, new HashSet<Psm> { psm });
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
                    if (peptideToPsmMatching.TryGetValue(peptide, out HashSet<Psm> psms))
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
                    if (sharedPepWithProteinGroups.TryGetValue(sharedPep.Sequence, out HashSet<ProteinGroup> v))
                        v.Add(proteinGroup);
                    else
                        sharedPepWithProteinGroups.Add(sharedPep.Sequence, new HashSet<ProteinGroup> { proteinGroup });

                    if (seqToPeptide.TryGetValue(sharedPep.Sequence, out HashSet<PeptideWithSetModifications> v1))
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

        private List<ProteinGroup> DoProteinFdr(List<ProteinGroup> proteinGroups)
        {
            Status("Calculating protein FDR...", nestedIds);

            if (noOneHitWonders)
            {
                proteinGroups = proteinGroups.Where(p => p.isDecoy || new HashSet<string>(p.AllPeptides.Select(x => x.Sequence)).Count > 1).ToList();
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

        #endregion Private Methods
    }
}