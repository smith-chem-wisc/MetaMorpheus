using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class ProteinScoringAndFdrEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly IEnumerable<Psm> newPsms;
        private readonly List<MassDiffAcceptor> searchModes;
        private readonly bool noOneHitWonders;
        private readonly bool treatModPeptidesAsDifferentPeptides;
        private List<ProteinGroup> proteinGroups;
        private readonly bool scoreAndMergeProteinGroups;

        #endregion Private Fields

        #region Public Constructors

        public ProteinScoringAndFdrEngine(List<ProteinGroup> proteinGroups, List<Psm> newPsms, List<MassDiffAcceptor> searchModes, bool noOneHitWonders, bool treatModPeptidesAsDifferentPeptides, bool scoreAndMergeProteinGroups, List<string> nestedIds) : base(nestedIds)
        {
            this.newPsms = newPsms;
            this.searchModes = searchModes;
            this.proteinGroups = proteinGroups;
            this.noOneHitWonders = noOneHitWonders;
            this.treatModPeptidesAsDifferentPeptides = treatModPeptidesAsDifferentPeptides;
            this.scoreAndMergeProteinGroups = scoreAndMergeProteinGroups;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            ProteinScoringAndFdrResults myAnalysisResults = new ProteinScoringAndFdrResults(this);
            Status("Running protein scoring and FDR engine!", nestedIds);

            if(scoreAndMergeProteinGroups)
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
        }

        private List<ProteinGroup> DoProteinFdr(List<ProteinGroup> proteinGroups)
        {
            Status("Calculating protein FDR...", nestedIds);

            if (noOneHitWonders)
            {
                if (treatModPeptidesAsDifferentPeptides)
                    proteinGroups = proteinGroups.Where(p => p.isDecoy || new HashSet<string>(p.AllPeptides.Select(x => x.Sequence)).Count > 1).ToList();
                else
                    proteinGroups = proteinGroups.Where(p => p.isDecoy || new HashSet<string>(p.AllPeptides.Select(x => x.BaseSequence)).Count > 1).ToList();
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