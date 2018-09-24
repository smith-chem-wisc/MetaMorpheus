using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class ProteinScoringAndFdrEngine : MetaMorpheusEngine
    {
        private readonly IEnumerable<PeptideSpectralMatch> NewPsms;
        private readonly bool NoOneHitWonders;
        private readonly bool TreatModPeptidesAsDifferentPeptides;
        private readonly bool MergeIndistinguishableProteinGroups;
        private readonly List<ProteinGroup> ProteinGroups;

        public ProteinScoringAndFdrEngine(List<ProteinGroup> proteinGroups, List<PeptideSpectralMatch> newPsms, bool noOneHitWonders, bool treatModPeptidesAsDifferentPeptides, bool mergeIndistinguishableProteinGroups, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            NewPsms = newPsms;
            ProteinGroups = proteinGroups;
            NoOneHitWonders = noOneHitWonders;
            TreatModPeptidesAsDifferentPeptides = treatModPeptidesAsDifferentPeptides;
            MergeIndistinguishableProteinGroups = mergeIndistinguishableProteinGroups;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            ProteinScoringAndFdrResults myAnalysisResults = new ProteinScoringAndFdrResults(this);

            ScoreProteinGroups(ProteinGroups, NewPsms);
            myAnalysisResults.SortedAndScoredProteinGroups = DoProteinFdr(ProteinGroups);

            return myAnalysisResults;
        }

        private static string StripDecoyIdentifier(string proteinGroupName) //we're keeping only the better scoring protein group for each target/decoy pair. to do that we need to strip decoy from the name temporarily. this is the "top-picked" method
        {
            return proteinGroupName.Contains("DECOY_") ? proteinGroupName.Replace("DECOY_", "") : proteinGroupName;
        }

        private void ScoreProteinGroups(List<ProteinGroup> proteinGroups, IEnumerable<PeptideSpectralMatch> psmList)
        {
            // add each protein groups PSMs
            var peptideToPsmMatching = new Dictionary<PeptideWithSetModifications, HashSet<PeptideSpectralMatch>>();
            foreach (var psm in psmList)
            {
                if (psm.FdrInfo.QValueNotch <= 0.01 && psm.FdrInfo.QValue <= 0.01)
                {
                    if ((TreatModPeptidesAsDifferentPeptides && psm.FullSequence != null) || (!TreatModPeptidesAsDifferentPeptides && psm.BaseSequence != null))
                    {
                        foreach (var pepWithSetMods in psm.BestMatchingPeptides.Select(p => p.Peptide))
                        {
                            if (!peptideToPsmMatching.TryGetValue(pepWithSetMods, out HashSet<PeptideSpectralMatch> psmsForThisPeptide))
                                peptideToPsmMatching.Add(pepWithSetMods, new HashSet<PeptideSpectralMatch> { psm });
                            else
                                psmsForThisPeptide.Add(psm);
                        }
                    }
                }
            }

            foreach (var proteinGroup in proteinGroups)
            {
                List<PeptideWithSetModifications> pepsToRemove = new List<PeptideWithSetModifications>();
                foreach (var peptide in proteinGroup.AllPeptides)
                {
                    // build PSM list for scoring
                    if (peptideToPsmMatching.TryGetValue(peptide, out HashSet<PeptideSpectralMatch> psms))
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
                for (int i = 0; i < (pg.Count - 1); i++)
                {
                    if (pg[i].ProteinGroupScore == pg[i + 1].ProteinGroupScore && pg[i].ProteinGroupScore != 0)
                    {
                        var pgsWithThisScore = pg.Where(p => p.ProteinGroupScore == pg[i].ProteinGroupScore).ToList();

                        // check to make sure they have the same peptides, then merge them
                        foreach (var p in pgsWithThisScore)
                        {
                            var seqs1 = new HashSet<string>(p.AllPeptides.Select(x => x.FullSequence + x.DigestionParams.Protease));
                            var seqs2 = new HashSet<string>(pg[i].AllPeptides.Select(x => x.FullSequence + x.DigestionParams.Protease));

                            if (p != pg[i] && seqs1.SetEquals(seqs2))
                            {
                                pg[i].MergeProteinGroupWith(p);
                            }
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
            if (NoOneHitWonders)
            {
                if (TreatModPeptidesAsDifferentPeptides)
                    proteinGroups = proteinGroups.Where(p => p.IsDecoy || new HashSet<string>(p.AllPeptides.Select(x => x.FullSequence)).Count > 1).ToList();
                else
                    proteinGroups = proteinGroups.Where(p => p.IsDecoy || new HashSet<string>(p.AllPeptides.Select(x => x.BaseSequence)).Count > 1).ToList();
            }

            // pair decoys and targets by accession
            // then use the best peptide notch-QValue as the score for the protein group
            Dictionary<string, List<ProteinGroup>> accessionToProteinGroup = new Dictionary<string, List<ProteinGroup>>();
            foreach (var pg in proteinGroups)
            {
                foreach (var protein in pg.Proteins)
                {
                    string stippedAccession = StripDecoyIdentifier(protein.Accession);

                    if (accessionToProteinGroup.TryGetValue(stippedAccession, out List<ProteinGroup> groups))
                        groups.Add(pg);
                    else
                        accessionToProteinGroup.Add(stippedAccession, new List<ProteinGroup> { pg });
                }

                pg.BestPeptideScore = pg.AllPsmsBelowOnePercentFDR.Max(psm => psm.Score);
                pg.BestPeptideQValue = pg.AllPsmsBelowOnePercentFDR.Min(psm => psm.FdrInfo.QValueNotch);
            }

            // pick the best notch-QValue for each paired accession
            foreach (var accession in accessionToProteinGroup)
            {
                if (accession.Value.Count > 1)
                {
                    var pgList = accession.Value.OrderBy(p => p.BestPeptideQValue).ThenByDescending(p => p.BestPeptideScore).ToList();
                    var pgToUse = pgList.First(); // pick lowest notch QValue and remove the rest
                    pgList.Remove(pgToUse);
                    proteinGroups = proteinGroups.Except(pgList).ToList();
                }
            }

            // order protein groups by notch-QValue
            var sortedProteinGroups = proteinGroups.OrderBy(b => b.BestPeptideQValue).ThenByDescending(p => p.BestPeptideScore).ToList();

            // calculate protein QValues
            int cumulativeTarget = 0;
            int cumulativeDecoy = 0;

            foreach (var proteinGroup in sortedProteinGroups)
            {
                if (proteinGroup.IsDecoy)
                    cumulativeDecoy++;
                else
                    cumulativeTarget++;

                proteinGroup.CumulativeTarget = cumulativeTarget;
                proteinGroup.CumulativeDecoy = cumulativeDecoy;
                proteinGroup.QValue = (double)cumulativeDecoy / cumulativeTarget;
            }

            return sortedProteinGroups;
        }
    }
}