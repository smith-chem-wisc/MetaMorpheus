using EngineLayer.ProteinParsimony;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class ProteinParsimonyEngine : MetaMorpheusEngine
    {
        /// <summary>
        /// All peptides meeting the prefiltering criteria for parsimony (e.g., peptides from non-ambiguous high-confidence PSMs)
        /// </summary>
        private readonly HashSet<PeptideWithSetModifications> _fdrFilteredPeptides;
        private readonly List<PeptideSpectralMatch> _fdrFilteredPsms;
        private readonly List<PeptideSpectralMatch> _allPsms;
        private const double FdrCutoffForParsimony = 0.01;

        /// <summary>
        /// User-selectable option that treats differently-modified forms of a peptide as different peptides for the purposes of parsimony
        /// </summary>
        private readonly bool _treatModPeptidesAsDifferentPeptides;

        public ProteinParsimonyEngine(List<PeptideSpectralMatch> allPsms, bool modPeptidesAreDifferent, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            _treatModPeptidesAsDifferentPeptides = modPeptidesAreDifferent;
            
            if (!allPsms.Any())
            {
                _fdrFilteredPsms = new List<PeptideSpectralMatch>();
            }

            // parsimony will only use non-ambiguous, high-confidence PSMs
            // KEEP decoys and contaminants for use in parsimony!
            if (modPeptidesAreDifferent)
            {
                _fdrFilteredPsms = allPsms.Where(p => p.FullSequence != null && p.FdrInfo.QValue <= FdrCutoffForParsimony && p.FdrInfo.QValueNotch <= FdrCutoffForParsimony).ToList();
            }
            else
            {
                _fdrFilteredPsms = allPsms.Where(p => p.BaseSequence != null && p.FdrInfo.QValue <= FdrCutoffForParsimony && p.FdrInfo.QValueNotch <= FdrCutoffForParsimony).ToList();
            }

            // if PSM is a decoy, add only decoy sequences; same for contaminants
            // peptides to use in parsimony = peptides observed in high-confidence PSMs
            _fdrFilteredPeptides = new HashSet<PeptideWithSetModifications>();

            foreach (var psm in _fdrFilteredPsms)
            {
                if (psm.IsDecoy)
                {
                    foreach (var peptide in psm.BestMatchingPeptides.Select(p => p.Peptide).Where(p => p.Protein.IsDecoy))
                    {
                        _fdrFilteredPeptides.Add(peptide);
                    }
                }
                else if (psm.IsContaminant)
                {
                    foreach (var peptide in psm.BestMatchingPeptides.Select(p => p.Peptide).Where(p => p.Protein.IsContaminant))
                    {
                        _fdrFilteredPeptides.Add(peptide);
                    }
                }
                else // PSM is target
                {
                    foreach (var peptide in psm.BestMatchingPeptides.Select(p => p.Peptide).Where(p => !p.Protein.IsDecoy && !p.Protein.IsContaminant))
                    {
                        _fdrFilteredPeptides.Add(peptide);
                    }
                }
            }

            // we're storing all PSMs (not just FDR-filtered ones) here because we will remove some protein associations 
            // from low-confidence PSMs if they can be explained by a parsimonious protein
            _allPsms = allPsms;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            ProteinParsimonyResults myAnalysisResults = new ProteinParsimonyResults(this);

            myAnalysisResults.ProteinGroups = RunProteinParsimonyEngine();

            return myAnalysisResults;
        }

        /// <summary>
        /// TODO: Summarize parsimony;
        /// Parsimony algorithm based on: https://www.ncbi.nlm.nih.gov/pubmed/14632076 Anal Chem. 2003 Sep 1;75(17):4646-58.
        /// TODO: Note describing that peptide objects with the same sequence are associated with different proteins
        /// </summary>
        private List<ProteinGroup> RunProteinParsimonyEngine()
        {
            // parsimonious list of proteins built by this protein parsimony engine
            HashSet<Protein> parsimoniousProteinList = new HashSet<Protein>();

            // list of peptides that can only be digestion products of one protein in the proteome (considering different protease digestion rules)
            HashSet<PeptideWithSetModifications> uniquePeptides = new HashSet<PeptideWithSetModifications>();

            // if there are no peptides observed, there are no proteins; return an empty list of protein groups
            if (_fdrFilteredPeptides.Count == 0)
            {
                return new List<ProteinGroup>();
            }

            // Parsimony stage 0: create peptide-protein associations if needed because the user wants a modification-agnostic parsimony
            if (!_treatModPeptidesAsDifferentPeptides)
            {
                foreach (var protease in _fdrFilteredPsms.GroupBy(p => p.DigestionParams.Protease))
                {
                    Dictionary<string, List<PeptideSpectralMatch>> sequenceWithPsms = new Dictionary<string, List<PeptideSpectralMatch>>();

                    // for each protease, match the base sequence of each peptide to its PSMs
                    foreach (PeptideSpectralMatch psm in protease)
                    {
                        if (sequenceWithPsms.TryGetValue(psm.BaseSequence, out List<PeptideSpectralMatch> peptidesForThisBaseSequence))
                        {
                            peptidesForThisBaseSequence.Add(psm);
                        }
                        else
                        {
                            sequenceWithPsms[psm.BaseSequence] = new List<PeptideSpectralMatch> { psm };
                        }
                    }

                    // create new peptide-protein associations
                    foreach (var baseSequence in sequenceWithPsms)
                    {
                        var peptidesWithNotchInfo = baseSequence.Value.SelectMany(p => p.BestMatchingPeptides).Distinct().ToList();

                        // if the base seq has >1 PeptideWithSetMods object and has >0 mods, it might need to be matched to new proteins
                        if (peptidesWithNotchInfo.Count > 1 && peptidesWithNotchInfo.Any(p => p.Peptide.NumMods > 0))
                        {
                            // list of proteins along with start/end residue in protein and the # missed cleavages
                            // this is needed to create new PeptideWithSetModification objects
                            var peptideInProteinInfo = new List<Tuple<Protein, DigestionParams, int, int, int, int>>();
                            foreach (var peptide in peptidesWithNotchInfo)
                            {
                                peptideInProteinInfo.Add(new Tuple<Protein, DigestionParams, int, int, int, int>(peptide.Peptide.Protein, peptide.Peptide.DigestionParams,
                                    peptide.Peptide.OneBasedStartResidueInProtein, peptide.Peptide.OneBasedEndResidueInProtein, peptide.Peptide.MissedCleavages, peptide.Notch));
                            }

                            // add the protein associations to the PSM
                            foreach (PeptideSpectralMatch psm in baseSequence.Value)
                            {
                                foreach (var proteinInfo in peptideInProteinInfo)
                                {
                                    var originalPep = psm.BestMatchingPeptides.First().Peptide;
                                    var pep = new PeptideWithSetModifications(proteinInfo.Item1, proteinInfo.Item2, proteinInfo.Item3, proteinInfo.Item4,
                                        originalPep.CleavageSpecificityForFdrCategory, originalPep.PeptideDescription, proteinInfo.Item5, originalPep.AllModsOneIsNterminus,
                                        originalPep.NumFixedMods);
                                    _fdrFilteredPeptides.Add(pep);
                                    psm.AddProteinMatch((proteinInfo.Item6, pep));
                                }
                            }
                        }
                    }
                }
            }

            // Parsimony stage 1: add proteins with unique peptides (for each protease)
            var peptidesGroupedByProtease = _fdrFilteredPeptides.GroupBy(p => p.DigestionParams.Protease);
            foreach (var peptidesForThisProtease in peptidesGroupedByProtease)
            {
                Dictionary<string, List<Protein>> peptideSequenceToProteinsForThisProtease = new Dictionary<string, List<Protein>>();
                Dictionary<string, List<PeptideWithSetModifications>> sequenceToPwsm = new Dictionary<string, List<PeptideWithSetModifications>>();

                foreach (PeptideWithSetModifications peptide in peptidesForThisProtease)
                {
                    string sequence = peptide.BaseSequence;
                    if (_treatModPeptidesAsDifferentPeptides)
                    {
                        //these and next set to full sequence but might be base sequence. treat modified as unique makes sense to use full
                        sequence = peptide.FullSequence;
                    }

                    if (peptideSequenceToProteinsForThisProtease.TryGetValue(sequence, out List<Protein> proteinsForThisPeptideSequence))
                    {
                        proteinsForThisPeptideSequence.Add(peptide.Protein);
                    }
                    else
                    {
                        peptideSequenceToProteinsForThisProtease.Add(sequence, new List<Protein> { peptide.Protein });
                    }

                    if (sequenceToPwsm.TryGetValue(sequence, out List<PeptideWithSetModifications> peptidesForThisSequence))
                    {
                        peptidesForThisSequence.Add(peptide);
                    }
                    else
                    {
                        sequenceToPwsm.Add(sequence, new List<PeptideWithSetModifications> { peptide });
                    }
                }

                foreach (var uniquePeptide in peptideSequenceToProteinsForThisProtease.Where(p => p.Value.Count == 1))
                {
                    // add the protein with the unique peptide to the parsimonious protein list
                    Protein proteinWithUniquePeptideSequence = uniquePeptide.Value.First();
                    parsimoniousProteinList.Add(proteinWithUniquePeptideSequence);

                    // add the unique peptide to the list of unique peptides
                    PeptideWithSetModifications uniquePwsm = sequenceToPwsm[uniquePeptide.Key].First();
                    uniquePeptides.Add(uniquePwsm);
                }
            }

            // Parsimony stage 2: build the peptide-protein matching structure for the parsimony greedy algorithm
            // and remove all peptides observed by proteins with unique peptides
            Dictionary<ParsimonySequence, List<Protein>> peptideSequenceToProteins = new Dictionary<ParsimonySequence, List<Protein>>();

            // this dictionary associates proteins w/ all peptide sequences (list will NOT shrink over time)
            // this is used in case of greedy algorithm ties to figure out which protein has more total peptides observed
            Dictionary<Protein, HashSet<ParsimonySequence>> proteinToPepSeqMatch = new Dictionary<Protein, HashSet<ParsimonySequence>>();

            foreach (var peptide in _fdrFilteredPeptides)
            {
                ParsimonySequence sequence = new ParsimonySequence(peptide, _treatModPeptidesAsDifferentPeptides);
                
                if (peptideSequenceToProteins.TryGetValue(sequence, out List<Protein> proteinsForThisPeptideSequence))
                {
                    proteinsForThisPeptideSequence.Add(peptide.Protein);
                }
                else
                {
                    peptideSequenceToProteins.Add(sequence, new List<Protein> { peptide.Protein });
                }

                if (proteinToPepSeqMatch.TryGetValue(peptide.Protein, out var peptideSequences))
                {
                    peptideSequences.Add(sequence);
                }
                else
                {
                    proteinToPepSeqMatch.Add(peptide.Protein, new HashSet<ParsimonySequence> { sequence });
                }
            }

            // remove the peptides observed by proteins with unique peptides
            HashSet<ParsimonySequence> toRemove = new HashSet<ParsimonySequence>();
            foreach (var seq in peptideSequenceToProteins)
            {
                bool observedAlready = seq.Value.Any(p => parsimoniousProteinList.Contains(p));

                if (observedAlready)
                {
                    toRemove.Add(seq.Key);
                }
            }
            foreach (var sequence in toRemove)
            {
                peptideSequenceToProteins.Remove(sequence);
            }

            if (peptideSequenceToProteins.Any())
            {
                // Parsimony stage 3: greedy algorithm

                // dictionary with proteins as keys and list of associated peptide sequences as the values.
                // this data structure makes parsimony easier because the algorithm can look up a protein's peptides 
                // to remove them from the list of available peptides. this list will shrink as the algorithm progresses
                var algDictionary = new Dictionary<Protein, HashSet<string>>(); 
                var algDictionaryProtease = new Dictionary<Protein, HashSet<ParsimonySequence>>();
                foreach (var kvp in peptideSequenceToProteins)
                {
                    foreach (var protein in kvp.Value)
                    {
                        if (algDictionaryProtease.TryGetValue(protein, out HashSet<ParsimonySequence> peptideSequencesWithProtease))
                        {
                            peptideSequencesWithProtease.Add(kvp.Key);
                        }
                        else
                        {
                            algDictionaryProtease.Add(protein, new HashSet<ParsimonySequence> { kvp.Key });
                        }

                        if (algDictionary.TryGetValue(protein, out HashSet<string> peptideSequences))
                        {
                            peptideSequences.Add(kvp.Key.Sequence);
                        }
                        else
                        {
                            algDictionary.Add(protein, new HashSet<string> { kvp.Key.Sequence });
                        }
                    }
                }

                // *** greedy algorithm loop
                int numNewSeqs = algDictionary.Max(p => p.Value.Count);
                while (numNewSeqs != 0)
                {
                    // gets list of proteins with the most unaccounted-for peptide sequences
                    var possibleBestProteinList = algDictionary.Where(p => p.Value.Count == numNewSeqs).ToList();

                    Protein bestProtein = possibleBestProteinList.First().Key;

                    // may need to select different protein in case of a greedy algorithm tie
                    // the protein with the most total peptide sequences wins in this case (doesn't matter if parsimony has grabbed them or not)
                    if (possibleBestProteinList.Count > 1)
                    {
                        int highestNumTotalPep = proteinToPepSeqMatch[bestProtein].Count;
                        foreach (var kvp in possibleBestProteinList)
                        {
                            if (proteinToPepSeqMatch[kvp.Key].Count > highestNumTotalPep)
                            {
                                highestNumTotalPep = proteinToPepSeqMatch[kvp.Key].Count;
                                bestProtein = kvp.Key;
                            }
                        }
                    }

                    parsimoniousProteinList.Add(bestProtein);

                    // remove observed peptide seqs
                    List<ParsimonySequence> temp = algDictionaryProtease[bestProtein].ToList();
                    foreach (ParsimonySequence peptideSequence in temp)
                    {
                        List<Protein> proteinsWithThisPeptide = peptideSequenceToProteins[peptideSequence];

                        foreach (var protein in proteinsWithThisPeptide)
                        {
                            algDictionary[protein].Remove(peptideSequence.Sequence);
                            algDictionaryProtease[protein].Remove(peptideSequence);
                        }
                    }

                    algDictionary.Remove(bestProtein);
                    algDictionaryProtease.Remove(bestProtein);
                    numNewSeqs = algDictionary.Any() ? algDictionary.Max(p => p.Value.Count) : 0;
                }

                // *** done with greedy algorithm

                // Parsimony stage 4: add back indistinguishable proteins (proteins that have identical peptide sets as parsimonious proteins)
                var allProteinsGroupedByNumPeptides = proteinToPepSeqMatch.GroupBy(p => p.Value.Count);
                var parsimonyProteinsGroupedByNumPeptides = parsimoniousProteinList.GroupBy(p => proteinToPepSeqMatch[p].Count);
                var indistinguishableProteins = new ConcurrentBag<Protein>();

                foreach (var group in allProteinsGroupedByNumPeptides)
                {
                    var parsimonyProteinsWithSameNumPeptides = parsimonyProteinsGroupedByNumPeptides.FirstOrDefault(p => p.Key == group.Key);
                    var list = group.ToList();

                    if (parsimonyProteinsWithSameNumPeptides != null)
                    {
                        Parallel.ForEach(Partitioner.Create(0, list.Count),
                            new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                            (range, loopState) =>
                            {
                                for (int i = range.Item1; i < range.Item2; i++)
                                {
                                    Protein otherProtein = list[i].Key;

                                    foreach (var parsimonyProtein in parsimonyProteinsWithSameNumPeptides)
                                    {
                                        // if the two proteins have the same set of peptide sequences, they're indistinguishable
                                        if (parsimonyProtein != otherProtein && proteinToPepSeqMatch[parsimonyProtein].SetEquals(proteinToPepSeqMatch[otherProtein]))
                                        {
                                            indistinguishableProteins.Add(otherProtein);
                                        }
                                    }
                                }
                            }
                        );
                    }
                }

                foreach (Protein protein in indistinguishableProteins)
                {
                    parsimoniousProteinList.Add(protein);
                }
            }

            // Parsimony stage 5: remove peptide objects that do not have proteins in the parsimonious list
            foreach (PeptideSpectralMatch psm in _allPsms)
            {
                // if this PSM has a protein in the parsimonious list, it removes the proteins NOT in the parsimonious list
                // otherwise, no proteins are removed (i.e., for PSMs that cannot be explained by a parsimonious protein, 
                // no protein associations are removed)
                if (psm.BestMatchingPeptides.Any(p => parsimoniousProteinList.Contains(p.Peptide.Protein)))
                {
                    psm.TrimProteinMatches(parsimoniousProteinList);
                }
            }

            // construct protein groups
            List<ProteinGroup> proteinGroups = ConstructProteinGroups(uniquePeptides);

            // finished with parsimony
            return proteinGroups;
        }

        /// <summary>
        /// Builds protein groups after parsimony. Each protein group will only have 1 protein at this point. 
        /// Indistinguishable protein groups will be merged during scoring for computational efficiency reasons 
        /// (it's easier to tell when two protein groups are identical after they're scored)
        /// </summary>
        private List<ProteinGroup> ConstructProteinGroups(HashSet<PeptideWithSetModifications> uniquePeptides)
        {
            List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<PeptideWithSetModifications>>();

            foreach (var peptide in _fdrFilteredPeptides)
            {
                if (proteinToPeptidesMatching.TryGetValue(peptide.Protein, out HashSet<PeptideWithSetModifications> peptidesHere))
                {
                    peptidesHere.Add(peptide);
                }
                else
                {
                    proteinToPeptidesMatching.Add(peptide.Protein, new HashSet<PeptideWithSetModifications> { peptide });
                }
            }

            foreach (var kvp in proteinToPeptidesMatching)
            {
                var allPeptidesHere = proteinToPeptidesMatching[kvp.Key];
                var uniquePeptidesHere = new HashSet<PeptideWithSetModifications>(allPeptidesHere.Where(p => uniquePeptides.Contains(p)));

                proteinGroups.Add(new ProteinGroup(new HashSet<Protein> { kvp.Key }, allPeptidesHere, uniquePeptidesHere));
            }

            foreach (var proteinGroup in proteinGroups)
            {
                proteinGroup.AllPeptides.RemoveWhere(p => !proteinGroup.Proteins.Contains(p.Protein));
                proteinGroup.DisplayModsOnPeptides = _treatModPeptidesAsDifferentPeptides;
            }

            return proteinGroups;
        }
    }
}
