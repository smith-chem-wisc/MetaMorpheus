using Chemistry;
using EngineLayer.ModernSearch;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.TruncationSearch
{
    public class TruncationSearchEngine : ModernSearchEngine
    {
        public TruncationSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
            MassDiffAcceptor massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled, List<string> nestedIds)
            : base(globalPsms, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, fileSpecificParameters, massDiffAcceptor,
                  maximumMassThatFragmentIonScoreIsDoubled, nestedIds)
        {

        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing truncation search... " +
                CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();

            Parallel.ForEach(threads, (scanIndex) =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>(PeptideIndex.Count);
                List<Product> peptideTheorProducts = new List<Product>();

                for (; scanIndex < ListOfSortedMs2Scans.Length; scanIndex += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        return;
                    }

                    Ms2ScanWithSpecificMass scan = ListOfSortedMs2Scans[scanIndex];

                    // do a fast rough first-pass scoring for this scan
                    IndexScoreScan(scan, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, CommonParameters.DissociationType);

                    // take indexed-scored peptides and re-score them using the more accurate but slower scoring algorithm
                    FineScorePeptides(idsOfPeptidesPossiblyObserved, scan, scanIndex, scoringTable, CommonParameters.DissociationType, peptideTheorProducts);

                    //report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing truncation search... " +
                            CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
                    }
                }
            });

            foreach (PeptideSpectralMatch psm in PeptideSpectralMatches.Where(p => p != null))
            {
                psm.ResolveAllAmbiguities();
            }

            DoTruncationAnalysis(PeptideSpectralMatches, ListOfSortedMs2Scans);

            return new MetaMorpheusEngineResults(this);
        }

        private void DoTruncationAnalysis(PeptideSpectralMatch[] psms, Ms2ScanWithSpecificMass[] ms2scans)
        {
            string truncationLabel = "Truncation";

            if (!ProteaseDictionary.Dictionary.ContainsKey(truncationLabel))
            {
                ProteaseDictionary.Dictionary.Add(truncationLabel, new Protease(truncationLabel, CleavageSpecificity.None, "", "", new List<DigestionMotif>()));
            }

            DigestionParams d = new DigestionParams(truncationLabel);
            double lowestAminoAcidMass = Proteomics.AminoAcidPolymer.Residue.ResidueMonoisotopicMass['G'];

            MassDiffAcceptor massDiffAcceptor = new DotMassDiffAcceptor("",
                new List<double> {
                    -1 * Constants.C13MinusC12,
                    0,
                    1 * Constants.C13MinusC12
                },
            CommonParameters.PrecursorMassTolerance);

            // mass difference must be negative and must be loss of at least one amino acid (G is smallest amino acid)
            List<PeptideSpectralMatch> possiblyTruncatedPsms = psms.Where(p => p != null && p.FullSequence != null && p.PeptideMonisotopicMass.HasValue &&
                p.PeptideMonisotopicMass.Value - p.ScanPrecursorMass > lowestAminoAcidMass).ToList();

            Parallel.ForEach(Partitioner.Create(0, possiblyTruncatedPsms.Count),
                new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                (range, loopState) =>
                {
                    List<Product> peptideTheorProducts = new List<Product>();

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops)
                        {
                            return;
                        }

                        PeptideSpectralMatch psm = possiblyTruncatedPsms[i];
                        Ms2ScanWithSpecificMass scan = ms2scans[psm.ScanIndex];
                        PeptideWithSetModifications pep = psm.BestMatchingPeptides.First().Peptide;

                        // check for c-terminal truncation
                        for (int r = pep.OneBasedStartResidueInProtein; r < pep.OneBasedEndResidueInProtein; r++)
                        {
                            var mods = pep.AllModsOneIsNterminus.Where(p => p.Key <= r + 1).ToDictionary(p => p.Key, p => p.Value);

                            PeptideWithSetModifications truncatedCTerm = new PeptideWithSetModifications(pep.Protein, d, pep.OneBasedStartResidueInProtein, r,
                                    CleavageSpecificity.Semi, truncationLabel, pep.MissedCleavages, mods, pep.NumFixedMods);

                            int notch = massDiffAcceptor.Accepts(truncatedCTerm.MonoisotopicMass, psm.ScanPrecursorMass);
                            if (notch >= 0)
                            {
                                truncatedCTerm.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, peptideTheorProducts);
                                List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, peptideTheorProducts, CommonParameters);
                                double thisScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan, matchedIons);
                                bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;

                                if (meetsScoreCutoff && thisScore > psm.Score)
                                {
                                    psm.AddOrReplace(truncatedCTerm, thisScore, 0, CommonParameters.ReportAllAmbiguity, matchedIons, 0);
                                    psm.ResolveAllAmbiguities();
                                }
                            }
                        }

                        // check for n-terminal truncation
                        for (int r = pep.OneBasedEndResidueInProtein; r > pep.OneBasedStartResidueInProtein; r--)
                        {
                            var mods = pep.AllModsOneIsNterminus.Where(p => p.Key >= r + 1).ToDictionary(p => p.Key - r + 1, p => p.Value);

                            PeptideWithSetModifications truncatedNTerm = new PeptideWithSetModifications(pep.Protein, d, r, pep.OneBasedEndResidueInProtein,
                                    CleavageSpecificity.Semi, truncationLabel, pep.MissedCleavages, mods, pep.NumFixedMods);

                            int notch = massDiffAcceptor.Accepts(truncatedNTerm.MonoisotopicMass, psm.ScanPrecursorMass);
                            if (notch >= 0)
                            {
                                truncatedNTerm.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, peptideTheorProducts);
                                List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, peptideTheorProducts, CommonParameters);
                                double thisScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan, matchedIons);
                                bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff && thisScore >= psm.Score;

                                if (meetsScoreCutoff && thisScore > psm.Score)
                                {
                                    psm.AddOrReplace(truncatedNTerm, thisScore, 0, CommonParameters.ReportAllAmbiguity, matchedIons, 0);
                                    psm.ResolveAllAmbiguities();
                                }
                            }
                        }
                    }
                });

            for (int i = 0; i < psms.Length; i++)
            {
                var psm = psms[i];

                if (psm == null || psm.FullSequence == null || psm.BestMatchingPeptides.First().Peptide.PeptideDescription != truncationLabel)
                {
                    psms[i] = null;
                }
            }
        }
    }
}
