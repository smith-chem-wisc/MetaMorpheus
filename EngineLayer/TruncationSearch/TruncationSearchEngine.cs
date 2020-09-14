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
        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

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
                    DoTruncationAnalysis(idsOfPeptidesPossiblyObserved, scan, scanIndex, lowestAminoAcidMass, peptideTheorProducts, truncationLabel,
                        d, massDiffAcceptor);

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

            return new MetaMorpheusEngineResults(this);
        }

        private void DoTruncationAnalysis(List<int> idsOfPeptidesPossiblyObserved, Ms2ScanWithSpecificMass scan, int scanIndex, double lowestAminoAcidMass,
            List<Product> peptideTheorProducts, string truncationLabel, DigestionParams digestionParams, MassDiffAcceptor massDiffAcceptor)
        {
            foreach (var id in idsOfPeptidesPossiblyObserved)
            {
                PeptideWithSetModifications pep = PeptideIndex[id];

                bool possiblyTruncated = pep.MonoisotopicMass - scan.PrecursorMass > lowestAminoAcidMass;

                if (!possiblyTruncated)
                {
                    continue;
                }

                // check for c-terminal truncation
                double truncatedCTermMass = WaterMonoisotopicMass;

                for (int r = pep.OneBasedStartResidueInProtein; r < pep.OneBasedEndResidueInProtein; r++)
                {
                    char aa = pep[r - pep.OneBasedStartResidueInProtein];
                    truncatedCTermMass += Proteomics.AminoAcidPolymer.Residue.ResidueMonoisotopicMass[aa];

                    int notch = massDiffAcceptor.Accepts(truncatedCTermMass, scan.PrecursorMass);
                    if (notch >= 0)
                    {
                        var mods = pep.AllModsOneIsNterminus.Where(p => p.Key <= r + 1).ToDictionary(p => p.Key, p => p.Value);

                        PeptideWithSetModifications truncatedCTerm = new PeptideWithSetModifications(pep.Protein, digestionParams, pep.OneBasedStartResidueInProtein, r,
                                CleavageSpecificity.Semi, truncationLabel, pep.MissedCleavages, mods, pep.NumFixedMods);

                        FineScoreTruncatedPeptide(truncatedCTerm, scan, scanIndex, peptideTheorProducts);
                    }
                }

                // check for n-terminal truncation
                double truncatedNTermMass = WaterMonoisotopicMass;

                for (int r = pep.OneBasedEndResidueInProtein; r > pep.OneBasedStartResidueInProtein; r--)
                {
                    char aa = pep[r - pep.OneBasedStartResidueInProtein];
                    truncatedNTermMass += Proteomics.AminoAcidPolymer.Residue.ResidueMonoisotopicMass[aa];

                    int notch = massDiffAcceptor.Accepts(truncatedNTermMass, scan.PrecursorMass);

                    if (notch >= 0)
                    {
                        //double modMass = pep.AllModsOneIsNterminus
                        var mods = pep.AllModsOneIsNterminus.Where(p => p.Key >= r + 1).ToDictionary(p => p.Key - r + 1, p => p.Value);

                        PeptideWithSetModifications truncatedNTerm = new PeptideWithSetModifications(pep.Protein, digestionParams, r, pep.OneBasedEndResidueInProtein,
                                    CleavageSpecificity.Semi, truncationLabel, pep.MissedCleavages, mods, pep.NumFixedMods);
                        double tempMass = truncatedNTerm.MonoisotopicMass;

                        if (truncatedNTerm.PreviousAminoAcid == 'M' && truncatedNTerm.OneBasedStartResidueInProtein == 2)
                        {
                            continue;
                        }

                        FineScoreTruncatedPeptide(truncatedNTerm, scan, scanIndex, peptideTheorProducts);
                    }
                }
            }
        }

        private PeptideSpectralMatch FineScoreTruncatedPeptide(PeptideWithSetModifications peptide, Ms2ScanWithSpecificMass scan, int scanIndex, List<Product> peptideTheorProducts)
        {
            peptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, peptideTheorProducts);

            List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan, peptideTheorProducts, CommonParameters);

            double thisScore = CalculatePeptideScore(scan.TheScan, matchedIons);
            int notch = MassDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMass);

            bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;
            bool scoreImprovement = PeptideSpectralMatches[scanIndex] == null || (thisScore - PeptideSpectralMatches[scanIndex].RunnerUpScore) > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

            if (meetsScoreCutoff && scoreImprovement)
            {
                if (PeptideSpectralMatches[scanIndex] == null)
                {
                    PeptideSpectralMatches[scanIndex] = new PeptideSpectralMatch(peptide, notch, thisScore, scanIndex, scan, CommonParameters, matchedIons);
                }
                else
                {
                    PeptideSpectralMatches[scanIndex].AddOrReplace(peptide, thisScore, notch, CommonParameters.ReportAllAmbiguity, matchedIons, 0);
                }
            }

            return PeptideSpectralMatches[scanIndex];
        }
    }
}
