using Chemistry;
using MassSpectrometry;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.ModernSearch
{
    public class ModernSearchEngine : MetaMorpheusEngine
    {
        protected const int FragmentBinsPerDalton = 1000;
        protected readonly List<int>[] FragmentIndex;
        protected readonly PeptideSpectralMatch[] PeptideSpectralMatches;
        protected readonly Ms2ScanWithSpecificMass[] ListOfSortedms2Scans;
        protected readonly List<CompactPeptide> PeptideIndex;
        protected readonly List<ProductType> ProductTypes;
        protected readonly int CurrentPartition;
        protected readonly MassDiffAcceptor MassDiffAcceptor;
        protected readonly List<DissociationType> DissociationTypes;
        protected readonly double MaximumMassThatFragmentIonScoreIsDoubled;

        public ModernSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans,
                List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition,
                CommonParameters commonParameters, MassDiffAcceptor massDiffAcceptor,
                double maximumMassThatFragmentIonScoreIsDoubled, List<string> nestedIds)
            : base(commonParameters, nestedIds)
        {
            PeptideSpectralMatches = globalPsms;
            ListOfSortedms2Scans = listOfSortedms2Scans;
            PeptideIndex = peptideIndex;
            FragmentIndex = fragmentIndex;
            ProductTypes = lp;
            CurrentPartition = currentPartition + 1;
            MassDiffAcceptor = massDiffAcceptor;
            DissociationTypes = DetermineDissociationType(lp);
            MaximumMassThatFragmentIonScoreIsDoubled = maximumMassThatFragmentIonScoreIsDoubled;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing modern search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;
            if (CommonParameters.CalculateEValue)
            {
                byteScoreCutoff = 1;
            }

            Parallel.ForEach(Partitioner.Create(0, ListOfSortedms2Scans.Length),
                new ParallelOptions { MaxDegreeOfParallelism = base.CommonParameters.MaxThreadsToUsePerFile },
                (range, loopState) =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        loopState.Stop();
                        return;
                    }

                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = ListOfSortedms2Scans[i];

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(scan);

                    // get allowed theoretical masses from the known experimental mass
                    // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)
                    // this is just PRELIMINARY precursor-mass filtering
                    // additional checks are made later to ensure that the theoretical precursor mass is acceptable
                    var notches = MassDiffAcceptor.GetAllowedPrecursorMassIntervals(scan.PrecursorMass);

                    double lowestMassPeptideToLookFor = double.NegativeInfinity;
                    double highestMassPeptideToLookFor = double.PositiveInfinity;

                    double largestMassDiff = notches.Max(p => p.AllowedInterval.Maximum);
                    double smallestMassDiff = notches.Min(p => p.AllowedInterval.Minimum);

                    if (!double.IsInfinity(largestMassDiff))
                    {
                        double largestOppositeMassDiff = -1 * (notches.Max(p => p.AllowedInterval.Maximum) - scan.PrecursorMass);
                        lowestMassPeptideToLookFor = scan.PrecursorMass + largestOppositeMassDiff;
                    }
                    if (!double.IsNegativeInfinity(smallestMassDiff))
                    {
                        double smallestOppositeMassDiff = -1 * (notches.Min(p => p.AllowedInterval.Minimum) - scan.PrecursorMass);
                        highestMassPeptideToLookFor = scan.PrecursorMass + smallestOppositeMassDiff;
                    }

                    // first-pass scoring
                    IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor);

                    // done with indexed scoring; refine scores and create PSMs
                    foreach (var id in idsOfPeptidesPossiblyObserved)
                    {
                        var compactPeptide = PeptideIndex[id];

                        var productMasses = compactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(ProductTypes);
                        Array.Sort(productMasses);

                        double scanPrecursorMass = scan.PrecursorMass;

                        var thisScore = CalculatePeptideScoreOld(scan.TheScan, CommonParameters.ProductMassTolerance, productMasses, scanPrecursorMass, DissociationTypes, CommonParameters.AddCompIons, 0);
                        int notch = MassDiffAcceptor.Accepts(scan.PrecursorMass, compactPeptide.MonoisotopicMassIncludingFixedMods);

                        bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;
                        bool scoreImprovement = PeptideSpectralMatches[i] == null || (thisScore - PeptideSpectralMatches[i].RunnerUpScore) > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

                        if (meetsScoreCutoff && scoreImprovement || base.CommonParameters.CalculateEValue)
                        {
                            if (PeptideSpectralMatches[i] == null)
                            {
                                PeptideSpectralMatches[i] = new PeptideSpectralMatch(compactPeptide, notch, thisScore, i, scan, CommonParameters.DigestionParams);
                            }
                            else
                            {
                                PeptideSpectralMatches[i].AddOrReplace(compactPeptide, thisScore, notch, CommonParameters.ReportAllAmbiguity);
                            }

                            if (base.CommonParameters.CalculateEValue)
                            {
                                PeptideSpectralMatches[i].AllScores.Add(thisScore);
                            }
                        }
                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedms2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing modern search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
                    }
                }
            });

            // remove peptides below the score cutoff that were stored to calculate expectation values
            if (base.CommonParameters.CalculateEValue)
            {
                for (int i = 0; i < PeptideSpectralMatches.Length; i++)
                {
                    if (PeptideSpectralMatches[i] != null && PeptideSpectralMatches[i].Score < CommonParameters.ScoreCutoff)
                    {
                        PeptideSpectralMatches[i] = null;
                    }
                }
            }

            return new MetaMorpheusEngineResults(this);
        }

        protected List<int> GetBinsToSearch(Ms2ScanWithSpecificMass scan)
        {
            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();
            foreach (var peakMz in scan.TheScan.MassSpectrum.XArray)
            {
                // assume charge state 1 to calculate mass tolerance
                double experimentalFragmentMass = ClassExtensions.ToMass(peakMz, 1);

                // get theoretical fragment bins within mass tolerance
                int obsFragmentFloorMass = (int)Math.Floor((CommonParameters.ProductMassTolerance.GetMinimumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);
                int obsFragmentCeilingMass = (int)Math.Ceiling((CommonParameters.ProductMassTolerance.GetMaximumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);

                // prevents double-counting peaks close in m/z and lower-bound out of range exceptions
                if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)
                {
                    obsFragmentFloorMass = obsPreviousFragmentCeilingMz;
                }
                obsPreviousFragmentCeilingMz = obsFragmentCeilingMass + 1;

                // prevent upper-bound index out of bounds errors;
                // lower-bound is handled by the previous "if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)" statement
                if (obsFragmentCeilingMass >= FragmentIndex.Length)
                {
                    obsFragmentCeilingMass = FragmentIndex.Length - 1;

                    if (obsFragmentFloorMass >= FragmentIndex.Length)
                    {
                        obsFragmentFloorMass = FragmentIndex.Length - 1;
                    }
                }

                // search mass bins within a tolerance
                for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                {
                    if (FragmentIndex[fragmentBin] != null)
                    {
                        binsToSearch.Add(fragmentBin);
                    }
                }

                // add complementary ions
                if (base.CommonParameters.AddCompIons)
                {
                    //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut and just straight up adding the bins assuming that they're z=1
                    foreach (DissociationType dissociationType in DissociationTypes)
                    {
                        if (ComplementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift))
                        {
                            protonMassShift = ClassExtensions.ToMass(protonMassShift, 1);
                            int compFragmentFloorMass = (int)Math.Round(((scan.PrecursorMass + protonMassShift) * FragmentBinsPerDalton)) - obsFragmentCeilingMass;
                            int compFragmentCeilingMass = (int)Math.Round(((scan.PrecursorMass + protonMassShift) * FragmentBinsPerDalton)) - obsFragmentFloorMass;

                            // prevent index out of bounds errors
                            if (compFragmentCeilingMass >= FragmentIndex.Length)
                            {
                                compFragmentCeilingMass = FragmentIndex.Length - 1;

                                if (compFragmentFloorMass >= FragmentIndex.Length)
                                {
                                    compFragmentFloorMass = FragmentIndex.Length - 1;
                                }
                            }
                            if (compFragmentFloorMass < 0)
                            {
                                compFragmentFloorMass = 0;
                            }

                            for (int fragmentBin = compFragmentFloorMass; fragmentBin <= compFragmentCeilingMass; fragmentBin++)
                            {
                                if (FragmentIndex[fragmentBin] != null)
                                {
                                    binsToSearch.Add(fragmentBin);
                                }
                            }
                        }
                        else
                        {
                            throw new NotImplementedException();
                        }
                    }
                }
            }
            return binsToSearch;
        }

        private int BinarySearchBinForPrecursorIndex(List<int> peptideIdsInThisBin, double peptideMassToLookFor)
        {
            int m = 0;
            int l = 0;
            int r = peptideIdsInThisBin.Count - 1;

            // binary search in the fragment bin for precursor mass
            while (l <= r)
            {
                m = l + ((r - l) / 2);

                if (r - l < 2)
                {
                    break;
                }
                if (PeptideIndex[peptideIdsInThisBin[m]].MonoisotopicMassIncludingFixedMods < peptideMassToLookFor)
                {
                    l = m + 1;
                }
                else
                {
                    r = m - 1;
                }
            }
            if (m > 0)
            {
                m--;
            }
            return m;
        }

        private void IndexedScoring(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = FragmentIndex[binsToSearch[i]];

                //get index for minimum monoisotopic allowed
                int lowestPeptideMassIndex = Double.IsInfinity(lowestMassPeptideToLookFor) ? 0 : BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, lowestMassPeptideToLookFor);

                // get index for highest mass allowed
                int highestPeptideMassIndex = peptideIdsInThisBin.Count - 1;

                if (!Double.IsInfinity(highestMassPeptideToLookFor))
                {
                    highestPeptideMassIndex = BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, highestMassPeptideToLookFor);

                    for (int j = highestPeptideMassIndex; j < peptideIdsInThisBin.Count; j++)
                    {
                        int nextId = peptideIdsInThisBin[j];
                        var nextPep = PeptideIndex[nextId];
                        if (nextPep.MonoisotopicMassIncludingFixedMods < highestMassPeptideToLookFor)
                        {
                            highestPeptideMassIndex = j;
                        }
                        else
                        {
                            break;
                        }
                    }
                }

                // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                {
                    int id = peptideIdsInThisBin[j];
                    scoringTable[id]++;

                    // add possible search results to the hashset of id's
                    if (scoringTable[id] == byteScoreCutoff && MassDiffAcceptor.Accepts(scanPrecursorMass, PeptideIndex[id].MonoisotopicMassIncludingFixedMods) >= 0)
                    {
                        idsOfPeptidesPossiblyObserved.Add(id);
                    }
                }

                if (MaximumMassThatFragmentIonScoreIsDoubled > 0)
                {
                    for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                    {
                        if (j < MaximumMassThatFragmentIonScoreIsDoubled * FragmentBinsPerDalton)
                        {
                            int id = peptideIdsInThisBin[j];
                            scoringTable[id]++;

                            // add possible search results to the hashset of id's
                            if (scoringTable[id] == byteScoreCutoff && MassDiffAcceptor.Accepts(scanPrecursorMass, PeptideIndex[id].MonoisotopicMassIncludingFixedMods) >= 0)
                            {
                                idsOfPeptidesPossiblyObserved.Add(id);
                            }
                        }
                    }
                }
            }
        }
    }
}