﻿using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.ModernSearch
{
    public class ModernSearchEngine : MetaMorpheusEngine
    {

        #region Protected Fields

        protected readonly List<int>[] fragmentIndex;
        protected readonly Psm[] globalPsms;
        protected readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;
        protected readonly List<CompactPeptide> peptideIndex;
        protected readonly List<ProductType> lp;
        protected readonly int currentPartition;
        protected readonly CommonParameters CommonParameters;
        protected readonly bool addCompIons;
        protected readonly MassDiffAcceptor massDiffAcceptor;
        protected readonly List<DissociationType> dissociationTypes;
        protected const int fragmentBinsPerDalton = 1000;

        #endregion Protected Fields

        #region Public Constructors

        public ModernSearchEngine(Psm[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition, CommonParameters CommonParameters, bool addCompIons, MassDiffAcceptor massDiffAcceptors, List<string> nestedIds) : base(nestedIds)
        {
            this.globalPsms = globalPsms;
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.fragmentIndex = fragmentIndex;
            this.lp = lp;
            this.currentPartition = currentPartition + 1;
            this.CommonParameters = CommonParameters;
            this.addCompIons = addCompIons;
            this.massDiffAcceptor = massDiffAcceptors;
            this.dissociationTypes = DetermineDissociationType(lp);
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing modern search... " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = Convert.ToByte((int)CommonParameters.ScoreCutoff);

            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUse }, range =>
            {
                byte[] scoringTable = new byte[peptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = listOfSortedms2Scans[i];

                    // filter ms2 fragment peaks by intensity
                    int numFragmentsToUse = scan.NumPeaks;
                    if (CommonParameters.TopNpeaks != null)
                        numFragmentsToUse = CommonParameters.TopNpeaks.Value;

                    var peaks = scan.TheScan.MassSpectrum.FilterByNumberOfMostIntense(numFragmentsToUse).ToList();
                    double largestIntensity = scan.TheScan.MassSpectrum.YofPeakWithHighestY;

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(peaks, largestIntensity, scan.PrecursorMass);

                    // get allowed theoretical masses from the known experimental mass
                    // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)
                    var notches = massDiffAcceptor.GetAllowedPrecursorMassIntervals(scan.PrecursorMass);

                    double lowestMassPeptideToLookFor = Double.NegativeInfinity;
                    double highestMassPeptideToLookFor = Double.PositiveInfinity;

                    double largestMassDiff = notches.Max(p => p.allowedInterval.Maximum);
                    double smallestMassDiff = notches.Min(p => p.allowedInterval.Minimum);

                    if (!Double.IsInfinity(largestMassDiff))
                    {
                        double largestOppositeMassDiff = -1 * (notches.Max(p => p.allowedInterval.Maximum) - scan.PrecursorMass);
                        lowestMassPeptideToLookFor = scan.PrecursorMass + largestOppositeMassDiff;
                    }
                    if (!Double.IsNegativeInfinity(smallestMassDiff))
                    {
                        double smallestOppositeMassDiff = -1 * (notches.Min(p => p.allowedInterval.Minimum) - scan.PrecursorMass);
                        highestMassPeptideToLookFor = scan.PrecursorMass + smallestOppositeMassDiff;
                    }

                    // first-pass scoring
                    IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor);

                    // done with indexed scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        int maxInitialScore = idsOfPeptidesPossiblyObserved.Max(id => scoringTable[id]);
                        var possiblyValidIds = idsOfPeptidesPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore);

                        foreach (var id in possiblyValidIds)
                        {
                            var peptide = peptideIndex[id];

                            var productMasses = peptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                            Array.Sort(productMasses);

                            double thePrecursorMass = scan.PrecursorMass;
                            double score = CalculatePeptideScore(scan.TheScan, CommonParameters.ProductMassTolerance, productMasses, thePrecursorMass, dissociationTypes, addCompIons);
                            int notch = massDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMassIncludingFixedMods);

                            if (score > CommonParameters.ScoreCutoff)
                            {
                                if (globalPsms[i] == null)
                                    globalPsms[i] = new Psm(peptide, notch, score, i, scan, CommonParameters.ExcelCompatible);
                                else
                                    globalPsms[i].AddOrReplace(peptide, score, notch, CommonParameters.ReportAllAmbiguity);
                            }
                        }
                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / listOfSortedms2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing modern search... " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));
                    }
                }
            });

            return new MetaMorpheusEngineResults(this);
        }

        protected List<int> IdentifyMostCommonBinsAll(List<int>[] fragmentIndex)
        {
            Tuple<int, int>[] mostCommonBins = new Tuple<int, int>[500];
            for (int i = 0; i < mostCommonBins.Length; i++)
                mostCommonBins[i] = new Tuple<int, int>(0, 0);
            for (int i = 0; i < fragmentIndex.Length; i++)
            {
                if (fragmentIndex[i] != null)
                {
                    if (mostCommonBins[499].Item2 < fragmentIndex[i].Count)
                    {
                        mostCommonBins[499] = new Tuple<int, int>(i, fragmentIndex[i].Count);
                        mostCommonBins = mostCommonBins.OrderByDescending(s => s.Item2).ToArray();
                    }
                }
            }
            return mostCommonBins.ToList().Select(s => s.Item1).ToList();
        }

        protected List<int> GetBinsToSearch(List<IMzPeak> peaks, double largestIntensity, double precursorMass)
        {
            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();
            foreach (IMzPeak peak in peaks)
            {
                if (CommonParameters.MinRatio == null || (peak.Intensity / largestIntensity) >= CommonParameters.MinRatio)
                {
                    // assume charge state 1 to calculate mz tolerance
                    int obsFragmentFloorMz = (int)Math.Floor((CommonParameters.ProductMassTolerance.GetMinimumValue(peak.Mz)) * fragmentBinsPerDalton);
                    if (obsFragmentFloorMz < obsPreviousFragmentCeilingMz)
                        obsFragmentFloorMz = obsPreviousFragmentCeilingMz;
                    int obsFragmentCeilingMz = (int)Math.Ceiling((CommonParameters.ProductMassTolerance.GetMaximumValue(peak.Mz)) * fragmentBinsPerDalton);
                    obsPreviousFragmentCeilingMz = obsFragmentCeilingMz + 1;
                    for (int fragmentBin = obsFragmentFloorMz; fragmentBin <= obsFragmentCeilingMz; fragmentBin++)
                        if (fragmentIndex[fragmentBin] != null)
                            binsToSearch.Add(fragmentBin);

                    if (addCompIons)
                    {
                        //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut and just straight up adding the bins assuming that they're z=1
                        for (int j = 0; j < dissociationTypes.Count; j++)
                        {
                            int compPrecursor = (int)((precursorMass + complementaryIonConversionDictionary[dissociationTypes[j]] + Constants.protonMass) * fragmentBinsPerDalton);
                            int compFragmentFloorMz = compPrecursor - obsFragmentCeilingMz;
                            int compFragmentCeilingMz = compPrecursor - obsFragmentFloorMz;
                            for (int fragmentBin = compFragmentFloorMz; fragmentBin <= compFragmentCeilingMz; fragmentBin++)
                                if (fragmentIndex[fragmentBin] != null)
                                    binsToSearch.Add(fragmentBin);
                        }
                    }
                }
            }
            return binsToSearch;
        }

        protected List<int> MostCommonBinsFound(List<int> binsToSearch, List<int> mostCommonBins, int intScoreCutoff, bool addCompIons)
        {
            int numChecksSkipped = 1;
            List<int> commonBinsFound = new List<int>();
            for (int commonBin = 0; commonBin < mostCommonBins.Count; commonBin++)
            {
                if (numChecksSkipped == intScoreCutoff)
                    break;
                if (binsToSearch.Contains(mostCommonBins[commonBin]))
                {
                    numChecksSkipped++;
                    binsToSearch.Remove(mostCommonBins[commonBin]);
                    commonBinsFound.Add(mostCommonBins[commonBin]);
                    if (addCompIons)
                        commonBin--;
                }
            }
            return commonBinsFound;
        }

        #endregion Protected Methods

        #region Private Methods

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
                    break;
                if (peptideIndex[peptideIdsInThisBin[m]].MonoisotopicMassIncludingFixedMods < peptideMassToLookFor)
                    l = m + 1;
                else
                    r = m - 1;
            }
            if (m > 0)
                m--;
            return m;
        }

        private void IndexedScoringForFrequentFragments(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = fragmentIndex[binsToSearch[i]];

                //get index for minimum monoisotopic allowed
                int m = Double.IsInfinity(lowestMassPeptideToLookFor) ? 0 : BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, lowestMassPeptideToLookFor);

                // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                if (!Double.IsInfinity(highestMassPeptideToLookFor))
                {
                    for (int h = m; h < peptideIdsInThisBin.Count; h++)
                    {
                        int id = peptideIdsInThisBin[h];
                        scoringTable[id]++;

                        if (peptideIndex[id].MonoisotopicMassIncludingFixedMods > highestMassPeptideToLookFor)
                            break;
                    }
                }
                else
                {
                    for (int h = m; h < peptideIdsInThisBin.Count; h++)
                    {
                        int id = peptideIdsInThisBin[h];
                        scoringTable[id]++;
                    }
                }
            }
        }

        private void IndexedScoring(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = fragmentIndex[binsToSearch[i]];

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
                        var nextPep = peptideIndex[nextId];
                        if (nextPep.MonoisotopicMassIncludingFixedMods < highestMassPeptideToLookFor)
                            highestPeptideMassIndex = j;
                        else
                            break;
                    }
                }

                // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                {
                    int id = peptideIdsInThisBin[j];
                    scoringTable[id]++;

                    // add possible search results to the hashset of id's
                    if (scoringTable[id] == byteScoreCutoff)
                        if (massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMassIncludingFixedMods) >= 0)
                            idsOfPeptidesPossiblyObserved.Add(id);
                }
            }
        }

        #endregion Private Methods

    }
}