using Chemistry;
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

            List<int> mostCommonBins = IdentifyMostCommonBinsAll(fragmentIndex);

            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUse }, range =>
            {
                byte[] scoringTable = new byte[peptideIndex.Count];
                HashSet<int> idsOfPeptidesPossiblyObserved = new HashSet<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = listOfSortedms2Scans[i];

                    // filter ms2 fragment peaks by intensity
                    int numFragmentsToUse = scan.NumPeaks;
                    if (CommonParameters.TopNpeaks != null)
                        numFragmentsToUse = (int)CommonParameters.TopNpeaks;

                    var peaks = scan.TheScan.MassSpectrum.FilterByNumberOfMostIntense(numFragmentsToUse).ToList();
                    double largestIntensity = scan.TheScan.MassSpectrum.YofPeakWithHighestY;

                    //get bins to add points to
                    List<int> allBinsToSearch = GetBinsToSearch(peaks, largestIntensity, scan.PrecursorMass);

                    //separate bins by common and uncommon fragments to improve search speed
                    List<int> commonBinsToSearch = MostCommonBinsFound(allBinsToSearch, mostCommonBins, byteScoreCutoff, addCompIons);

                    // get allowed precursor masses
                    AllowedIntervalWithNotch[] notches = massDiffAcceptor.GetAllowedPrecursorMassIntervals(scan.PrecursorMass).ToArray();
                    double[] lowestMassPeptideToLookForArray = new double[notches.Count()];
                    double[] highestMassPeptideToLookForArray = new double[notches.Count()];
                    for (int notch = 0; notch < notches.Length; notch++)
                    {
                        lowestMassPeptideToLookForArray[notch] = notches[notch].allowedInterval.Minimum;
                        highestMassPeptideToLookForArray[notch] = notches[notch].allowedInterval.Maximum;
                    }

                    //First pass scoring commonly viewed fragments but not searching for candidates
                    IndexedScoringCommon(commonBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookForArray, highestMassPeptideToLookForArray);
                    //Second pass scoring all other peaks and collecting theoretical peptides with high enough scores
                    IndexedScoringUncommon(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookForArray, highestMassPeptideToLookForArray);

                    // done with indexed scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        int maxInitialScore = idsOfPeptidesPossiblyObserved.Max(id => scoringTable[id]);

                        foreach (var id in idsOfPeptidesPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore))
                        {
                            var candidatePeptide = peptideIndex[id];
                            double[] fragmentMasses = candidatePeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp).Distinct().Where(p => !Double.IsNaN(p)).OrderBy(p => p).ToArray();
                            double peptideScore = CalculatePeptideScore(scan.TheScan, CommonParameters.ProductMassTolerance, fragmentMasses, scan.PrecursorMass, dissociationTypes, addCompIons);
                            int notch = massDiffAcceptor.Accepts(scan.PrecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);

                            if (globalPsms[i] == null)
                                globalPsms[i] = new Psm(candidatePeptide, notch, peptideScore, i, scan);
                            else
                                globalPsms[i].AddOrReplace(candidatePeptide, peptideScore, notch, CommonParameters.ReportAllAmbiguity);
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

        protected static List<int> IdentifyMostCommonBinsAll(List<int>[] fragmentIndex)
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
                if (peak.Mz > precursorMass)
                    break;
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

        protected static List<int> MostCommonBinsFound(List<int> binsToSearch, List<int> mostCommonBins, int intScoreCutoff, bool addCompIons)
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

        private static int FindLowestMassIndexAllowed(List<CompactPeptide> peptideIndex, List<int> peptideIdsInThisBin, double lowestMassPeptideToLookFor)
        {
            int m = 0;
            int l = 0;
            int r = peptideIdsInThisBin.Count - 1;

            // binary search in the fragment bin for lowest acceptable precursor mass
            while (l <= r)
            {
                m = l + ((r - l) / 2);

                if (r - l < 2)
                    break;
                if (peptideIndex[peptideIdsInThisBin[m]].MonoisotopicMassIncludingFixedMods < lowestMassPeptideToLookFor)
                    l = m + 1;
                else
                    r = m - 1;
            }
            if (m > 0)
                m--;
            return m;
        }

        private void IndexedScoringCommon(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, HashSet<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double[] lowestMassPeptideToLookForArray, double[] highestMassPeptideToLookForArray)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = fragmentIndex[binsToSearch[i]];
                for (int notch = 0; notch < lowestMassPeptideToLookForArray.Length; notch++)
                {
                    double lowestMassPeptideToLookFor = lowestMassPeptideToLookForArray[notch];
                    double highestMassPeptideToLookFor = highestMassPeptideToLookForArray[notch];
                    //get index for minimum monoisotopic allowed
                    int m = Double.IsInfinity(lowestMassPeptideToLookFor) ? 0 : FindLowestMassIndexAllowed(peptideIndex, peptideIdsInThisBin, lowestMassPeptideToLookFor);

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
        }

        private void IndexedScoringUncommon(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, HashSet<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double[] lowestMassPeptideToLookForArray, double[] highestMassPeptideToLookForArray)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = fragmentIndex[binsToSearch[i]];
                for (int notch = 0; notch < lowestMassPeptideToLookForArray.Length; notch++)
                {
                    double lowestMassPeptideToLookFor = lowestMassPeptideToLookForArray[notch];
                    double highestMassPeptideToLookFor = highestMassPeptideToLookForArray[notch];
                    //get index for minimum monoisotopic allowed
                    int m = Double.IsInfinity(lowestMassPeptideToLookFor) ? 0 : FindLowestMassIndexAllowed(peptideIndex, peptideIdsInThisBin, lowestMassPeptideToLookFor);

                    // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                    if (!Double.IsInfinity(highestMassPeptideToLookFor))
                    {
                        for (int h = m; h < peptideIdsInThisBin.Count; h++)
                        {
                            int id = peptideIdsInThisBin[h];
                            scoringTable[id]++;

                            // add possible search results to the hashset of id's, should be only difference between IndexedScoringUncommon and IndexedScoringCommon
                            if (scoringTable[id] == byteScoreCutoff)
                                idsOfPeptidesPossiblyObserved.Add(id);

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

                            // add possible search results to the hashset of id's, should be only difference between IndexedScoringUncommon and IndexedScoringCommon
                            if (scoringTable[id] == byteScoreCutoff)
                                idsOfPeptidesPossiblyObserved.Add(id);
                        }
                    }
                }
            }
        }

        #endregion Private Methods

    }
}