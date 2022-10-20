using Chemistry;
using MassSpectrometry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
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
        protected readonly Ms2ScanWithSpecificMass[] ListOfSortedMs2Scans;
        protected readonly List<PeptideWithSetModifications> PeptideIndex;
        protected readonly int CurrentPartition;
        protected readonly MassDiffAcceptor MassDiffAcceptor;
        protected readonly DissociationType DissociationType;
        protected readonly double MaxMassThatFragmentIonScoreIsDoubled;

        public ModernSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, MassDiffAcceptor massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled,
            List<string> nestedIds) : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            PeptideSpectralMatches = globalPsms;
            ListOfSortedMs2Scans = listOfSortedms2Scans;
            PeptideIndex = peptideIndex;
            FragmentIndex = fragmentIndex;
            CurrentPartition = currentPartition + 1;
            MassDiffAcceptor = massDiffAcceptor;
            DissociationType = commonParameters.DissociationType;
            MaxMassThatFragmentIonScoreIsDoubled = maximumMassThatFragmentIonScoreIsDoubled;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing modern search... " +
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
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing modern search... " +
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

        /// <summary>
        /// This is a first-pass scoring method which is supposed to be *very fast* but may sometimes miscount the number of truly matched fragments. The number
        /// of matched fragments should always be overestimated, never underestimated.
        /// </summary>
        protected void IndexScoreScan(Ms2ScanWithSpecificMass scan, byte[] scoringTable, byte byteScoreCutoff, List<int> peptidesPossiblyObserved, DissociationType dissociationType)
        {
            // get allowed theoretical masses from the known experimental mass
            // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)	
            // this is just PRELIMINARY precursor-mass filtering	
            // additional checks are made later to ensure that the theoretical precursor mass is acceptable
            List<AllowedIntervalWithNotch> notches = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(scan.PrecursorMass).ToList();
            double lowestMassPeptideToLookFor = notches.Min(p => p.AllowedInterval.Minimum);
            double highestMassPeptideToLookFor = notches.Max(p => p.AllowedInterval.Maximum);

            // clear the scoring table to score the new scan (conserves memory compared to allocating a new array)
            Array.Clear(scoringTable, 0, scoringTable.Length);
            peptidesPossiblyObserved.Clear();

            if (dissociationType == DissociationType.LowCID)
            {
                double[] masses = scan.TheScan.MassSpectrum.XArray;
                double[] intensities = scan.TheScan.MassSpectrum.YArray;

                for (int i = 0; i < masses.Length; i++)
                {
                    //convert to an int since we're in discrete 1.0005...
                    int fragmentBin = (int)(Math.Round(masses[i].ToMass(1) / 1.0005079) * 1.0005079 * FragmentBinsPerDalton);

                    List<int> bin = FragmentIndex[fragmentBin];

                    if (bin != null)
                    {
                        // filter bin by peptide mass
                        var (start, end) = GetFirstAndLastIndexesInBinToIncrement(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, scan.PrecursorMass);

                        // add +1 to each peptide score
                        IncrementPeptideScoresInBin(start, end, bin, scoringTable, scan, byteScoreCutoff, peptidesPossiblyObserved, CommonParameters.DissociationType);
                    }

                    // add complementary ions
                    if (CommonParameters.AddCompIons)
                    {
                        if (complementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift)) //TODO: this is broken for EThcD because that method needs two conversions
                        {
                            protonMassShift = ClassExtensions.ToMass(protonMassShift, 1);
                            fragmentBin = (int)Math.Round((scan.PrecursorMass + protonMassShift - masses[i]) / 1.0005079);

                            bin = FragmentIndex[fragmentBin];

                            if (bin != null)
                            {
                                // filter bin by peptide mass
                                var (start, end) = GetFirstAndLastIndexesInBinToIncrement(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, scan.PrecursorMass);

                                // add +1 to each peptide score
                                IncrementPeptideScoresInBin(start, end, bin, scoringTable, scan, byteScoreCutoff, peptidesPossiblyObserved, CommonParameters.DissociationType);
                            }
                        }
                        else
                        {
                            throw new NotImplementedException();
                        }
                    }
                }
            }
            else
            {
                for (int i = 0; i < scan.ExperimentalFragments.Length; i++)
                {
                    double mass = scan.ExperimentalFragments[i].MonoisotopicMass;

                    // get theoretical fragment bins within mass tolerance
                    int obsFragmentFloorMass = Math.Max(0,
                        (int)Math.Floor((CommonParameters.ProductMassTolerance.GetMinimumValue(mass)) * FragmentBinsPerDalton));
                    int obsFragmentCeilingMass = Math.Min(FragmentIndex.Length - 1,
                        (int)Math.Ceiling((CommonParameters.ProductMassTolerance.GetMaximumValue(mass)) * FragmentBinsPerDalton));

                    for (int b = obsFragmentFloorMass; b <= obsFragmentCeilingMass; b++)
                    {
                        List<int> bin = FragmentIndex[b];

                        if (bin == null)
                        {
                            continue;
                        }

                        // filter bin by peptide mass
                        var (start, end) = GetFirstAndLastIndexesInBinToIncrement(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, scan.PrecursorMass);

                        // add +1 to each peptide score
                        IncrementPeptideScoresInBin(start, end, bin, scoringTable, scan, byteScoreCutoff, peptidesPossiblyObserved, CommonParameters.DissociationType);
                    }

                    if (CommonParameters.AddCompIons)
                    {
                        //TODO: this is broken for EThcD because that method needs two conversions
                        if (complementaryIonConversionDictionary.TryGetValue(CommonParameters.DissociationType, out double protonMassShift))
                        {
                            protonMassShift = ClassExtensions.ToMass(protonMassShift, 1);

                            int compFragmentFloorMass = Math.Max(0,
                                (int)Math.Round(((scan.PrecursorMass + protonMassShift) * FragmentBinsPerDalton)) - obsFragmentCeilingMass);
                            int compFragmentCeilingMass = Math.Min(FragmentIndex.Length - 1,
                                (int)Math.Round(((scan.PrecursorMass + protonMassShift) * FragmentBinsPerDalton)) - obsFragmentFloorMass);

                            for (int b = compFragmentFloorMass; b <= compFragmentCeilingMass; b++)
                            {
                                List<int> bin = FragmentIndex[b];

                                if (bin == null)
                                {
                                    continue;
                                }

                                // filter bin by peptide mass
                                var (start, end) = GetFirstAndLastIndexesInBinToIncrement(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, scan.PrecursorMass);

                                // add +1 to each peptide score
                                IncrementPeptideScoresInBin(start, end, bin, scoringTable, scan, byteScoreCutoff, peptidesPossiblyObserved, CommonParameters.DissociationType);
                            }
                        }
                        else
                        {
                            throw new NotImplementedException();
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Finds the first and last bin-indexes of the peptides to add a +1 score to, based on the precursor mass and precursor mass tolerance.
        /// </summary>
        protected (int start, int end) GetFirstAndLastIndexesInBinToIncrement(double lowestPeptideMassToLookFor, double highestPeptideMassToLookFor, List<int> bin, double precursorMass)
        {
            int start = 0;
            int end = bin.Count - 1;

            if (!double.IsPositiveInfinity(highestPeptideMassToLookFor))
            {
                end = BinarySearchBinForPrecursorIndex(bin, highestPeptideMassToLookFor, PeptideIndex);
            }

            if (!double.IsNegativeInfinity(lowestPeptideMassToLookFor))
            {
                start = BinarySearchBinForPrecursorIndex(bin, lowestPeptideMassToLookFor, PeptideIndex);
            }

            return (start, end);
        }

        /// <summary>
        /// Returns the bin-index of the first peptide with a mass less than or equal to the specified mass. Returns 0 if there 
        /// are no peptides with masses smaller than the specified mass.
        /// </summary>
        protected static int BinarySearchBinForPrecursorIndex(List<int> bin, double peptideMassToLookFor, List<PeptideWithSetModifications> peptideIndex)
        {
            int m = 0;
            int l = 0;
            int r = bin.Count - 1;

            // binary search in the fragment bin for precursor mass
            while (l <= r)
            {
                m = (l + r) / 2;

                if (r - l < 2)
                {
                    for (; r >= 0; r--)
                    {
                        if (peptideIndex[bin[r]].MonoisotopicMass <= peptideMassToLookFor)
                        {
                            return r;
                        }
                    }
                }

                if (peptideMassToLookFor > peptideIndex[bin[m]].MonoisotopicMass)
                {
                    l = m + 1;
                }
                else
                {
                    r = m - 1;
                }
            }

            // this happens only if there are no peptides in the bin less than the looked-for mass
            return 0;
        }

        /// <summary>
        /// Adds a +1 score to all the peptides in the fragment mass bin that meet the precursor mass tolerance.
        /// </summary>
        protected void IncrementPeptideScoresInBin(int start, int end, List<int> bin, byte[] scoringTable, Ms2ScanWithSpecificMass scan, byte byteScoreCutoff,
            List<int> peptidesPossiblyObserved, DissociationType dissociationType)
        {
            if (dissociationType == DissociationType.LowCID)
            {
                // add score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                for (int j = start; j <= end; j++)
                {
                    int peptideId = bin[j];

                    // TODO: revisit this at some point.. this will be mega slow. every peptide with 
                    // at least one fragment in the scan will get fine-scored.

                    // add possible search results to the list of id's (only once)
                    if (scoringTable[peptideId] == 0 && MassDiffAcceptor.Accepts(scan.PrecursorMass, PeptideIndex[peptideId].MonoisotopicMass) >= 0)
                    {
                        peptidesPossiblyObserved.Add(peptideId);
                    }

                    // mark the peptide as potentially observed so it doesn't get added more than once
                    scoringTable[peptideId] = 1;
                }
            }
            else
            {
                // add +1 to each peptide score
                for (int p = start; p <= end; p++)
                {
                    int peptideId = bin[p];
                    byte score = ++scoringTable[peptideId];

                    // if the peptide has met the score cutoff, add it to the list of peptides 
                    // possibly observed so it can be re-scored with the "fine scoring" algorithm
                    if (score == byteScoreCutoff && MassDiffAcceptor.Accepts(scan.PrecursorMass, PeptideIndex[peptideId].MonoisotopicMass) >= 0)
                    {
                        peptidesPossiblyObserved.Add(peptideId);
                    }
                }
            }
        }

        /// <summary>
        /// This is a second-pass scoring method which is costly (in terms of computational time and RAM) but calculates the "normal" MetaMorpheus score instead
        /// of the approximation computed by the IndexScoreScan method.
        /// </summary>
        protected PeptideSpectralMatch FineScorePeptide(int id, Ms2ScanWithSpecificMass scan, int scanIndex, List<Product> peptideTheorProducts)
        {
            PeptideWithSetModifications peptide = PeptideIndex[id];

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

        protected void FineScorePeptides(List<int> peptideIds, Ms2ScanWithSpecificMass scan, int scanIndex, byte[] scoringTable, 
            DissociationType dissociationType, List<Product> peptideTheorProducts)
        {
            // this method re-scores the top-scoring peptides until no peptide in the rough-scored list can out-score
            // the best-scoring peptide. this guarantees that peptides will be scored accurately, according to metamorpheus score,
            // while high speed is maintained. however, this means that the delta score is only an approximation. since PEP 
            // analysis trains on delta score, the modern search output is not guaranteed to be the same as classic search 
            // output, though it will be very close.
            byte bestScore = 0;

            foreach (int id in peptideIds.OrderByDescending(p => scoringTable[p]))
            {
                if (scoringTable[id] < bestScore && dissociationType != DissociationType.LowCID)
                {
                    FineScorePeptide(id, scan, scanIndex, peptideTheorProducts);

                    break;
                }

                PeptideSpectralMatch psm = FineScorePeptide(id, scan, scanIndex, peptideTheorProducts);

                if (psm != null && psm.Score > bestScore)
                {
                    bestScore = (byte)Math.Floor(psm.Score);
                }
            }
        }

        /// <summary>
        /// Deprecated.
        /// </summary>
        protected void IndexedScoring(List<int>[] FragmentIndex, List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor,
            double highestMassPeptideToLookFor, List<PeptideWithSetModifications> peptideIndex, MassDiffAcceptor massDiffAcceptor, double maxMassThatFragmentIonScoreIsDoubled, DissociationType dissociationType)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = FragmentIndex[binsToSearch[i]];

                //get index for minimum monoisotopic allowed
                int lowestPeptideMassIndex = Double.IsInfinity(lowestMassPeptideToLookFor) ? 0 : BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, lowestMassPeptideToLookFor, peptideIndex);

                // get index for highest mass allowed
                int highestPeptideMassIndex = peptideIdsInThisBin.Count - 1;

                if (!Double.IsInfinity(highestMassPeptideToLookFor))
                {
                    highestPeptideMassIndex = BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, highestMassPeptideToLookFor, peptideIndex);

                    for (int j = highestPeptideMassIndex; j < peptideIdsInThisBin.Count; j++)
                    {
                        int nextId = peptideIdsInThisBin[j];
                        var nextPep = peptideIndex[nextId];
                        if (nextPep.MonoisotopicMass < highestMassPeptideToLookFor)
                        {
                            highestPeptideMassIndex = j;
                        }
                        else
                        {
                            break;
                        }
                    }
                }

                if (dissociationType == DissociationType.LowCID)
                {
                    // add intensity for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                    for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                    {
                        int id = peptideIdsInThisBin[j];

                        // add possible search results to the hashset of id's (only once)
                        if (scoringTable[id] == 0 && massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMass) >= 0)
                        {
                            idsOfPeptidesPossiblyObserved.Add(id);
                        }

                        // mark the peptide as potentially observed so it doesn't get added more than once
                        scoringTable[id] = 1;
                    }
                }
                else
                {
                    // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                    for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                    {
                        int id = peptideIdsInThisBin[j];
                        scoringTable[id]++;

                        // add possible search results to the hashset of id's
                        if (scoringTable[id] == byteScoreCutoff && massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMass) >= 0)
                        {
                            idsOfPeptidesPossiblyObserved.Add(id);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Deprecated.
        /// </summary>
        protected List<int> GetBinsToSearch(Ms2ScanWithSpecificMass scan, List<int>[] FragmentIndex, DissociationType dissociationType)
        {
            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();

            if (dissociationType == DissociationType.LowCID)
            {
                double[] masses = scan.TheScan.MassSpectrum.XArray;
                double[] intensities = scan.TheScan.MassSpectrum.YArray;

                for (int i = 0; i < masses.Length; i++)
                {
                    //convert to an int since we're in discrete 1.0005...
                    int fragmentBin = (int)(Math.Round(masses[i].ToMass(1) / 1.0005079) * 1.0005079 * FragmentBinsPerDalton);

                    if (FragmentIndex[fragmentBin] != null)
                    {
                        binsToSearch.Add(fragmentBin);
                    }

                    // add complementary ions
                    if (CommonParameters.AddCompIons)
                    {
                        if (complementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift)) //TODO: this is broken for EThcD because that method needs two conversions
                        {
                            protonMassShift = ClassExtensions.ToMass(protonMassShift, 1);
                            fragmentBin = (int)Math.Round((scan.PrecursorMass + protonMassShift - masses[i]) / 1.0005079);

                            if (FragmentIndex[fragmentBin] != null)
                            {
                                binsToSearch.Add(fragmentBin);
                            }
                        }
                        else
                        {
                            throw new NotImplementedException();
                        }
                    }
                }
            }
            else
            {
                foreach (var envelope in scan.ExperimentalFragments)
                {
                    // assume charge state 1 to calculate mass tolerance
                    double experimentalFragmentMass = envelope.MonoisotopicMass;

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
                    if (CommonParameters.AddCompIons)
                    {
                        //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut and just straight up adding the bins assuming that they're z=1

                        if (complementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift)) //TODO: this is broken for EThcD because that method needs two conversions
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
    }
}