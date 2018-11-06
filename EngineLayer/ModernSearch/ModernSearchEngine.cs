using Chemistry;
using MassSpectrometry;
using Proteomics.Fragmentation;
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
        protected readonly Ms2ScanWithSpecificMass[] ListOfSortedMs2Scans;
        protected readonly List<PeptideWithSetModifications> PeptideIndex;
        protected readonly int CurrentPartition;
        protected readonly MassDiffAcceptor MassDiffAcceptor;
        protected readonly DissociationType DissociationType;
        protected readonly double MaxMassThatFragmentIonScoreIsDoubled;

        public ModernSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex, List<int>[] fragmentIndex, int currentPartition, CommonParameters commonParameters, MassDiffAcceptor massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled, List<string> nestedIds) : base(commonParameters, nestedIds)
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
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing modern search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = (byte)commonParameters.ScoreCutoff;
            if (commonParameters.CalculateEValue)
                byteScoreCutoff = 1;

            Parallel.ForEach(Partitioner.Create(0, ListOfSortedMs2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (range, loopState) =>
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
                    Ms2ScanWithSpecificMass scan = ListOfSortedMs2Scans[i];

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(scan);

                    // get allowed theoretical masses from the known experimental mass
                    // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)
                    // this is just PRELIMINARY precursor-mass filtering
                    // additional checks are made later to ensure that the theoretical precursor mass is acceptable
                    IEnumerable<AllowedIntervalWithNotch> notches = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(scan.PrecursorMass);

                    double lowestMassPeptideToLookFor = notches.Min(p => p.AllowedInterval.Minimum);
                    double highestMassPeptideToLookFor = notches.Max(p => p.AllowedInterval.Maximum);

                    // first-pass scoring
                    IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor, PeptideIndex, MassDiffAcceptor, MaxMassThatFragmentIonScoreIsDoubled);

                    // done with indexed scoring; refine scores and create PSMs
                    foreach (int id in idsOfPeptidesPossiblyObserved)
                    {
                        PeptideWithSetModifications peptide = PeptideIndex[id];

                        List<Product> peptideTheorProducts = peptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both).ToList();

                        List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan, peptideTheorProducts, commonParameters);
                        
                        double thisScore = CalculatePeptideScore(scan.TheScan, matchedIons, 0);
                        int notch = MassDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMass);

                        bool meetsScoreCutoff = thisScore >= commonParameters.ScoreCutoff;
                        bool scoreImprovement = PeptideSpectralMatches[i] == null || (thisScore - PeptideSpectralMatches[i].RunnerUpScore) > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

                        if (meetsScoreCutoff && scoreImprovement || commonParameters.CalculateEValue)
                        {
                            if (PeptideSpectralMatches[i] == null)
                            {
                                PeptideSpectralMatches[i] = new PeptideSpectralMatch(peptide, notch, thisScore, i, scan, commonParameters.DigestionParams, matchedIons);
                            }
                            else
                            {
                                PeptideSpectralMatches[i].AddOrReplace(peptide, thisScore, notch, commonParameters.ReportAllAmbiguity, matchedIons);
                            }

                            if (commonParameters.CalculateEValue)
                            {
                                PeptideSpectralMatches[i].AllScores.Add(thisScore);
                            }
                        }
                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing modern search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));
                    }
                }
            });

            // remove peptides below the score cutoff that were stored to calculate expectation values
            if (commonParameters.CalculateEValue)
            {
                for (int i = 0; i < PeptideSpectralMatches.Length; i++)
                {
                    if (PeptideSpectralMatches[i] != null && PeptideSpectralMatches[i].Score < commonParameters.ScoreCutoff)
                    {
                        PeptideSpectralMatches[i] = null;
                    }
                }
            }

            foreach (PeptideSpectralMatch psm in PeptideSpectralMatches.Where(p => p != null))
            {
                psm.ResolveAllAmbiguities();
            }

            return new MetaMorpheusEngineResults(this);
        }

        protected List<int> GetBinsToSearch(Ms2ScanWithSpecificMass scan)
        {
            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();
            
            foreach (var envelope in scan.ExperimentalFragments)
            {
                // assume charge state 1 to calculate mass tolerance
                double experimentalFragmentMass = envelope.monoisotopicMass;

                // get theoretical fragment bins within mass tolerance
                int obsFragmentFloorMass = (int)Math.Floor((commonParameters.ProductMassTolerance.GetMinimumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);
                int obsFragmentCeilingMass = (int)Math.Ceiling((commonParameters.ProductMassTolerance.GetMaximumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);

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
                if (commonParameters.AddCompIons)
                {
                    //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut and just straight up adding the bins assuming that they're z=1

                    if (complementaryIonConversionDictionary.TryGetValue(commonParameters.DissociationType, out double protonMassShift)) //TODO: this is broken for EThcD because that method needs two conversions
                    {
                        protonMassShift = ClassExtensions.ToMass(protonMassShift, 1);
                        int compFragmentFloorMass = (int)Math.Round(((scan.PrecursorMass + protonMassShift) * FragmentBinsPerDalton)) - obsFragmentCeilingMass;
                        int compFragmentCeilingMass = (int)Math.Round(((scan.PrecursorMass + protonMassShift) * FragmentBinsPerDalton)) - obsFragmentFloorMass;

                        // prevent index out of bounds errors
                        if (compFragmentCeilingMass >= FragmentIndex.Length)
                        {
                            compFragmentCeilingMass = FragmentIndex.Length - 1;

                            if (compFragmentFloorMass >= FragmentIndex.Length)
                                compFragmentFloorMass = FragmentIndex.Length - 1;
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
            return binsToSearch;
        }

        protected static int BinarySearchBinForPrecursorIndex(List<int> peptideIdsInThisBin, double peptideMassToLookFor, List<PeptideWithSetModifications> peptideIndex)
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
                if (peptideIndex[peptideIdsInThisBin[m]].MonoisotopicMass < peptideMassToLookFor)
                    l = m + 1;
                else
                    r = m - 1;
            }
            if (m > 0)
                m--;
            return m;
        }

        protected void IndexedScoring(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor,
            double highestMassPeptideToLookFor, List<PeptideWithSetModifications> peptideIndex, MassDiffAcceptor massDiffAcceptor, double maxMassThatFragmentIonScoreIsDoubled)
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
                    if (scoringTable[id] == byteScoreCutoff && massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMass) >= 0)
                        idsOfPeptidesPossiblyObserved.Add(id);
                }

                if (maxMassThatFragmentIonScoreIsDoubled > 0)
                {
                    for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                    {
                        if (j < maxMassThatFragmentIonScoreIsDoubled * FragmentBinsPerDalton)
                        {
                            int id = peptideIdsInThisBin[j];
                            scoringTable[id]++;

                            // add possible search results to the hashset of id's
                            if (scoringTable[id] == byteScoreCutoff && massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMass) >= 0)
                                idsOfPeptidesPossiblyObserved.Add(id);
                        }
                    }
                }
            }
        }
    }
}