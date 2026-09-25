using Chemistry;
using MassSpectrometry;
using Omics;
using Omics.Fragmentation;
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
        protected Indexing.FragmentIndex FragmentIndex { get; private set; }
        protected readonly SpectralMatch[] PeptideSpectralMatches;
        protected readonly Ms2ScanWithSpecificMass[] ListOfSortedMs2Scans;
        protected readonly List<IBioPolymerWithSetMods> PeptideIndex;
        protected readonly int CurrentPartition;
        protected readonly MassDiffAcceptor MassDiffAcceptor;
        protected readonly DissociationType DissociationType;
        protected readonly double MaxMassThatFragmentIonScoreIsDoubled;

        public ModernSearchEngine(SpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, IEnumerable<IBioPolymerWithSetMods> peptideIndex,
            Indexing.FragmentIndex fragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, MassDiffAcceptor massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled,
            List<string> nestedIds) : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            PeptideSpectralMatches = globalPsms;
            ListOfSortedMs2Scans = listOfSortedms2Scans;
            PeptideIndex = peptideIndex as List<IBioPolymerWithSetMods> ?? peptideIndex?.ToList();
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

            foreach (SpectralMatch psm in PeptideSpectralMatches.Where(p => p != null))
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
            List<AllowedIntervalWithNotch> notches = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(scan.GetPrecursorMassForSearch(CommonParameters)).ToList();
            double lowestMassPeptideToLookFor = notches.Min(p => p.Minimum);
            double highestMassPeptideToLookFor = notches.Max(p => p.Maximum);

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

                    ReadOnlySpan<int> bin = FragmentIndex[fragmentBin];

                    if (!bin.IsEmpty)
                    {
                        // filter bin by peptide mass
                        var (start, end) = GetFirstAndLastIndexesInBinToIncrement(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, scan.PrecursorMass);

                        // add +1 to each peptide score
                        IncrementPeptideScoresInBin(start, end, bin, scoringTable, scan, byteScoreCutoff, peptidesPossiblyObserved, CommonParameters.DissociationType);
                    }

                    // add complementary ions
                    if (CommonParameters.AddCompIons)
                    {
                        if (complementaryIonConversionDictionary.ContainsKey(CommonParameters.DissociationType))
                        {
                            foreach (double massshift in complementaryIonConversionDictionary[CommonParameters.DissociationType])
                            {
                                double protonMassShift = massshift.ToMass(1);
                                fragmentBin = (int)Math.Round((scan.PrecursorMass + protonMassShift - masses[i]) / 1.0005079);

                                bin = FragmentIndex[fragmentBin];

                                if (!bin.IsEmpty)
                                {
                                    // filter bin by peptide mass
                                    var (start, end) = GetFirstAndLastIndexesInBinToIncrement(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, scan.PrecursorMass);

                                    // add +1 to each peptide score
                                    IncrementPeptideScoresInBin(start, end, bin, scoringTable, scan, byteScoreCutoff, peptidesPossiblyObserved, CommonParameters.DissociationType);
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
                        ReadOnlySpan<int> bin = FragmentIndex[b];

                        if (bin.IsEmpty)
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
                        if (complementaryIonConversionDictionary.ContainsKey(CommonParameters.DissociationType))
                        {
                            foreach (double massShift in complementaryIonConversionDictionary[CommonParameters.DissociationType])
                            {
                                double protonMassShift = massShift.ToMass(1);

                                int compFragmentFloorMass = Math.Max(0,
                                    (int)Math.Round(((scan.PrecursorMass + protonMassShift) * FragmentBinsPerDalton)) - obsFragmentCeilingMass);
                                int compFragmentCeilingMass = Math.Min(FragmentIndex.Length - 1,
                                    (int)Math.Round(((scan.PrecursorMass + protonMassShift) * FragmentBinsPerDalton)) - obsFragmentFloorMass);

                                for (int b = compFragmentFloorMass; b <= compFragmentCeilingMass; b++)
                                {
                                    ReadOnlySpan<int> bin = FragmentIndex[b];

                                    if (bin.IsEmpty)
                                    {
                                        continue;
                                    }

                                    // filter bin by peptide mass
                                    var (start, end) = GetFirstAndLastIndexesInBinToIncrement(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, scan.PrecursorMass);

                                    // add +1 to each peptide score
                                    IncrementPeptideScoresInBin(start, end, bin, scoringTable, scan, byteScoreCutoff, peptidesPossiblyObserved, CommonParameters.DissociationType);
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
        }

        /// <summary>
        /// Finds the first and last bin-indexes of the peptides to add a +1 score to, based on the precursor mass and precursor mass tolerance.
        /// </summary>
        protected (int start, int end) GetFirstAndLastIndexesInBinToIncrement(double lowestPeptideMassToLookFor, double highestPeptideMassToLookFor, ReadOnlySpan<int> bin, double precursorMass)
        {
            int start = 0;
            int end = bin.Length - 1;

            if (!double.IsPositiveInfinity(highestPeptideMassToLookFor))
            {
                end = BinarySearchBinForPrecursorIndex(bin, highestPeptideMassToLookFor, PeptideIndex);

                // every peptide in this bin is heavier than the window allows
                if (end < 0)
                {
                    return (0, -1);
                }
            }

            if (!double.IsNegativeInfinity(lowestPeptideMassToLookFor))
            {
                start = BinarySearchBinForFirstAtOrAbove(bin, lowestPeptideMassToLookFor, PeptideIndex);
            }

            return (start, end);
        }

        /// <summary>
        /// A peptide's mass as the bin searches must see it: an undefined mass reads as lower than everything.
        ///
        /// Both searches need the bin's masses to be monotone in their predicate, and every comparison against
        /// NaN is false, so a raw NaN reads as "too heavy" to one search and "too light" to the other no matter
        /// where it sits. Sorting the peptide index by mass puts NaNs at the front of every bin they occupy, so
        /// reading them as negative infinity makes both predicates monotone: "at or below the upper bound" is
        /// true for a prefix, and "at or above the lower bound" is false for a prefix.
        ///
        /// They are placed below the window rather than removed from the index. A tolerance-based acceptor can
        /// never match one, but OpenSearchMode.Accepts returns 0 for anything, and an open search leaves both
        /// bounds infinite so neither search runs and the whole bin is scored. Peptides with an unknown residue
        /// are findable that way, and dropping them from the index would silently stop open and glyco searches
        /// reporting them.
        /// </summary>
        private static double MassForBinSearch(double monoisotopicMass)
        {
            return double.IsNaN(monoisotopicMass) ? double.NegativeInfinity : monoisotopicMass;
        }

        /// <summary>
        /// The index of the first peptide in the bin with a mass at or above the specified mass, or the bin's
        /// length if there is none, which makes the returned range empty.
        ///
        /// The start of the window needs a lower bound; the end needs an upper bound. Both used to be taken
        /// from the upper-bound search below, which is the wrong question to ask for the start: it answers
        /// with the LAST entry at or below the bound, so a run of equal masses sitting on the window's lower
        /// edge was clipped down to its final member. Equal masses are not a corner case here - the same
        /// peptide sequence occurs in several proteins, and a reversed decoy carries its target's mass - and
        /// notch intervals are built off a peptide mass, so the lower edge lands exactly on such a run
        /// routinely. Measured on the mouse proteome: five copies of QQAQNIEKMSK share one bin set for scan
        /// 27831 and four of them scored nothing.
        /// </summary>
        protected static int BinarySearchBinForFirstAtOrAbove<T>(ReadOnlySpan<int> bin, double peptideMassToLookFor, List<T> peptideIndex) where T : IBioPolymerWithSetMods
        {
            int low = 0;
            int high = bin.Length - 1;
            int result = bin.Length;

            while (low <= high)
            {
                int mid = low + ((high - low) / 2);

                if (MassForBinSearch(peptideIndex[bin[mid]].MonoisotopicMass) >= peptideMassToLookFor)
                {
                    result = mid;
                    high = mid - 1;
                }
                else
                {
                    low = mid + 1;
                }
            }

            return result;
        }

        /// <summary>
        /// The index of the last peptide in the bin with a mass at or below the specified mass, or -1 if
        /// there is none.
        /// </summary>
        protected static int BinarySearchBinForPrecursorIndex<T>(ReadOnlySpan<int> bin, double peptideMassToLookFor, List<T> peptideIndex) where T : IBioPolymerWithSetMods
        {
            // Plain upper-bound search: find the last index whose mass is <= the target.
            //
            // The previous implementation narrowed the window and then linear-scanned downwards from the
            // window's r, which gave a path-dependent answer. On an exact tie between the target and a stored
            // mass it took r = m - 1, discarding the very index it was looking for, and returned the first
            // index of an equal-mass run rather than the last. Because the path depends on bin.Count, the same
            // target could return different answers for bins of different length. Since this value is the
            // inclusive end of the range that receives coarse score increments, too small an index means
            // candidates are silently never scored: probed against a linear-scan ground truth over 104,472
            // (bin length, target) pairs it was wrong 1,390 times, once returning 0 where 719 was correct.
            int low = 0;
            int high = bin.Length - 1;
            int result = -1;

            while (low <= high)
            {
                int mid = low + ((high - low) / 2);

                if (MassForBinSearch(peptideIndex[bin[mid]].MonoisotopicMass) <= peptideMassToLookFor)
                {
                    result = mid;
                    low = mid + 1;
                }
                else
                {
                    high = mid - 1;
                }
            }

            // -1 rather than 0 when nothing is at or below the looked-for mass, so callers can tell that apart
            // from index 0 being the answer. Conflating them made a bin whose peptides are all too heavy
            // still score its first entry; that fired on a third of window lookups on the mouse proteome.
            return result;
        }

        /// <summary>
        /// Adds a +1 score to all the peptides in the fragment mass bin that meet the precursor mass tolerance.
        /// </summary>
        protected void IncrementPeptideScoresInBin(int start, int end, ReadOnlySpan<int> bin, byte[] scoringTable, Ms2ScanWithSpecificMass scan, byte byteScoreCutoff,
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
                    if (scoringTable[peptideId] == 0 && MassDiffAcceptor.Accepts(scan.GetPrecursorMassForSearch(CommonParameters), PeptideIndex[peptideId].MonoisotopicMass) >= 0)
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
                    if (score == byteScoreCutoff && MassDiffAcceptor.Accepts(scan.GetPrecursorMassForSearch(CommonParameters), PeptideIndex[peptideId].MonoisotopicMass) >= 0)
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
        protected SpectralMatch FineScorePeptide(int id, Ms2ScanWithSpecificMass scan, int scanIndex, List<Product> peptideTheorProducts)
        {
            IBioPolymerWithSetMods peptide = PeptideIndex[id];

            peptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, peptideTheorProducts, CommonParameters.FragmentationParameters);

            List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan, peptideTheorProducts, CommonParameters);

            double thisScore = CalculatePeptideScore(scan.TheScan, matchedIons);
            // Use the mass this search selects on, so the final notch is consistent with the coarse
            // selection above.
            int notch = MassDiffAcceptor.Accepts(scan.GetPrecursorMassForSearch(CommonParameters), peptide.MonoisotopicMass);

            bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;
            bool scoreImprovement = PeptideSpectralMatches[scanIndex] == null || (thisScore - PeptideSpectralMatches[scanIndex].RunnerUpScore) > -SpectralMatch.ToleranceForScoreDifferentiation;

            if (meetsScoreCutoff && scoreImprovement)
            {
                if (PeptideSpectralMatches[scanIndex] == null)
                {
                    // Same match-type switch ClassicSearchEngine makes, so a modern search over a
                    // nucleic acid database reports oligos as OSMs instead of mislabelling them PSMs.
                    PeptideSpectralMatches[scanIndex] = GlobalVariables.AnalyteType == AnalyteType.Oligo
                        ? new OligoSpectralMatch(peptide, notch, thisScore, scanIndex, scan, CommonParameters, matchedIons)
                        : new PeptideSpectralMatch(peptide, notch, thisScore, scanIndex, scan, CommonParameters, matchedIons);
                }
                else
                {
                    PeptideSpectralMatches[scanIndex].AddOrReplace(peptide, thisScore, notch, CommonParameters.ReportAllAmbiguity, matchedIons);
                }
            }

            return PeptideSpectralMatches[scanIndex];
        }

        protected void FineScorePeptides(List<int> peptideIds, Ms2ScanWithSpecificMass scan, int scanIndex, byte[] scoringTable, 
            DissociationType dissociationType, List<Product> peptideTheorProducts)
        {
            // Every rough-scored candidate is fine-scored.
            //
            // This used to stop early: candidates were taken in descending rough-score order and the loop broke
            // once no remaining candidate's rough score could beat the best fine score found so far. That is a
            // sound bound — the rough score over-counts matched fragments, never under-counts — but the bound
            // was accumulated per partition. Each partition ran the loop over only its own candidates starting
            // from zero, so the number of candidates fine-scored, and with it the runner-up and the Delta Score
            // derived from it, depended on how the database had been split. Splitting is a memory decision, so
            // reported Delta Score, PEP and q-values moved with the amount of RAM available. Scoring everything
            // makes the result a function of the data alone: 1-, 2- and 4-partition runs are byte-identical.
            //
            // Measured on the mouse proteome this costs nothing outside run-to-run variance, partly because the
            // ordering below materialises the whole sort either way, and partly because partitioning itself
            // shrinks each loop. The ordering is kept so that equal-scoring matches are still recorded in a
            // stable order.
            foreach (int id in peptideIds.OrderByDescending(p => scoringTable[p]))
            {
                FineScorePeptide(id, scan, scanIndex, peptideTheorProducts);
            }
        }

        /// <summary>
        /// Whether every peptide's mass, as the bin searches see it (undefined reads as negative infinity), is at or above the one
        /// before it. The indexing engine sorts the index this way, which is what lets a mass bound become a peptide id bound.
        /// </summary>
        protected static bool IsSortedForBinSearch(List<IBioPolymerWithSetMods> peptideIndex)
        {
            for (int id = 1; id < peptideIndex.Count; id++)
            {
                if (MassForBinSearch(peptideIndex[id].MonoisotopicMass) < MassForBinSearch(peptideIndex[id - 1].MonoisotopicMass))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// The last peptide id whose mass, as the bin searches see it, is at or below <paramref name="peptideMassToLookFor"/>, or -1
        /// if there is none. Only meaningful on an index for which <see cref="IsSortedForBinSearch"/> is true.
        /// </summary>
        protected static int LastPeptideIdAtOrBelow(List<IBioPolymerWithSetMods> peptideIndex, double peptideMassToLookFor)
        {
            int low = 0;
            int high = peptideIndex.Count - 1;
            int result = -1;

            while (low <= high)
            {
                int mid = low + ((high - low) / 2);

                if (MassForBinSearch(peptideIndex[mid].MonoisotopicMass) <= peptideMassToLookFor)
                {
                    result = mid;
                    low = mid + 1;
                }
                else
                {
                    high = mid - 1;
                }
            }

            return result;
        }

        /// <summary>
        /// The position of the last id in the bin that is at or below <paramref name="maxPeptideId"/>, or -1 if there is none. Ids
        /// ascend within a bin, so this is an integer search over the bin itself, with no peptide read.
        /// </summary>
        protected static int LastBinPositionAtOrBelowId(ReadOnlySpan<int> bin, int maxPeptideId)
        {
            int low = 0;
            int high = bin.Length - 1;
            int result = -1;

            while (low <= high)
            {
                int mid = low + ((high - low) / 2);

                if (bin[mid] <= maxPeptideId)
                {
                    result = mid;
                    low = mid + 1;
                }
                else
                {
                    high = mid - 1;
                }
            }

            return result;
        }

        /// <summary>
        /// The last peptide index checked with <see cref="IsSortedForBinSearch"/> and the answer, so the check runs once per index
        /// rather than once per scan. Held by reference: an index is built once and not changed while it is searched. A search
        /// starts all its threads at once, so the check is made under a lock: it reads every peptide, and each thread repeating it
        /// cost more than the scoring it was guarding.
        /// </summary>
        private sealed record BinSearchOrder(List<IBioPolymerWithSetMods> PeptideIndex, bool Sorted);

        private volatile BinSearchOrder _binSearchOrder;

        private readonly object _binSearchOrderLock = new object();

        private int _binSearchOrderChecks;

        /// <summary>
        /// How many times this engine has checked a peptide index's order, so a test can see the check is not repeated.
        /// </summary>
        internal int BinSearchOrderChecks => System.Threading.Volatile.Read(ref _binSearchOrderChecks);

        private bool PeptideIndexIsSortedForBinSearch(List<IBioPolymerWithSetMods> peptideIndex)
        {
            BinSearchOrder known = _binSearchOrder;
            if (known == null || !ReferenceEquals(known.PeptideIndex, peptideIndex))
            {
                lock (_binSearchOrderLock)
                {
                    known = _binSearchOrder;
                    if (known == null || !ReferenceEquals(known.PeptideIndex, peptideIndex))
                    {
                        System.Threading.Interlocked.Increment(ref _binSearchOrderChecks);
                        known = new BinSearchOrder(peptideIndex, IsSortedForBinSearch(peptideIndex));
                        _binSearchOrder = known;
                    }
                }
            }
            return known.Sorted;
        }

        /// <summary>
        /// Deprecated.
        /// </summary>
        protected void IndexedScoring(Indexing.FragmentIndex FragmentIndex, List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor,
            double highestMassPeptideToLookFor, List<IBioPolymerWithSetMods> peptideIndex, MassDiffAcceptor massDiffAcceptor, double maxMassThatFragmentIonScoreIsDoubled, DissociationType dissociationType)
        {
            // The window ends at the last peptide no heavier than the upper bound. On an index sorted by mass, as the indexing engine
            // builds it, ids are in mass order, so that is one peptide id for the whole scan, and each bin's end is the last id at or
            // below it: the same position the per-bin mass search finds, without reading a peptide in every bin.
            bool windowEndsAtAPeptideId = !Double.IsInfinity(highestMassPeptideToLookFor) && PeptideIndexIsSortedForBinSearch(peptideIndex);
            int lastPeptideIdInWindow = windowEndsAtAPeptideId ? LastPeptideIdAtOrBelow(peptideIndex, highestMassPeptideToLookFor) : -1;

            // OpenSearchMode accepts every mass, so a candidate's mass need not be read to ask it; the glyco and crosslink searches
            // use it, and that read was a trip to a peptide object for every candidate reaching the cutoff. Exactly that type only:
            // a subclass may override Accepts.
            bool acceptorAcceptsEveryMass = massDiffAcceptor.GetType() == typeof(OpenSearchMode);

            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++) //binsToSearch is the list of fragment in Spectra
            {
                ReadOnlySpan<int> peptideIdsInThisBin = FragmentIndex[binsToSearch[i]];

                // An empty bin used to be a null list here, which would have thrown; callers only ever pass
                // populated bins. Skipping keeps that contract from depending on the caller getting it right.
                if (peptideIdsInThisBin.IsEmpty)
                {
                    continue;
                }

                //get index for minimum monoisotopic allowed
                int lowestPeptideMassIndex = Double.IsInfinity(lowestMassPeptideToLookFor) ? 0 : BinarySearchBinForFirstAtOrAbove(peptideIdsInThisBin, lowestMassPeptideToLookFor, peptideIndex);

                // get index for highest mass allowed
                int highestPeptideMassIndex = peptideIdsInThisBin.Length - 1;

                if (windowEndsAtAPeptideId)
                {
                    highestPeptideMassIndex = LastBinPositionAtOrBelowId(peptideIdsInThisBin, lastPeptideIdInWindow);

                    // nothing in this bin is light enough for the window
                    if (highestPeptideMassIndex < 0)
                    {
                        continue;
                    }
                }
                else if (!Double.IsInfinity(highestMassPeptideToLookFor)) //check if the highest mass is infinity
                {
                    // An index out of mass order cannot bound the window by id, so each bin is searched by mass. The walk below never
                    // moves the end: the search already returns the last entry at or below the bound.
                    highestPeptideMassIndex = BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, highestMassPeptideToLookFor, peptideIndex); //get index for maximum monoisotopic allowed

                    // nothing in this bin is light enough for the window
                    if (highestPeptideMassIndex < 0)
                    {
                        continue;
                    }

                    for (int j = highestPeptideMassIndex; j < peptideIdsInThisBin.Length; j++) //find the highest peptide mass index 
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
                        if (scoringTable[id] == 0 && (acceptorAcceptsEveryMass || massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMass) >= 0))
                        {
                            idsOfPeptidesPossiblyObserved.Add(id);
                        }

                        // mark the peptide as potentially observed so it doesn't get added more than once
                        scoringTable[id] = 1;
                    }
                }
                else
                {   
                    // account the peptide index shown in the bin
                    for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++) // iterate through the peptide index in the bin
                    {
                        int id = peptideIdsInThisBin[j];
                        scoringTable[id]++;

                        // if the score of the peptide >3 (counts > 3 times), and the mass difference is accepted, add the peptide to the list of peptides possibly observed
                        if (scoringTable[id] == byteScoreCutoff && (acceptorAcceptsEveryMass || massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMass) >= 0))
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
        protected List<int> GetBinsToSearch(Ms2ScanWithSpecificMass scan, Indexing.FragmentIndex FragmentIndex, DissociationType dissociationType)
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

                    if (!FragmentIndex[fragmentBin].IsEmpty)
                    {
                        binsToSearch.Add(fragmentBin);
                    }

                    // add complementary ions
                    if (CommonParameters.AddCompIons)
                    {
                        if (complementaryIonConversionDictionary.ContainsKey(CommonParameters.DissociationType))
                        {
                            foreach (double massshift in complementaryIonConversionDictionary[CommonParameters.DissociationType])
                            {
                                double protonMassShift = massshift.ToMass(1);
                                fragmentBin = (int)Math.Round((scan.PrecursorMass + protonMassShift - masses[i]) / 1.0005079);

                                if (!FragmentIndex[fragmentBin].IsEmpty)
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
                        if (!FragmentIndex[fragmentBin].IsEmpty)
                        {
                            binsToSearch.Add(fragmentBin);
                        }
                    }

                    // add complementary ions
                    if (CommonParameters.AddCompIons)
                    {
                        //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut and just straight up adding the bins assuming that they're z=1

                        if (complementaryIonConversionDictionary.ContainsKey(CommonParameters.DissociationType)) 
                        {
                            foreach (double massShift in complementaryIonConversionDictionary[CommonParameters.DissociationType])
                            {
                                double protonMassShift = massShift.ToMass(1);
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
                                    if (!FragmentIndex[fragmentBin].IsEmpty)
                                    {
                                        binsToSearch.Add(fragmentBin);
                                    }
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
