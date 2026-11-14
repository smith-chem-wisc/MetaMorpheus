using Chemistry;
using EngineLayer.Util;
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

        /// <summary>
        /// Whether the per-scan score table stamps cells instead of clearing them. Worth it when a scan
        /// touches a small slice of the peptide index; see <see cref="ScanScoringTable.IsWorthStamping"/> and
        /// <see cref="ScanScoringTable.MaxWindowShareWorthStamping"/>.
        /// </summary>
        protected readonly bool UseStampedScoringTable;

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
            UseStampedScoringTable = ScanScoringTable.IsWorthStamping(massDiffAcceptor) && PrecursorWindowsAreNarrowEnoughToStamp();
        }

        /// <summary>
        /// Whether this engine's scans, on average, put few enough of the index in their precursor window for stamping to win;
        /// see <see cref="ScanScoringTable.MaxWindowShareWorthStamping"/>. When that cannot be measured (no scans, no peptides,
        /// or an index not sorted by mass) the acceptor alone decides, as it did before. Either table gives the same scores.
        /// </summary>
        private bool PrecursorWindowsAreNarrowEnoughToStamp()
        {
            if (ListOfSortedMs2Scans == null || PeptideIndex == null || !PeptideIndexIsSortedForBinSearch(PeptideIndex))
            {
                return true;
            }

            var precursorMasses = ListOfSortedMs2Scans.Select(scan => scan.GetPrecursorMassForSearch(CommonParameters));
            return ScanScoringTable.IsWindowNarrowEnoughToStamp(MeanWindowShareOfIndex(PeptideIndex, precursorMasses, MassDiffAcceptor));
        }

        /// <summary>
        /// The share of <paramref name="peptideIndex"/> inside the precursor window <see cref="IndexScoreScan"/> scores for each mass,
        /// averaged over <paramref name="precursorMasses"/>, or zero when there are none. The coarse scoring touches only peptides
        /// in that window, so this bounds how much of the scoring table a scan touches. Two binary searches per scan, and only
        /// meaningful on an index for which <see cref="IsSortedForBinSearch"/> is true.
        /// </summary>
        internal static double MeanWindowShareOfIndex(List<IBioPolymerWithSetMods> peptideIndex, IEnumerable<double> precursorMasses,
            MassDiffAcceptor massDiffAcceptor)
        {
            if (peptideIndex.Count == 0)
            {
                return 0;
            }

            long peptidesInWindows = 0;
            int scans = 0;
            foreach (double precursorMass in precursorMasses)
            {
                // The same min/max collapse IndexScoreScan makes: one window spanning every allowed notch.
                List<AllowedIntervalWithNotch> notches = massDiffAcceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(precursorMass).ToList();
                int lastBelowWindow = LastPeptideIdAtOrBelow(peptideIndex, notches.Min(p => p.Minimum));
                int lastInWindow = LastPeptideIdAtOrBelow(peptideIndex, notches.Max(p => p.Maximum));
                peptidesInWindows += Math.Max(0, lastInWindow - lastBelowWindow);
                scans++;
            }

            return scans == 0 ? 0 : (double)peptidesInWindows / scans / peptideIndex.Count;
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
                var scoringTable = new ScanScoringTable(PeptideIndex.Count, UseStampedScoringTable);
                List<int> idsOfPeptidesPossiblyObserved = new List<int>(PeptideIndex.Count);
                List<Product> peptideTheorProducts = new List<Product>();
                var scoreSorter = new DescendingScoreSorter();

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
                    FineScorePeptides(idsOfPeptidesPossiblyObserved, scan, scanIndex, scoringTable, CommonParameters.DissociationType, peptideTheorProducts, scoreSorter);

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
        protected void IndexScoreScan(Ms2ScanWithSpecificMass scan, ScanScoringTable scoringTable, byte byteScoreCutoff, List<int> peptidesPossiblyObserved, DissociationType dissociationType)
        {
            // get allowed theoretical masses from the known experimental mass
            // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)	
            // this is just PRELIMINARY precursor-mass filtering	
            // additional checks are made later to ensure that the theoretical precursor mass is acceptable
            List<AllowedIntervalWithNotch> notches = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(scan.GetPrecursorMassForSearch(CommonParameters)).ToList();
            double lowestMassPeptideToLookFor = notches.Min(p => p.Minimum);
            double highestMassPeptideToLookFor = notches.Max(p => p.Maximum);

            // retire the previous scan's scores (conserves memory compared to allocating a new array)
            scoringTable.BeginScan();
            peptidesPossiblyObserved.Clear();

            // The window as one range of peptide ids for the whole scan, so each bin below is cut to it by searching the bin's own
            // ids instead of reading a peptide's mass at every probe; null on an index out of mass order. See PeptideIdWindowOf.
            // A bin holding one of the few peptides out of order in an index read back from its cache is still cut by mass.
            (int FirstId, int LastId)? peptideIdWindow = PeptideIdWindowOf(lowestMassPeptideToLookFor, highestMassPeptideToLookFor);
            ulong[] binsHoldingOutOfOrder = peptideIdWindow == null ? null : BinsHoldingOutOfOrderPeptides();

            // no peptide is heavy enough and light enough at once, so no bin has anything to score; a bin cut by mass might
            if (peptideIdWindow is (int firstIdInWindow, int lastIdInWindow) && firstIdInWindow > lastIdInWindow && binsHoldingOutOfOrder == null)
            {
                return;
            }

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
                        var (start, end) = BinRangeToIncrement(peptideIdWindow, lowestMassPeptideToLookFor, highestMassPeptideToLookFor,
                            HoldsOutOfOrderPeptide(binsHoldingOutOfOrder, fragmentBin), bin, scan.PrecursorMass);

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
                                    var (start, end) = BinRangeToIncrement(peptideIdWindow, lowestMassPeptideToLookFor, highestMassPeptideToLookFor,
                                        HoldsOutOfOrderPeptide(binsHoldingOutOfOrder, fragmentBin), bin, scan.PrecursorMass);

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
                        var (start, end) = BinRangeToIncrement(peptideIdWindow, lowestMassPeptideToLookFor, highestMassPeptideToLookFor,
                            HoldsOutOfOrderPeptide(binsHoldingOutOfOrder, b), bin, scan.PrecursorMass);

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
                                    var (start, end) = BinRangeToIncrement(peptideIdWindow, lowestMassPeptideToLookFor, highestMassPeptideToLookFor,
                                        HoldsOutOfOrderPeptide(binsHoldingOutOfOrder, b), bin, scan.PrecursorMass);

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
        /// The precursor window [lowest, highest] as the ids of the first and last peptide inside it, or null when the index is
        /// not sorted by mass, or not nearly enough (see <see cref="AnalyzeBinSearchOrder"/>). On a nearly sorted index the ids
        /// bound the peptides in mass order, and a bin holding a peptide out of order is cut by mass instead. First above last means
        /// no peptide in order is inside. An infinite bound leaves that end of the index open, as it leaves that end of every bin
        /// open in <see cref="GetFirstAndLastIndexesInBinToIncrement"/>.
        ///
        /// The indexing engine sorts the index by mass (undefined masses first, read as negative infinity) and numbers peptides
        /// by position, and ids ascend within every bin. So "mass at or below the bound" is "id at or below the last id at or below
        /// the bound", and likewise at the lower end, and a bin's part of the window is two searches over the bin's own ids. The
        /// mass searches this replaces read a peptide object, and rounded its mass, at every probe of every bin: on the mouse
        /// proteome (9.25 M peptides) that was ~18,500 peptide reads per scan against ~2,500 scores incremented, and 86% of the
        /// search's CPU.
        /// </summary>
        protected (int FirstId, int LastId)? PeptideIdWindowOf(double lowestPeptideMassToLookFor, double highestPeptideMassToLookFor)
        {
            if (PeptideIndex == null)
            {
                return null;
            }

            BinSearchOrder order = BinSearchOrderOf(PeptideIndex);
            if (!order.IdWindowUsable)
            {
                return null;
            }

            int firstId = double.IsNegativeInfinity(lowestPeptideMassToLookFor) ? 0
                : FirstPeptideIdAtOrAbove(PeptideIndex, lowestPeptideMassToLookFor, order.BoundMassOfOutOfOrder);
            int lastId = double.IsPositiveInfinity(highestPeptideMassToLookFor) ? PeptideIndex.Count - 1
                : LastPeptideIdAtOrBelow(PeptideIndex, highestPeptideMassToLookFor, order.BoundMassOfOutOfOrder);
            return (firstId, lastId);
        }

        /// <summary>
        /// The first and last positions in the bin of the peptides to add a +1 score to: with the scan's window as peptide ids
        /// (<see cref="PeptideIdWindowOf"/>), two searches over the bin's ids; without, or for a bin holding a peptide out of mass
        /// order, <see cref="GetFirstAndLastIndexesInBinToIncrement"/>. Either way the positions that receive a score are the same:
        /// first above last when none do.
        /// </summary>
        protected (int start, int end) BinRangeToIncrement((int FirstId, int LastId)? peptideIdWindow, double lowestPeptideMassToLookFor,
            double highestPeptideMassToLookFor, bool binHoldsOutOfOrderPeptide, ReadOnlySpan<int> bin, double precursorMass)
        {
            if (peptideIdWindow is (int firstId, int lastId) && !binHoldsOutOfOrderPeptide)
            {
                return (FirstBinPositionAtOrAboveId(bin, firstId), LastBinPositionAtOrBelowId(bin, lastId));
            }

            return GetFirstAndLastIndexesInBinToIncrement(lowestPeptideMassToLookFor, highestPeptideMassToLookFor, bin, precursorMass);
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
        protected void IncrementPeptideScoresInBin(int start, int end, ReadOnlySpan<int> bin, ScanScoringTable scoringTable, Ms2ScanWithSpecificMass scan, byte byteScoreCutoff,
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
                    scoringTable.Set(peptideId, 1);
                }
            }
            else
            {
                // add +1 to each peptide score
                for (int p = start; p <= end; p++)
                {
                    int peptideId = bin[p];
                    byte score = scoringTable.Increment(peptideId);

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

        protected void FineScorePeptides(List<int> peptideIds, Ms2ScanWithSpecificMass scan, int scanIndex, ScanScoringTable scoringTable, 
            DissociationType dissociationType, List<Product> peptideTheorProducts, DescendingScoreSorter scoreSorter)
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
            foreach (int id in scoreSorter.Sort(peptideIds, scoringTable))
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
        /// if there is none. Only meaningful on an index for which <see cref="IsSortedForBinSearch"/> is true, or on a nearly sorted
        /// one given the stand-in masses <see cref="AnalyzeBinSearchOrder"/> found for it: the masses searched are then in order,
        /// and the answer is exact for every peptide in order.
        /// </summary>
        protected static int LastPeptideIdAtOrBelow(List<IBioPolymerWithSetMods> peptideIndex, double peptideMassToLookFor,
            Dictionary<int, double> boundMassOfOutOfOrder = null)
        {
            int low = 0;
            int high = peptideIndex.Count - 1;
            int result = -1;

            while (low <= high)
            {
                int mid = low + ((high - low) / 2);

                if (BoundMass(peptideIndex, mid, boundMassOfOutOfOrder) <= peptideMassToLookFor)
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
        /// The first peptide id whose mass, as the bin searches see it, is at or above <paramref name="peptideMassToLookFor"/>, or the
        /// index's count if there is none. The lower-bound twin of <see cref="LastPeptideIdAtOrBelow"/>, and like it only meaningful on
        /// an index for which <see cref="IsSortedForBinSearch"/> is true, or on a nearly sorted one given its stand-in masses.
        /// </summary>
        protected static int FirstPeptideIdAtOrAbove(List<IBioPolymerWithSetMods> peptideIndex, double peptideMassToLookFor,
            Dictionary<int, double> boundMassOfOutOfOrder = null)
        {
            int low = 0;
            int high = peptideIndex.Count - 1;
            int result = peptideIndex.Count;

            while (low <= high)
            {
                int mid = low + ((high - low) / 2);

                if (BoundMass(peptideIndex, mid, boundMassOfOutOfOrder) >= peptideMassToLookFor)
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
        /// The mass the id-bound searches read for a peptide: its own, as the bin searches see it, or the stand-in
        /// <see cref="AnalyzeBinSearchOrder"/> gave a peptide out of order.
        /// </summary>
        private static double BoundMass(List<IBioPolymerWithSetMods> peptideIndex, int id, Dictionary<int, double> boundMassOfOutOfOrder)
        {
            return boundMassOfOutOfOrder != null && boundMassOfOutOfOrder.TryGetValue(id, out double standIn)
                ? standIn
                : MassForBinSearch(peptideIndex[id].MonoisotopicMass);
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
        /// The position of the first id in the bin that is at or above <paramref name="minPeptideId"/>, or the bin's length if there
        /// is none. The lower-bound twin of <see cref="LastBinPositionAtOrBelowId"/>: an integer search over the bin, no peptide read.
        /// </summary>
        protected static int FirstBinPositionAtOrAboveId(ReadOnlySpan<int> bin, int minPeptideId)
        {
            int low = 0;
            int high = bin.Length - 1;
            int result = bin.Length;

            while (low <= high)
            {
                int mid = low + ((high - low) / 2);

                if (bin[mid] >= minPeptideId)
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
        /// The last peptide index checked with <see cref="AnalyzeBinSearchOrder"/> and the answer, so the check runs once per index
        /// rather than once per scan. Held by reference: an index is built once and not changed while it is searched. A search
        /// starts all its threads at once, so the check is made under a lock: it reads every peptide, and each thread repeating it
        /// cost more than the scoring it was guarding.
        /// </summary>
        /// <param name="Sorted">Every peptide is in mass order, as the bin searches see it.</param>
        /// <param name="BoundMassOfOutOfOrder">Null when sorted, or when too many peptides are out of order for the id window. Else,
        /// for each peptide out of order, the mass the id-bound searches read for it instead of its own.</param>
        private sealed record BinSearchOrder(List<IBioPolymerWithSetMods> PeptideIndex, bool Sorted, Dictionary<int, double> BoundMassOfOutOfOrder)
        {
            public bool IdWindowUsable => Sorted || BoundMassOfOutOfOrder != null;
        }

        private volatile BinSearchOrder _binSearchOrder;

        private readonly object _binSearchOrderLock = new object();

        private int _binSearchOrderChecks;

        /// <summary>
        /// How many times this engine has checked a peptide index's order, so a test can see the check is not repeated.
        /// </summary>
        internal int BinSearchOrderChecks => System.Threading.Volatile.Read(ref _binSearchOrderChecks);

        private bool PeptideIndexIsSortedForBinSearch(List<IBioPolymerWithSetMods> peptideIndex) => BinSearchOrderOf(peptideIndex).Sorted;

        private BinSearchOrder BinSearchOrderOf(List<IBioPolymerWithSetMods> peptideIndex)
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
                        var (sorted, boundMassOfOutOfOrder) = AnalyzeBinSearchOrder(peptideIndex);
                        known = new BinSearchOrder(peptideIndex, sorted, boundMassOfOutOfOrder);
                        _binSearchOrder = known;
                    }
                }
            }
            return known;
        }

        /// <summary>
        /// Whether the index is in mass order, as the bin searches see it, and if not, whether it is nearly enough for the id window.
        ///
        /// An index read back from its cache is not quite in order. Its peptides' masses are recomputed on reading, summing the
        /// modifications in position order rather than the order digestion added them, so a few come out a unit in the last place
        /// apart from the masses the index was sorted by: 9 of the bottom-up nightly's 9.8 M peptides, out of order at 5 places.
        /// Without this, every spectra file searched with a read-back index, which is every file but the one that built it, and
        /// every file once the index is cached, lost the id window.
        ///
        /// Walking the ids in order, a peptide is in order when its mass is at or above the heaviest mass in order before it; the
        /// rest are out of order. The masses of the peptides in order then never decrease, so two peptides in order are never out
        /// of order with each other, and any two that are include a peptide out of order. So a bin holding no peptide out of order
        /// has its masses in order, the per-bin mass search cuts it at exactly the peptides inside the window, and so does the id
        /// window as long as its bounds are right for the peptides in order. They are, if the id-bound searches read, for each
        /// peptide out of order, the heaviest mass in order before it: the masses searched then never decrease, and each peptide in
        /// order is read at its own mass. A bin that does hold a peptide out of order is cut by mass, as before.
        ///
        /// More than one peptide in a thousand out of order (a shuffled index, say) is not nearly sorted: the id window is not used.
        /// </summary>
        protected static (bool Sorted, Dictionary<int, double> BoundMassOfOutOfOrder) AnalyzeBinSearchOrder(List<IBioPolymerWithSetMods> peptideIndex)
        {
            int allowedOutOfOrder = peptideIndex.Count / 1000;
            Dictionary<int, double> boundMassOfOutOfOrder = null;
            double heaviestInOrder = double.NegativeInfinity;

            for (int id = 0; id < peptideIndex.Count; id++)
            {
                double mass = MassForBinSearch(peptideIndex[id].MonoisotopicMass);
                if (mass >= heaviestInOrder)
                {
                    heaviestInOrder = mass;
                    continue;
                }

                boundMassOfOutOfOrder ??= new Dictionary<int, double>();
                if (boundMassOfOutOfOrder.Count == allowedOutOfOrder)
                {
                    return (false, null);
                }
                boundMassOfOutOfOrder.Add(id, heaviestInOrder);
            }

            return (boundMassOfOutOfOrder == null, boundMassOfOutOfOrder);
        }

        /// <summary>
        /// For the bins of <see cref="FragmentIndex"/> holding a peptide out of order in a nearly sorted <see cref="PeptideIndex"/>, so
        /// <see cref="BinRangeToIncrement"/> cuts them by mass. Computed once per index pair, under a lock, as the order is.
        /// </summary>
        private sealed record OutOfOrderBins(BinSearchOrder Order, Indexing.FragmentIndex FragmentIndex, ulong[] Bins);

        private volatile OutOfOrderBins _outOfOrderBins;

        private readonly object _outOfOrderBinsLock = new object();

        /// <summary>
        /// One bit per fragment bin, set when the bin holds a peptide out of mass order, or null when no bin does (a sorted index).
        /// </summary>
        private ulong[] BinsHoldingOutOfOrderPeptides()
        {
            BinSearchOrder order = BinSearchOrderOf(PeptideIndex);
            if (order.BoundMassOfOutOfOrder == null || FragmentIndex == null)
            {
                return null;
            }

            OutOfOrderBins known = _outOfOrderBins;
            if (known == null || !ReferenceEquals(known.Order, order) || !ReferenceEquals(known.FragmentIndex, FragmentIndex))
            {
                lock (_outOfOrderBinsLock)
                {
                    known = _outOfOrderBins;
                    if (known == null || !ReferenceEquals(known.Order, order) || !ReferenceEquals(known.FragmentIndex, FragmentIndex))
                    {
                        known = new OutOfOrderBins(order, FragmentIndex, FindBinsHolding(FragmentIndex, order.BoundMassOfOutOfOrder.Keys,
                            PeptideIndex.Count, CommonParameters.MaxThreadsToUsePerFile));
                        _outOfOrderBins = known;
                    }
                }
            }
            return known.Bins;
        }

        /// <summary>
        /// One bit per bin of <paramref name="fragmentIndex"/>, set when the bin holds any of <paramref name="peptideIds"/>.
        /// </summary>
        protected static ulong[] FindBinsHolding(Indexing.FragmentIndex fragmentIndex, IEnumerable<int> peptideIds, int peptideCount, int maxThreads)
        {
            var wanted = new bool[peptideCount];
            foreach (int id in peptideIds)
            {
                wanted[id] = true;
            }

            var bins = new ulong[(fragmentIndex.Length + 63) / 64];

            // ranges start at multiples of 2^16, so no two share a word of the bit set
            System.Threading.Tasks.Parallel.ForEach(System.Collections.Concurrent.Partitioner.Create(0, fragmentIndex.Length, 1 << 16),
                new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, maxThreads) }, range =>
                {
                    for (int b = range.Item1; b < range.Item2; b++)
                    {
                        foreach (int id in fragmentIndex[b])
                        {
                            if (wanted[id])
                            {
                                bins[b >> 6] |= 1UL << (b & 63);
                                break;
                            }
                        }
                    }
                });

            return bins;
        }

        protected static bool HoldsOutOfOrderPeptide(ulong[] binsHoldingOutOfOrder, int bin)
            => binsHoldingOutOfOrder != null && (binsHoldingOutOfOrder[bin >> 6] & (1UL << (bin & 63))) != 0;

        /// <summary>
        /// Deprecated.
        /// </summary>
        protected void IndexedScoring(Indexing.FragmentIndex FragmentIndex, List<int> binsToSearch, ScanScoringTable scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor,
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
                        scoringTable.Set(id, 1);
                    }
                }
                else
                {   
                    // account the peptide index shown in the bin
                    for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++) // iterate through the peptide index in the bin
                    {
                        int id = peptideIdsInThisBin[j];

                        // if the score of the peptide >3 (counts > 3 times), and the mass difference is accepted, add the peptide to the list of peptides possibly observed
                        if (scoringTable.Increment(id) == byteScoreCutoff && (acceptorAcceptsEveryMass || massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMass) >= 0))
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
