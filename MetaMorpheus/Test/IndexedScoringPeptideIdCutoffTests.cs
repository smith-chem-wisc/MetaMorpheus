using EngineLayer;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// Locks in the upper end of the scoring window in <c>ModernSearchEngine.IndexedScoring</c>, the path the glyco and crosslink
    /// engines score through.
    /// </summary>
    /// <remarks>
    /// Glyco and crosslink searches count, in every fragment bin a scan hits, every peptide no heavier than the precursor plus
    /// 1 Da. Finding that end of each bin was a binary search reading <c>PeptideWithSetModifications.MonoisotopicMass</c> off an
    /// object at every probe, repeated for every bin of every scan; on a semi-tryptic N+O search the scoring was over half of all
    /// CPU. The peptide index is sorted by mass (undefined masses first) and ids are positions in it, so "mass at or below the
    /// bound" is the same as "id at or below the last id at or below the bound". That id is found once per scan, and each bin's
    /// end becomes a search over the bin's own ids. These tests pin the helpers against linear scans, and pin IndexedScoring
    /// against a verbatim copy of the per-bin mass search it replaces: the same ids, in the same order, with the same scores.
    /// </remarks>
    [TestFixture]
    public static class IndexedScoringPeptideIdCutoffTests
    {
        /// <summary>
        /// Anagram peptides give runs of equal masses, and X residues give undefined masses, the two cases the window has got
        /// wrong before. The rest are ordinary proteins for distinct masses and gaps.
        /// </summary>
        private static List<IBioPolymerWithSetMods> SortedPeptideIndex()
        {
            var sequences = new List<string>
            {
                "ACDEFKCADEFKDACEFKEDCAFKFEDCAKGHILMRHGILMRIGHLMRLIGHMRMLIGHR",
                "PEPTIDEKPETPIDEKTPEPIDEKXAAGKAXAGKAAXGKQQAQNIEKMSKAQQNQIEKMSK",
                "WYVSTRNPQGHKLLVVEKAGSTPRDDEEKKXLLPPGGKSSTTNNQQRWWYYFFK",
            };
            var random = new Random(20260915);
            const string residues = "ACDEFGHIKLMNPQRSTVWY";
            for (int p = 0; p < 40; p++)
            {
                sequences.Add(new string(Enumerable.Range(0, 120).Select(_ => residues[random.Next(residues.Length)]).ToArray()));
            }

            var digestion = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 1);
            var peptides = sequences.Select((s, i) => new Protein(s, "P" + i))
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<IBioPolymerWithSetMods>()
                .ToList();
            return IndexingEngine.SortByMonoisotopicMass(peptides);
        }

        private static double MassForBinSearch(IBioPolymerWithSetMods peptide)
            => double.IsNaN(peptide.MonoisotopicMass) ? double.NegativeInfinity : peptide.MonoisotopicMass;

        /// <summary>
        /// Targets on every distinct mass, just either side of it, and beyond both ends of the index.
        /// </summary>
        private static List<double> Targets(List<IBioPolymerWithSetMods> peptideIndex)
        {
            var defined = peptideIndex.Select(p => p.MonoisotopicMass).Where(m => !double.IsNaN(m)).Distinct().OrderBy(m => m).ToList();
            var targets = new List<double> { defined[0] - 100, defined[^1] + 100, double.MaxValue, double.MinValue };
            foreach (double m in defined)
            {
                targets.Add(m);
                targets.Add(m + 1e-6);
                targets.Add(m - 1e-6);
            }
            return targets;
        }

        [Test]
        public static void TheFixtureHasEqualMassRunsAndUndefinedMasses()
        {
            var peptideIndex = SortedPeptideIndex();
            var masses = peptideIndex.Select(p => p.MonoisotopicMass).ToList();
            Assert.That(masses.Count(double.IsNaN), Is.GreaterThan(0), "undefined masses");
            Assert.That(masses.Where(m => !double.IsNaN(m)).GroupBy(m => m).Any(g => g.Count() > 2), Is.True, "runs of equal masses");
            Assert.That(Probe.IsSorted(peptideIndex), Is.True, "SortByMonoisotopicMass output must count as sorted");
        }

        [Test]
        public static void LastPeptideIdAtOrBelowMatchesALinearScan()
        {
            var peptideIndex = SortedPeptideIndex();
            var failures = new List<string>();
            foreach (double target in Targets(peptideIndex))
            {
                int expected = -1;
                for (int id = peptideIndex.Count - 1; id >= 0; id--)
                {
                    if (MassForBinSearch(peptideIndex[id]) <= target) { expected = id; break; }
                }
                int actual = Probe.LastIdAtOrBelow(peptideIndex, target);
                if (actual != expected && failures.Count < 10)
                {
                    failures.Add($"target {target:F6}: expected {expected}, got {actual}");
                }
            }
            Assert.That(failures, Is.Empty, string.Join("; ", failures));
            Assert.That(Probe.LastIdAtOrBelow(new List<IBioPolymerWithSetMods>(), 1000), Is.EqualTo(-1), "an empty index has no such id");
        }

        [Test]
        public static void LastPositionAtOrBelowMatchesALinearScan()
        {
            var random = new Random(7);
            var failures = new List<string>();
            for (int trial = 0; trial < 2000 && failures.Count < 10; trial++)
            {
                int length = random.Next(0, 40);
                int[] bin = Enumerable.Range(0, 500).OrderBy(_ => random.Next()).Take(length).OrderBy(id => id).ToArray();
                int maxId = random.Next(-2, 502);

                int expected = -1;
                for (int k = bin.Length - 1; k >= 0; k--)
                {
                    if (bin[k] <= maxId) { expected = k; break; }
                }
                int actual = Probe.LastPositionAtOrBelow(bin, maxId);
                if (actual != expected)
                {
                    failures.Add($"bin [{string.Join(",", bin)}] maxId {maxId}: expected {expected}, got {actual}");
                }
            }
            Assert.That(failures, Is.Empty, string.Join("; ", failures));
        }

        [Test]
        public static void AnIndexOutOfMassOrderIsNotSorted()
        {
            var peptideIndex = SortedPeptideIndex();
            var swapped = peptideIndex.ToList();
            int last = swapped.Count - 1;
            (swapped[last], swapped[last - 5]) = (swapped[last - 5], swapped[last]);
            Assume.That(MassForBinSearch(swapped[last - 5]) > MassForBinSearch(swapped[last]), "the swap must break mass order");

            Assert.That(Probe.IsSorted(swapped), Is.False);
            Assert.That(Probe.IsSorted(new List<IBioPolymerWithSetMods>()), Is.True, "an empty index is sorted");
        }

        /// <summary>
        /// On a sorted index, every bin and window scores exactly what the per-bin mass search scored: the same ids added in the
        /// same order, and the same scoring table.
        /// </summary>
        [Test]
        public static void IndexedScoringOnASortedIndexScoresExactlyWhatThePerBinMassSearchScored(
            [Values(DissociationType.HCD, DissociationType.LowCID)] DissociationType dissociationType)
        {
            AssertScoringMatchesReference(SortedPeptideIndex(), dissociationType);
        }

        /// <summary>
        /// An index that is not in mass order cannot use the id cutoff, so it must still score what the per-bin mass search scored.
        /// </summary>
        [Test]
        public static void IndexedScoringOnAnUnsortedIndexStillScoresWhatThePerBinMassSearchScored()
        {
            var random = new Random(99);
            var shuffled = SortedPeptideIndex().OrderBy(_ => random.Next()).ToList();
            Assume.That(Probe.IsSorted(shuffled), Is.False);
            AssertScoringMatchesReference(shuffled, DissociationType.HCD);
        }

        /// <summary>
        /// A tolerance acceptor rejects most candidates, so reading each candidate's mass for it still decides what is scored.
        /// </summary>
        [Test]
        public static void IndexedScoringWithAToleranceAcceptorScoresWhatThePerBinMassSearchScored()
        {
            AssertScoringMatchesReference(SortedPeptideIndex(), DissociationType.HCD, new SinglePpmAroundZeroSearchMode(5));
        }

        /// <summary>
        /// Only an acceptor that is exactly <see cref="OpenSearchMode"/> is known to accept every mass. A subclass may override
        /// Accepts, so it must still be asked.
        /// </summary>
        [Test]
        public static void IndexedScoringStillAsksAnAcceptorDerivedFromOpenSearchMode()
        {
            AssertScoringMatchesReference(SortedPeptideIndex(), DissociationType.HCD, new RejectEverythingOpenSearchMode());
        }

        private sealed class RejectEverythingOpenSearchMode : OpenSearchMode
        {
            public override int Accepts(double scanPrecursorMass, double peptideMass) => -1;
        }

        /// <summary>
        /// Every scan's search asks whether the index is in mass order, and a search starts all its threads at once. The index is
        /// checked once, not once per thread: the check reads every peptide.
        /// </summary>
        [Test]
        public static void ScansStartingTogetherCheckTheIndexOrderOnce()
        {
            // Many references to a few peptides: sorted, and long enough that an unguarded check overlaps between threads.
            var few = SortedPeptideIndex();
            var peptideIndex = Enumerable.Range(0, 400_000).Select(i => few[i * few.Count / 400_000]).ToList();
            var fragmentIndex = new FragmentIndex(new[] { 0, 3 }, new[] { 0, 1, 2 });
            var probe = new Probe();

            const int threads = 16;
            using var barrier = new System.Threading.Barrier(threads);
            System.Threading.Tasks.Parallel.For(0, threads, new System.Threading.Tasks.ParallelOptions { MaxDegreeOfParallelism = threads }, _ =>
            {
                barrier.SignalAndWait();
                probe.ScoreOn(fragmentIndex, new List<int> { 0 }, new byte[peptideIndex.Count], 1, new List<int>(), 1000,
                    double.NegativeInfinity, 2000, peptideIndex, new OpenSearchMode(), DissociationType.HCD);
            });

            Assert.That(probe.IndexOrderChecks, Is.EqualTo(1));
        }

        private static void AssertScoringMatchesReference(List<IBioPolymerWithSetMods> peptideIndex, DissociationType dissociationType)
            => AssertScoringMatchesReference(peptideIndex, dissociationType, new OpenSearchMode());

        private static void AssertScoringMatchesReference(List<IBioPolymerWithSetMods> peptideIndex, DissociationType dissociationType, MassDiffAcceptor acceptor)
        {
            var random = new Random(31337);
            var ids = Enumerable.Range(0, peptideIndex.Count).ToArray();

            // 300 bins, each an ascending random subset of the ids, including empty bins and bins holding every id.
            var bins = new List<int[]>();
            for (int b = 0; b < 300; b++)
            {
                int length = b % 50 == 0 ? 0 : b % 50 == 1 ? ids.Length : random.Next(1, Math.Min(ids.Length, 400));
                bins.Add(ids.OrderBy(_ => random.Next()).Take(length).OrderBy(id => id).ToArray());
            }
            var binStart = new int[bins.Count + 1];
            for (int b = 0; b < bins.Count; b++)
            {
                binStart[b + 1] = binStart[b] + bins[b].Length;
            }
            var fragmentIndex = new FragmentIndex(binStart, bins.SelectMany(b => b).ToArray());

            var targets = Targets(peptideIndex);
            var windows = new List<(double lowest, double highest)>
            {
                (double.NegativeInfinity, double.PositiveInfinity),
            };
            for (int w = 0; w < 400; w++)
            {
                double highest = targets[random.Next(targets.Count)];
                windows.Add((double.NegativeInfinity, highest)); // the glyco and crosslink window
                windows.Add((highest - random.NextDouble() * 500, highest));
            }

            var failures = new List<string>();
            foreach (var (lowest, highest) in windows)
            {
                var binsToSearch = Enumerable.Range(0, bins.Count).OrderBy(_ => random.Next()).Take(random.Next(1, 80)).ToList();
                byte cutoff = (byte)random.Next(1, 4);
                double precursor = highest - 1;

                var expectedTable = new byte[peptideIndex.Count];
                var expectedIds = new List<int>();
                ReferenceIndexedScoring(fragmentIndex, binsToSearch, expectedTable, cutoff, expectedIds, precursor, lowest, highest, peptideIndex, acceptor, dissociationType);

                var actualTable = new byte[peptideIndex.Count];
                var actualIds = new List<int>();
                Probe.Score(fragmentIndex, binsToSearch, actualTable, cutoff, actualIds, precursor, lowest, highest, peptideIndex, acceptor, dissociationType);

                if (!actualIds.SequenceEqual(expectedIds) || !actualTable.SequenceEqual(expectedTable))
                {
                    failures.Add($"window [{lowest:F4}, {highest:F4}] cutoff {cutoff}: ids {actualIds.Count} vs {expectedIds.Count}");
                    if (failures.Count >= 5) break;
                }
            }
            Assert.That(failures, Is.Empty, string.Join("; ", failures));
        }

        /// <summary>
        /// ModernSearchEngine.IndexedScoring as it was before the id cutoff, copied verbatim apart from reading the helpers through
        /// <see cref="Probe"/>: the reference the new scoring must reproduce exactly.
        /// </summary>
        private static void ReferenceIndexedScoring(FragmentIndex FragmentIndex, List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor,
            double highestMassPeptideToLookFor, List<IBioPolymerWithSetMods> peptideIndex, MassDiffAcceptor massDiffAcceptor, DissociationType dissociationType)
        {
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                ReadOnlySpan<int> peptideIdsInThisBin = FragmentIndex[binsToSearch[i]];
                if (peptideIdsInThisBin.IsEmpty)
                {
                    continue;
                }

                int lowestPeptideMassIndex = Double.IsInfinity(lowestMassPeptideToLookFor) ? 0 : Probe.FirstAtOrAbove(peptideIdsInThisBin, lowestMassPeptideToLookFor, peptideIndex);
                int highestPeptideMassIndex = peptideIdsInThisBin.Length - 1;

                if (!Double.IsInfinity(highestMassPeptideToLookFor))
                {
                    highestPeptideMassIndex = Probe.PrecursorIndex(peptideIdsInThisBin, highestMassPeptideToLookFor, peptideIndex);
                    if (highestPeptideMassIndex < 0)
                    {
                        continue;
                    }

                    for (int j = highestPeptideMassIndex; j < peptideIdsInThisBin.Length; j++)
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
                    for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                    {
                        int id = peptideIdsInThisBin[j];
                        if (scoringTable[id] == 0 && massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMass) >= 0)
                        {
                            idsOfPeptidesPossiblyObserved.Add(id);
                        }
                        scoringTable[id] = 1;
                    }
                }
                else
                {
                    for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                    {
                        int id = peptideIdsInThisBin[j];
                        scoringTable[id]++;
                        if (scoringTable[id] == byteScoreCutoff && massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMass) >= 0)
                        {
                            idsOfPeptidesPossiblyObserved.Add(id);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Reaches ModernSearchEngine's protected members; the base constructor only assigns fields.
        /// </summary>
        private sealed class Probe : ModernSearchEngine
        {
            internal Probe() : base(null, null, null, null, 0, new CommonParameters(), null, new OpenSearchMode(), 0, new List<string>()) { }

            internal int IndexOrderChecks => BinSearchOrderChecks;

            internal void ScoreOn(FragmentIndex fragmentIndex, List<int> binsToSearch, byte[] scoringTable, byte cutoff, List<int> observed, double precursor,
                double lowest, double highest, List<IBioPolymerWithSetMods> peptideIndex, MassDiffAcceptor acceptor, DissociationType dissociationType)
                => IndexedScoring(fragmentIndex, binsToSearch, scoringTable, cutoff, observed, precursor, lowest, highest, peptideIndex, acceptor, 0, dissociationType);

            internal static bool IsSorted(List<IBioPolymerWithSetMods> peptideIndex) => IsSortedForBinSearch(peptideIndex);

            internal static int LastIdAtOrBelow(List<IBioPolymerWithSetMods> peptideIndex, double mass) => LastPeptideIdAtOrBelow(peptideIndex, mass);

            internal static int LastPositionAtOrBelow(int[] bin, int maxId) => LastBinPositionAtOrBelowId(bin, maxId);

            internal static int FirstAtOrAbove(ReadOnlySpan<int> bin, double mass, List<IBioPolymerWithSetMods> peptideIndex)
                => BinarySearchBinForFirstAtOrAbove(bin, mass, peptideIndex);

            internal static int PrecursorIndex(ReadOnlySpan<int> bin, double mass, List<IBioPolymerWithSetMods> peptideIndex)
                => BinarySearchBinForPrecursorIndex(bin, mass, peptideIndex);

            internal static void Score(FragmentIndex fragmentIndex, List<int> binsToSearch, byte[] scoringTable, byte cutoff, List<int> observed, double precursor,
                double lowest, double highest, List<IBioPolymerWithSetMods> peptideIndex, MassDiffAcceptor acceptor, DissociationType dissociationType)
                => new Probe().IndexedScoring(fragmentIndex, binsToSearch, scoringTable, cutoff, observed, precursor, lowest, highest, peptideIndex, acceptor, 0, dissociationType);
        }
    }
}
