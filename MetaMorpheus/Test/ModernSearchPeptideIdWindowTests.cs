using Chemistry;
using EngineLayer;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// The modern search's coarse scoring cuts every fragment bin a scan hits to the scan's precursor window. It used to find
    /// both ends of that cut with a binary search over peptide masses in every bin, reading a peptide object at every probe. On
    /// an index sorted by mass, with ids ascending in every bin, the window is one range of peptide ids for the whole scan, and
    /// each bin is cut by searching its own ids. These tests hold the new cut to the old one, kept verbatim below, bin by bin and
    /// window by window, and hold a whole search to what it found before, match for match.
    /// </summary>
    [TestFixture]
    public static class ModernSearchPeptideIdWindowTests
    {
        /// <summary>
        /// Anagram peptides give runs of equal masses and X residues give undefined masses, the two cases a mass window has got
        /// wrong before; the rest are random proteins for distinct masses and gaps.
        /// </summary>
        private static List<IBioPolymerWithSetMods> SortedPeptideIndex()
        {
            var sequences = new List<string>
            {
                "ACDEFKCADEFKDACEFKEDCAFKFEDCAKGHILMRHGILMRIGHLMRLIGHMRMLIGHR",
                "PEPTIDEKPETPIDEKTPEPIDEKXAAGKAXAGKAAXGKQQAQNIEKMSKAQQNQIEKMSK",
                "WYVSTRNPQGHKLLVVEKAGSTPRDDEEKKXLLPPGGKSSTTNNQQRWWYYFFK",
            };
            var random = new Random(20261015);
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
        /// Bounds on every distinct mass, just either side of it, beyond both ends of the index, both infinities and NaN.
        /// </summary>
        private static List<double> Bounds(List<IBioPolymerWithSetMods> peptideIndex)
        {
            var defined = peptideIndex.Select(p => p.MonoisotopicMass).Where(m => !double.IsNaN(m)).Distinct().OrderBy(m => m).ToList();
            var bounds = new List<double> { defined[0] - 100, defined[^1] + 100, double.MaxValue, double.MinValue,
                double.NegativeInfinity, double.PositiveInfinity, double.NaN };
            foreach (double m in defined)
            {
                bounds.Add(m);
                bounds.Add(m + 1e-6);
                bounds.Add(m - 1e-6);
            }
            return bounds;
        }

        [Test]
        public static void TheFixtureHasEqualMassRunsAndUndefinedMasses()
        {
            var masses = SortedPeptideIndex().Select(p => p.MonoisotopicMass).ToList();
            Assert.That(masses.Count(double.IsNaN), Is.GreaterThan(0), "undefined masses");
            Assert.That(masses.Where(m => !double.IsNaN(m)).GroupBy(m => m).Any(g => g.Count() > 2), Is.True, "runs of equal masses");
            Assert.That(Probe.IsSorted(SortedPeptideIndex()), Is.True, "SortByMonoisotopicMass output must count as sorted");
        }

        [Test]
        public static void FirstPeptideIdAtOrAboveMatchesALinearScan()
        {
            var peptideIndex = SortedPeptideIndex();
            var failures = new List<string>();
            foreach (double bound in Bounds(peptideIndex))
            {
                int expected = peptideIndex.Count;
                for (int id = 0; id < peptideIndex.Count; id++)
                {
                    if (MassForBinSearch(peptideIndex[id]) >= bound) { expected = id; break; }
                }
                int actual = Probe.FirstIdAtOrAbove(peptideIndex, bound);
                if (actual != expected && failures.Count < 10)
                {
                    failures.Add($"bound {bound:R}: expected {expected}, got {actual}");
                }
            }
            Assert.That(failures, Is.Empty, string.Join("; ", failures));
            Assert.That(Probe.FirstIdAtOrAbove(new List<IBioPolymerWithSetMods>(), 1000), Is.EqualTo(0), "an empty index has none");
        }

        [Test]
        public static void FirstBinPositionAtOrAboveIdMatchesALinearScan()
        {
            var random = new Random(11);
            var failures = new List<string>();
            for (int trial = 0; trial < 2000 && failures.Count < 10; trial++)
            {
                int length = random.Next(0, 40);
                int[] bin = Enumerable.Range(0, 500).OrderBy(_ => random.Next()).Take(length).OrderBy(id => id).ToArray();
                int minId = random.Next(-2, 502);
                int expected = bin.Length;
                for (int k = 0; k < bin.Length; k++)
                {
                    if (bin[k] >= minId) { expected = k; break; }
                }
                int actual = Probe.FirstPositionAtOrAbove(bin, minId);
                if (actual != expected)
                {
                    failures.Add($"bin [{string.Join(",", bin)}] minId {minId}: expected {expected}, got {actual}");
                }
            }
            Assert.That(failures, Is.Empty, string.Join("; ", failures));
        }

        /// <summary>
        /// Every bin, for every window: the positions the id window gives receive a score exactly when the per-bin mass search
        /// gave them one. Windows run over every combination of finite, infinite and undefined bounds, including windows with
        /// nothing inside, where the search now stops before the first bin; the reference shows no bin had anything to score there.
        /// </summary>
        [Test]
        public static void OnASortedIndexTheIdWindowCutsEveryBinWhereThePerBinMassSearchDid()
        {
            var peptideIndex = SortedPeptideIndex();
            var probe = new Probe(peptideIndex);
            var random = new Random(31337);
            List<int[]> bins = RandomBins(peptideIndex.Count, random);
            List<double> bounds = Bounds(peptideIndex);

            var windows = new List<(double Lowest, double Highest)>();
            foreach (double a in new[] { double.NegativeInfinity, double.PositiveInfinity, double.NaN })
            {
                foreach (double b in new[] { double.NegativeInfinity, double.PositiveInfinity, double.NaN })
                {
                    windows.Add((a, b));
                }
            }
            for (int w = 0; w < 600; w++)
            {
                double highest = bounds[random.Next(bounds.Count)];
                windows.Add((double.NegativeInfinity, highest));
                windows.Add((bounds[random.Next(bounds.Count)], double.PositiveInfinity));
                windows.Add((highest - random.NextDouble() * 50, highest));
                windows.Add((bounds[random.Next(bounds.Count)], bounds[random.Next(bounds.Count)])); // often inverted
            }

            int compared = 0, nonEmpty = 0, emptyWindows = 0;
            var failures = new List<string>();
            foreach (var (lowest, highest) in windows)
            {
                (int FirstId, int LastId)? idWindow = probe.Window(lowest, highest);
                Assert.That(idWindow, Is.Not.Null, "a sorted index must give an id window");
                bool windowIsEmpty = idWindow.Value.FirstId > idWindow.Value.LastId;
                emptyWindows += windowIsEmpty ? 1 : 0;

                foreach (int[] bin in bins)
                {
                    var (refStart, refEnd) = ReferenceRange(lowest, highest, bin, peptideIndex);
                    var (start, end) = probe.Range(idWindow, lowest, highest, bin);
                    bool refEmpty = refStart > refEnd;
                    bool empty = start > end;
                    compared++;
                    nonEmpty += refEmpty ? 0 : 1;

                    if (refEmpty != empty || (!refEmpty && (refStart != start || refEnd != end)) || (windowIsEmpty && !refEmpty))
                    {
                        if (failures.Count < 10)
                        {
                            failures.Add($"window [{lowest:R}, {highest:R}] bin of {bin.Length}: reference ({refStart}, {refEnd}), id window ({start}, {end})");
                        }
                    }
                }
            }

            Assert.That(failures, Is.Empty, string.Join("; ", failures));
            Assert.That(nonEmpty, Is.GreaterThan(compared / 10), "most comparisons should score something, or the test proves little");
            Assert.That(emptyWindows, Is.GreaterThan(0), "windows with nothing inside must be covered");
        }

        /// <summary>
        /// An index out of mass order cannot be cut by id, so the search keeps the per-bin mass search.
        /// </summary>
        [Test]
        public static void AnIndexOutOfMassOrderHasNoIdWindowAndKeepsThePerBinMassSearch()
        {
            var random = new Random(99);
            var shuffled = SortedPeptideIndex().OrderBy(_ => random.Next()).ToList();
            Assume.That(Probe.IsSorted(shuffled), Is.False);
            var probe = new Probe(shuffled);
            List<int[]> bins = RandomBins(shuffled.Count, random);

            foreach (var (lowest, highest) in new[] { (500.0, 1500.0), (double.NegativeInfinity, 900.0), (700.0, double.PositiveInfinity) })
            {
                var idWindow = probe.Window(lowest, highest);
                Assert.That(idWindow, Is.Null);
                foreach (int[] bin in bins)
                {
                    Assert.That(probe.Range(idWindow, lowest, highest, bin), Is.EqualTo(ReferenceRange(lowest, highest, bin, shuffled)));
                }
            }
        }

        /// <summary>
        /// A whole modern search over real spectra and a real database finds exactly what it found before: the same peptides, by
        /// reference, in the same order, with the same notches, scores and runner-up scores, for each fragmentation mode the
        /// coarse scoring branches on, with and without complementary ions, and for windows bounded on both sides, on one side and
        /// on neither. Low-resolution CID with complementary ions is left out: it throws on these spectra before and after this
        /// change, because that branch computes a complementary bin without scaling it and can index below zero.
        /// </summary>
        [Test]
        [TestCase(DissociationType.HCD, false, "OneMM 10 ppm")]
        [TestCase(DissociationType.HCD, true, "OneMM 10 ppm")]
        [TestCase(DissociationType.HCD, false, "Open")]
        [TestCase(DissociationType.HCD, false, "[-187, +inf) Da")]
        [TestCase(DissociationType.HCD, false, "[-500, +5] Da")]
        [TestCase(DissociationType.LowCID, false, "OneMM 10 ppm")]
        public static void ModernSearchFindsExactlyWhatThePerBinMassSearchFound(DissociationType dissociationType, bool addCompIons, string acceptorName)
        {
            var (parameters, scans, index) = RealSearchInputs(dissociationType, addCompIons);
            MassDiffAcceptor acceptor = Acceptor(acceptorName);

            AssertSameMatches(scans, index.PeptideIndex, index.FragmentIndex, parameters, acceptor, expectIdWindow: true, minimumMatches: 50);
        }

        /// <summary>
        /// The same, on an index whose ids were shuffled (and each bin re-sorted so ids still ascend in it): out of mass order, so
        /// the search falls back to the per-bin mass search and must still find what it found before. A mass search over an index
        /// out of mass order misses most candidates, so this finds few matches; it is the agreement that is being tested.
        /// </summary>
        [Test]
        public static void ModernSearchOnAnIndexOutOfMassOrderFindsWhatItFoundBefore()
        {
            var (parameters, scans, index) = RealSearchInputs(DissociationType.HCD, false);
            var random = new Random(5);
            int[] newIdOf = Enumerable.Range(0, index.PeptideIndex.Count).OrderBy(_ => random.Next()).ToArray();
            var shuffledPeptides = new IBioPolymerWithSetMods[index.PeptideIndex.Count];
            for (int id = 0; id < newIdOf.Length; id++)
            {
                shuffledPeptides[newIdOf[id]] = index.PeptideIndex[id];
            }

            var binStart = (int[])index.FragmentIndex.BinStart.Clone();
            var peptideIds = new int[index.FragmentIndex.PeptideIds.Length];
            for (int b = 0; b < index.FragmentIndex.Length; b++)
            {
                for (int k = binStart[b]; k < binStart[b + 1]; k++)
                {
                    peptideIds[k] = newIdOf[index.FragmentIndex.PeptideIds[k]];
                }
                Array.Sort(peptideIds, binStart[b], binStart[b + 1] - binStart[b]);
            }

            var peptideIndex = shuffledPeptides.ToList();
            Assume.That(Probe.IsSorted(peptideIndex), Is.False);
            AssertSameMatches(scans, peptideIndex, new FragmentIndex(binStart, peptideIds), parameters, Acceptor("OneMM 10 ppm"), expectIdWindow: false, minimumMatches: 1);
        }

        private static MassDiffAcceptor Acceptor(string name) => name switch
        {
            "OneMM 10 ppm" => SearchTask.GetMassDiffAcceptor(new PpmTolerance(10), MassDiffAcceptorType.OneMM, null),
            "Open" => new OpenSearchMode(),
            "[-187, +inf) Da" => new IntervalMassDiffAcceptor("modopen", new[] { new DoubleRange(-187, double.PositiveInfinity) }),
            "[-500, +5] Da" => new IntervalMassDiffAcceptor("wide", new[] { new DoubleRange(-500, 5) }),
            _ => throw new ArgumentException(name),
        };

        private static void AssertSameMatches(Ms2ScanWithSpecificMass[] scans, List<IBioPolymerWithSetMods> peptideIndex, FragmentIndex fragmentIndex,
            CommonParameters parameters, MassDiffAcceptor acceptor, bool expectIdWindow, int minimumMatches)
        {
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            var expected = new SpectralMatch[scans.Length];
            new ReferenceEngine(expected, scans, peptideIndex, fragmentIndex, parameters, fsp, acceptor).Run();

            var actual = new SpectralMatch[scans.Length];
            var engine = new ModernSearchEngine(actual, scans, peptideIndex, fragmentIndex, 0, parameters, fsp, acceptor, 0, new List<string>());
            engine.Run();
            Assert.That(new Probe(peptideIndex).Window(400, 500) != null, Is.EqualTo(expectIdWindow), "premise: which cut the search takes");

            int psms = 0;
            var failures = new List<string>();
            for (int i = 0; i < scans.Length; i++)
            {
                string difference = Difference(expected[i], actual[i]);
                psms += expected[i] == null ? 0 : 1;
                if (difference != null && failures.Count < 10)
                {
                    failures.Add($"scan {i}: {difference}");
                }
            }
            Assert.That(failures, Is.Empty, string.Join("; ", failures));
            Assert.That(psms, Is.GreaterThanOrEqualTo(minimumMatches), "the search must find enough to compare");
        }

        private static string Difference(SpectralMatch x, SpectralMatch y)
        {
            if (x == null || y == null)
            {
                return x == null && y == null ? null : $"one search matched and the other did not ({x != null} vs {y != null})";
            }
            if (BitConverter.DoubleToInt64Bits(x.Score) != BitConverter.DoubleToInt64Bits(y.Score)
                || BitConverter.DoubleToInt64Bits(x.RunnerUpScore) != BitConverter.DoubleToInt64Bits(y.RunnerUpScore))
            {
                return $"score {x.Score:R}/{x.RunnerUpScore:R} vs {y.Score:R}/{y.RunnerUpScore:R}";
            }
            if (x.Notch != y.Notch || x.FullSequence != y.FullSequence)
            {
                return $"{x.FullSequence} notch {x.Notch} vs {y.FullSequence} notch {y.Notch}";
            }
            var hx = x.BestMatchingBioPolymersWithSetMods.ToList();
            var hy = y.BestMatchingBioPolymersWithSetMods.ToList();
            if (hx.Count != hy.Count)
            {
                return $"{hx.Count} vs {hy.Count} hypotheses";
            }
            for (int k = 0; k < hx.Count; k++)
            {
                if (!ReferenceEquals(hx[k].SpecificBioPolymer, hy[k].SpecificBioPolymer) || hx[k].Notch != hy[k].Notch
                    || !hx[k].MatchedIons.Select(m => (m.NeutralTheoreticalProduct.Annotation, BitConverter.DoubleToInt64Bits(m.Mz)))
                        .SequenceEqual(hy[k].MatchedIons.Select(m => (m.NeutralTheoreticalProduct.Annotation, BitConverter.DoubleToInt64Bits(m.Mz)))))
                {
                    return $"hypothesis {k} differs";
                }
            }
            return null;
        }

        private static readonly Dictionary<(DissociationType, bool), (CommonParameters, Ms2ScanWithSpecificMass[], IndexingResults)> RealInputs = new();

        /// <summary>
        /// The small yeast run against its database with decoys, oxidized methionine variable and carbamidomethyl fixed: about
        /// 1,500 precursors against tens of thousands of peptides. Built once per fragmentation setting.
        /// </summary>
        private static (CommonParameters, Ms2ScanWithSpecificMass[], IndexingResults) RealSearchInputs(DissociationType dissociationType, bool addCompIons)
        {
            lock (RealInputs)
            {
                if (RealInputs.TryGetValue((dissociationType, addCompIons), out var cached))
                {
                    return cached;
                }

                var parameters = new CommonParameters(dissociationType: dissociationType, addCompIons: addCompIons, maxThreadsToUsePerFile: 4,
                    digestionParams: new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 7));
                string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
                var index = BuildRealIndex(parameters);

                string spectra = Path.Combine(folder, "SmallCalibratible_Yeast.mzML");
                var file = MsDataFileReader.GetDataFile(spectra).LoadAllStaticData();
                var scans = MetaMorpheusTask.GetMs2Scans(file, spectra, parameters).OrderBy(s => s.GetPrecursorMassForSearch(parameters)).ToArray();

                var inputs = (parameters, scans, index);
                RealInputs[(dissociationType, addCompIons)] = inputs;
                return inputs;
            }
        }

        /// <summary>
        /// The small yeast database with decoys, oxidized methionine variable and carbamidomethyl fixed, indexed afresh each call.
        /// </summary>
        private static IndexingResults BuildRealIndex(CommonParameters parameters)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
            var proteins = ProteinDbLoader.LoadProteinFasta(Path.Combine(folder, "smalldb.fasta"), true, DecoyType.Reverse, false, out _);
            var variable = GlobalVariables.AllModsKnown.Where(m => m.IdWithMotif == "Oxidation on M").ToList();
            var fixedMods = GlobalVariables.AllModsKnown.Where(m => m.IdWithMotif == "Carbamidomethyl on C").ToList();
            return (IndexingResults)new IndexingEngine(proteins, variable, fixedMods, null, null, null, 0, DecoyType.Reverse, parameters,
                new List<(string, CommonParameters)> { ("", parameters) }, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();
        }

        /// <summary>
        /// 300 bins, each an ascending random subset of the ids, including empty bins and bins holding every id.
        /// </summary>
        private static List<int[]> RandomBins(int peptideCount, Random random)
        {
            var ids = Enumerable.Range(0, peptideCount).ToArray();
            var bins = new List<int[]>();
            for (int b = 0; b < 300; b++)
            {
                int length = b % 50 == 0 ? 0 : b % 50 == 1 ? ids.Length : random.Next(1, Math.Min(ids.Length, 400));
                bins.Add(ids.OrderBy(_ => random.Next()).Take(length).OrderBy(id => id).ToArray());
            }
            return bins;
        }

        /// <summary>
        /// A sorted index made nearly sorted the way reading it back from its cache does: a few peptides' masses moved by a unit
        /// in the ninth decimal place, the precision masses are rounded to. Each is the middle of a run of three equal masses, alternately
        /// moved down (so it sits below the peptide before it) and up (so the peptide after it sits below it). Returns the moved ids.
        /// </summary>
        private static List<int> MakeNearlySorted(List<IBioPolymerWithSetMods> sortedIndex, int howMany)
        {
            var massField = typeof(PeptideWithSetModifications).GetField("_monoisotopicMass",
                System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            var candidates = new List<int>();
            for (int id = 1; id + 1 < sortedIndex.Count; id++)
            {
                double mass = sortedIndex[id].MonoisotopicMass;
                // a run of exactly three, so moving the middle one disturbs one other peptide at most
                if (!double.IsNaN(mass) && mass == sortedIndex[id - 1].MonoisotopicMass && mass == sortedIndex[id + 1].MonoisotopicMass
                    && (id < 2 || sortedIndex[id - 2].MonoisotopicMass != mass) && (id + 2 >= sortedIndex.Count || sortedIndex[id + 2].MonoisotopicMass != mass))
                {
                    candidates.Add(id);
                }
            }
            Assert.That(candidates.Count, Is.GreaterThanOrEqualTo(howMany), "the fixture needs runs of equal masses to disturb");

            var moved = new List<int>();
            for (int k = 0; k < howMany; k++)
            {
                int id = candidates[(int)((long)k * candidates.Count / howMany)];
                double mass = sortedIndex[id].MonoisotopicMass;
                massField.SetValue(sortedIndex[id], (double?)(k % 2 == 0 ? mass - 1e-9 : mass + 1e-9));
                moved.Add(id);
            }
            return moved;
        }

        /// <summary>
        /// About 15,000 peptides: the anagram and X-residue proteins of <see cref="SortedPeptideIndex"/> and 400 random ones, so a
        /// few peptides out of order stay within the one in a thousand a nearly sorted index may have.
        /// </summary>
        private static List<IBioPolymerWithSetMods> LargeSortedPeptideIndex()
        {
            var random = new Random(4242);
            const string residues = "ACDEFGHIKLMNPQRSTVWY";
            var sequences = new List<string>
            {
                "ACDEFKCADEFKDACEFKEDCAFKFEDCAKGHILMRHGILMRIGHLMRLIGHMRMLIGHR",
                "PEPTIDEKPETPIDEKTPEPIDEKXAAGKAXAGKAAXGKQQAQNIEKMSKAQQNQIEKMSK",
            };
            for (int p = 0; p < 400; p++)
            {
                sequences.Add(new string(Enumerable.Range(0, 120).Select(_ => residues[random.Next(residues.Length)]).ToArray()));
            }
            var digestion = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 1);
            var peptides = sequences.Select((s, i) => new Protein(s, "L" + i))
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<IBioPolymerWithSetMods>()
                .ToList();
            return IndexingEngine.SortByMonoisotopicMass(peptides);
        }

        [Test]
        public static void TheOrderCheckTellsSortedNearlySortedAndShuffledApart()
        {
            var sorted = LargeSortedPeptideIndex();
            Assert.That(Probe.Analyze(sorted), Is.EqualTo((true, (Dictionary<int, double>)null)));

            List<int> moved = MakeNearlySorted(sorted, 6);
            var (isSorted, standIns) = Probe.Analyze(sorted);
            Assert.That(isSorted, Is.False);
            Assert.That(standIns, Is.Not.Null);
            // a peptide moved down is out of order itself; one moved up puts the rest of its equal-mass run out of order
            for (int k = 0; k < moved.Count; k++)
            {
                Assert.That(standIns.ContainsKey(k % 2 == 0 ? moved[k] : moved[k] + 1), Is.True, $"moved peptide {k}");
            }
            Assert.That(standIns.Keys.All(id => moved.Any(m => Math.Abs(id - m) <= 20)), Is.True, "only the disturbed runs are out of order");

            // the peptides left in order are in order, and each stand-in is the heaviest mass in order before its peptide
            double heaviest = double.NegativeInfinity;
            for (int id = 0; id < sorted.Count; id++)
            {
                double mass = double.IsNaN(sorted[id].MonoisotopicMass) ? double.NegativeInfinity : sorted[id].MonoisotopicMass;
                if (standIns.TryGetValue(id, out double standIn))
                {
                    Assert.That(standIn, Is.EqualTo(heaviest), $"stand-in of {id}");
                    Assert.That(mass, Is.LessThan(heaviest), $"{id} is out of order");
                }
                else
                {
                    Assert.That(mass, Is.GreaterThanOrEqualTo(heaviest), $"{id} is in order");
                    heaviest = mass;
                }
            }

            var random = new Random(3);
            var shuffled = sorted.OrderBy(_ => random.Next()).ToList();
            Assert.That(Probe.Analyze(shuffled), Is.EqualTo((false, (Dictionary<int, double>)null)), "far out of order is not nearly sorted");
        }

        /// <summary>
        /// The nearly sorted case of the bin-by-bin comparison: every bin, for every window, including windows bounded exactly at,
        /// and a unit in the last place either side of, every mass the moved peptides and their neighbours have. Bins holding a moved
        /// peptide go to the mass search, the rest are cut by id, and every one receives scores exactly where the old per-bin mass
        /// search gave them.
        /// </summary>
        [Test]
        public static void OnANearlySortedIndexEveryBinIsCutWhereThePerBinMassSearchCutIt()
        {
            var peptideIndex = LargeSortedPeptideIndex();
            List<int> moved = MakeNearlySorted(peptideIndex, 6);
            var (_, standIns) = Probe.Analyze(peptideIndex);
            Assert.That(standIns, Is.Not.Null, "premise: nearly sorted");
            var probe = new Probe(peptideIndex);
            var random = new Random(8080);
            List<int[]> bins = RandomBins(peptideIndex.Count, random);
            // and bins built around each moved peptide and its neighbours, so every disturbed run is cut many ways
            foreach (int id in moved)
            {
                for (int width = 1; width <= 4; width++)
                {
                    bins.Add(Enumerable.Range(Math.Max(0, id - width), 2 * width + 1).Where(i => i < peptideIndex.Count).ToArray());
                    bins.Add(Enumerable.Range(Math.Max(0, id - width), 2 * width + 1).Where(i => i < peptideIndex.Count && i != id).ToArray());
                }
            }

            var bounds = new List<double>();
            foreach (int id in moved)
            {
                for (int i = Math.Max(0, id - 3); i <= Math.Min(peptideIndex.Count - 1, id + 3); i++)
                {
                    double m = peptideIndex[i].MonoisotopicMass;
                    bounds.AddRange(new[] { m, m - 1e-9, m + 1e-9, m - 5e-10, m + 5e-10, m - 2e-9, m + 2e-9 });
                }
            }
            var windows = new List<(double Lowest, double Highest)>();
            foreach (double a in bounds)
            {
                windows.Add((a, double.PositiveInfinity));
                windows.Add((double.NegativeInfinity, a));
                foreach (double b in bounds)
                {
                    windows.Add((a, b));
                }
            }
            List<double> general = Bounds(peptideIndex);
            for (int w = 0; w < 400; w++)
            {
                windows.Add((general[random.Next(general.Count)], general[random.Next(general.Count)]));
            }

            int compared = 0, nonEmpty = 0, byMass = 0;
            var failures = new List<string>();
            foreach (var (lowest, highest) in windows)
            {
                (int FirstId, int LastId)? idWindow = probe.Window(lowest, highest);
                Assert.That(idWindow, Is.Not.Null, "a nearly sorted index must give an id window");

                foreach (int[] bin in bins)
                {
                    bool holdsOutOfOrder = bin.Any(standIns.ContainsKey);
                    var (refStart, refEnd) = ReferenceRange(lowest, highest, bin, peptideIndex);
                    var (start, end) = probe.Range(idWindow, lowest, highest, bin, holdsOutOfOrder);
                    bool refEmpty = refStart > refEnd;
                    bool empty = start > end;
                    compared++;
                    nonEmpty += refEmpty ? 0 : 1;
                    byMass += holdsOutOfOrder ? 1 : 0;

                    if (refEmpty != empty || (!refEmpty && (refStart != start || refEnd != end)))
                    {
                        if (failures.Count < 10)
                        {
                            failures.Add($"window [{lowest:R}, {highest:R}] bin of {bin.Length} (out of order: {holdsOutOfOrder}): reference ({refStart}, {refEnd}), now ({start}, {end})");
                        }
                    }
                }
            }

            Assert.That(failures, Is.Empty, string.Join("; ", failures));
            Assert.That(nonEmpty, Is.GreaterThan(compared / 10), "most comparisons should score something");
            Assert.That(byMass, Is.GreaterThan(compared / 20), "bins holding a peptide out of order must be covered");
        }

        /// <summary>
        /// The bins marked for the mass search are exactly the bins of a real fragment index holding a peptide out of order.
        /// </summary>
        [Test]
        public static void TheBinsMarkedAreExactlyTheBinsHoldingAPeptideOutOfOrder()
        {
            var parameters = new CommonParameters(maxThreadsToUsePerFile: 4,
                digestionParams: new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 7));
            IndexingResults index = BuildRealIndex(parameters);
            MakeNearlySorted(index.PeptideIndex, 8);
            var (_, standIns) = Probe.Analyze(index.PeptideIndex);
            Assert.That(standIns, Is.Not.Null, "premise: nearly sorted");

            ulong[] marked = Probe.BinsHolding(index.FragmentIndex, standIns.Keys, index.PeptideIndex.Count);
            int holding = 0, wrong = 0;
            for (int b = 0; b < index.FragmentIndex.Length; b++)
            {
                bool expected = false;
                foreach (int id in index.FragmentIndex[b])
                {
                    expected |= standIns.ContainsKey(id);
                }
                holding += expected ? 1 : 0;
                wrong += expected != Probe.Holds(marked, b) ? 1 : 0;
            }
            Assert.That(wrong, Is.EqualTo(0));
            Assert.That(holding, Is.GreaterThan(0));
            Assert.That(Probe.Holds(null, 0), Is.False, "no bit set means no bin holds one");
        }

        /// <summary>
        /// A whole modern search over a nearly sorted index (a fresh copy of the real one with peptides moved as reading back does)
        /// finds exactly what the per-bin mass search found, for every fragmentation mode and window shape the sorted case covers.
        /// </summary>
        [Test]
        [TestCase(DissociationType.HCD, false, "OneMM 10 ppm")]
        [TestCase(DissociationType.HCD, true, "OneMM 10 ppm")]
        [TestCase(DissociationType.HCD, false, "Open")]
        [TestCase(DissociationType.HCD, false, "[-187, +inf) Da")]
        [TestCase(DissociationType.HCD, false, "[-500, +5] Da")]
        [TestCase(DissociationType.LowCID, false, "OneMM 10 ppm")]
        public static void ModernSearchOnANearlySortedIndexFindsExactlyWhatThePerBinMassSearchFound(DissociationType dissociationType, bool addCompIons, string acceptorName)
        {
            var (parameters, scans, _) = RealSearchInputs(dissociationType, addCompIons);
            IndexingResults index = BuildRealIndex(parameters);
            MakeNearlySorted(index.PeptideIndex, 8);
            Assert.That(Probe.Analyze(index.PeptideIndex).Sorted, Is.False, "premise: not sorted");

            AssertSameMatches(scans, index.PeptideIndex, index.FragmentIndex, parameters, Acceptor(acceptorName), expectIdWindow: true, minimumMatches: 50);
        }

        /// <summary>
        /// A scan whose window holds a peptide out of order and no peptide in order: the id window is empty, but the bins holding
        /// that peptide are cut by mass, and the search still finds it, as the per-bin mass search did. The scan is made from the
        /// peptide's own fragments, and its window is half a unit in the last place either side of the peptide's moved mass, so the
        /// rest of its run, a unit above, is outside.
        /// </summary>
        [Test]
        public static void AScanWhoseWindowHoldsOnlyAPeptideOutOfOrderStillFindsIt()
        {
            var (parameters, _, _) = RealSearchInputs(DissociationType.HCD, false);
            IndexingResults index = BuildRealIndex(parameters);
            List<int> moved = MakeNearlySorted(index.PeptideIndex, 8);
            IBioPolymerWithSetMods peptide = index.PeptideIndex[moved[0]]; // the first is moved down, below the rest of its run

            var products = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            double[] mz = products.Select(x => x.NeutralMass.ToMz(1)).OrderBy(x => x).ToArray();
            double[] intensities = mz.Select(_ => 1000.0).ToArray();
            double precursorMz = peptide.MonoisotopicMass.ToMz(1);
            var dataScan = new MsDataScan(new MzSpectrum(mz, intensities, false), 2, 2, true, Polarity.Positive, 1.0, new MzRange(50, 5000), null,
                MZAnalyzerType.Orbitrap, intensities.Sum(), null, null, "scan=2", precursorMz, 1, 1000, precursorMz, 2, DissociationType.HCD, 1, precursorMz);
            var scan = new Ms2ScanWithSpecificMass(dataScan, precursorMz, 1, "synthetic", parameters);
            var acceptor = new IntervalMassDiffAcceptor("half a unit", new[] { new DoubleRange(-5e-10, 5e-10) });

            var probe = new Probe(index.PeptideIndex);
            var window = acceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(scan.PrecursorMass).Single();
            var idWindow = probe.Window(window.Minimum, window.Maximum);
            Assert.That(idWindow.Value.FirstId, Is.GreaterThan(idWindow.Value.LastId), "premise: no peptide in order is inside the window");

            AssertSameMatches(new[] { scan }, index.PeptideIndex, index.FragmentIndex, parameters, acceptor, expectIdWindow: true, minimumMatches: 1);
        }

        #region the per-bin mass search, verbatim from ModernSearchEngine before the id window

        private static (int start, int end) ReferenceRange(double lowestPeptideMassToLookFor, double highestPeptideMassToLookFor, ReadOnlySpan<int> bin,
            List<IBioPolymerWithSetMods> peptideIndex)
        {
            int start = 0;
            int end = bin.Length - 1;

            if (!double.IsPositiveInfinity(highestPeptideMassToLookFor))
            {
                end = ReferenceLastAtOrBelow(bin, highestPeptideMassToLookFor, peptideIndex);

                // every peptide in this bin is heavier than the window allows
                if (end < 0)
                {
                    return (0, -1);
                }
            }

            if (!double.IsNegativeInfinity(lowestPeptideMassToLookFor))
            {
                start = ReferenceFirstAtOrAbove(bin, lowestPeptideMassToLookFor, peptideIndex);
            }

            return (start, end);
        }

        private static double ReferenceMass(double monoisotopicMass) => double.IsNaN(monoisotopicMass) ? double.NegativeInfinity : monoisotopicMass;

        private static int ReferenceFirstAtOrAbove(ReadOnlySpan<int> bin, double peptideMassToLookFor, List<IBioPolymerWithSetMods> peptideIndex)
        {
            int low = 0;
            int high = bin.Length - 1;
            int result = bin.Length;

            while (low <= high)
            {
                int mid = low + ((high - low) / 2);

                if (ReferenceMass(peptideIndex[bin[mid]].MonoisotopicMass) >= peptideMassToLookFor)
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

        private static int ReferenceLastAtOrBelow(ReadOnlySpan<int> bin, double peptideMassToLookFor, List<IBioPolymerWithSetMods> peptideIndex)
        {
            int low = 0;
            int high = bin.Length - 1;
            int result = -1;

            while (low <= high)
            {
                int mid = low + ((high - low) / 2);

                if (ReferenceMass(peptideIndex[bin[mid]].MonoisotopicMass) <= peptideMassToLookFor)
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
        /// ModernSearchEngine as it was before the id window: RunSpecific (less progress reporting) and IndexScoreScan verbatim,
        /// cutting each bin with <see cref="ReferenceRange"/>. Increments and fine scoring are the engine's own, which this change
        /// does not touch.
        /// </summary>
        private sealed class ReferenceEngine : ModernSearchEngine
        {
            internal ReferenceEngine(SpectralMatch[] psms, Ms2ScanWithSpecificMass[] scans, List<IBioPolymerWithSetMods> peptideIndex, FragmentIndex fragmentIndex,
                CommonParameters parameters, List<(string, CommonParameters)> fsp, MassDiffAcceptor acceptor)
                : base(psms, scans, peptideIndex, fragmentIndex, 0, parameters, fsp, acceptor, 0, new List<string>())
            {
            }

            protected override MetaMorpheusEngineResults RunSpecific()
            {
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
                        Ms2ScanWithSpecificMass scan = ListOfSortedMs2Scans[scanIndex];
                        ReferenceIndexScoreScan(scan, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, CommonParameters.DissociationType);
                        FineScorePeptides(idsOfPeptidesPossiblyObserved, scan, scanIndex, scoringTable, CommonParameters.DissociationType, peptideTheorProducts, scoreSorter);
                    }
                });

                foreach (SpectralMatch psm in PeptideSpectralMatches.Where(p => p != null))
                {
                    psm.ResolveAllAmbiguities();
                }

                return new MetaMorpheusEngineResults(this);
            }

            private void ReferenceIndexScoreScan(Ms2ScanWithSpecificMass scan, ScanScoringTable scoringTable, byte byteScoreCutoff, List<int> peptidesPossiblyObserved, DissociationType dissociationType)
            {
                List<AllowedIntervalWithNotch> notches = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(scan.GetPrecursorMassForSearch(CommonParameters)).ToList();
                double lowestMassPeptideToLookFor = notches.Min(p => p.Minimum);
                double highestMassPeptideToLookFor = notches.Max(p => p.Maximum);

                // retire the previous scan's scores (conserves memory compared to allocating a new array)
                scoringTable.BeginScan();
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
                            var (start, end) = ReferenceRange(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, PeptideIndex);

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
                                        var (start, end) = ReferenceRange(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, PeptideIndex);

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
                            var (start, end) = ReferenceRange(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, PeptideIndex);

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
                                        var (start, end) = ReferenceRange(lowestMassPeptideToLookFor, highestMassPeptideToLookFor, bin, PeptideIndex);

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
        }

        #endregion

        /// <summary>
        /// Reaches ModernSearchEngine's protected members; the base constructor only assigns fields.
        /// </summary>
        private sealed class Probe : ModernSearchEngine
        {
            internal Probe(List<IBioPolymerWithSetMods> peptideIndex)
                : base(null, null, peptideIndex, null, 0, new CommonParameters(), null, new OpenSearchMode(), 0, new List<string>())
            {
            }

            internal (int FirstId, int LastId)? Window(double lowest, double highest) => PeptideIdWindowOf(lowest, highest);

            internal (int start, int end) Range((int FirstId, int LastId)? idWindow, double lowest, double highest, int[] bin, bool holdsOutOfOrder = false)
                => BinRangeToIncrement(idWindow, lowest, highest, holdsOutOfOrder, bin, 0);

            internal static (bool Sorted, Dictionary<int, double> StandIns) Analyze(List<IBioPolymerWithSetMods> peptideIndex) => AnalyzeBinSearchOrder(peptideIndex);

            internal static ulong[] BinsHolding(FragmentIndex fragmentIndex, IEnumerable<int> ids, int peptideCount) => FindBinsHolding(fragmentIndex, ids, peptideCount, 4);

            internal static bool Holds(ulong[] bins, int bin) => HoldsOutOfOrderPeptide(bins, bin);

            internal static bool IsSorted(List<IBioPolymerWithSetMods> peptideIndex) => IsSortedForBinSearch(peptideIndex);

            internal static int FirstIdAtOrAbove(List<IBioPolymerWithSetMods> peptideIndex, double mass) => FirstPeptideIdAtOrAbove(peptideIndex, mass);

            internal static int FirstPositionAtOrAbove(int[] bin, int minId) => FirstBinPositionAtOrAboveId(bin, minId);
        }
    }
}
