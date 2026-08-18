using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Covers the memory-aware partition guard, the mass sort, the deterministic digestion merge, the flat
    /// fragment-index format and the index-folder naming.
    /// </summary>
    [TestFixture]
    public static class IndexPartitioningTest
    {
        private const BindingFlags PrivateStatic = BindingFlags.NonPublic | BindingFlags.Static;

        private static List<Protein> MakeProteins(int count)
        {
            var proteins = new List<Protein>();
            for (int i = 0; i < count; i++)
            {
                // varied but deterministic sequences, long enough to yield several tryptic peptides each
                proteins.Add(new Protein("MKAVLDPRTGEKVIRSNYFQWLMHDACEGK" + new string('A', i % 7) + "KPEPTIDERKQQ", "P" + i));
            }
            return proteins;
        }

        /// <summary>
        /// Proteins whose tryptic peptides are anagrams of one another, so a large fraction of the peptide
        /// list shares an exact monoisotopic mass. This is the situation that matters in real searches —
        /// a reversed decoy has its target's residue composition, so ~41 % of a real peptide index has an
        /// equal-mass neighbour — and it is what makes tie ordering observable.
        /// </summary>
        private static List<Protein> MakeAnagramProteins()
        {
            // several distinct residue multisets, so the list is a handful of large equal-mass groups rather
            // than one degenerate group; 6 letters gives 720 permutations each
            var permutations = new List<string>();
            foreach (string multiset in new[] { "AGSTDNV", "AGSTDQV", "AGSTDEV" })
            {
                Permute(multiset.ToCharArray(), 0, permutations);
            }

            var proteins = new List<Protein>();
            for (int i = 0; i < permutations.Count; i += 30)
            {
                // "...K" boundaries so trypsin yields one peptide per permutation
                string sequence = "M" + string.Join("", permutations.Skip(i).Take(30).Select(p => p + "K"));
                proteins.Add(new Protein(sequence, "ANA" + i));
            }

            // A handful of peptides with an undefined residue, so some masses are NaN. Real databases contain
            // these (159 of 1.19 M in the mouse proteome), and they are what makes the two sorts diverge:
            // Array.Sort over double[] keys hoists NaN to the front in a pre-pass, which reorders everything
            // else, whereas a CompareTo-based comparison does not.
            proteins.Add(new Protein("MXAGSTDKXXVLDPRTGEKXQWERTYK", "NAN1"));
            proteins.Add(new Protein("MAGXSTDKVLXDPRTGEKQWXERTYK", "NAN2"));
            return proteins;
        }

        private static void Permute(char[] letters, int start, List<string> into)
        {
            if (start == letters.Length - 1)
            {
                into.Add(new string(letters));
                return;
            }
            for (int i = start; i < letters.Length; i++)
            {
                (letters[start], letters[i]) = (letters[i], letters[start]);
                Permute(letters, start + 1, into);
                (letters[start], letters[i]) = (letters[i], letters[start]);
            }
        }

        private static int Suggest(List<Protein> proteins, CommonParameters parameters, int requested,
            out long estimatedBytes, out long budgetBytes, out bool cappedByLimit)
        {
            return IndexPartitioning.SuggestTotalPartitions(proteins, parameters, new List<Modification>(),
                new List<Modification>(), null, null, null, 30000.0, requested,
                out estimatedBytes, out budgetBytes, out cappedByLimit);
        }

        // ------------------------------------------------------------------ partition guard

        /// <summary>
        /// A small database on a machine with normal memory must not be partitioned at all — otherwise the
        /// guard would change results for searches that were fine as configured.
        /// </summary>
        [Test]
        public static void SuggestTotalPartitions_SmallDatabase_LeavesRequestedCountAlone()
        {
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin"));
            int suggested = Suggest(MakeProteins(20), parameters, 1, out long estimated, out long budget, out bool capped);

            Assert.That(suggested, Is.EqualTo(1));
            Assert.That(capped, Is.False);
            Assert.That(estimated, Is.GreaterThan(0), "a non-empty database should produce a non-zero estimate");
            Assert.That(budget, Is.GreaterThan(0));
        }

        /// <summary>Never reduce a partition count the user (or XLSearchTask) deliberately set higher.</summary>
        [Test]
        public static void SuggestTotalPartitions_NeverLowersRequestedCount()
        {
            var parameters = new CommonParameters(totalPartitions: 12,
                digestionParams: new DigestionParams(protease: "trypsin"));
            int suggested = Suggest(MakeProteins(20), parameters, 12, out _, out _, out _);

            Assert.That(suggested, Is.EqualTo(12));
        }

        /// <summary>
        /// The estimate must scale with the database: twice the residues must not need fewer partitions.
        /// Monotonicity is the property the guard actually relies on, and it holds regardless of how much
        /// memory the test machine happens to have.
        /// </summary>
        [Test]
        public static void SuggestTotalPartitions_EstimateGrowsWithDatabaseSize()
        {
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin"));

            Suggest(MakeProteins(50), parameters, 1, out long smallEstimate, out _, out _);
            Suggest(MakeProteins(500), parameters, 1, out long largeEstimate, out _, out _);

            Assert.That(largeEstimate, Is.GreaterThan(smallEstimate));
            // 10x the proteins, so within a factor of two of 10x the bytes
            Assert.That(largeEstimate / (double)smallEstimate, Is.GreaterThan(5.0));
        }

        /// <summary>An empty or null database must be a no-op rather than a divide-by-zero.</summary>
        [Test]
        public static void SuggestTotalPartitions_NoProteins_ReturnsRequestedCount()
        {
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin"));

            Assert.That(Suggest(new List<Protein>(), parameters, 3, out long estimatedEmpty, out _, out _), Is.EqualTo(3));
            Assert.That(estimatedEmpty, Is.EqualTo(0));
            Assert.That(Suggest(null, parameters, 3, out long estimatedNull, out _, out _), Is.EqualTo(3));
            Assert.That(estimatedNull, Is.EqualTo(0));
        }

        /// <summary>
        /// RNA digestion parameters are not DigestionParams, and the indexing engine refuses them anyway;
        /// the guard must decline rather than throw.
        /// </summary>
        [Test]
        public static void SuggestTotalPartitions_NonProteinDigestion_ReturnsRequestedCount()
        {
            var parameters = new CommonParameters(digestionParams: new global::Transcriptomics.Digestion.RnaDigestionParams());

            Assert.That(Suggest(MakeProteins(5), parameters, 2, out long estimated, out _, out _), Is.EqualTo(2));
            Assert.That(estimated, Is.EqualTo(0));
        }

        /// <summary>
        /// The decision itself, with memory supplied rather than measured, so the assertions are about the
        /// code and not about how much RAM the test machine happens to have. 55 % of 10 GB is 5.5 GB budgeted
        /// per partition.
        /// </summary>
        [Test]
        public static void PartitionsForBudget_RaisesOnlyAsFarAsNeeded()
        {
            const long tenGb = 10L * 1024 * 1024 * 1024;

            // comfortably inside the budget -> leave it alone
            Assert.That(IndexPartitioning.PartitionsForBudget(1L << 30, tenGb, 1, 5000, out bool capped1), Is.EqualTo(1));
            Assert.That(capped1, Is.False);

            // 11 GB needed against a 5.5 GB budget -> 2 partitions
            Assert.That(IndexPartitioning.PartitionsForBudget(11L * 1024 * 1024 * 1024, tenGb, 1, 5000, out _), Is.EqualTo(2));

            // 22 GB needed -> 4 partitions
            Assert.That(IndexPartitioning.PartitionsForBudget(22L * 1024 * 1024 * 1024, tenGb, 1, 5000, out _), Is.EqualTo(4));

            // a request above the computed need always wins
            Assert.That(IndexPartitioning.PartitionsForBudget(11L * 1024 * 1024 * 1024, tenGb, 16, 5000, out _), Is.EqualTo(16));
        }

        /// <summary>
        /// Absurd databases must flag that even the maximum will not fit, rather than silently returning a
        /// number that cannot work — and must never ask for more partitions than there are proteins.
        /// </summary>
        [Test]
        public static void PartitionsForBudget_ReportsCapAndRespectsProteinCount()
        {
            const long tenGb = 10L * 1024 * 1024 * 1024;

            int suggested = IndexPartitioning.PartitionsForBudget(long.MaxValue / 2, tenGb, 1, 100_000, out bool capped);
            Assert.That(capped, Is.True, "must report that the partition cap was hit");
            Assert.That(suggested, Is.EqualTo(256), "capped at MaxPartitions");

            // 22 GB would want 4 partitions, but there are only 3 proteins to split
            Assert.That(IndexPartitioning.PartitionsForBudget(22L * 1024 * 1024 * 1024, tenGb, 1, 3, out _), Is.EqualTo(3));

            // an unknown/empty protein count must not produce zero partitions
            Assert.That(IndexPartitioning.PartitionsForBudget(22L * 1024 * 1024 * 1024, tenGb, 1, 0, out _), Is.EqualTo(1));
        }

        /// <summary>A machine reporting no usable memory must be a no-op, not a divide-by-zero.</summary>
        [Test]
        public static void PartitionsForBudget_NoBudget_ReturnsRequestedCount()
        {
            Assert.That(IndexPartitioning.PartitionsForBudget(1L << 40, 0, 2, 500, out bool capped), Is.EqualTo(2));
            Assert.That(capped, Is.False);
        }

        /// <summary>The warning has to name both numbers and the new count, since it is the only signal a user gets.</summary>
        [Test]
        public static void PartitionIncreaseWarning_NamesCountsAndFlagsOutputDifference()
        {
            string warning = IndexPartitioning.PartitionIncreaseWarning(1, 4, 20_000_000_000, 12_000_000_000);

            Assert.That(warning, Does.Contain("1 partition"));
            Assert.That(warning, Does.Contain("4"));
            Assert.That(warning, Does.Contain("18.6 GB").Or.Contain("18.63 GB"), "estimate in GB");
            Assert.That(warning, Does.Contain("Delta Score"), "must say results are not bit-identical");
        }

        // ------------------------------------------------------------------ mass sort

        /// <summary>
        /// The cached-key sort must reproduce the permutation the previous in-place delegate sort produced,
        /// not merely a correct mass ordering. Equal masses are pervasive (a reversed decoy has its target's
        /// exact mass) and both sorts are unstable, so a different tie order would shift every peptide id.
        /// </summary>
        [Test]
        public static void SortByMonoisotopicMass_ReproducesDelegateSortPermutation()
        {
            var digestionParams = new DigestionParams(protease: "trypsin", minPeptideLength: 5);
            var peptides = MakeAnagramProteins()
                .SelectMany(p => p.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .ToList();

            Assert.That(peptides.Count, Is.GreaterThan(100), "need a non-trivial list for the sort to matter");
            var sortedMasses = peptides.Select(p => p.MonoisotopicMass).OrderBy(m => m).ToList();
            int tiedPairs = sortedMasses.Zip(sortedMasses.Skip(1), (a, b) => a == b ? 1 : 0).Sum();
            Assert.That(tiedPairs, Is.GreaterThan(peptides.Count / 2),
                "this test is only meaningful when most peptides share a mass with a neighbour");
            Assert.That(peptides.Count(p => double.IsNaN(p.MonoisotopicMass)), Is.GreaterThan(0),
                "this test is only meaningful when some masses are NaN; see MakeAnagramProteins");

            var expected = new List<PeptideWithSetModifications>(peptides);
            expected.Sort((x, y) => x.MonoisotopicMass.CompareTo(y.MonoisotopicMass));

            var actual = IndexingEngine.SortByMonoisotopicMass(peptides);

            Assert.That(actual.Count, Is.EqualTo(expected.Count));
            for (int i = 0; i < expected.Count; i++)
            {
                Assert.That(ReferenceEquals(actual[i], expected[i]), Is.True,
                    $"peptide at index {i} differs from the delegate sort's permutation");
            }
        }

        /// <summary>NaN masses occur in real databases; they must land where the old comparison put them.</summary>
        [Test]
        public static void SortByMonoisotopicMass_EmptyAndSingleInputsAreHandled()
        {
            Assert.That(IndexingEngine.SortByMonoisotopicMass(new List<PeptideWithSetModifications>()), Is.Empty);

            var one = new Protein("MKAVLDPRTGEK", "P0")
                .Digest(new DigestionParams(protease: "trypsin", minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>().Take(1).ToList();
            Assert.That(IndexingEngine.SortByMonoisotopicMass(one).Count, Is.EqualTo(1));
        }

        // ------------------------------------------------------- deterministic digestion

        /// <summary>
        /// The peptide index used to depend on which digestion thread finished first, which changed reported
        /// Delta Scores between otherwise identical runs. Two runs must now agree element for element,
        /// including order.
        /// </summary>
        [Test]
        public static void IndexingEngine_PeptideIndexOrderIsReproducible()
        {
            var proteins = MakeProteins(300);
            var parameters = new CommonParameters(maxThreadsToUsePerFile: 4, scoreCutoff: 1,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            List<string> RunOnce()
            {
                var engine = new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(), null, null, null,
                    0, DecoyType.Reverse, parameters, fsp, 30000, false, new List<FileInfo>(),
                    TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
                var results = (IndexingResults)engine.Run();
                return results.PeptideIndex.Select(p => p.FullSequence + "@" + p.Protein.Accession + ":" + p.OneBasedStartResidue).ToList();
            }

            var first = RunOnce();
            var second = RunOnce();

            Assert.That(first.Count, Is.GreaterThan(100), "need enough peptides that thread interleaving could differ");
            Assert.That(second, Is.EqualTo(first), "peptide index order must not depend on thread completion order");
        }

        // ------------------------------------------------------ CommonParameters clone

        /// <summary>
        /// Every setting except TotalPartitions must survive the clone. Compared by reflection so that a
        /// parameter added to the constructor later, and forgotten here, fails this test rather than silently
        /// resetting to its default mid-search.
        /// </summary>
        [Test]
        public static void CloneWithNewTotalPartitions_PreservesEveryOtherSetting()
        {
            var original = new CommonParameters(
                taskDescriptor: "descriptor",
                dissociationType: DissociationType.ETD,
                ms2childScanDissociationType: DissociationType.EThcD,
                ms3childScanDissociationType: DissociationType.HCD,
                separationType: "CZE",
                doPrecursorDeconvolution: false,
                useProvidedPrecursorInfo: false,
                deconvolutionIntensityRatio: 7,
                deconvolutionMaxAssumedChargeState: 9,
                reportAllAmbiguity: false,
                addCompIons: true,
                totalPartitions: 3,
                qValueThreshold: 0.02,
                pepQValueThreshold: 0.5,
                qValueCutoffForPepCalculation: 0.004,
                scoreCutoff: 4,
                numberOfPeaksToKeepPerWindow: 111,
                minimumAllowedIntensityRatioToBasePeak: 0.02,
                windowWidthThomsons: 12,
                numberOfWindows: 6,
                normalizePeaksAccrossAllWindows: true,
                trimMs1Peaks: true,
                trimMsMsPeaks: false,
                productMassTolerance: new PpmTolerance(11),
                precursorMassTolerance: new PpmTolerance(6),
                maxThreadsToUsePerFile: 2,
                digestionParams: new DigestionParams(protease: "Asp-N", maxMissedCleavages: 3, minPeptideLength: 6),
                assumeOrphanPeaksAreZ1Fragments: false,
                maxHeterozygousVariants: 2,
                minVariantDepth: 5,
                addTruncations: true,
                useMostAbundantPrecursorIntensity: false,
                precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant,
                rtPredictorName: "SSRCalc3");

            // CustomIons is not a constructor parameter — the constructor reads it from a global dictionary,
            // so give the original a distinctive value. Without this, a clone that simply re-read the global
            // would look correct and the assertion below would prove nothing.
            typeof(CommonParameters).GetProperty(nameof(CommonParameters.CustomIons))
                .SetValue(original, new List<ProductType> { ProductType.c, ProductType.zDot });

            CommonParameters clone = original.CloneWithNewTotalPartitions(9);

            Assert.That(clone.TotalPartitions, Is.EqualTo(9));
            Assert.That(original.TotalPartitions, Is.EqualTo(3), "the original must not be mutated");

            var mismatches = new List<string>();
            foreach (PropertyInfo property in typeof(CommonParameters).GetProperties(BindingFlags.Public | BindingFlags.Instance))
            {
                if (property.Name == nameof(CommonParameters.TotalPartitions) || property.GetIndexParameters().Length > 0)
                {
                    continue;
                }

                object a = property.GetValue(original);
                object b = property.GetValue(clone);

                bool equal = a is System.Collections.IEnumerable && a is not string
                    ? string.Join("|", ((System.Collections.IEnumerable)a).Cast<object>()) ==
                      string.Join("|", ((System.Collections.IEnumerable)(b ?? Array.Empty<object>())).Cast<object>())
                    : Equals(a?.ToString(), b?.ToString());

                if (!equal)
                {
                    mismatches.Add($"{property.Name}: '{a}' -> '{b}'");
                }
            }

            Assert.That(mismatches, Is.Empty, "clone lost settings: " + string.Join("; ", mismatches));
        }

        // ------------------------------------------- flat fragment index serialization

        private static List<int>[] RoundTrip(List<int>[] fragmentIndex, out long fileLength)
        {
            var write = typeof(MetaMorpheusTask).GetMethod("WriteFragmentIndex", PrivateStatic);
            var read = typeof(MetaMorpheusTask).GetMethod("ReadFragmentIndex", PrivateStatic);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "roundtrip-" + Guid.NewGuid().ToString("N") + ".bin");
            try
            {
                write.Invoke(null, new object[] { fragmentIndex, path });
                fileLength = new FileInfo(path).Length;
                return (List<int>[])read.Invoke(null, new object[] { path });
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// The search binary-searches inside a bin, so contents and order both have to survive, and empty
        /// bins have to come back null rather than as empty lists — the search engines test for null.
        /// </summary>
        [Test]
        public static void FragmentIndex_RoundTripPreservesBinsContentsAndOrder()
        {
            var index = new List<int>[50];
            index[0] = new List<int> { 0 };                     // first bin populated
            index[7] = new List<int> { 3, 3, 9, 41 };           // duplicates allowed
            index[8] = new List<int>();                         // empty list, must come back as null
            index[23] = Enumerable.Range(0, 1000).ToList();     // long bin
            index[49] = new List<int> { int.MaxValue, 0 };      // last bin, extreme value, descending
            var expectedNulls = Enumerable.Range(0, 50).Where(i => index[i] == null || index[i].Count == 0).ToArray();

            List<int>[] actual = RoundTrip(index, out long fileLength);

            Assert.That(actual.Length, Is.EqualTo(index.Length));
            foreach (int i in expectedNulls)
            {
                Assert.That(actual[i], Is.Null, $"bin {i} should be null");
            }
            Assert.That(actual[0], Is.EqualTo(index[0]));
            Assert.That(actual[7], Is.EqualTo(index[7]));
            Assert.That(actual[23], Is.EqualTo(index[23]));
            Assert.That(actual[49], Is.EqualTo(index[49]));
            // 3 header ints + 50 count ints + 1007 entries
            Assert.That(fileLength, Is.EqualTo((3 + 50 + 1007) * sizeof(int)));
        }

        /// <summary>
        /// Entries are staged through a fixed buffer, so a bin larger than the buffer exercises the flush path.
        /// Real fragment indexes are hundreds of millions of entries, so this must be right.
        /// </summary>
        [Test]
        public static void FragmentIndex_RoundTripBinLargerThanWriteBuffer()
        {
            const int entries = (1 << 20) + 12345; // one buffer flush plus a remainder
            var index = new List<int>[3];
            index[1] = new List<int>(entries);
            for (int i = 0; i < entries; i++)
            {
                index[1].Add(entries - i); // descending, so order is checkable
            }
            index[2] = new List<int> { 7 };

            List<int>[] actual = RoundTrip(index, out long fileLength);

            Assert.That(actual[0], Is.Null);
            Assert.That(actual[1].Count, Is.EqualTo(entries));
            Assert.That(actual[1][0], Is.EqualTo(entries), "first entry, and order preserved");
            Assert.That(actual[1][entries - 1], Is.EqualTo(1), "last entry");
            Assert.That(actual[2], Is.EqualTo(index[2]));
            Assert.That(fileLength, Is.EqualTo((3L + 3 + entries + 1) * sizeof(int)));
        }

        /// <summary>An index with no fragments at all still has to round-trip; the precursor index can be sparse.</summary>
        [Test]
        public static void FragmentIndex_RoundTripEmptyIndex()
        {
            List<int>[] actual = RoundTrip(new List<int>[10], out long fileLength);

            Assert.That(actual.Length, Is.EqualTo(10));
            Assert.That(actual.All(b => b == null), Is.True);
            Assert.That(fileLength, Is.EqualTo((3 + 10) * sizeof(int)));
        }

        /// <summary>
        /// An index written by an older MetaMorpheus must be rejected, not misinterpreted. The reader is the
        /// last line of defence: GenerateSecondIndexes reads without a try/catch.
        /// </summary>
        [Test]
        public static void FragmentIndex_ReadRejectsForeignFile()
        {
            var read = typeof(MetaMorpheusTask).GetMethod("ReadFragmentIndex", PrivateStatic);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "foreign-" + Guid.NewGuid().ToString("N") + ".bin");
            File.WriteAllBytes(path, new byte[] { 9, 9, 9, 9, 1, 0, 0, 0, 2, 0, 0, 0 });
            try
            {
                var exception = Assert.Throws<TargetInvocationException>(() => read.Invoke(null, new object[] { path }));
                Assert.That(exception.InnerException, Is.TypeOf<MetaMorpheusException>());
            }
            finally
            {
                File.Delete(path);
            }
        }

        // ------------------------------------------------- task-level guard plumbing

        /// <summary>Minimal task so the protected guard helper can be exercised without running a search.</summary>
        private class GuardProbeTask : MetaMorpheusTask
        {
            public GuardProbeTask() : base(MyTask.Search) { CommonParameters = new CommonParameters(); }

            protected override MyTaskResults RunSpecific(string outputFolder, List<DbForTask> dbFilenameList,
                List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
                => throw new NotSupportedException("probe only");

            public CommonParameters Probe(List<Protein> proteins, CommonParameters parameters, ref bool warned)
                => RaisePartitionsToFitMemory(proteins, parameters, new List<Modification>(),
                    new List<Modification>(), null, null, null, 30000.0, ref warned);
        }

        /// <summary>
        /// When the database fits, the guard must hand back the very same instance — SetAllFileSpecificCommonParams
        /// often returns the task's own CommonParameters, and allocating a copy for no reason would make the
        /// object the search runs with differ from the one the task reports.
        /// </summary>
        [Test]
        public static void RaisePartitionsToFitMemory_FitsInMemory_ReturnsSameInstanceAndDoesNotWarn()
        {
            var task = new GuardProbeTask();
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin"));
            bool warned = false;

            var warnings = new List<string>();
            void Collect(object sender, EngineLayer.StringEventArgs e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += Collect;
            try
            {
                CommonParameters result = task.Probe(MakeProteins(20), parameters, ref warned);

                Assert.That(ReferenceEquals(result, parameters), Is.True, "must not copy when no change is needed");
                Assert.That(warned, Is.False);
                Assert.That(warnings, Is.Empty);
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= Collect;
            }
        }

        /// <summary>
        /// A database far too large for any machine must come back with more partitions, on a copy, having
        /// warned exactly once even across repeated calls (the guard runs per spectra file).
        /// </summary>
        [Test]
        public static void RaisePartitionsToFitMemory_TooLargeForMemory_ReturnsRaisedCopyAndWarnsOnce()
        {
            var task = new GuardProbeTask();
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            // One long protein repeated: only the sample is digested, but the residue count it is scaled by is
            // ~10^9, so the estimate is hundreds of GB and exceeds any machine's budget. Keeps the test
            // independent of how much RAM the runner has, and costs a list of references rather than sequences.
            var oneProtein = new Protein("M" + string.Concat(Enumerable.Repeat("AGSTDVK", 700)), "HUGE");
            var huge = Enumerable.Repeat(oneProtein, 200_000).ToList();
            bool warned = false;

            var warnings = new List<string>();
            void Collect(object sender, EngineLayer.StringEventArgs e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += Collect;
            try
            {
                CommonParameters first = task.Probe(huge, parameters, ref warned);
                CommonParameters second = task.Probe(huge, parameters, ref warned);

                Assert.That(first.TotalPartitions, Is.GreaterThan(1), "an impossible database must be partitioned");
                Assert.That(ReferenceEquals(first, parameters), Is.False, "must return a copy, not mutate the input");
                Assert.That(parameters.TotalPartitions, Is.EqualTo(1), "the caller's parameters must be untouched");
                Assert.That(second.TotalPartitions, Is.EqualTo(first.TotalPartitions), "decision must be stable");
                Assert.That(warned, Is.True);
                Assert.That(warnings.Count, Is.EqualTo(1), "must warn once, not once per spectra file");
                Assert.That(warnings[0], Does.Contain("TotalPartitions"));
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= Collect;
            }
        }

        /// <summary>
        /// If even the maximum partition count will not fit, say so rather than leaving the user to wonder why
        /// the search is paging. The first warning alone would imply the problem had been solved.
        /// </summary>
        [Test]
        public static void RaisePartitionsToFitMemory_BeyondPartitionCap_WarnsThatItStillWillNotFit()
        {
            var task = new GuardProbeTask();
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            var oneProtein = new Protein("M" + string.Concat(Enumerable.Repeat("AGSTDVK", 700)), "HUGE");
            var absurd = Enumerable.Repeat(oneProtein, 4_000_000).ToList();
            bool warned = false;

            var warnings = new List<string>();
            void Collect(object sender, EngineLayer.StringEventArgs e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += Collect;
            try
            {
                CommonParameters result = task.Probe(absurd, parameters, ref warned);

                Assert.That(result.TotalPartitions, Is.EqualTo(256), "capped at MaxPartitions");
                Assert.That(warnings.Count, Is.EqualTo(2), "the raise, plus a warning that it still will not fit");
                Assert.That(warnings[1], Does.Contain("may page heavily"));
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= Collect;
            }
        }

        /// <summary>
        /// The sample is capped by peptide count, not protein count: a non-specific protease yields orders of
        /// magnitude more peptides per protein, and every sampled peptide is also fragmented. Without the cap
        /// the estimate becomes slower than the index build it is sizing.
        /// </summary>
        [Test, Timeout(60_000)]
        public static void SuggestTotalPartitions_ManyPeptidesPerProtein_StillReturnsQuickly()
        {
            // ~250 cleavage sites each, so >50,000 peptides are reached well inside the 200-protein sample
            var proteins = Enumerable.Range(0, 300)
                .Select(i => new Protein("M" + string.Concat(Enumerable.Repeat("AGSTDVK", 250)), "BIG" + i))
                .ToList();
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));

            var stopwatch = System.Diagnostics.Stopwatch.StartNew();
            int suggested = Suggest(proteins, parameters, 1, out long estimated, out _, out _);
            stopwatch.Stop();

            Assert.That(estimated, Is.GreaterThan(0));
            Assert.That(suggested, Is.GreaterThanOrEqualTo(1));
            Assert.That(stopwatch.Elapsed.TotalSeconds, Is.LessThan(30),
                "the estimate must stay cheap however many peptides a protein yields");
        }

        // ------------------------------------------------------------- index folder naming

        /// <summary>
        /// Folders are named by a one-second-resolution timestamp. Two partitions indexed inside the same
        /// second used to share a folder and overwrite each other's params and peptide index, which left the
        /// earlier partition unfindable — fatal for XL search, which re-reads each partition's index.
        /// </summary>
        [Test]
        public static void GenerateOutputFolderForIndices_SameSecondCallsGetDistinctFolders()
        {
            var method = typeof(MetaMorpheusTask).GetMethod("GenerateOutputFolderForIndices", PrivateStatic);
            string directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "FolderNaming-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(directory);
            string databasePath = Path.Combine(directory, "db.fasta");
            File.WriteAllText(databasePath, ">P0\nMKAVLDPRTGEK\n");
            var dbList = new List<DbForTask> { new DbForTask(databasePath, false) };

            try
            {
                var folders = new List<string>();
                for (int i = 0; i < 5; i++)
                {
                    folders.Add((string)method.Invoke(null, new object[] { dbList }));
                }

                Assert.That(folders.Distinct().Count(), Is.EqualTo(5), "each call must get its own folder");
                Assert.That(folders.All(Directory.Exists), Is.True, "each folder must actually be created");
            }
            finally
            {
                Directory.Delete(directory, true);
            }
        }
    }
}
