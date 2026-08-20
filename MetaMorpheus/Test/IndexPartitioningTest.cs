using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
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
using System.Runtime.InteropServices;
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

        /// <summary>
        /// A DIAparameters instance that is merely distinct, for identity comparison. Its real constructor
        /// wants an XIC constructor and a grouping engine; the clone test never dereferences the object, it
        /// only asserts the same reference comes across, so an uninitialised instance is sufficient and keeps
        /// the fixture from depending on DIA internals.
        /// </summary>
        private static EngineLayer.DIA.DIAparameters DistinctDiaParameters()
            => (EngineLayer.DIA.DIAparameters)System.Runtime.CompilerServices.RuntimeHelpers
                .GetUninitializedObject(typeof(EngineLayer.DIA.DIAparameters));

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
        /// DissociationType.Custom resolves its product types from a process-global list that the
        /// CommonParameters constructor clears, and only MetaMorpheusEngine.Run() restores it. The estimator
        /// fragments outside any engine, so without restoring it the fragment term would silently be zero and
        /// a Custom search would under-partition — the failure the guard exists to prevent. Compared against
        /// HCD with the same ion types, which should give a similar estimate.
        /// </summary>
        [Test]
        public static void SuggestTotalPartitions_CustomDissociation_StillCountsFragments()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 5);
            List<ProductType> globalBefore = digestion.ProductsFromDissociationType()[DissociationType.Custom];
            try
            {
                // b/y is what HCD produces for peptides, so the two estimates should be comparable
                digestion.ProductsFromDissociationType()[DissociationType.Custom] =
                    new List<ProductType> { ProductType.b, ProductType.y };
                var customParams = new CommonParameters(dissociationType: DissociationType.Custom, digestionParams: digestion);

                // premises: the constructor captured the ion list and then emptied the global
                Assert.That(customParams.CustomIons, Is.Not.Empty, "CustomIons must be captured by the constructor");
                Assert.That(digestion.ProductsFromDissociationType()[DissociationType.Custom], Is.Empty,
                    "the constructor is expected to have cleared the global list");

                var proteins = MakeProteins(60);
                Suggest(proteins, customParams, 1, out long customEstimate, out _, out _);

                var hcdParams = new CommonParameters(dissociationType: DissociationType.HCD,
                    digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
                Suggest(proteins, hcdParams, 1, out long hcdEstimate, out _, out _);

                Assert.That(customEstimate, Is.GreaterThan(0), "a Custom search must not estimate zero fragments");
                Assert.That(customEstimate, Is.EqualTo(hcdEstimate).Within(25).Percent,
                    "Custom with b/y ions should size like HCD");
            }
            finally
            {
                digestion.ProductsFromDissociationType()[DissociationType.Custom] = globalBefore;
            }
        }

        /// <summary>
        /// Under a batch scheduler the job's allocation is the ceiling, and sites that do not enforce it with
        /// cgroups leave the process seeing the whole node — so a job entitled to a few GB of a large node
        /// would otherwise size its index for the node and be killed for exceeding its allocation.
        /// </summary>
        [Test]
        public static void SchedulerMemoryLimitBytes_ReadsSlurmAllocation()
        {
            string perNodeBefore = Environment.GetEnvironmentVariable("SLURM_MEM_PER_NODE");
            string perCpuBefore = Environment.GetEnvironmentVariable("SLURM_MEM_PER_CPU");
            string cpusBefore = Environment.GetEnvironmentVariable("SLURM_CPUS_ON_NODE");
            try
            {
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", null);
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_CPU", null);
                Environment.SetEnvironmentVariable("SLURM_CPUS_ON_NODE", null);
                Assert.That(IndexPartitioning.SchedulerMemoryLimitBytes(), Is.EqualTo(0),
                    "no scheduler variables means no limit to honour");

                // --mem=16G
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", "16384");
                Assert.That(IndexPartitioning.SchedulerMemoryLimitBytes(), Is.EqualTo(16L * 1024 * 1024 * 1024));

                // --mem=0 means unlimited, so it must not be read as a zero-byte allocation
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", "0");
                Assert.That(IndexPartitioning.SchedulerMemoryLimitBytes(), Is.EqualTo(0));

                // --mem-per-cpu=4G on 8 cores
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", null);
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_CPU", "4096");
                Environment.SetEnvironmentVariable("SLURM_CPUS_ON_NODE", "8");
                Assert.That(IndexPartitioning.SchedulerMemoryLimitBytes(), Is.EqualTo(32L * 1024 * 1024 * 1024));

                // per-node wins when both are present, matching how Slurm sets them
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", "8192");
                Assert.That(IndexPartitioning.SchedulerMemoryLimitBytes(), Is.EqualTo(8L * 1024 * 1024 * 1024));

                // garbage must not be taken as a limit
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", "not-a-number");
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_CPU", null);
                Assert.That(IndexPartitioning.SchedulerMemoryLimitBytes(), Is.EqualTo(0));
            }
            finally
            {
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", perNodeBefore);
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_CPU", perCpuBefore);
                Environment.SetEnvironmentVariable("SLURM_CPUS_ON_NODE", cpusBefore);
            }
        }

        /// <summary>
        /// A small Slurm allocation must actually shrink the budget the guard works from, not merely be read.
        /// </summary>
        [Test]
        public static void SuggestTotalPartitions_SmallSlurmAllocation_RaisesPartitions()
        {
            string before = Environment.GetEnvironmentVariable("SLURM_MEM_PER_NODE");
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            // only ~200 are digested for the sample; the rest scale the residue count, so this stays fast
            var proteins = MakeProteins(40_000);
            try
            {
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", null);
                int unconstrained = Suggest(proteins, parameters, 1, out long estimate, out long wideBudget, out _);

                // an allocation far smaller than the estimate must force more partitions than the same
                // database needs when nothing is constraining it
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", "512");
                int constrained = Suggest(proteins, parameters, 1, out _, out long narrowBudget, out _);

                Assert.That(narrowBudget, Is.LessThan(wideBudget), "the allocation must shrink the budget");
                Assert.That(narrowBudget, Is.LessThanOrEqualTo(512L * 1024 * 1024), "budget cannot exceed the allocation");
                Assert.That(constrained, Is.GreaterThan(unconstrained), "a tight allocation must raise partitions");
                Assert.That(estimate, Is.GreaterThan(0));
            }
            finally
            {
                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", before);
            }
        }

        /// <summary>
        /// The budgeted share of the supplied ceiling, read from the class rather than hard-coded, so
        /// recalibrating the fraction cannot silently invalidate the expectations below.
        /// </summary>
        private static double BudgetFraction => (double)typeof(IndexPartitioning)
            .GetField("MemoryBudgetFraction", BindingFlags.NonPublic | BindingFlags.Static)
            .GetRawConstantValue();

        /// <summary>
        /// The decision itself, with the free-memory figure supplied rather than measured, so the assertions
        /// are about the code and not about how much RAM the test machine happens to have.
        /// </summary>
        [Test]
        public static void PartitionsForBudget_RaisesOnlyAsFarAsNeeded()
        {
            const long tenGb = 10L * 1024 * 1024 * 1024;
            long budget = (long)(tenGb * BudgetFraction);

            // comfortably inside the budget -> leave it alone
            Assert.That(IndexPartitioning.PartitionsForBudget(budget / 4, 0, tenGb, 1, 5000, out bool capped1), Is.EqualTo(1));
            Assert.That(capped1, Is.False);

            // just over one budget -> two partitions
            Assert.That(IndexPartitioning.PartitionsForBudget(budget + 1, 0, tenGb, 1, 5000, out _), Is.EqualTo(2));

            // four budgets exactly -> four partitions
            Assert.That(IndexPartitioning.PartitionsForBudget(budget * 4, 0, tenGb, 1, 5000, out _), Is.EqualTo(4));

            // a request above the computed need always wins
            Assert.That(IndexPartitioning.PartitionsForBudget(budget + 1, 0, tenGb, 16, 5000, out _), Is.EqualTo(16));
        }

        /// <summary>
        /// Absurd databases must flag that even the maximum will not fit, rather than silently returning a
        /// number that cannot work — and must never ask for more partitions than there are proteins.
        /// </summary>
        [Test]
        public static void PartitionsForBudget_ReportsCapAndRespectsProteinCount()
        {
            const long tenGb = 10L * 1024 * 1024 * 1024;

            int suggested = IndexPartitioning.PartitionsForBudget(long.MaxValue / 2, 0, tenGb, 1, 100_000, out bool capped);
            Assert.That(capped, Is.True, "must report that the partition cap was hit");
            Assert.That(suggested, Is.EqualTo(256), "capped at MaxPartitions");

            // 22 GB would want 4 partitions, but there are only 3 proteins to split
            Assert.That(IndexPartitioning.PartitionsForBudget(22L * 1024 * 1024 * 1024, 0, tenGb, 1, 3, out _), Is.EqualTo(3));

            // an unknown/empty protein count must not produce zero partitions
            Assert.That(IndexPartitioning.PartitionsForBudget(22L * 1024 * 1024 * 1024, 0, tenGb, 1, 0, out _), Is.EqualTo(1));
        }

        /// <summary>
        /// Bin space does not shrink when the database is split, so it has to come off the budget rather
        /// than be divided into it. Halving what is left must double the partitions, and a fixed cost that
        /// swallows the whole budget must be a no-op rather than a negative-divisor crash.
        /// </summary>
        [Test]
        public static void PartitionsForBudget_FixedBytesComeOffTheBudget()
        {
            const long tenGb = 10L * 1024 * 1024 * 1024;
            long budget = (long)(tenGb * BudgetFraction);

            // with half the budget spent on bin space, the same index needs twice as many partitions
            Assert.That(IndexPartitioning.PartitionsForBudget(budget, 0, tenGb, 1, 5000, out _), Is.EqualTo(1));
            Assert.That(IndexPartitioning.PartitionsForBudget(budget, budget / 2, tenGb, 1, 5000, out _), Is.EqualTo(2));

            // bin space alone exceeding the budget leaves nothing to divide into
            Assert.That(IndexPartitioning.PartitionsForBudget(budget, budget, tenGb, 3, 5000, out bool capped), Is.EqualTo(3));
            Assert.That(capped, Is.False);
        }

        /// <summary>A machine reporting no usable memory must be a no-op, not a divide-by-zero.</summary>
        [Test]
        public static void PartitionsForBudget_NoBudget_ReturnsRequestedCount()
        {
            Assert.That(IndexPartitioning.PartitionsForBudget(1L << 40, 0, 0, 2, 500, out bool capped), Is.EqualTo(2));
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
            Assert.That(warning, Does.Contain("free"), "must say what the budget was derived from");
            Assert.That(warning, Does.Contain("Identifications and scores are unaffected"),
                "must say what does not move");
            Assert.That(warning, Does.Contain("q-values can shift"),
                "and must not claim the partition count is entirely result-neutral, because it is not");
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

        /// <summary>
        /// The compressed layout must hold exactly what the old list-per-bin layout held. Rather than assert
        /// hand-picked bins, this rebuilds the reference form from the same peptides with the arithmetic the
        /// engine used to inline, and compares every bin — an offset that is off by one shifts contents into
        /// the neighbouring bin instead of failing loudly, so a sample would not find it.
        ///
        /// Ascending order within a bin is asserted separately because the search binary-searches inside a
        /// bin by peptide mass, and the peptide index is sorted by mass, so ascending id is what makes that
        /// search valid. A fill pass that visited peptides in any other order would still round-trip.
        /// </summary>
        [Test]
        public static void IndexingEngine_CompressedFragmentIndexMatchesBinPerListLayout()
        {
            const double maxFragmentSize = 3000;
            const int binsPerDalton = 1000;

            var proteins = MakeAnagramProteins();
            var parameters = new CommonParameters(scoreCutoff: 1,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            var engine = new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(), null, null, null,
                0, DecoyType.Reverse, parameters, fsp, maxFragmentSize, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var results = (IndexingResults)engine.Run();

            var reference = new List<int>[(int)Math.Ceiling(maxFragmentSize) * binsPerDalton + 1];
            var fragments = new List<Product>();
            for (int peptideId = 0; peptideId < results.PeptideIndex.Count; peptideId++)
            {
                results.PeptideIndex[peptideId].Fragment(parameters.DissociationType,
                    parameters.DigestionParams.FragmentationTerminus, fragments, parameters.FragmentationParameters);

                foreach (var fragment in fragments)
                {
                    double mass = fragment.NeutralMass;
                    if (mass < maxFragmentSize && mass > 0)
                    {
                        int bin = (int)Math.Round(mass * binsPerDalton);
                        (reference[bin] ??= new List<int>()).Add(peptideId);
                    }
                }
            }

            Assert.That(results.FragmentIndex.Length, Is.EqualTo(reference.Length), "bin count");
            Assert.That(results.FragmentIndex.EntryCount, Is.EqualTo(reference.Sum(b => b?.Count ?? 0)), "total entries");
            Assert.That(results.FragmentIndex.EntryCount, Is.GreaterThan(10_000),
                "the comparison is only worth anything on a well-populated index");

            var mismatches = new List<string>();
            int populated = 0;
            for (int bin = 0; bin < reference.Length && mismatches.Count < 10; bin++)
            {
                ReadOnlySpan<int> actual = results.FragmentIndex[bin];
                int[] expected = reference[bin]?.ToArray() ?? Array.Empty<int>();

                if (expected.Length > 0)
                {
                    populated++;
                }

                if (!actual.SequenceEqual(expected))
                {
                    mismatches.Add($"bin {bin}: expected [{string.Join(",", expected)}] got [{string.Join(",", actual.ToArray())}]");
                }
                else if (actual.Length > 1)
                {
                    for (int i = 1; i < actual.Length; i++)
                    {
                        if (actual[i] < actual[i - 1])
                        {
                            mismatches.Add($"bin {bin} is not ascending at {i}");
                            break;
                        }
                    }
                }
            }

            Assert.That(populated, Is.GreaterThan(1000), $"need many populated bins for this to mean anything; got {populated}");
            Assert.That(mismatches, Is.Empty, string.Join("; ", mismatches));
        }

        // ------------------------------------------------------ CommonParameters clone

        /// <summary>
        /// Every setting except TotalPartitions must survive the clone.
        ///
        /// Two things are needed for that claim to hold. Every constructor parameter must be given a value
        /// distinguishable from its default — a parameter left at its default would reproduce the same default
        /// if the clone dropped it, and the sweep would see no mismatch — which is asserted as a premise below.
        /// And reference-typed settings are compared by identity, not by ToString: the clone is supposed to hand
        /// the very same instance across, and several of these types (deconvolution and fragmentation
        /// parameters especially) do not override ToString, so a distinctively configured object replaced by a
        /// freshly defaulted one of the same type would compare equal.
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
                productMassTolerance_LowRes: new AbsoluteTolerance(0.44),
                deconvolutionMassTolerance: new PpmTolerance(7),
                maxThreadsToUsePerFile: 2,
                digestionParams: new DigestionParams(protease: "Asp-N", maxMissedCleavages: 3, minPeptideLength: 6),
                listOfModsVariable: new List<(string, string)> { ("Common Biological", "Phosphorylation on S") },
                listOfModsFixed: new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C") },
                assumeOrphanPeaksAreZ1Fragments: false,
                maxHeterozygousVariants: 2,
                minVariantDepth: 5,
                addTruncations: true,
                precursorDeconParams: new ClassicDeconvolutionParameters(2, 9, 5, 4),
                productDeconParams: new ClassicDeconvolutionParameters(1, 7, 6, 3),
                useMostAbundantPrecursorIntensity: false,
                diaParameters: DistinctDiaParameters(),
                fragmentationParams: new FragmentationParams(),
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

            var defaults = new CommonParameters();
            var indistinguishable = new List<string>();
            var mismatches = new List<string>();

            foreach (PropertyInfo property in typeof(CommonParameters).GetProperties(BindingFlags.Public | BindingFlags.Instance))
            {
                if (property.Name == nameof(CommonParameters.TotalPartitions) || property.GetIndexParameters().Length > 0)
                {
                    continue;
                }

                object a = property.GetValue(original);
                object b = property.GetValue(clone);
                object def = property.GetValue(defaults);

                // premise: the fixture must actually differ from a default instance, or dropping this
                // parameter from the clone would be invisible
                bool differsFromDefault = property.PropertyType.IsValueType || property.PropertyType == typeof(string)
                    ? !Equals(a, def)
                    : !ReferenceEquals(a, def);
                if (!differsFromDefault)
                {
                    indistinguishable.Add(property.Name);
                }

                // reference types must come across as the same instance; value types and strings by equality
                bool equal = a == null || property.PropertyType.IsValueType || property.PropertyType == typeof(string)
                    ? Equals(a, b)
                    : ReferenceEquals(a, b);

                if (!equal)
                {
                    mismatches.Add($"{property.Name}: '{a}' -> '{b}'");
                }
            }

            Assert.That(mismatches, Is.Empty, "clone lost settings: " + string.Join("; ", mismatches));
            Assert.That(indistinguishable, Is.Empty,
                "these properties match a default-constructed CommonParameters, so this test cannot detect the "
                + "clone dropping them: " + string.Join(", ", indistinguishable));
        }

        // ------------------------------------------------ fragment bin boundary search

        /// <summary>
        /// The bin boundary search returns the inclusive end of the range that receives coarse score
        /// increments, so an index that is too small means candidates are never scored at all. It is
        /// therefore checked against a linear-scan ground truth rather than against itself, across bin
        /// <summary>
        /// Exposes the protected static boundary search. It takes a ReadOnlySpan now, which reflection cannot
        /// pass, so the test reaches it by derivation instead. The constructor exists only to satisfy the base
        /// class; nothing ever instantiates this.
        /// </summary>
        private sealed class BinSearchProbe : ModernSearchEngine
        {
            private BinSearchProbe() : base(null, null, null, null, 0, new CommonParameters(), null, new OpenSearchMode(), 0, new List<string>()) { }

            internal static int Search(ReadOnlySpan<int> bin, double mass, List<PeptideWithSetModifications> peptideIndex)
                => BinarySearchBinForPrecursorIndex(bin, mass, peptideIndex);

            /// <summary>Reaches the protected instance scoring method; the base constructor only assigns fields.</summary>
            internal static List<int> Score(FragmentIndex fragmentIndex, List<int> binsToSearch, double highestMass,
                List<PeptideWithSetModifications> peptideIndex)
            {
                var observed = new List<int>();
                new BinSearchProbe().IndexedScoring(fragmentIndex, binsToSearch, new byte[peptideIndex.Count], 1,
                    observed, peptideIndex[0].MonoisotopicMass, double.NegativeInfinity, highestMass, peptideIndex,
                    new OpenSearchMode(), 0, DissociationType.HCD);
                return observed;
            }
        }

        /// lengths — the previous implementation was path-dependent, so the same target could give different
        /// answers for bins of different length, which is how database partitioning leaked into scoring.
        /// Exact ties between the target and a stored mass are the case that used to fail, and equal masses
        /// are pervasive because a reversed decoy carries its target's exact mass.
        /// </summary>
        [Test]
        public static void BinarySearchBinForPrecursorIndex_MatchesLinearScanAtEveryBinLength()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 5);
            var peptideIndex = MakeAnagramProteins()
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .Where(p => !double.IsNaN(p.MonoisotopicMass))
                .OrderBy(p => p.MonoisotopicMass)
                .ToList();

            var masses = peptideIndex.Select(p => p.MonoisotopicMass).ToList();
            var distinctMasses = masses.Distinct().OrderBy(m => m).ToList();
            Assert.That(masses.Count - distinctMasses.Count, Is.GreaterThan(0),
                "the interesting case is exact ties; this index must contain some");


            var targets = new List<double>();
            foreach (double m in distinctMasses)
            {
                targets.Add(m);              // exact tie - the case that used to fail
                targets.Add(m + 1e-6);
                targets.Add(m - 1e-6);
            }
            targets.Add(masses[0] - 100);    // below everything
            targets.Add(masses[^1] + 100);   // above everything

            var failures = new List<string>();
            foreach (int binLength in new[] { masses.Count, masses.Count / 2, masses.Count / 3, 101, 17, 8, 3, 1 })
            {
                if (binLength < 1) continue;
                var bin = Enumerable.Range(0, binLength).ToList();

                foreach (double target in targets)
                {
                    int actual = BinSearchProbe.Search(CollectionsMarshal.AsSpan(bin), target, peptideIndex);

                    int expected = 0;
                    for (int k = binLength - 1; k >= 0; k--)
                    {
                        if (masses[bin[k]] <= target) { expected = k; break; }
                    }

                    if (actual != expected && failures.Count < 10)
                    {
                        failures.Add($"binLength={binLength} target={target:F5} expected={expected} actual={actual}");
                    }
                }
            }

            Assert.That(failures, Is.Empty, "boundary disagrees with a linear scan: " + string.Join("; ", failures));
        }

        /// <summary>
        /// An empty bin must be skipped, not scored. Under the old list-per-bin layout an unpopulated bin was
        /// null and callers filtered on that; a span is never null, so `bin != null` silently became "always
        /// true" and empty bins reached the scoring loop, which indexed element 0 of nothing. Real crosslink
        /// and glyco searches died on this, so the check is on the observable behaviour, not just on no-throw.
        /// Both the finite and infinite upper bound are covered: they take different paths to the bin length.
        /// </summary>
        [Test]
        [TestCase(double.PositiveInfinity)]
        [TestCase(1e6)]
        public static void IndexedScoring_EmptyBinsAreSkippedNotScored(double highestMass)
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 5);
            var peptideIndex = MakeAnagramProteins()
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .Where(p => !double.IsNaN(p.MonoisotopicMass))
                .OrderBy(p => p.MonoisotopicMass)
                .Take(40)
                .ToList();

            // bins 0, 2 and 4 are empty; the populated ones hold ascending ids as the real index does
            var populated = new Dictionary<int, int[]> { [1] = new[] { 0, 1, 2 }, [3] = new[] { 3 } };
            var binStart = new List<int> { 0 };
            var peptideIds = new List<int>();
            for (int bin = 0; bin < 5; bin++)
            {
                if (populated.TryGetValue(bin, out int[] ids))
                {
                    peptideIds.AddRange(ids);
                }
                binStart.Add(peptideIds.Count);
            }
            var fragmentIndex = new FragmentIndex(binStart.ToArray(), peptideIds.ToArray());

            List<int> observed = BinSearchProbe.Score(fragmentIndex, new List<int> { 0, 1, 2, 3, 4 }, highestMass, peptideIndex);

            Assert.That(observed, Is.EquivalentTo(new[] { 0, 1, 2, 3 }),
                "exactly the ids in the populated bins, and nothing conjured out of the empty ones");

            // and searching only empty bins must be a no-op rather than a crash
            Assert.That(BinSearchProbe.Score(fragmentIndex, new List<int> { 0, 2, 4 }, highestMass, peptideIndex), Is.Empty);
        }

        /// <summary>
        /// Every ambiguity-joined psmtsv column reads BestMatchingBioPolymersWithSetMods, which sorts by a
        /// total order, so the printed order does not depend on which candidate reached AddOrReplace first.
        /// Three properties read the backing list instead and so printed in arrival order: Mass Diff (Da),
        /// Mass Diff (ppm), and the most-abundant-mode error. That made those columns move with the index
        /// partition count, and SpectralMatch.CompareTo tie-breaks on PrecursorMassErrorPpm.First(), so it
        /// reached psm ordering and the cumulative target/decoy counts too.
        ///
        /// The peptides here are added worst-score-first and best-score-first; both orders must print the
        /// same, and must agree with the sorted view.
        /// </summary>
        [Test]
        public static void PrecursorMassError_IsOrderedLikeEveryOtherAmbiguityColumn()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 5);
            var peptides = MakeAnagramProteins()
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .Where(p => !double.IsNaN(p.MonoisotopicMass))
                .GroupBy(p => p.MonoisotopicMass)
                .Select(g => g.First())
                .OrderBy(p => p.FullSequence, StringComparer.Ordinal)
                .Take(4)
                .ToList();
            Assert.That(peptides.Select(p => p.MonoisotopicMass).Distinct().Count(), Is.EqualTo(4),
                "distinct masses, so the mass-error list distinguishes the ordering");

            static List<double> ErrorsFor(List<PeptideWithSetModifications> inAdditionOrder)
            {
                var spectrum = new MzSpectrum(new double[,] { });
                var dataScan = new MsDataScan(spectrum, 1, 2, true, Polarity.Positive, 0, new MzRange(0, 2000), "",
                    MZAnalyzerType.Orbitrap, 0, null, null, "");
                var scan = new Ms2ScanWithSpecificMass(dataScan, 500, 1, "", new CommonParameters());

                SpectralMatch match = null;
                foreach (var peptide in inAdditionOrder)
                {
                    // one identical score for all of them, so every candidate is kept as ambiguous and the
                    // only thing that can differ between the two calls is insertion order
                    if (match == null)
                    {
                        match = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
                    }
                    else
                    {
                        match.AddOrReplace(peptide, 10, 0, true, new List<MatchedFragmentIon>());
                    }
                }

                Assert.That(match.NumDifferentMatchingPeptides, Is.EqualTo(inAdditionOrder.Count),
                    "all candidates must be retained, otherwise the ordering is not being exercised");
                return match.PrecursorMassErrorDa;
            }

            List<double> forward = ErrorsFor(peptides);
            List<double> reversed = ErrorsFor(Enumerable.Reverse(peptides).ToList());

            Assert.That(reversed, Is.EqualTo(forward), "insertion order must not reach the printed column");

            // and the order must be the sorted one the sibling columns use, not merely self-consistent
            var spectrumForSorted = new MzSpectrum(new double[,] { });
            var scanForSorted = new Ms2ScanWithSpecificMass(
                new MsDataScan(spectrumForSorted, 1, 2, true, Polarity.Positive, 0, new MzRange(0, 2000), "",
                    MZAnalyzerType.Orbitrap, 0, null, null, ""), 500, 1, "", new CommonParameters());
            SpectralMatch sorted = new PeptideSpectralMatch(peptides[0], 0, 10, 0, scanForSorted, new CommonParameters(), new List<MatchedFragmentIon>());
            foreach (var peptide in peptides.Skip(1))
            {
                sorted.AddOrReplace(peptide, 10, 0, true, new List<MatchedFragmentIon>());
            }
            var expected = sorted.BestMatchingBioPolymersWithSetMods
                .Select(p => Math.Round(sorted.ScanPrecursorMass - p.SpecificBioPolymer.MonoisotopicMass, 5))
                .ToList();
            Assert.That(forward, Is.EqualTo(expected));
        }

        // ------------------------------------------- flat index serialization

        private static List<int>[] RoundTrip(List<int>[] precursorIndex, out long fileLength)
        {
            var write = typeof(MetaMorpheusTask).GetMethod("WritePrecursorIndex", PrivateStatic);
            var read = typeof(MetaMorpheusTask).GetMethod("ReadPrecursorIndex", PrivateStatic);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "roundtrip-" + Guid.NewGuid().ToString("N") + ".bin");
            try
            {
                write.Invoke(null, new object[] { precursorIndex, path });
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
        public static void PrecursorIndex_RoundTripPreservesBinsContentsAndOrder()
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
        public static void PrecursorIndex_RoundTripBinLargerThanWriteBuffer()
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
        public static void PrecursorIndex_RoundTripEmptyIndex()
        {
            List<int>[] actual = RoundTrip(new List<int>[10], out long fileLength);

            Assert.That(actual.Length, Is.EqualTo(10));
            Assert.That(actual.All(b => b == null), Is.True);
            Assert.That(fileLength, Is.EqualTo((3 + 10) * sizeof(int)));
        }

        private static FragmentIndex RoundTripFragments(FragmentIndex index, out long fileLength)
        {
            var write = typeof(MetaMorpheusTask).GetMethod("WriteFragmentIndex", PrivateStatic);
            var read = typeof(MetaMorpheusTask).GetMethod("ReadFragmentIndex", PrivateStatic);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "roundtrip-frag-" + Guid.NewGuid().ToString("N") + ".bin");
            try
            {
                write.Invoke(null, new object[] { index, path });
                fileLength = new FileInfo(path).Length;
                return (FragmentIndex)read.Invoke(null, new object[] { path });
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// The search binary-searches inside a bin, so bin boundaries, contents and order all have to survive.
        /// The compressed layout puts every bin in one array, so an off-by-one in the offsets would shift a
        /// bin's contents into its neighbour rather than fail loudly — hence checking every bin, not a sample.
        /// </summary>
        [Test]
        public static void FragmentIndex_RoundTripPreservesBinBoundariesAndOrder()
        {
            // bins deliberately include the first, the last, empties either side of a populated one, and a
            // bin longer than the 16 Mi-int I/O chunk boundary is covered separately below
            var binStart = new List<int> { 0 };
            var peptideIds = new List<int>();
            var contents = new List<int>[6];
            contents[0] = new List<int> { 0 };
            contents[1] = new List<int>();
            contents[2] = new List<int> { 3, 3, 9, 41 };          // duplicates allowed
            contents[3] = new List<int>();
            contents[4] = Enumerable.Range(0, 1000).ToList();
            contents[5] = new List<int> { int.MaxValue, 0 };      // last bin, extreme value, descending
            foreach (var bin in contents)
            {
                peptideIds.AddRange(bin);
                binStart.Add(peptideIds.Count);
            }

            FragmentIndex actual = RoundTripFragments(new FragmentIndex(binStart.ToArray(), peptideIds.ToArray()), out long fileLength);

            Assert.That(actual.Length, Is.EqualTo(6));
            Assert.That(actual.EntryCount, Is.EqualTo(peptideIds.Count));
            for (int i = 0; i < contents.Length; i++)
            {
                Assert.That(actual[i].ToArray(), Is.EqualTo(contents[i].ToArray()), $"bin {i}");
            }
            // 4 header ints + 7 offsets + 1007 entries
            Assert.That(fileLength, Is.EqualTo((4 + 7 + 1007) * sizeof(int)));
        }

        /// <summary>
        /// Both arrays are written through a fixed int chunk, so an index whose entries exceed one chunk
        /// exercises the loop rather than the single-shot path. Real indexes are hundreds of millions of
        /// entries, so the chunked path is the one that actually runs in production.
        /// </summary>
        [Test]
        public static void FragmentIndex_RoundTripSurvivesMoreEntriesThanOneIoChunk()
        {
            const int entries = (16 * 1024 * 1024) + 12345; // one chunk plus a remainder
            var peptideIds = new int[entries];
            for (int i = 0; i < entries; i++)
            {
                peptideIds[i] = entries - i; // descending, so order is checkable
            }

            FragmentIndex actual = RoundTripFragments(new FragmentIndex(new[] { 0, 1, entries }, peptideIds), out long fileLength);

            Assert.That(actual.Length, Is.EqualTo(2));
            Assert.That(actual[0].Length, Is.EqualTo(1));
            Assert.That(actual[1].Length, Is.EqualTo(entries - 1));
            Assert.That(actual[1][0], Is.EqualTo(entries - 1), "first entry of the long bin, and order preserved");
            Assert.That(actual[1][entries - 2], Is.EqualTo(1), "last entry");
            Assert.That(fileLength, Is.EqualTo((4L + 3 + entries) * sizeof(int)));
        }

        /// <summary>An index with no fragments at all still has to round-trip.</summary>
        [Test]
        public static void FragmentIndex_RoundTripEmptyIndex()
        {
            FragmentIndex actual = RoundTripFragments(new FragmentIndex(new int[11], Array.Empty<int>()), out long fileLength);

            Assert.That(actual.Length, Is.EqualTo(10));
            Assert.That(actual.EntryCount, Is.Zero);
            Assert.That(Enumerable.Range(0, 10).All(i => actual[i].IsEmpty), Is.True);
            Assert.That(fileLength, Is.EqualTo((4 + 11) * sizeof(int)));
        }

        /// <summary>
        /// An index written by an older MetaMorpheus must be rejected, not misinterpreted. The reader is the
        /// last line of defence: GenerateSecondIndexes reads without a try/catch. The stale-version case is
        /// the one that actually happens — the fragment index kept its magic number when the layout changed
        /// from bin counts to offsets, so only the version distinguishes them.
        /// </summary>
        [Test]
        [TestCase(new byte[] { 9, 9, 9, 9, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 }, TestName = "foreign magic")]
        [TestCase(new byte[] { 0x4D, 0x4D, 0x46, 0x49, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 }, TestName = "previous format version")]
        [TestCase(new byte[] { 0x4D, 0x4D, 0x50, 0x49, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 }, TestName = "precursor index in the fragment index slot")]
        public static void FragmentIndex_ReadRejectsForeignFile(byte[] header)
        {
            var read = typeof(MetaMorpheusTask).GetMethod("ReadFragmentIndex", PrivateStatic);
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "foreign-" + Guid.NewGuid().ToString("N") + ".bin");
            File.WriteAllBytes(path, header);
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

            public CommonParameters Probe(List<Protein> proteins, CommonParameters parameters, ref int? decided)
                => RaisePartitionsToFitMemory(proteins, parameters, new List<Modification>(),
                    new List<Modification>(), null, null, null, 30000.0, ref decided);
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
            int? decided = null;

            var warnings = new List<string>();
            void Collect(object sender, EngineLayer.StringEventArgs e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += Collect;
            try
            {
                CommonParameters result = task.Probe(MakeProteins(20), parameters, ref decided);

                Assert.That(ReferenceEquals(result, parameters), Is.True, "must not copy when no change is needed");
                Assert.That(decided, Is.EqualTo(1), "decision recorded so later files reuse it");
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
            int? decided = null;

            var warnings = new List<string>();
            void Collect(object sender, EngineLayer.StringEventArgs e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += Collect;
            try
            {
                CommonParameters first = task.Probe(huge, parameters, ref decided);
                CommonParameters second = task.Probe(huge, parameters, ref decided);

                Assert.That(first.TotalPartitions, Is.GreaterThan(1), "an impossible database must be partitioned");
                Assert.That(ReferenceEquals(first, parameters), Is.False, "must return a copy, not mutate the input");
                Assert.That(parameters.TotalPartitions, Is.EqualTo(1), "the caller's parameters must be untouched");
                Assert.That(second.TotalPartitions, Is.EqualTo(first.TotalPartitions), "decision must be stable");
                Assert.That(decided, Is.EqualTo(first.TotalPartitions), "decision cached for later files");
                Assert.That(warnings.Count, Is.EqualTo(1), "must warn once, not once per spectra file");
                Assert.That(warnings[0], Does.Contain("TotalPartitions"));
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= Collect;
            }
        }

        /// <summary>
        /// The guard runs once per spectra file, but the decision must be taken once per task. Available memory
        /// shrinks as PSMs accumulate, so re-deriving per file could index file 1 in one partition and file 5 in
        /// four — which invalidates the disk cache for every partition, since the count is part of the cache
        /// key, and leaves files in one run searched under different partitionings.
        /// </summary>
        [Test]
        public static void RaisePartitionsToFitMemory_DecisionAlreadyMade_IsReusedWithoutReDeriving()
        {
            var task = new GuardProbeTask();
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin"));
            // as though an earlier spectra file had already settled on 7, on a database that on its own
            // would need only 1
            int? decided = 7;

            var warnings = new List<string>();
            void Collect(object sender, EngineLayer.StringEventArgs e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += Collect;
            try
            {
                CommonParameters result = task.Probe(MakeProteins(20), parameters, ref decided);

                Assert.That(decided, Is.EqualTo(7), "an existing decision must not be re-derived");
                Assert.That(result.TotalPartitions, Is.EqualTo(7), "later files must reuse the first file's count");
                Assert.That(warnings, Is.Empty, "a decision already reported must not warn again");
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= Collect;
            }
        }

        /// <summary>
        /// If even the maximum partition count will not fit, say so rather than leaving the user to wonder why
        /// the search is paging — the first warning alone would imply the problem had been solved. Driven
        /// through the pure message builder with the estimate and budget supplied, because asserting an exact
        /// partition count from an end-to-end call would depend on how much memory the runner happens to have.
        /// </summary>
        [Test]
        public static void PartitionWarnings_ReportTheRaiseAndWhetherItStillWillNotFit()
        {
            const long tenGb = 10L * 1024 * 1024 * 1024;

            var fits = IndexPartitioning.PartitionWarnings(1, 4, 22L * 1024 * 1024 * 1024, tenGb, cappedByLimit: false)
                .Where(w => w.Contains("TotalPartitions") || w.Contains("may page heavily")).ToList();
            Assert.That(fits.Count, Is.EqualTo(1), "a raise that fits needs one message");
            Assert.That(fits[0], Does.Contain("TotalPartitions to 4"));
            Assert.That(fits[0], Does.Not.Contain("may page heavily"));

            var capped = IndexPartitioning.PartitionWarnings(1, 256, long.MaxValue / 2, tenGb, cappedByLimit: true)
                .Where(w => w.Contains("TotalPartitions") || w.Contains("may page heavily")).ToList();
            Assert.That(capped.Count, Is.EqualTo(2), "the raise, plus a warning that it still will not fit");
            Assert.That(capped[0], Does.Contain("TotalPartitions to 256"));
            Assert.That(capped[1], Does.Contain("may page heavily"));
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
