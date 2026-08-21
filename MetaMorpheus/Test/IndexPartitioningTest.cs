using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.CrosslinkSearch;
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
                Assert.That(estimate, Is.GreaterThan(0));

                Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", "512");
                Suggest(proteins, parameters, 1, out _, out long narrowBudget, out _);
                Assert.That(narrowBudget, Is.LessThan(wideBudget), "the allocation must shrink the budget");
                Assert.That(narrowBudget, Is.LessThanOrEqualTo(512L * 1024 * 1024), "budget cannot exceed the allocation");

                // Tightening the allocation must never ask for fewer partitions, and at some point must ask
                // for more. Which allocation crosses that line depends on the byte-cost constants and on how
                // much memory the machine has; hard-coding one couples the test to both, which is how this
                // test broke when the fragment-entry cost was corrected. Walk them instead.
                long[] allocationsMb = { 65536, 16384, 4096, 1024, 512, 384, 320, 288, 272, 264, 224, 200 };
                int previous = unconstrained;
                int everRaised = 0;

                foreach (long mb in allocationsMb)
                {
                    Environment.SetEnvironmentVariable("SLURM_MEM_PER_NODE", mb.ToString());
                    int here = Suggest(proteins, parameters, 1, out _, out _, out _);

                    Assert.That(here, Is.GreaterThanOrEqualTo(previous),
                        $"a tighter allocation ({mb} MB) asked for fewer partitions than a looser one");
                    previous = here;

                    if (here > unconstrained)
                    {
                        everRaised = here;
                    }
                }

                Assert.That(everRaised, Is.GreaterThan(unconstrained),
                    "some Slurm allocation in this range must force more partitions than an unconstrained machine");
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

            // bin space alone exceeding the budget is a measurement, not a missing one: split as far as the
            // protein count allows and report that it still will not fit, rather than returning the request
            Assert.That(IndexPartitioning.PartitionsForBudget(budget, budget, tenGb, 3, 5000, out bool capped), Is.EqualTo(256));
            Assert.That(capped, Is.True);

            // and never more partitions than there are proteins to put in them
            Assert.That(IndexPartitioning.PartitionsForBudget(budget, budget, tenGb, 1, 9, out _), Is.EqualTo(9));

            // an unmeasurable machine is still a no-op, which is a different thing entirely
            Assert.That(IndexPartitioning.PartitionsForBudget(budget, budget, 0, 3, 5000, out bool unmeasured), Is.EqualTo(3));
            Assert.That(unmeasured, Is.False);
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

        /// <summary>
        /// Digestion is parallel and each thread keeps its own list, so the peptide index is concatenated in
        /// thread order and its ids depend on the thread count. Reproducibility across thread counts is the
        /// property that has to hold, not a particular id order - asserting the order would pin an arbitrary
        /// choice. What must not change is the answer: the same scan searched against an index built with one
        /// thread and with four must produce the same match with the same score.
        ///
        /// This is the thread-count analogue of the partition-count invariance the rest of this branch is
        /// about, and it is checked end to end for the same reason: the peptide ids genuinely do differ
        /// between the two, so only the output can show that it does not matter.
        /// </summary>
        [Test]
        public static void SearchResultIsIndependentOfTheDigestionThreadCount()
        {
            var (scan, proteins) = SyntheticScanFor("MKPEPTIDERTIDEK", "MKAAAAAAAAKGGGGGGGGKSSSSSSSSK");

            static (string Sequence, double Score, string Ids) SearchWith(int threads, List<Protein> proteins, Ms2ScanWithSpecificMass scan)
            {
                var parameters = new CommonParameters(maxThreadsToUsePerFile: threads,
                    digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
                var fsp = new List<(string, CommonParameters)> { ("", parameters) };

                var index = (IndexingResults)new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(),
                    null, null, null, 0, DecoyType.Reverse, parameters, fsp, 30000, false, new List<FileInfo>(),
                    TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();

                var matches = new SpectralMatch[1];
                new ModernSearchEngine(matches, new[] { scan }, index.PeptideIndex, index.FragmentIndex, 0, parameters,
                    fsp, new OpenSearchMode(), 0, new List<string>()).Run();
                matches[0]?.ResolveAllAmbiguities();

                string ids = string.Join(",", index.PeptideIndex.Select(p => p.FullSequence));
                return (matches[0]?.BaseSequence, matches[0]?.Score ?? 0, ids);
            }

            var one = SearchWith(1, proteins, scan);
            var four = SearchWith(4, proteins, scan);

            Assert.That(one.Sequence, Is.Not.Null, "the fixture must produce a match at all");
            Assert.That(four.Sequence, Is.EqualTo(one.Sequence), "same peptide");
            Assert.That(four.Score, Is.EqualTo(one.Score).Within(1e-9), "same score");
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

            private BinSearchProbe(List<PeptideWithSetModifications> peptideIndex)
                : base(null, null, peptideIndex, null, 0, new CommonParameters(), null, new OpenSearchMode(), 0, new List<string>()) { }

            private BinSearchProbe(CommonParameters parameters)
                : base(null, null, null, null, 0, parameters, null, new OpenSearchMode(), 0, new List<string>()) { }

            internal static int Search(ReadOnlySpan<int> bin, double mass, List<PeptideWithSetModifications> peptideIndex)
                => BinarySearchBinForPrecursorIndex(bin, mass, peptideIndex);

            internal static int FirstAtOrAbove(ReadOnlySpan<int> bin, double mass, List<PeptideWithSetModifications> peptideIndex)
                => BinarySearchBinForFirstAtOrAbove(bin, mass, peptideIndex);

            internal static List<int> BinsToSearch(Ms2ScanWithSpecificMass scan, FragmentIndex fragmentIndex, CommonParameters parameters)
                => new BinSearchProbe(parameters).GetBinsToSearch(scan, fragmentIndex, parameters.DissociationType);

            /// <summary>The window is an instance method because it reads PeptideIndex off the engine.</summary>
            internal static (int start, int end) Window(double lowest, double highest, ReadOnlySpan<int> bin,
                List<PeptideWithSetModifications> peptideIndex)
                => new BinSearchProbe(peptideIndex).GetFirstAndLastIndexesInBinToIncrement(lowest, highest, bin, 0);

            /// <summary>Reaches the protected instance scoring method; the base constructor only assigns fields.</summary>
            internal static List<int> Score(FragmentIndex fragmentIndex, List<int> binsToSearch, double highestMass,
                List<PeptideWithSetModifications> peptideIndex)
                => Score(fragmentIndex, binsToSearch, double.NegativeInfinity, highestMass, peptideIndex);

            internal static List<int> Score(FragmentIndex fragmentIndex, List<int> binsToSearch, double lowestMass,
                double highestMass, List<PeptideWithSetModifications> peptideIndex)
            {
                var observed = new List<int>();
                new BinSearchProbe().IndexedScoring(fragmentIndex, binsToSearch, new byte[peptideIndex.Count], 1,
                    observed, peptideIndex.First(p => !double.IsNaN(p.MonoisotopicMass)).MonoisotopicMass,
                    lowestMass, highestMass, peptideIndex, new OpenSearchMode(), 0, DissociationType.HCD);
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

                    int expected = -1;
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

        /// <summary>
        /// The search binary-searches a bin by peptide mass, which needs the bin to be monotone in whichever
        /// predicate is being searched. A peptide with an undefined mass breaks that if taken literally: every
        /// comparison against NaN is false, so it reads as "too heavy" to one search and "too light" to the
        /// other, wherever it sits. Sorting the peptide index by mass puts NaNs at the front of every bin they
        /// occupy, and the searches read them as negative infinity, which makes both predicates monotone.
        ///
        /// The property that matters is that the window never excludes an entry whose mass is inside it, and
        /// never includes an undefined one. Such peptides stay in the index: a tolerance acceptor can never
        /// match one, but OpenSearchMode accepts anything and an open search leaves both bounds infinite, so
        /// neither search runs and the whole bin is scored. Dropping them would silently stop open and glyco
        /// searches reporting peptides with an unknown residue. The database here includes proteins with X
        /// residues, so the index really does contain undefined masses.
        /// </summary>
        [Test]
        public static void FragmentIndex_UndefinedMassPeptidesCannotClipTheSearchWindow()
        {
            const double maxFragmentSize = 3000;
            var parameters = new CommonParameters(scoreCutoff: 1,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            var engine = new IndexingEngine(MakeAnagramProteins(), new List<Modification>(), new List<Modification>(),
                null, null, null, 0, DecoyType.Reverse, parameters, fsp, maxFragmentSize, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var results = (IndexingResults)engine.Run();

            int undefined = results.PeptideIndex.Count(p => double.IsNaN(p.MonoisotopicMass));
            Assert.That(undefined, Is.GreaterThan(0), "the fixture must actually contain undefined masses");

            var indexed = new HashSet<int>();
            for (int bin = 0; bin < results.FragmentIndex.Length; bin++)
            {
                foreach (int id in results.FragmentIndex[bin])
                {
                    indexed.Add(id);
                }
            }
            Assert.That(indexed.Any(id => double.IsNaN(results.PeptideIndex[id].MonoisotopicMass)), Is.True,
                "undefined-mass peptides must still be indexed, or an open search could never find them");

            // the invariant the search depends on, checked against a linear scan over real windows
            var masses = results.PeptideIndex.Select(p => p.MonoisotopicMass).ToList();
            var clipped = new List<string>();
            int windowsChecked = 0;

            for (int bin = 0; bin < results.FragmentIndex.Length && clipped.Count < 5; bin++)
            {
                ReadOnlySpan<int> entries = results.FragmentIndex[bin];
                if (entries.Length < 3)
                {
                    continue;
                }

                foreach (int anchorIndex in new[] { 0, entries.Length / 2, entries.Length - 1 })
                {
                    double centre = masses[entries[anchorIndex]];
                    if (double.IsNaN(centre))
                    {
                        continue;
                    }

                    double lo = centre - 1.5;
                    double hi = centre + 1.5;
                    windowsChecked++;

                    var (start, end) = BinSearchProbe.Window(lo, hi, entries, results.PeptideIndex);

                    for (int j = 0; j < entries.Length; j++)
                    {
                        double mass = masses[entries[j]];
                        if (mass >= lo && mass <= hi && (j < start || j > end))
                        {
                            clipped.Add($"bin {bin} entry {j} (mass {mass:F4}) is inside [{lo:F4},{hi:F4}] but outside [{start},{end}]");
                            break;
                        }

                        // and an undefined mass must never be counted as inside a finite window
                        if (double.IsNaN(mass) && j >= start && j <= end)
                        {
                            clipped.Add($"bin {bin} entry {j} has an undefined mass but is inside [{start},{end}]");
                            break;
                        }
                    }
                }
            }

            Assert.That(windowsChecked, Is.GreaterThan(100), "need a decent number of windows for this to mean anything");
            Assert.That(clipped, Is.Empty, string.Join("; ", clipped));
        }

        /// <summary>
        /// The two ends of the window need different searches. The end is an upper bound, the start a lower
        /// bound, and taking the start from the upper-bound search answers the wrong question: it gives the
        /// LAST entry at or below the bound, so a run of equal masses sitting exactly on the window's lower
        /// edge collapses to its final member and the rest score nothing.
        ///
        /// Equal masses are the normal case, not a corner case: the same sequence occurs in several proteins
        /// and a reversed decoy carries its target's mass, and notch intervals are derived from a peptide
        /// mass so the edge lands on such a run routinely. The fixture is built from anagram peptides for
        /// exactly this reason.
        /// </summary>
        [Test]
        public static void GetFirstAndLastIndexesInBinToIncrement_KeepsEveryEqualMassOnTheLowerEdge()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 5);
            var peptideIndex = MakeAnagramProteins()
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .Where(p => !double.IsNaN(p.MonoisotopicMass))
                .OrderBy(p => p.MonoisotopicMass)
                .ToList();
            var masses = peptideIndex.Select(p => p.MonoisotopicMass).ToList();

            // the longest run of exactly equal masses
            int runStart = 0, runLength = 1, bestStart = 0, bestLength = 1;
            for (int i = 1; i < masses.Count; i++)
            {
                if (masses[i] == masses[i - 1]) { runLength++; }
                else { runStart = i; runLength = 1; }
                if (runLength > bestLength) { bestLength = runLength; bestStart = runStart; }
            }
            Assert.That(bestLength, Is.GreaterThan(2), "the fixture must contain a run of equal masses");

            var bin = Enumerable.Range(0, masses.Count).ToList();
            double edge = masses[bestStart];

            var (start, end) = BinSearchProbe.Window(edge, masses[^1] + 1, CollectionsMarshal.AsSpan(bin), peptideIndex);

            Assert.That(start, Is.EqualTo(bestStart),
                $"the window starts at the first of the {bestLength} entries with mass {edge:F5}, not the last");
            Assert.That(end, Is.GreaterThanOrEqualTo(bestStart + bestLength - 1), "and covers the whole run");

            // stated directly on the search itself
            Assert.That(BinSearchProbe.FirstAtOrAbove(CollectionsMarshal.AsSpan(bin), edge, peptideIndex), Is.EqualTo(bestStart));
            Assert.That(BinSearchProbe.Search(CollectionsMarshal.AsSpan(bin), edge, peptideIndex),
                Is.EqualTo(bestStart + bestLength - 1), "the upper bound still answers with the last of the run");
        }

        /// <summary>
        /// A bin whose peptides are all heavier than the window must contribute nothing. Returning 0 for both
        /// "nothing found" and "index 0 is the answer" made it score its first entry instead.
        /// </summary>
        [Test]
        public static void GetFirstAndLastIndexesInBinToIncrement_WindowBelowTheWholeBinIsEmpty()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 5);
            var peptideIndex = MakeAnagramProteins()
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .Where(p => !double.IsNaN(p.MonoisotopicMass))
                .OrderBy(p => p.MonoisotopicMass)
                .Take(10)
                .ToList();
            var bin = Enumerable.Range(0, peptideIndex.Count).ToList();
            double lightest = peptideIndex[0].MonoisotopicMass;

            var (start, end) = BinSearchProbe.Window(lightest - 200, lightest - 100, CollectionsMarshal.AsSpan(bin), peptideIndex);
            Assert.That(end, Is.LessThan(start), "an empty range, so the scoring loop does not execute");

            (start, end) = BinSearchProbe.Window(lightest - 200, peptideIndex[^1].MonoisotopicMass + 100, CollectionsMarshal.AsSpan(bin), peptideIndex);
            Assert.That(start, Is.Zero);
            Assert.That(end, Is.EqualTo(bin.Count - 1));
        }

        /// <summary>
        /// The bin an indexed fragment lands in is derived in one place, shared by the counting and the
        /// filling pass. Low-resolution CID rounds the mass to the nearest 1.0005079 first, and fragments
        /// outside the index are dropped. Both are exercised here because both passes have to agree: if only
        /// one applied the rounding the two passes would disagree about how many entries a bin holds, and the
        /// fill would write past its slice.
        /// </summary>
        [Test]
        [TestCase(DissociationType.HCD)]
        [TestCase(DissociationType.LowCID)]
        public static void IndexingEngine_FragmentBinningAgreesBetweenBothPasses(DissociationType dissociationType)
        {
            // small enough that a good number of fragments fall outside it and take the discard branch
            const double maxFragmentSize = 900;
            var parameters = new CommonParameters(dissociationType: dissociationType, scoreCutoff: 1,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            var engine = new IndexingEngine(MakeAnagramProteins(), new List<Modification>(), new List<Modification>(),
                null, null, null, 0, DecoyType.Reverse, parameters, fsp, maxFragmentSize, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var results = (IndexingResults)engine.Run();

            // the count pass sized the arrays and the fill pass filled them; a disagreement shows up as a
            // slice that was not completely written, which reads as a stray zero where an id should be
            Assert.That(results.FragmentIndex.EntryCount, Is.GreaterThan(1000), "the fixture must populate the index");
            Assert.That(results.FragmentIndex.BinStart[^1], Is.EqualTo(results.FragmentIndex.EntryCount),
                "the last offset must be the total, so every counted entry has a slot");

            int occupied = 0;
            for (int bin = 0; bin < results.FragmentIndex.Length; bin++)
            {
                ReadOnlySpan<int> entries = results.FragmentIndex[bin];
                if (entries.IsEmpty)
                {
                    continue;
                }
                occupied++;
                for (int i = 1; i < entries.Length; i++)
                {
                    Assert.That(entries[i], Is.GreaterThanOrEqualTo(entries[i - 1]), $"bin {bin} is not ascending");
                }
            }
            Assert.That(occupied, Is.GreaterThan(100));

            // and nothing was indexed above the maximum fragment mass
            Assert.That(results.FragmentIndex.Length, Is.EqualTo((int)Math.Ceiling(maxFragmentSize) * 1000 + 1));
        }

        /// <summary>The index type must reject arrays it cannot describe rather than fail later at a lookup.</summary>
        [Test]
        public static void FragmentIndex_RejectsArraysItCannotDescribe()
        {
            Assert.That(() => new FragmentIndex(null, new int[1]), Throws.ArgumentNullException);
            Assert.That(() => new FragmentIndex(new int[2], null), Throws.ArgumentNullException);
            Assert.That(() => new FragmentIndex(Array.Empty<int>(), Array.Empty<int>()),
                Throws.ArgumentException, "no offsets at all cannot describe even an empty index");

            // the smallest well-formed index: one offset, so zero bins
            var empty = new FragmentIndex(new int[1], Array.Empty<int>());
            Assert.That(empty.Length, Is.Zero);
            Assert.That(empty.EntryCount, Is.Zero);
        }

        /// <summary>A database whose proteins have no residues must not divide by zero while sampling.</summary>
        [Test]
        public static void SuggestTotalPartitions_ProteinsWithNoResidues_ReturnsRequestedCount()
        {
            var parameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin"));
            var empty = Enumerable.Range(0, 10).Select(i => new Protein("", "EMPTY" + i)).ToList();

            Assert.That(Suggest(empty, parameters, 2, out long estimate, out _, out bool capped), Is.EqualTo(2));
            Assert.That(estimate, Is.Zero);
            Assert.That(capped, Is.False);
        }

        // ------------------------------------------------------ end-to-end modern search

        /// <summary>
        /// Builds a scan out of a peptide's own theoretical fragments, so the expected answer is known
        /// without a fixture file. Returns the scan and the peptide index it should be searched against.
        /// </summary>
        private static (Ms2ScanWithSpecificMass Scan, List<Protein> Proteins) SyntheticScanFor(string target, string decoyFiller,
            CommonParameters parameters = null)
        {
            parameters ??= new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            var protein = new Protein(target, "TARGET");
            var peptide = protein.Digest(parameters.DigestionParams, new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>().First();

            var products = new List<Product>();
            peptide.Fragment(parameters.DissociationType, FragmentationTerminus.Both, products, parameters.FragmentationParameters);
            var neutral = products.Select(p => p.NeutralMass).Where(m => m > 0).Distinct().OrderBy(m => m).ToList();
            Assert.That(neutral.Count, Is.GreaterThan(6), "the target must produce enough fragments to clear the score cutoff");

            double[] mz = neutral.Select(m => m.ToMz(1)).ToArray();
            double[] intensity = neutral.Select(_ => 1000.0).ToArray();
            var dataScan = new MsDataScan(new MzSpectrum(mz, intensity, false), 1, 2, true, Polarity.Positive, 1,
                new MzRange(0, 5000), "", MZAnalyzerType.Orbitrap, intensity.Sum(), null, null, "");

            var envelopes = neutral
                .Select(m => new IsotopicEnvelope(new List<(double, double)> { (m.ToMz(1), 1000.0) }, m, 1, 1000.0, 0))
                .ToArray();

            var scan = new Ms2ScanWithSpecificMass(dataScan, peptide.MonoisotopicMass.ToMz(1), 1, "", parameters, envelopes);
            return (scan, new List<Protein> { protein, new Protein(decoyFiller, "FILLER") });
        }

        /// <summary>
        /// A whole modern search over a synthetic scan. This is the only test here that walks RunSpecific,
        /// IndexScoreScan, IncrementPeptideScoresInBin and FineScorePeptide, which the span rewrite touched
        /// in every one of them and which no unit test reached before.
        /// </summary>
        [Test]
        [TestCase(DissociationType.HCD, false, TestName = "HCD")]
        [TestCase(DissociationType.HCD, true, TestName = "HCD with complementary ions")]
        public static void ModernSearchEngine_FindsThePeptideItsOwnFragmentsCameFrom(DissociationType dissociationType, bool addCompIons)
        {
            // Complementary ions add a second bin window per peak, which the span rewrite touched separately.
            // Low-resolution CID is deliberately not asserted here: it does not find the peptide, because the
            // indexing and scoring sides round the bin differently - Math.Round against a cast to int - and
            // disagree on exactly half of all bins. That predates this branch and is not fixed here, so
            // asserting either outcome would be asserting current behaviour rather than correct behaviour.
            var parameters = new CommonParameters(dissociationType: dissociationType, addCompIons: addCompIons,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            var (scan, proteins) = SyntheticScanFor("MKPEPTIDERTIDEK", "MKAAAAAAAAKGGGGGGGGKSSSSSSSSK", parameters);
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            var indexEngine = new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(), null, null, null,
                0, DecoyType.Reverse, parameters, fsp, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var index = (IndexingResults)indexEngine.Run();

            var matches = new SpectralMatch[1];
            new ModernSearchEngine(matches, new[] { scan }, index.PeptideIndex, index.FragmentIndex, 0, parameters,
                fsp, new OpenSearchMode(), 0, new List<string>()).Run();

            Assert.That(matches[0], Is.Not.Null, "the peptide its own fragments came from must be found");
            matches[0].ResolveAllAmbiguities();
            Assert.That(matches[0].BaseSequence, Is.EqualTo("PEPTIDER"));
            Assert.That(matches[0].Score, Is.GreaterThan(5));
        }

        /// <summary>
        /// Fine scoring used to stop early once no remaining candidate's coarse score could beat the best
        /// fine score, which made the runner-up depend on how the database had been split. Every candidate
        /// is now scored, so a second, worse candidate still reaches the psm and sets RunnerUpScore. An
        /// unset RunnerUpScore reads as the score cutoff, so that is what this distinguishes.
        /// </summary>
        [Test]
        public static void FineScorePeptides_ScoresTheRunnerUpAndNotJustTheWinner()
        {
            // the filler shares a tryptic peptide mass region with the target, so it is coarse-scored too
            var (scan, proteins) = SyntheticScanFor("MKPEPTIDERTIDEK", "MKPEPTLDERSSSSSSSSK");
            var parameters = new CommonParameters(scoreCutoff: 3,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            var indexEngine = new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(), null, null, null,
                0, DecoyType.Reverse, parameters, fsp, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var index = (IndexingResults)indexEngine.Run();

            var matches = new SpectralMatch[1];
            new ModernSearchEngine(matches, new[] { scan }, index.PeptideIndex, index.FragmentIndex, 0, parameters,
                fsp, new OpenSearchMode(), 0, new List<string>()).Run();

            Assert.That(matches[0], Is.Not.Null);
            Assert.That(matches[0].RunnerUpScore, Is.GreaterThan(parameters.ScoreCutoff),
                "a runner-up was scored, rather than the loop stopping at the winner and leaving the cutoff behind");
            Assert.That(matches[0].DeltaScore, Is.LessThan(matches[0].Score),
                "delta score is measured against a real runner-up");
        }

        /// <summary>
        /// GetBinsToSearch must return only populated bins. Under the old list-per-bin layout an unpopulated
        /// bin was null and this filtered on that; a span is a struct and is never null, so the same test
        /// compiled and silently became "always true". Every empty bin in the window then reached
        /// IndexedScoring, which indexed element 0 of nothing — crosslink and glyco searches died on it.
        /// </summary>
        [Test]
        [TestCase(DissociationType.HCD, false, TestName = "HCD")]
        [TestCase(DissociationType.HCD, true, TestName = "HCD with complementary ions")]
        [TestCase(DissociationType.LowCID, false, TestName = "low-resolution CID")]
        public static void GetBinsToSearch_ReturnsOnlyPopulatedBins(DissociationType dissociationType, bool addCompIons)
        {
            var parameters = new CommonParameters(dissociationType: dissociationType, addCompIons: addCompIons,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 5));
            var (scan, proteins) = SyntheticScanFor("MKPEPTIDERTIDEK", "MKAAAAAAAAKGGGGGGGGKSSSSSSSSK", parameters);
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            var indexEngine = new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(), null, null, null,
                0, DecoyType.Reverse, parameters, fsp, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var index = (IndexingResults)indexEngine.Run();

            List<int> bins = BinSearchProbe.BinsToSearch(scan, index.FragmentIndex, parameters);

            Assert.That(bins, Is.Not.Empty, "the scan's own fragments must land in populated bins");
            Assert.That(bins.All(b => !index.FragmentIndex[b].IsEmpty), Is.True,
                "an empty bin must never be returned; downstream scoring indexes into it without checking");
        }

        /// <summary>
        /// An open search accepts any precursor mass, so a peptide whose mass is undefined is a legitimate
        /// answer. This is the case that regressed when such peptides were dropped from the index, and it is
        /// why the binary-search problem is fixed in the searches rather than by removing them.
        /// </summary>
        [Test]
        public static void FragmentIndex_UndefinedMassPeptideIsStillReachable()
        {
            var parameters = new CommonParameters(scoreCutoff: 1,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            // trypsin cuts this into MNNNK and QXQ; QXQ has an undefined mass
            var proteins = new List<Protein> { new Protein("MNNNKQXQ", "X-CONTAINING") };
            var engine = new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(), null, null, null,
                1, DecoyType.Reverse, parameters, fsp, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var results = (IndexingResults)engine.Run();

            Assert.That(results.PeptideIndex.Any(p => double.IsNaN(p.MonoisotopicMass)), Is.True,
                "the fixture must produce a peptide with an undefined mass");

            var indexed = new HashSet<int>();
            for (int bin = 0; bin < results.FragmentIndex.Length; bin++)
            {
                foreach (int id in results.FragmentIndex[bin])
                {
                    indexed.Add(id);
                }
            }

            Assert.That(indexed.Any(id => double.IsNaN(results.PeptideIndex[id].MonoisotopicMass)), Is.True,
                "it must be reachable through the fragment index, or an open search cannot score it");
        }

        /// <summary>
        /// Both bin searches read an undefined mass as negative infinity, which is only monotone if the sort
        /// actually puts those peptides first. It does, because Array.Sort here compares with double.CompareTo
        /// and CompareTo orders NaN below everything — unlike the &lt; operator, which is false against NaN in
        /// both directions. Swapping the comparison back to operators, or to Array.Sort over a double[] key
        /// array, would move them and silently break both searches, so the invariant is pinned here.
        /// </summary>
        [Test]
        public static void SortByMonoisotopicMass_PutsUndefinedMassesFirst()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 1);
            var peptides = MakeAnagramProteins()
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .ToList();

            int undefined = peptides.Count(p => double.IsNaN(p.MonoisotopicMass));
            Assert.That(undefined, Is.GreaterThan(0), "the fixture must contain undefined masses");

            var sorted = IndexingEngine.SortByMonoisotopicMass(peptides);

            Assert.That(sorted.Take(undefined).All(p => double.IsNaN(p.MonoisotopicMass)), Is.True,
                "every undefined mass must come first");
            Assert.That(sorted.Skip(undefined).Any(p => double.IsNaN(p.MonoisotopicMass)), Is.False,
                "and none may appear after a defined one");

            var defined = sorted.Skip(undefined).Select(p => p.MonoisotopicMass).ToList();
            Assert.That(defined, Is.Ordered, "the rest ascend, which is what the searches binary-search over");
        }

        /// <summary>
        /// The precursor index bins on the peptide's own mass, so an undefined one has no bin to go in. It is
        /// already skipped; this pins that, because the same peptide is deliberately kept in the fragment
        /// index and the two indexes have opposite answers for it.
        /// </summary>
        [Test]
        public static void CreateNewPrecursorIndex_SkipsUndefinedMasses()
        {
            var parameters = new CommonParameters(scoreCutoff: 1,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };

            var engine = new IndexingEngine(new List<Protein> { new Protein("MNNNKQXQ", "X-CONTAINING") },
                new List<Modification>(), new List<Modification>(), null, null, null, 0, DecoyType.None,
                parameters, fsp, 30000, true, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var results = (IndexingResults)engine.Run();

            var undefinedIds = Enumerable.Range(0, results.PeptideIndex.Count)
                .Where(i => double.IsNaN(results.PeptideIndex[i].MonoisotopicMass))
                .ToHashSet();
            Assert.That(undefinedIds, Is.Not.Empty, "the fixture must produce an undefined mass");

            Assert.That(results.PrecursorIndex, Is.Not.Null);
            var inPrecursorIndex = results.PrecursorIndex.Where(b => b != null).SelectMany(b => b).ToHashSet();
            Assert.That(inPrecursorIndex.Overlaps(undefinedIds), Is.False,
                "an undefined mass has no precursor bin");

            var inFragmentIndex = new HashSet<int>();
            for (int bin = 0; bin < results.FragmentIndex.Length; bin++)
            {
                foreach (int id in results.FragmentIndex[bin])
                {
                    inFragmentIndex.Add(id);
                }
            }
            Assert.That(inFragmentIndex.Overlaps(undefinedIds), Is.True,
                "but it stays in the fragment index, where an open search can still reach it");
        }

        /// <summary>
        /// The contrast that defines correct handling of an undefined mass, asserted end to end on one
        /// fixture: a tolerance acceptor must never match it, an open acceptor must. Dropping such peptides
        /// from the index satisfied the first and broke the second, which is what the Windows suite caught.
        /// </summary>
        [Test]
        public static void ModernSearchEngine_UndefinedMassMatchesOnlyUnderAnOpenAcceptor()
        {
            var parameters = new CommonParameters(scoreCutoff: 1,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };
            var proteins = new List<Protein> { new Protein("MNNNKQXQ", "X-CONTAINING") };

            var indexEngine = new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(), null, null, null,
                0, DecoyType.None, parameters, fsp, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var index = (IndexingResults)indexEngine.Run();

            int undefinedId = Enumerable.Range(0, index.PeptideIndex.Count)
                .First(i => double.IsNaN(index.PeptideIndex[i].MonoisotopicMass));

            // a scan made from the undefined-mass peptide's own defined fragments
            var products = new List<Product>();
            index.PeptideIndex[undefinedId].Fragment(parameters.DissociationType, FragmentationTerminus.Both,
                products, parameters.FragmentationParameters);
            var neutral = products.Select(pr => pr.NeutralMass).Where(m => m > 0 && !double.IsNaN(m)).Distinct().ToList();
            Assert.That(neutral, Is.Not.Empty, "an X-containing peptide still has some defined fragments");

            double[] mz = neutral.Select(m => m.ToMz(1)).ToArray();
            var dataScan = new MsDataScan(new MzSpectrum(mz, neutral.Select(_ => 1000.0).ToArray(), false), 1, 2, true,
                Polarity.Positive, 1, new MzRange(0, 5000), "", MZAnalyzerType.Orbitrap, 1000.0 * mz.Length, null, null, "");
            var envelopes = neutral.Select(m => new IsotopicEnvelope(new List<(double, double)> { (m.ToMz(1), 1000.0) }, m, 1, 1000.0, 0)).ToArray();
            var scan = new Ms2ScanWithSpecificMass(dataScan, 500.0, 1, "", parameters, envelopes);

            static bool Scored(IndexingResults index, Ms2ScanWithSpecificMass scan, CommonParameters parameters,
                List<(string, CommonParameters)> fsp, MassDiffAcceptor acceptor, int undefinedId)
            {
                var matches = new SpectralMatch[1];
                new ModernSearchEngine(matches, new[] { scan }, index.PeptideIndex, index.FragmentIndex, 0,
                    parameters, fsp, acceptor, 0, new List<string>()).Run();
                return matches[0] != null
                    && matches[0].BestMatchingBioPolymersWithSetMods.Any(b => ReferenceEquals(b.SpecificBioPolymer, index.PeptideIndex[undefinedId]));
            }

            Assert.That(Scored(index, scan, parameters, fsp, new OpenSearchMode(), undefinedId), Is.True,
                "an open acceptor takes any precursor mass, so the undefined-mass peptide is a valid answer");
            Assert.That(Scored(index, scan, parameters, fsp, new SinglePpmAroundZeroSearchMode(5), undefinedId), Is.False,
                "a tolerance acceptor can never match an undefined mass");
        }

        /// <summary>
        /// The crosslink and glyco engines score through IndexedScoring rather than through IndexScoreScan, so
        /// the two window bugs have to be checked on that path too. Both are exercised here on one bin: an
        /// undefined mass must not be scored under a finite window, and a run of equal masses on the window's
        /// lower edge must be scored in full.
        /// </summary>
        [Test]
        public static void IndexedScoring_HandlesUndefinedMassesAndEqualMassRuns()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 1);
            var peptideIndex = IndexingEngine.SortByMonoisotopicMass(MakeAnagramProteins()
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .ToList());
            var masses = peptideIndex.Select(p => p.MonoisotopicMass).ToList();

            int undefined = masses.Count(double.IsNaN);
            Assert.That(undefined, Is.GreaterThan(0), "the fixture must contain undefined masses");

            // longest run of equal defined masses
            int runStart = undefined, runLength = 1, bestStart = undefined, bestLength = 1;
            for (int i = undefined + 1; i < masses.Count; i++)
            {
                if (masses[i] == masses[i - 1]) { runLength++; } else { runStart = i; runLength = 1; }
                if (runLength > bestLength) { bestLength = runLength; bestStart = runStart; }
            }
            Assert.That(bestLength, Is.GreaterThan(2), "and a run of equal masses");

            // one bin holding every peptide, so the undefined masses lead it exactly as they do in a real index
            var bin = Enumerable.Range(0, masses.Count).ToList();
            double edge = masses[bestStart];

            var observed = BinSearchProbe.Score(new FragmentIndex(new[] { 0, bin.Count }, bin.ToArray()),
                new List<int> { 0 }, edge, masses[^1] + 1, peptideIndex);

            Assert.That(observed.Any(id => double.IsNaN(masses[id])), Is.False,
                "no undefined mass may be scored under a finite window");
            for (int i = bestStart; i < bestStart + bestLength; i++)
            {
                Assert.That(observed, Does.Contain(i),
                    $"entry {i} sits exactly on the window's lower edge and must be scored");
            }
        }

        /// <summary>
        /// The crosslink and glyco engines already had a correct lower-bound search — CrosslinkSearchEngine
        /// and GlycoPeptides both take the bitwise complement for the insertion point and then walk back over
        /// a run of equal masses. Modern search was the odd one out, asking the upper-bound search for its
        /// window start. This pins the two against each other so they cannot drift apart again: for every
        /// target, the modern window's start must be the index that the crosslink engine's independent
        /// implementation already returns.
        /// </summary>
        [Test]
        public static void WindowStart_AgreesWithTheCrosslinkEnginesLowerBoundSearch()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 5);
            // anagrams for the equal-mass runs that distinguish the two searches, plus varied proteins so the
            // comparison also covers ordinary distinct masses and the gaps between them
            var peptideIndex = MakeAnagramProteins().Concat(MakeProteins(300))
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .Where(p => !double.IsNaN(p.MonoisotopicMass))
                .OrderBy(p => p.MonoisotopicMass)
                .ToList();
            var masses = peptideIndex.Select(p => p.MonoisotopicMass).ToArray();
            var bin = Enumerable.Range(0, masses.Length).ToList();

            var targets = new List<double>();
            foreach (double m in masses.Distinct())
            {
                targets.Add(m);            // exactly on a run of equal masses, the case that differed
                targets.Add(m + 1e-4);
                targets.Add(m - 1e-4);
            }
            targets.Add(masses[0] - 500);
            targets.Add(masses[^1] + 500);

            var disagreements = new List<string>();
            foreach (double target in targets)
            {
                int modern = BinSearchProbe.FirstAtOrAbove(CollectionsMarshal.AsSpan(bin), target, peptideIndex);
                int crosslink = CrosslinkSearchEngine.BinarySearchGetIndex(masses, target);

                if (modern != crosslink && disagreements.Count < 10)
                {
                    disagreements.Add($"target {target:F5}: modern {modern}, crosslink {crosslink}");
                }
            }

            Assert.That(targets.Count, Is.GreaterThan(100));
            Assert.That(disagreements, Is.Empty,
                "the two lower-bound searches must agree: " + string.Join("; ", disagreements));
        }

        /// <summary>
        /// Classic search does not use the fragment index at all, and it handles an undefined mass correctly
        /// by accident of comparison semantics: the allowed interval it derives is NaN-bounded, and every
        /// "scanMass &lt;= maximum" test against NaN is false, so no scan is considered. Modern search must
        /// reach the same answer through completely different machinery. Asserting the two agree is what
        /// makes "open matches, tolerance does not" a property of MetaMorpheus rather than of one engine.
        /// </summary>
        [Test]
        public static void UndefinedMassContract_IsTheSameInClassicAndModernSearch()
        {
            var parameters = new CommonParameters(scoreCutoff: 1,
                digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));
            var proteins = new List<Protein> { new Protein("MNNNKQXQ", "X-CONTAINING") };
            var undefined = proteins[0].Digest(parameters.DigestionParams, new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>()
                .First(p => double.IsNaN(p.MonoisotopicMass));

            // classic search asks the acceptor for the intervals it should look in
            var openIntervals = new OpenSearchMode()
                .GetAllowedPrecursorMassIntervalsFromTheoreticalMass(undefined.MonoisotopicMass).ToList();
            var toleranceIntervals = new SinglePpmAroundZeroSearchMode(5)
                .GetAllowedPrecursorMassIntervalsFromTheoreticalMass(undefined.MonoisotopicMass).ToList();

            // a scan can only be considered when its mass is inside an interval; against a NaN bound every
            // comparison is false, so a tolerance acceptor considers nothing while an open one considers all
            const double someScanMass = 500.0;
            Assert.That(openIntervals.Any(i => someScanMass >= i.Minimum && someScanMass <= i.Maximum), Is.True,
                "classic search would consider scans for this peptide under an open acceptor");
            Assert.That(toleranceIntervals.Any(i => someScanMass >= i.Minimum && someScanMass <= i.Maximum), Is.False,
                "and none under a tolerance acceptor");

            // modern search reaches the same two answers through the bin window instead
            var fsp = new List<(string, CommonParameters)> { ("", parameters) };
            var engine = new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(), null, null, null,
                0, DecoyType.None, parameters, fsp, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var index = (IndexingResults)engine.Run();
            int undefinedId = Enumerable.Range(0, index.PeptideIndex.Count)
                .First(i => double.IsNaN(index.PeptideIndex[i].MonoisotopicMass));

            int binWithIt = Enumerable.Range(0, index.FragmentIndex.Length)
                .First(b => index.FragmentIndex[b].ToArray().Contains(undefinedId));
            ReadOnlySpan<int> entries = index.FragmentIndex[binWithIt];
            int position = entries.ToArray().ToList().IndexOf(undefinedId);

            var (openStart, openEnd) = BinSearchProbe.Window(double.NegativeInfinity, double.PositiveInfinity, entries, index.PeptideIndex);
            Assert.That(position, Is.InRange(openStart, openEnd),
                "an open search leaves both bounds infinite, so the whole bin including this peptide is scored");

            var (tolStart, tolEnd) = BinSearchProbe.Window(499.0, 501.0, entries, index.PeptideIndex);
            Assert.That(position < tolStart || position > tolEnd, Is.True,
                "a finite window must place an undefined mass outside it");
        }

        /// <summary>
        /// The case that actually distinguishes reading an undefined mass as negative infinity from taking it
        /// literally: a bin whose entries are mostly undefined. With a handful of NaNs leading a large bin the
        /// binary search's first probe never lands among them, so the bug hides; when they are the majority it
        /// probes one immediately, every comparison is false, the search walks left and reports an empty
        /// window — losing the one real candidate in the bin.
        ///
        /// This is not a contrived shape. A bin holds the peptides sharing one fragment mass, so it is often
        /// only a handful of entries, and a database with many unknown residues puts several of them in bins
        /// that small.
        /// </summary>
        [Test]
        public static void WindowSurvivesABinThatIsMostlyUndefinedMasses()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 1);

            // three peptides with an unknown residue and one ordinary one
            var proteins = new List<Protein>
            {
                new Protein("QXQK", "X1"), new Protein("NXNK", "X2"), new Protein("AXAK", "X3"),
                new Protein("PEPTIDEK", "REAL"),
            };
            var peptideIndex = IndexingEngine.SortByMonoisotopicMass(proteins
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .ToList());

            var undefinedIds = Enumerable.Range(0, peptideIndex.Count)
                .Where(i => double.IsNaN(peptideIndex[i].MonoisotopicMass)).ToList();
            var realIds = Enumerable.Range(0, peptideIndex.Count)
                .Where(i => !double.IsNaN(peptideIndex[i].MonoisotopicMass)).ToList();
            Assert.That(undefinedIds.Count, Is.GreaterThanOrEqualTo(3), "need the undefined masses to dominate");
            Assert.That(realIds, Is.Not.Empty);

            // the bin the sort would produce: undefined masses first, then the real one
            int realId = realIds[0];
            var bin = undefinedIds.Concat(new[] { realId }).ToList();
            double realMass = peptideIndex[realId].MonoisotopicMass;

            var (start, end) = BinSearchProbe.Window(realMass - 1, realMass + 1,
                CollectionsMarshal.AsSpan(bin), peptideIndex);

            int position = bin.IndexOf(realId);
            Assert.That(position, Is.InRange(start, end),
                "the one real candidate in the bin must be inside the window, not lost behind the undefined ones");
            Assert.That(start, Is.EqualTo(position),
                "and the window must start at it, excluding every undefined mass");
        }

        /// <summary>
        /// The case reported in issue #1542 in 2018: a bin of twelve entries with ascending masses, looking
        /// for the third one. The search returned 0 instead of 2.
        ///
        /// The reply at the time was that the imprecision was deliberate - the search "errs on the side of
        /// caution in that it goes back a space", so as not to skip a valid hit. That reasoning holds for the
        /// start of the precursor window, where starting early only widens it. The same function also
        /// computed the window's end, though, and there going back a space drops the last valid candidate
        /// instead of including an extra one. Probed against a linear scan over 104,472 (bin length, target)
        /// pairs, the old implementation disagreed on 1,390 of them, once answering 0 where 719 was correct.
        ///
        /// Both ends are now exact, and each end uses the search that answers its own question.
        /// </summary>
        [Test]
        public static void BinarySearchBinForPrecursorIndex_ReportedCaseFromIssue1542()
        {
            var digestion = new DigestionParams(protease: "trypsin", minPeptideLength: 5);
            var peptideIndex = MakeProteins(400)
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<PeptideWithSetModifications>()
                .Where(p => !double.IsNaN(p.MonoisotopicMass))
                .GroupBy(p => p.MonoisotopicMass)
                .Select(g => g.First())
                .OrderBy(p => p.MonoisotopicMass)
                .Take(12)
                .ToList();
            Assert.That(peptideIndex.Count, Is.EqualTo(12), "the reported case is a bin of twelve");

            var bin = Enumerable.Range(0, 12).ToList();
            double third = peptideIndex[2].MonoisotopicMass;

            Assert.That(BinSearchProbe.Search(CollectionsMarshal.AsSpan(bin), third, peptideIndex),
                Is.EqualTo(2), "the last entry at or below the third mass is the third entry, not the first");
            Assert.That(BinSearchProbe.FirstAtOrAbove(CollectionsMarshal.AsSpan(bin), third, peptideIndex),
                Is.EqualTo(2), "and so is the first entry at or above it");

            // every position in the bin, not just the reported one
            var wrong = new List<string>();
            for (int i = 0; i < 12; i++)
            {
                double m = peptideIndex[i].MonoisotopicMass;
                int upper = BinSearchProbe.Search(CollectionsMarshal.AsSpan(bin), m, peptideIndex);
                int lower = BinSearchProbe.FirstAtOrAbove(CollectionsMarshal.AsSpan(bin), m, peptideIndex);
                if (upper != i || lower != i)
                {
                    wrong.Add($"index {i}: upper {upper}, lower {lower}");
                }
            }
            Assert.That(wrong, Is.Empty, string.Join("; ", wrong));
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
