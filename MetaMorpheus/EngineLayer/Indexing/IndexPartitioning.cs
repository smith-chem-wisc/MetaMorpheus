using System;
using System.Collections.Generic;
using MassSpectrometry;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.Indexing
{
    /// <summary>
    /// Chooses a TotalPartitions that keeps one index build inside the memory the process actually has.
    /// One partition of a mammalian database needs several GB; at TotalPartitions = 1 a 16 GB machine
    /// pages instead of failing, which looks like a hang rather than a configuration problem.
    /// </summary>
    public static class IndexPartitioning
    {
        /// <summary>
        /// Share of currently-free physical memory one index build may occupy; the rest is left for the
        /// spectra, the PSMs and quantification.
        ///
        /// Measured against free memory rather than installed memory: a fixed fraction of installed is too
        /// generous on a busy laptop and far too mean on an idle cluster node with hundreds of free GB.
        /// The cost is that the partition count is not reproducible across runs, which matters only to the
        /// extent that results still depend on it — see ModernSearchEngine.FineScorePeptides.
        /// </summary>
        private const double MemoryBudgetFraction = 0.55;

        /// <summary>
        /// Heap bytes per peptide in the peptide index. Measured on mouse at 9.21 M peptides: 6,926 MB of
        /// live heap after digestion, i.e. ~750 B/peptide. Server GC keeps far more garbage before
        /// collecting, so the same work occupies about twice the footprint and the estimate has to know
        /// which mode it is running in.
        /// </summary>
        private static long BytesPerPeptide => System.Runtime.GCSettings.IsServerGC ? 1500 : 750;

        /// <summary>
        /// Bytes per fragment-index entry at the build's peak, which is more than the finished index costs.
        /// The finished array is one int per entry; while FragmentIndexBuilder.Build runs, every block's
        /// emission run is alive alongside it, and the runs grow by doubling, so their capacity overshoots.
        /// Twelve covers all three. Under-counting is the direction that makes the guard useless.
        ///
        /// This used to be measured (~11 B under workstation GC, ~20 B under server GC at 4,000 targets /
        /// 125.09 M binned fragments) and used to depend on the GC mode, because the entry carried a share
        /// of per-bin List&lt;int&gt; overhead and of the garbage repeated doubling left behind.
        ///
        /// Counting fragments separately rather than folding them into a per-peptide constant is what makes
        /// the estimate hold for searches whose peptides are long: fragments per peptide grows with peptide
        /// length, so a top-down or unbounded-length search has far more index entries per peptide than the
        /// tryptic case the peptide constant was calibrated on.
        /// </summary>
        private const long BytesPerFragmentEntry = 3 * sizeof(int);

        /// <summary>
        /// Bytes of bin offsets: one int per bin plus a terminator. Set by the maximum fragment mass, not by
        /// the database, so unlike everything else here it does not shrink when the database is partitioned
        /// and has to be taken off the budget instead of divided into it. 120 MB at the default 30,000 Da.
        /// </summary>
        private static long BinSpaceBytes(double maxFragmentSize) =>
            sizeof(int) * ((long)Math.Ceiling(maxFragmentSize) * FragmentBinsPerDalton + 2);

        /// <summary>Bins per dalton, matching IndexingEngine.</summary>
        private const int FragmentBinsPerDalton = 1000;

        /// <summary>
        /// Fragment entries one partition's index is sized to hold. The fragment index is one flat int[], so
        /// <see cref="Array.MaxLength"/> is a ceiling that no amount of free memory lifts, and on a large node the
        /// memory budget alone would leave a database that exceeds it in one partition.
        ///
        /// Half the ceiling, because partitions are cut by protein count rather than by fragments. Measured on
        /// reviewed human UniProt (78,120 entries), the largest contiguous slice holds 1.21x the mean residue
        /// count at 2 partitions, 1.60x at 4, 1.75x at 64 and 1.95x at 256. The fragment count itself is not
        /// where the risk is: on the same database with decoys, EstimateIndexSize gave 584.1 M entries against
        /// 579.3 M built at 2 missed cleavages, and 1,461.7 M against 1,457.5 M at 4.
        /// </summary>
        private static long MaxFragmentEntriesPerPartition => Array.MaxLength / 2;

        /// <summary>Proteins digested to learn the peptide yield of the configured digestion settings.</summary>
        private const int SampleSize = 200;

        /// <summary>Upper bound on peptides examined, so the estimate cannot outcost the thing it sizes.</summary>
        private const int MaxSampledPeptides = 50_000;

        private const int MaxPartitions = 256;

        /// <summary>
        /// Returns the smallest partition count &gt;= <paramref name="requestedPartitions"/> whose estimated
        /// peak index footprint fits the memory budget and keeps each partition under the fragment entry limit.
        /// Never returns less than what was requested.
        /// </summary>
        /// <param name="cappedByMemory">
        /// Even the maximum partition count does not fit the memory budget. Memory only: the entry limit's own cap
        /// is a different problem with a different fix, and <see cref="PartitionWarnings"/> works it out from
        /// <paramref name="estimatedFragmentEntries"/>.
        /// </param>
        public static int SuggestTotalPartitions(IReadOnlyList<IBioPolymer> proteins, CommonParameters commonParameters,
            List<Modification> fixedModifications, List<Modification> variableModifications,
            List<SilacLabel> silacLabels, SilacLabel startLabel, SilacLabel endLabel, double maxFragmentSize,
            int requestedPartitions, out long estimatedBytes, out long budgetBytes, out bool cappedByMemory,
            out long estimatedFragmentEntries)
        {
            estimatedBytes = 0;
            estimatedFragmentEntries = 0;
            budgetBytes = AvailableBytes();
            cappedByMemory = false;

            if (proteins == null || proteins.Count == 0 || commonParameters.DigestionParams is not DigestionParams digestionParams)
            {
                return requestedPartitions;
            }

            // mirrors how IndexingEngine derives its turnover labels from the two search parameters
            (SilacLabel, SilacLabel)? turnoverLabels = startLabel != null || endLabel != null
                ? (startLabel, endLabel)
                : null;

            (long peptideCount, long fragmentCount) = EstimateIndexSize(proteins, commonParameters, digestionParams,
                fixedModifications, variableModifications, silacLabels, turnoverLabels, maxFragmentSize);
            estimatedBytes = peptideCount * BytesPerPeptide + fragmentCount * BytesPerFragmentEntry;
            estimatedFragmentEntries = fragmentCount;

            int forMemory = PartitionsForBudget(estimatedBytes, BinSpaceBytes(maxFragmentSize), budgetBytes,
                requestedPartitions, proteins.Count, out cappedByMemory);
            int forEntries = PartitionsForEntryLimit(fragmentCount, requestedPartitions, proteins.Count, out _);

            return Math.Max(forMemory, forEntries);
        }

        /// <summary>
        /// How many partitions keep each partition's fragment index under <see cref="MaxFragmentEntriesPerPartition"/>,
        /// never going below <paramref name="requestedPartitions"/> and never exceeding one partition per protein.
        /// Independent of memory: a node with room for a 3 B-entry index still cannot put it in one int[].
        /// </summary>
        internal static int PartitionsForEntryLimit(long estimatedFragmentEntries, int requestedPartitions, int proteinCount,
            out bool cappedByLimit)
        {
            cappedByLimit = false;
            if (estimatedFragmentEntries <= 0)
            {
                return requestedPartitions;
            }

            double exact = Math.Ceiling(estimatedFragmentEntries / (double)MaxFragmentEntriesPerPartition);
            cappedByLimit = exact > MaxPartitions;
            int needed = (int)Math.Min(MaxPartitions, exact);
            return Math.Max(requestedPartitions, Math.Min(needed, Math.Max(1, proteinCount)));
        }

        /// <summary>
        /// The decision itself, separated from measuring the database so it can be exercised directly:
        /// how many partitions does <paramref name="estimatedBytes"/> need to fit the share of
        /// <paramref name="availableBytes"/> budgeted for an index, never going below
        /// <paramref name="requestedPartitions"/> and never exceeding one partition per protein.
        /// <paramref name="fixedBytes"/> is the part of an index that partitioning does not shrink, so it
        /// comes off the budget rather than being divided into it.
        /// </summary>
        internal static int PartitionsForBudget(long estimatedBytes, long fixedBytes, long availableBytes,
            int requestedPartitions, int proteinCount, out bool cappedByLimit)
        {
            cappedByLimit = false;

            // nothing measurable to size against, so leave the setting alone
            if (availableBytes <= 0)
            {
                return requestedPartitions;
            }

            long budget = (long)(availableBytes * MemoryBudgetFraction) - fixedBytes;

            // the part partitioning cannot shrink already exceeds the budget. This is a measurement, not a
            // missing one, so split as far as possible and say so rather than returning the request unchanged
            if (budget <= 0)
            {
                cappedByLimit = true;
                return Math.Max(requestedPartitions, Math.Min(MaxPartitions, Math.Max(1, proteinCount)));
            }

            double exact = Math.Ceiling(estimatedBytes / (double)budget);
            cappedByLimit = exact > MaxPartitions;
            int needed = (int)Math.Min(MaxPartitions, exact);
            return Math.Max(requestedPartitions, Math.Min(needed, Math.Max(1, proteinCount)));
        }

        /// <summary>
        /// Digests a stride-spaced sample to get peptides per residue for the configured settings, then
        /// scales by the whole database's residue count. Sampling rather than assuming a fixed yield is
        /// what makes this hold up across proteases, missed-cleavage counts and variable-mod loads.
        /// </summary>
        private static (long Peptides, long FragmentEntries) EstimateIndexSize(IReadOnlyList<IBioPolymer> proteins,
            CommonParameters commonParameters, DigestionParams digestionParams,
            List<Modification> fixedModifications, List<Modification> variableModifications,
            List<SilacLabel> silacLabels, (SilacLabel, SilacLabel)? turnoverLabels, double maxFragmentSize)
        {
            long totalResidues = 0;
            foreach (IBioPolymer protein in proteins)
            {
                totalResidues += protein.BaseSequence.Length;
            }

            // DissociationType.Custom resolves its product types from a process-global list that the
            // CommonParameters constructor clears, and only MetaMorpheusEngine.Run() restores it. This
            // samples outside any engine, so without restoring it here a Custom search would fragment into
            // an empty product list, count zero fragments, and under-partition — silently, and precisely in
            // the case the guard exists for. Put the global back afterwards so the estimate leaks no state.
            bool isCustom = commonParameters.DissociationType == DissociationType.Custom;
            List<ProductType> previousCustomProducts = null;
            if (isCustom)
            {
                previousCustomProducts = digestionParams.ProductsFromDissociationType()[DissociationType.Custom];
                commonParameters.SetCustomProductTypes();
            }

            try
            {

            int stride = Math.Max(1, proteins.Count / SampleSize);
            long sampledResidues = 0;
            long sampledPeptides = 0;
            long sampledFragments = 0;
            var products = new List<Product>();
            for (int i = 0; i < proteins.Count; i += stride)
            {
                // Stop at a protein boundary once enough peptides have been seen, so the estimate itself
                // cannot become the slow part. A non-specific ("single") protease yields orders of magnitude
                // more peptides per protein than a tryptic one, and each is fragmented here.
                if (sampledPeptides >= MaxSampledPeptides)
                {
                    break;
                }

                // Digesting a zero-length sequence throws, and this runs while the task is still being
                // configured, so it would surface as a crash report rather than as a search that ran.
                if (proteins[i].BaseSequence.Length == 0)
                {
                    continue;
                }

                sampledResidues += proteins[i].BaseSequence.Length;
                foreach (var peptide in proteins[i].Digest(digestionParams, fixedModifications, variableModifications, silacLabels, turnoverLabels))
                {
                    sampledPeptides++;

                    // count the fragments this peptide would actually contribute to the index, applying the
                    // same mass window the indexing loop applies
                    peptide.Fragment(commonParameters.DissociationType, digestionParams.FragmentationTerminus,
                        products, commonParameters.FragmentationParameters);
                    foreach (Product product in products)
                    {
                        if (product.NeutralMass > 0 && product.NeutralMass < maxFragmentSize)
                        {
                            sampledFragments++;
                        }
                    }
                }
            }

            if (sampledResidues == 0)
            {
                return (0, 0);
            }

            double scale = totalResidues / (double)sampledResidues;
            return ((long)(sampledPeptides * scale), (long)(sampledFragments * scale));

            }
            finally
            {
                if (isCustom)
                {
                    digestionParams.ProductsFromDissociationType()[DissociationType.Custom] = previousCustomProducts;
                }
            }
        }

        /// <summary>
        /// Physical memory that is genuinely free, clamped by a batch scheduler's allocation when there is one.
        ///
        /// TotalAvailableMemoryBytes is installed (or container-limited) memory rather than headroom, and
        /// subtracting only this process's managed heap would ignore every other process, the file cache and
        /// our own unmanaged allocations. MemoryLoadBytes is the system-wide in-use figure, which is the right
        /// basis — but it reads 0 until the first collection has run, so an unset value is forced rather than
        /// taken to mean the whole machine is free.
        ///
        /// A scheduler's allocation is a separate ceiling and is not always enforced with cgroups; where it is
        /// not, this process sees the whole node, so a job entitled to 16 GB of a 512 GB node would otherwise
        /// size its index against the node's free memory and be killed for exceeding its allocation.
        /// </summary>
        private static long AvailableBytes()
        {
            GCMemoryInfo info = GC.GetGCMemoryInfo();
            long total = info.TotalAvailableMemoryBytes;
            if (total <= 0)
            {
                return 0;
            }

            long inUse = info.MemoryLoadBytes;
            if (inUse <= 0)
            {
                GC.Collect(0, GCCollectionMode.Forced, blocking: true);
                inUse = GC.GetGCMemoryInfo().MemoryLoadBytes;
            }

            // No usable in-use figure: fall back to a deliberately pessimistic slice rather than reporting the
            // whole machine, which would under-partition exactly when the guard matters most.
            long free = inUse > 0 && inUse < total ? total - inUse : total / 4;

            long schedulerLimit = SchedulerMemoryLimitBytes();
            if (schedulerLimit > 0)
            {
                // our own footprint already counts against the allocation
                free = Math.Min(free, Math.Max(schedulerLimit / 4, schedulerLimit - Environment.WorkingSet));
            }

            return free;
        }

        /// <summary>
        /// The memory this job has been allocated by a batch scheduler, in bytes, or 0 if it is not running
        /// under one. Slurm advertises <c>SLURM_MEM_PER_NODE</c> for <c>--mem</c> and
        /// <c>SLURM_MEM_PER_CPU</c> (times <c>SLURM_CPUS_ON_NODE</c>) for <c>--mem-per-cpu</c>, both in MB;
        /// 0 or absent means no limit was requested.
        /// </summary>
        internal static long SchedulerMemoryLimitBytes()
        {
            const long megabyte = 1024 * 1024;

            long perNode = ReadPositiveLong("SLURM_MEM_PER_NODE");
            if (perNode > 0)
            {
                return perNode * megabyte;
            }

            long perCpu = ReadPositiveLong("SLURM_MEM_PER_CPU");
            long cpus = ReadPositiveLong("SLURM_CPUS_ON_NODE");
            return perCpu > 0 && cpus > 0 ? perCpu * cpus * megabyte : 0;
        }

        private static long ReadPositiveLong(string variable)
        {
            return long.TryParse(Environment.GetEnvironmentVariable(variable), out long value) && value > 0
                ? value
                : 0;
        }

        /// <summary>
        /// Everything the user should be told about a raised partition count: that it was raised and why, and --
        /// when even the maximum is not enough -- what is still wrong. Returned rather than emitted so the set of
        /// messages can be asserted without depending on the test machine's memory.
        /// </summary>
        public static IEnumerable<string> PartitionWarnings(int requestedPartitions, int suggestedPartitions,
            long estimatedBytes, long budgetBytes, bool cappedByMemory, long estimatedFragmentEntries = 0)
        {
            // The count is the larger of what memory needs and what the entry limit needs, so explain whichever set
            // it. Asking only whether the entry limit would raise the count at all gets this wrong wherever memory is
            // tighter -- a 16 GB machine, long before the entry limit matters -- and then hides the memory figures
            // that explain the number. The entry rule depends on nothing but the estimate, so it is recomputed here;
            // without a protein-count bound, which caps both rules at the same number and so decides nothing. When
            // both rules need the same count the entry message is the one given.
            int forEntries = PartitionsForEntryLimit(estimatedFragmentEntries, requestedPartitions, int.MaxValue, out bool cappedByEntries);
            if (forEntries > requestedPartitions && forEntries >= suggestedPartitions)
            {
                yield return PartitionEntryLimitWarning(requestedPartitions, suggestedPartitions, estimatedFragmentEntries);
            }
            else
            {
                yield return PartitionIncreaseWarning(requestedPartitions, suggestedPartitions, estimatedBytes, budgetBytes);
            }

            // Two different ways for the maximum to fall short, with different fixes. Kept apart: a cap from the
            // entry limit reported as a memory problem would suggest paging when the risk is the build failing.
            if (cappedByEntries)
            {
                yield return $"Even {suggestedPartitions} partitions leaves more than {MaxFragmentEntriesPerPartition / 1e9:N2} " +
                             $"billion fragment index entries in each, so building the fragment index may fail. Consider a " +
                             $"smaller database, fewer missed cleavages, fewer variable modifications, or a lower maximum " +
                             $"fragment mass.";
            }

            if (cappedByMemory)
            {
                yield return $"Even {suggestedPartitions} partitions is estimated not to fit in memory; this search " +
                             $"may page heavily. Consider a smaller database, a shorter maximum peptide length, " +
                             $"or fewer variable modifications.";
            }
        }

        /// <summary>Message for the user when the entry limit is what raised the partition count.</summary>
        public static string PartitionEntryLimitWarning(int requestedPartitions, int suggestedPartitions, long estimatedFragmentEntries)
        {
            return $"Indexing this database in {requestedPartitions} partition(s) is estimated to need " +
                   $"{estimatedFragmentEntries / 1e9:N2} billion fragment index entries. Increasing TotalPartitions to " +
                   $"{suggestedPartitions} keeps each partition under {MaxFragmentEntriesPerPartition / 1e9:N2} billion, half " +
                   $"of what one partition's index can hold, because partitions are split by protein count and come out " +
                   $"uneven. Identifications and scores are unaffected; reported PEP and q-values can shift slightly.";
        }

        /// <summary>Message for the user when the requested partition count cannot fit.</summary>
        public static string PartitionIncreaseWarning(int requestedPartitions, int suggestedPartitions,
            long estimatedBytes, long budgetBytes)
        {
            double freeGb = budgetBytes / 1073741824.0;
            return $"Indexing this database in {requestedPartitions} partition(s) is estimated to need " +
                   $"{estimatedBytes / 1073741824.0:N1} GB, which does not fit the " +
                   $"{freeGb * MemoryBudgetFraction:N1} GB budgeted for the index " +
                   $"({MemoryBudgetFraction:P0} of {freeGb:N1} GB free). Increasing TotalPartitions to {suggestedPartitions} " +
                   $"so the index stays in memory. Identifications and scores are unaffected; reported PEP and " +
                   $"q-values can shift slightly, because the partition count is derived from free memory.";
        }
    }
}
