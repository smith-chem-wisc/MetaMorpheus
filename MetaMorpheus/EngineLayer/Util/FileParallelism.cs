using System;
using System.Collections.Concurrent;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Util
{
    /// <summary>
    /// How many spectra files a task searches at once, and with how many threads each.
    /// </summary>
    /// <param name="FilesInParallel"> Files searched at the same time; 1 means one after another. </param>
    /// <param name="ThreadsPerFile"> MaxThreadsToUsePerFile handed to each file's loading, scan extraction and search. </param>
    /// <param name="LimitedBy"> What set <paramref name="FilesInParallel"/>, in words fit for the manuscript prose. </param>
    public sealed record FileParallelismPlan(int FilesInParallel, int ThreadsPerFile, string LimitedBy);

    /// <summary>
    /// Divides a task's thread budget across spectra files searched at once, bounded by the number of files, a minimum number of
    /// threads for each, and the memory the machine has free. The decision is kept apart from measuring the machine so it can be
    /// tested without spectra or a particular amount of RAM.
    /// </summary>
    public static class FileParallelism
    {
        /// <summary>
        /// Fewest threads a file is given before searching another file at the same time stops paying: below this the search of
        /// each file, which is almost all CPU, slows more than running files side by side gains.
        /// </summary>
        public const int MinimumThreadsPerFile = 4;

        /// <summary>
        /// Share of the free memory the additional files may take up, leaving the rest for the results they accumulate and for
        /// everything else on the machine.
        /// </summary>
        public const double MemoryBudgetFraction = 0.8;

        /// <summary>
        /// Memory the search's scoring tables take: each thread searching holds two tables of one byte per peptide in the index.
        /// The threads are one budget shared by every file searched at once, so this is paid once per task, not once per file.
        /// </summary>
        public static long ScoringTableBytes(int threadBudget, int peptideCount)
        {
            return 2L * Math.Max(1, threadBudget) * peptideCount;
        }

        /// <param name="fileCount"> Spectra files the task will search. </param>
        /// <param name="threadBudget"> The task's MaxThreadsToUsePerFile, treated as the budget for all files together. </param>
        /// <param name="availableBytes"> Free physical memory, measured with the index and one file already in memory. </param>
        /// <param name="bytesPerFile"> Estimated memory one more file adds while it is searched; 0 when unknown. </param>
        /// <param name="maximumFilesInParallel"> A user cap: 0 for no cap, 1 to search files one after another. </param>
        /// <param name="fixedBytes"> Memory the search needs however many files run at once (see <see cref="ScoringTableBytes"/>), so it
        /// comes off the budget rather than being charged to each file, as <c>IndexPartitioning.PartitionsForBudget</c> does. </param>
        public static FileParallelismPlan Decide(int fileCount, int threadBudget, long availableBytes, long bytesPerFile, int maximumFilesInParallel = 0,
            long fixedBytes = 0)
        {
            int budget = Math.Max(1, threadBudget);
            if (fileCount <= 1)
            {
                return new FileParallelismPlan(1, budget, "there is one spectra file");
            }
            if (maximumFilesInParallel == 1)
            {
                return new FileParallelismPlan(1, budget, "the task setting MaximumSpectraFilesInParallel = 1");
            }

            int files = fileCount;
            string limitedBy = "the number of spectra files";

            int byThreads = Math.Max(1, budget / MinimumThreadsPerFile);
            if (byThreads < files)
            {
                files = byThreads;
                limitedBy = $"the thread budget of {budget} at no fewer than {MinimumThreadsPerFile} threads per file";
            }

            if (bytesPerFile > 0)
            {
                long budgetBytes = Math.Max(0, (long)(Math.Max(0, availableBytes) * MemoryBudgetFraction) - Math.Max(0, fixedBytes));
                long byMemoryLong = 1 + budgetBytes / bytesPerFile;
                int byMemory = (int)Math.Min(int.MaxValue, byMemoryLong);
                if (byMemory < files)
                {
                    files = byMemory;
                    limitedBy = "free memory";
                }
            }

            if (maximumFilesInParallel > 1 && maximumFilesInParallel < files)
            {
                files = maximumFilesInParallel;
                limitedBy = $"the task setting MaximumSpectraFilesInParallel = {maximumFilesInParallel}";
            }

            return new FileParallelismPlan(files, Math.Max(1, budget / files), limitedBy);
        }

        /// <summary>
        /// Runs <paramref name="searchFile"/> for every file index from 0 to <paramref name="fileCount"/> - 1, no more than
        /// <paramref name="filesInParallel"/> of them at a time, and hands the next index out only when a worker is free.
        /// </summary>
        /// <remarks>
        /// NoBuffering, rather than Parallel.For: Parallel.For claims indices in chunks that double in size (1, then 2, then
        /// 4...) and a worker keeps the whole chunk it claimed, so a file behind a slow one in the same chunk cannot be picked
        /// up by a worker that is free. With enough files - more than twice <paramref name="filesInParallel"/> - that brings
        /// back the idle file slot this class is here to remove.
        /// </remarks>
        public static void ForEachFile(int fileCount, int filesInParallel, Action<int> searchFile)
        {
            Parallel.ForEach(
                Partitioner.Create(Enumerable.Range(0, fileCount), EnumerablePartitionerOptions.NoBuffering),
                new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, filesInParallel) },
                searchFile);
        }
    }
}
