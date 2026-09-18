using System;

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

        /// <param name="fileCount"> Spectra files the task will search. </param>
        /// <param name="threadBudget"> The task's MaxThreadsToUsePerFile, treated as the budget for all files together. </param>
        /// <param name="availableBytes"> Free physical memory, measured with the index and one file already in memory. </param>
        /// <param name="bytesPerFile"> Estimated memory one more file adds while it is searched; 0 when unknown. </param>
        /// <param name="maximumFilesInParallel"> A user cap: 0 for no cap, 1 to search files one after another. </param>
        public static FileParallelismPlan Decide(int fileCount, int threadBudget, long availableBytes, long bytesPerFile, int maximumFilesInParallel = 0)
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
                long byMemoryLong = 1 + (long)(Math.Max(0, availableBytes) * MemoryBudgetFraction) / bytesPerFile;
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
    }
}
