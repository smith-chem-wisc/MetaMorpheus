using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Runtime.ExceptionServices;
using System.Threading.Tasks;

namespace EngineLayer.Util
{
    /// <summary>
    /// Writes one line per item, in the items' order, building the lines on several threads.
    ///
    /// The rows of a results table do not depend on each other, and building them (formatting every column of every row) is
    /// most of the cost of writing one: on the aging reproducer's four-file search, about 85% of writing a .psmtsv on one
    /// thread, against 11-15% for the write itself on an SSD. So the lines of the next batch are built, on up to
    /// <c>maxThreads</c> threads, while the lines of the current batch are written.
    ///
    /// The text is what writing <c>line(item)</c> for each item on the calling thread writes: every line is built by the
    /// same function from the same item, and the lines are written in item order. Several columns are formatted in the
    /// current culture, so each thread that builds lines takes the calling thread's culture while it does, and gets its own
    /// back afterwards. A failure is thrown as itself, as the one-thread loop would throw it.
    /// </summary>
    public static class OrderedLines
    {
        /// <summary>
        /// Lines built per batch. Two batches are held at once, the one being written and the one being built: for
        /// .psmtsv rows of 1-3 kB, a few tens of megabytes.
        /// </summary>
        public const int BatchSize = 4096;

        /// <summary>
        /// Fewer items than this are written on the calling thread: the tasks would cost more than they save.
        /// </summary>
        public const int MinimumItemsForThreads = 256;

        public static void Write<T>(TextWriter output, IReadOnlyList<T> items, Func<T, string> line, int maxThreads)
        {
            if (maxThreads <= 1 || items.Count < MinimumItemsForThreads)
            {
                for (int i = 0; i < items.Count; i++)
                {
                    output.WriteLine(line(items[i]));
                }
                return;
            }

            CultureInfo culture = CultureInfo.CurrentCulture;
            CultureInfo uiCulture = CultureInfo.CurrentUICulture;
            Task<string[]> next = Build(items, 0, line, maxThreads, culture, uiCulture);
            try
            {
                for (int start = 0; start < items.Count; start += BatchSize)
                {
                    string[] lines = Result(next);
                    next = start + BatchSize < items.Count
                        ? Build(items, start + BatchSize, line, maxThreads, culture, uiCulture)
                        : null;
                    foreach (string text in lines)
                    {
                        output.WriteLine(text);
                    }
                }
            }
            catch
            {
                // the batch being built is not left running behind a failure; its own failure is not the one to report
                if (next != null)
                {
                    try { next.Wait(); } catch { }
                }
                throw;
            }
        }

        private static Task<string[]> Build<T>(IReadOnlyList<T> items, int start, Func<T, string> line, int maxThreads,
            CultureInfo culture, CultureInfo uiCulture)
        {
            int count = Math.Min(BatchSize, items.Count - start);
            return Task.Run(() =>
            {
                var lines = new string[count];
                Parallel.For(0, count, new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                    () =>
                    {
                        var own = (CultureInfo.CurrentCulture, CultureInfo.CurrentUICulture);
                        CultureInfo.CurrentCulture = culture;
                        CultureInfo.CurrentUICulture = uiCulture;
                        return own;
                    },
                    (i, _, own) =>
                    {
                        lines[i] = line(items[start + i]);
                        return own;
                    },
                    own =>
                    {
                        CultureInfo.CurrentCulture = own.Item1;
                        CultureInfo.CurrentUICulture = own.Item2;
                    });
                return lines;
            });
        }

        /// <summary>
        /// The batch's lines; a line that failed to build is thrown as its own exception, not wrapped by the task or the loop.
        /// </summary>
        private static string[] Result(Task<string[]> batch)
        {
            try
            {
                return batch.GetAwaiter().GetResult();
            }
            catch (AggregateException e) when (e.InnerExceptions.Count == 1)
            {
                ExceptionDispatchInfo.Capture(e.InnerException).Throw();
                throw; // not reached
            }
        }
    }
}
