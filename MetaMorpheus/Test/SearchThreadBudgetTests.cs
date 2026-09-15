using EngineLayer;
using EngineLayer.GlycoSearch;
using EngineLayer.Util;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Test
{
    /// <summary>
    /// Locks in how a task's search threads are shared by the spectra files it searches at once.
    /// </summary>
    /// <remarks>
    /// Splitting the budget evenly once, before the search, left most of it idle at the end of a run: on six files the smaller
    /// files finished minutes before the largest, which kept searching on its sixth of the budget. <see cref="SearchThreadBudget"/>
    /// hands a finished file's threads to the files still searching. The share rules are tested one call at a time, so they do
    /// not depend on timing; <see cref="SearchThreadBudget.Run"/> is tested for handing out every item exactly once, never
    /// running more workers than the budget, and moving threads to a run that is still working.
    /// </remarks>
    [TestFixture]
    public static class SearchThreadBudgetTests
    {
        [Test]
        public static void OneFileTakesTheWholeBudgetAndNoMore()
        {
            var budget = new SearchThreadBudget(4);
            var share = budget.Join();

            for (int i = 0; i < 4; i++)
            {
                Assert.That(share.TryTakeThread(), Is.True, $"thread {i + 1} of 4");
            }
            Assert.That(share.TryTakeThread(), Is.False, "the budget is spent");
            Assert.That(budget.ThreadsInUse, Is.EqualTo(4));
            Assert.That(share.ShouldReturnThread(), Is.False, "nobody else wants a thread");
        }

        [Test]
        public static void ABudgetBelowOneIsOneThread()
        {
            Assert.That(new SearchThreadBudget(0).TotalThreads, Is.EqualTo(1));
            Assert.That(new SearchThreadBudget(-3).TotalThreads, Is.EqualTo(1));
        }

        [Test]
        public static void AFileThatStartsLaterGetsItsFairShareFromAFileAlreadySearching()
        {
            var budget = new SearchThreadBudget(8);
            var first = budget.Join();
            while (first.TryTakeThread()) { }
            Assert.That(first.Threads, Is.EqualTo(8));

            var second = budget.Join();
            Assert.That(second.TryTakeThread(), Is.False, "every thread is in use");
            Assert.That(first.ShouldReturnThread(), Is.True, "the first file is over its share of 4 while the second has none");

            // The first file's workers give threads back one at a time, after their current scan, and the waiting second file takes
            // each one. A thread is given back only while every thread is in use, so a free thread nobody has taken yet stops the
            // first file from giving back more than the second can use.
            int given = 0;
            while (first.ShouldReturnThread())
            {
                first.ReturnThread();
                Assert.That(first.ShouldReturnThread(), Is.False, "with a thread free, the first file waits for the second to take it");
                Assert.That(second.TryTakeThread(), Is.True, "the second file takes the thread given back");
                Assert.That(++given, Is.LessThanOrEqualTo(4), "gave back more than the second file's share");
            }
            Assert.That(first.Threads, Is.EqualTo(4));
            Assert.That(second.Threads, Is.EqualTo(4));
            Assert.That(second.TryTakeThread(), Is.False, "the budget is spent");
            Assert.That(first.ShouldReturnThread() || second.ShouldReturnThread(), Is.False, "both files are at their share");
        }

        [Test]
        public static void AFileThatHasHandedOutAllItsWorkGivesItsThreadsToTheFilesStillSearching()
        {
            var budget = new SearchThreadBudget(8);
            var finished = budget.Join();
            var searching = budget.Join();
            while (finished.TryTakeThread()) { }
            while (searching.TryTakeThread()) { }
            Assert.That((finished.Threads, searching.Threads), Is.EqualTo((4, 4)));

            finished.FinishHandingOutWork();
            Assert.That(finished.TryTakeThread(), Is.False, "a file with no work left takes no more threads");
            Assert.That(searching.ShouldReturnThread(), Is.False, "a file with no work left does not make others give threads back");

            for (int i = 0; i < 4; i++)
            {
                finished.ReturnThread();
            }
            Assert.That(budget.SearchingFiles, Is.EqualTo(1));
            while (searching.TryTakeThread()) { }
            Assert.That(searching.Threads, Is.EqualTo(8), "the file still searching gets the whole budget");
        }

        [Test]
        public static void ThreadsLeftOverFromDividingTheBudgetAreNotLeftIdle()
        {
            // 8 threads for 3 files is a share of 2 each, with 2 left over.
            var budget = new SearchThreadBudget(8);
            var a = budget.Join();
            var b = budget.Join();
            var c = budget.Join();
            foreach (var share in new[] { a, b, c })
            {
                Assert.That(share.TryTakeThread() && share.TryTakeThread(), Is.True);
            }

            Assert.That(a.TryTakeThread(), Is.True, "no file is below its share, so a spare thread goes to whoever asks");
            Assert.That(a.TryTakeThread(), Is.True);
            Assert.That(b.TryTakeThread(), Is.False, "the budget is spent");
            Assert.That(a.ShouldReturnThread(), Is.False, "a is over its share but nobody is below theirs");
            Assert.That(budget.ThreadsInUse, Is.EqualTo(8));
        }

        [Test]
        public static void WithMoreFilesThanThreadsEveryFileThatHasAThreadKeepsIt()
        {
            var budget = new SearchThreadBudget(2);
            var a = budget.Join();
            var b = budget.Join();
            var c = budget.Join();

            Assert.That(a.TryTakeThread(), Is.True);
            Assert.That(b.TryTakeThread(), Is.True);
            Assert.That(c.TryTakeThread(), Is.False);
            Assert.That(a.ShouldReturnThread() || b.ShouldReturnThread(), Is.False, "a file never gives up its only thread");

            a.FinishHandingOutWork();
            a.ReturnThread();
            Assert.That(c.TryTakeThread(), Is.True, "c gets the thread a finished with");
        }

        [Test]
        public static void WaitingForAThreadWakesWhenOneIsReturned()
        {
            var budget = new SearchThreadBudget(1);
            var holder = budget.Join();
            var waiter = budget.Join();
            Assert.That(holder.TryTakeThread(), Is.True);

            var returner = Task.Run(() =>
            {
                Thread.Sleep(100);
                holder.FinishHandingOutWork();
                holder.ReturnThread();
            });

            var timer = System.Diagnostics.Stopwatch.StartNew();
            while (!waiter.TryTakeThread())
            {
                waiter.WaitForChange(TimeSpan.FromSeconds(10));
                Assert.That(timer.Elapsed, Is.LessThan(TimeSpan.FromSeconds(10)), "the waiter was never woken");
            }
            returner.Wait();
            Assert.That(waiter.Threads, Is.EqualTo(1));
        }

        [Test]
        public static void RunHandsOutEveryItemExactlyOnce([Values(1, 3, 8)] int threads)
        {
            const int items = 5000;
            var visits = new int[items];
            new SearchThreadBudget(threads).Run(items, nextItem =>
            {
                for (int i = nextItem(); i >= 0; i = nextItem())
                {
                    Interlocked.Increment(ref visits[i]);
                }
            });

            Assert.That(visits.All(v => v == 1), Is.True, $"items visited other than once: {visits.Count(v => v != 1)}");
        }

        [Test]
        public static void RunWithNoItemsReturnsAndLeavesNoThreadInUse()
        {
            var budget = new SearchThreadBudget(4);
            bool workerRan = false;
            budget.Run(0, nextItem =>
            {
                workerRan = true;
                Assert.That(nextItem(), Is.EqualTo(-1));
            });
            Assert.That(budget.ThreadsInUse, Is.EqualTo(0));
            Assert.That(budget.SearchingFiles, Is.EqualTo(0));
            Assert.That(workerRan, Is.False, "no worker is started when there is nothing to hand out");
        }

        [Test]
        public static void RunNeverUsesMoreThreadsThanTheBudgetAndReturnsThemAll()
        {
            var budget = new SearchThreadBudget(3);
            int running = 0;
            int peak = 0;

            Parallel.For(0, 4, _ => budget.Run(300, nextItem =>
            {
                for (int i = nextItem(); i >= 0; i = nextItem())
                {
                    int now = Interlocked.Increment(ref running);
                    int seen;
                    while (now > (seen = Volatile.Read(ref peak)) && Interlocked.CompareExchange(ref peak, now, seen) != seen) { }
                    Thread.SpinWait(2000);
                    Interlocked.Decrement(ref running);
                }
            }));

            Assert.That(peak, Is.InRange(1, 3), "items were worked on by more threads at once than the budget");
            Assert.That(budget.ThreadsInUse, Is.EqualTo(0));
            Assert.That(budget.SearchingFiles, Is.EqualTo(0));
        }

        [Test]
        public static void ARunStillWorkingGetsTheThreadsOfARunThatFinished()
        {
            var budget = new SearchThreadBudget(8);
            int longRunWorkers = 0;
            int longRunPeakAfterShortRunFinished = 0;
            var shortRunDone = new ManualResetEventSlim();

            var shortRun = Task.Run(() =>
            {
                budget.Run(40, nextItem =>
                {
                    for (int i = nextItem(); i >= 0; i = nextItem())
                    {
                        Thread.Sleep(5);
                    }
                });
                shortRunDone.Set();
            });

            var longRun = Task.Run(() => budget.Run(1500, nextItem =>
            {
                int now = Interlocked.Increment(ref longRunWorkers);
                try
                {
                    for (int i = nextItem(); i >= 0; i = nextItem())
                    {
                        if (shortRunDone.IsSet)
                        {
                            int workers = Volatile.Read(ref longRunWorkers);
                            int seen;
                            while (workers > (seen = Volatile.Read(ref longRunPeakAfterShortRunFinished))
                                && Interlocked.CompareExchange(ref longRunPeakAfterShortRunFinished, workers, seen) != seen) { }
                        }
                        Thread.Sleep(2);
                    }
                }
                finally
                {
                    Interlocked.Decrement(ref longRunWorkers);
                }
            }));

            Assert.That(Task.WaitAll(new[] { shortRun, longRun }, TimeSpan.FromMinutes(2)), Is.True, "the runs did not finish");
            Assert.That(longRunPeakAfterShortRunFinished, Is.EqualTo(8), "after the short run finished, the long run should be searching on the whole budget");
        }

        /// <summary>
        /// Spectra files finish loading at different times, so a file usually starts searching after others hold every thread. Its
        /// workers must get threads while the others are still searching, not only after they finish.
        /// </summary>
        [Test]
        public static void ARunThatStartsLaterGetsItsShareFromARunAlreadyWorking()
        {
            var budget = new SearchThreadBudget(8);
            var laterRunStarted = new ManualResetEventSlim();
            var laterRunDone = new ManualResetEventSlim();
            int earlierWorkers = 0;
            int laterWorkers = 0;
            int laterPeakWhileEarlierWorking = 0;

            var earlier = Task.Run(() => budget.Run(100_000, nextItem =>
            {
                Interlocked.Increment(ref earlierWorkers);
                try
                {
                    // Keeps working until the later run has finished, so every thread the later run had came from this one.
                    for (int i = nextItem(); i >= 0 && !laterRunDone.IsSet; i = nextItem())
                    {
                        Thread.Sleep(1);
                    }
                }
                finally
                {
                    Interlocked.Decrement(ref earlierWorkers);
                }
            }));

            SpinWait.SpinUntil(() => budget.ThreadsInUse == 8, TimeSpan.FromSeconds(30));
            Assert.That(budget.ThreadsInUse, Is.EqualTo(8), "the earlier run should have taken the whole budget");

            var later = Task.Run(() =>
            {
                laterRunStarted.Set();
                budget.Run(400, nextItem =>
                {
                    int now = Interlocked.Increment(ref laterWorkers);
                    try
                    {
                        for (int i = nextItem(); i >= 0; i = nextItem())
                        {
                            int workers = Volatile.Read(ref laterWorkers);
                            int seen;
                            while (workers > (seen = Volatile.Read(ref laterPeakWhileEarlierWorking))
                                && Interlocked.CompareExchange(ref laterPeakWhileEarlierWorking, workers, seen) != seen) { }
                            Thread.Sleep(1);
                        }
                    }
                    finally
                    {
                        Interlocked.Decrement(ref laterWorkers);
                    }
                });
                laterRunDone.Set();
            });

            Assert.That(later.Wait(TimeSpan.FromMinutes(2)), Is.True, "the later run never finished: the earlier run kept every thread");
            Assert.That(earlier.Wait(TimeSpan.FromMinutes(2)), Is.True);
            Assert.That(laterPeakWhileEarlierWorking, Is.EqualTo(4), "the later run should reach its share of 4 while the earlier run is still working");
            Assert.That(budget.ThreadsInUse, Is.EqualTo(0));
        }

        [Test]
        public static void AWorkerExceptionEndsTheRunAndIsRethrown()
        {
            var budget = new SearchThreadBudget(4);
            var thrown = Assert.Throws<AggregateException>(() => budget.Run(1000, nextItem =>
            {
                for (int i = nextItem(); i >= 0; i = nextItem())
                {
                    if (i == 10)
                    {
                        throw new InvalidOperationException("scan 10 failed");
                    }
                }
            }));
            Assert.That(thrown.Flatten().InnerExceptions.OfType<InvalidOperationException>().Any(), Is.True);
            Assert.That(budget.ThreadsInUse, Is.EqualTo(0), "a failed worker still returns its thread");
        }

        /// <summary>
        /// The glyco engine searches through the budget it is given, and returns every thread when it is done.
        /// </summary>
        [Test]
        [NonParallelizable] // the glyco engine constructor writes process-wide glycan state (GlycanBox statics)
        public static void GlycoSearchEngineSearchesThroughTheBudgetItIsGiven()
        {
            var budget = new SearchThreadBudget(3);
            var commonParameters = new CommonParameters(dissociationType: MassSpectrometry.DissociationType.HCD, maxThreadsToUsePerFile: 3);
            var engine = new GlycoSearchEngine(new List<GlycoSpectralMatch>[0], new Ms2ScanWithSpecificMass[0],
                new List<PeptideWithSetModifications>(), null, null, 0, commonParameters, null, "OGlycan.gdb", null,
                glycoSearchType: GlycoSearchType.OGlycanSearch, 30, 3, false, new List<string>()) // not null: an engine-started handler another test left subscribed joins the ids
            {
                ThreadBudget = budget
            };

            engine.Run();

            Assert.That(budget.RunsStarted, Is.EqualTo(1), "the engine must search through the budget it was given");
            Assert.That(budget.ThreadsInUse, Is.EqualTo(0));
        }
    }
}
