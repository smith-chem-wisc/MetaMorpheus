using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;

namespace EngineLayer.Util
{
    /// <summary>
    /// A task's search threads, shared by the spectra files it searches at the same time.
    /// </summary>
    /// <remarks>
    /// Splitting the budget evenly once, before the search (<see cref="FileParallelism.Decide"/>), leaves most of it idle at the end of
    /// a run: files finish at different times, and the last one keeps searching on its own fraction of the budget. Here every file
    /// that is searching joins the budget and takes threads as it can use them. A file is entitled to an equal share of the budget
    /// among the files still handing out work; threads no file below its share is asking for go to whoever asks, and a file over its
    /// share gives a thread back, after the scan it is on, only when another file is below its share and every thread is in use. So
    /// a file that finishes hands its threads to the files still searching, and a file that starts later gets its share from the
    /// files already searching.
    /// <para>
    /// Which thread searches which scan does not change the results: each scan index is handed out exactly once and its results go
    /// to that scan's own slot.
    /// </para>
    /// </remarks>
    public sealed class SearchThreadBudget
    {
        /// <summary>
        /// How long a run waits for a thread before looking again, in case a change was missed. Returning a thread and finishing a
        /// file both wake waiting runs at once, so this is only a safety net.
        /// </summary>
        private static readonly TimeSpan WaitBeforeLookingAgain = TimeSpan.FromMilliseconds(50);

        private readonly object _lock = new object();
        private readonly List<Share> _shares = new List<Share>();
        private int _threadsInUse;
        private int _runsStarted;

        /// <summary>
        /// Starts one worker on its own thread. A test replaces this to make starting a thread fail, which is otherwise only
        /// reachable when the machine is out of threads or memory.
        /// </summary>
        internal Func<Action, Task> StartWorkerThread = body =>
            Task.Factory.StartNew(body, CancellationToken.None, TaskCreationOptions.LongRunning, TaskScheduler.Default);

        /// <param name="totalThreads"> The threads every file together may use; below 1 means 1. </param>
        public SearchThreadBudget(int totalThreads)
        {
            TotalThreads = Math.Max(1, totalThreads);
        }

        public int TotalThreads { get; }

        public int ThreadsInUse
        {
            get { lock (_lock) { return _threadsInUse; } }
        }

        /// <summary>
        /// Files that have joined and still have work to hand out.
        /// </summary>
        public int SearchingFiles
        {
            get { lock (_lock) { return CountSearching(); } }
        }

        /// <summary>
        /// Calls to <see cref="Run"/> so far, so a test can see that a search went through this budget.
        /// </summary>
        internal int RunsStarted => Volatile.Read(ref _runsStarted);

        /// <summary>
        /// A file starts handing out work.
        /// </summary>
        public Share Join()
        {
            var share = new Share(this);
            lock (_lock)
            {
                _shares.Add(share);
            }
            return share;
        }

        /// <summary>
        /// Hands the items 0 to <paramref name="itemCount"/> - 1 out, each exactly once, to workers started while this budget grants
        /// threads, and returns when every worker has finished.
        /// </summary>
        /// <param name="worker">
        /// Runs on its own thread. It calls the function it is given for its next item until that returns -1, which means the items
        /// are all handed out or this worker's thread is wanted by another file. Anything a worker allocates for itself is allocated
        /// once per worker, before its first item.
        /// </param>
        /// <exception cref="AggregateException"> A worker threw; no further workers are started, and every thread is returned. </exception>
        public void Run(int itemCount, Action<Func<int>> worker)
        {
            Interlocked.Increment(ref _runsStarted);
            Share share = Join();
            int lastHandedOut = -1;
            int failed = 0;

            int NextItem()
            {
                if (Volatile.Read(ref failed) != 0 || share.ShouldReturnThread())
                {
                    return -1;
                }
                int item = Interlocked.Increment(ref lastHandedOut);
                if (item >= itemCount)
                {
                    // Keep the counter from creeping further past the end as idle workers keep asking.
                    Interlocked.Exchange(ref lastHandedOut, itemCount);
                    return -1;
                }
                return item;
            }

            var workers = new List<Task>();
            try
            {
                while (Volatile.Read(ref lastHandedOut) < itemCount - 1 && Volatile.Read(ref failed) == 0 && !GlobalVariables.StopLoops)
                {
                    if (!share.TryTakeThread())
                    {
                        share.WaitForChange(WaitBeforeLookingAgain);
                        continue;
                    }

                    try
                    {
                        workers.Add(StartWorkerThread(() =>
                        {
                            try
                            {
                                worker(NextItem);
                            }
                            catch
                            {
                                Interlocked.Exchange(ref failed, 1);
                                throw;
                            }
                            finally
                            {
                                share.ReturnThread();
                            }
                        }));
                    }
                    catch
                    {
                        // The thread was taken above and the worker whose finally would return it never ran, so return it here.
                        // Leaving it out would count it as in use for the rest of the task, shrinking the budget every file
                        // still to be searched draws on.
                        share.ReturnThread();
                        throw;
                    }
                }
            }
            finally
            {
                share.FinishHandingOutWork();
                try
                {
                    Task.WaitAll(workers.ToArray());
                }
                finally
                {
                    share.Leave();
                }
            }
        }

        /// <summary>
        /// Must be called under the lock.
        /// </summary>
        private int CountSearching()
        {
            int searching = 0;
            foreach (Share share in _shares)
            {
                if (share.HasWork)
                {
                    searching++;
                }
            }
            return searching;
        }

        /// <summary>
        /// A file's equal part of the budget among the files still handing out work; never below one thread. Must be called under the lock.
        /// </summary>
        private int FairShare()
        {
            return Math.Max(1, TotalThreads / Math.Max(1, CountSearching()));
        }

        /// <summary>
        /// Whether a file other than <paramref name="asking"/> has work and fewer threads than its share. Must be called under the lock.
        /// </summary>
        private bool AnotherFileIsBelowItsShare(Share asking, int fairShare)
        {
            foreach (Share share in _shares)
            {
                if (share != asking && share.HasWork && share.HeldThreads < fairShare)
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// One file's hold on the budget.
        /// </summary>
        public sealed class Share
        {
            private readonly SearchThreadBudget _budget;

            internal Share(SearchThreadBudget budget)
            {
                _budget = budget;
                HasWork = true;
            }

            /// <summary> Under the budget's lock. </summary>
            internal int HeldThreads { get; private set; }

            /// <summary> Under the budget's lock. </summary>
            internal bool HasWork { get; private set; }

            public int Threads
            {
                get { lock (_budget._lock) { return HeldThreads; } }
            }

            /// <summary>
            /// Takes a thread if one is free and this file is below its share, or if no other file below its share wants it.
            /// </summary>
            public bool TryTakeThread()
            {
                lock (_budget._lock)
                {
                    if (!HasWork || _budget._threadsInUse >= _budget.TotalThreads)
                    {
                        return false;
                    }
                    int fairShare = _budget.FairShare();
                    if (HeldThreads >= fairShare && _budget.AnotherFileIsBelowItsShare(this, fairShare))
                    {
                        return false;
                    }
                    HeldThreads++;
                    _budget._threadsInUse++;
                    return true;
                }
            }

            /// <summary>
            /// True when this file is over its share, every thread is in use, and another file is below its share: one of this file's
            /// workers should stop after its current item and return its thread.
            /// </summary>
            public bool ShouldReturnThread()
            {
                lock (_budget._lock)
                {
                    int fairShare = _budget.FairShare();
                    return HeldThreads > fairShare
                        && _budget._threadsInUse >= _budget.TotalThreads
                        && _budget.AnotherFileIsBelowItsShare(this, fairShare);
                }
            }

            public void ReturnThread()
            {
                lock (_budget._lock)
                {
                    if (HeldThreads == 0)
                    {
                        throw new InvalidOperationException("Returned a thread this file does not hold.");
                    }
                    HeldThreads--;
                    _budget._threadsInUse--;
                    Monitor.PulseAll(_budget._lock);
                }
            }

            /// <summary>
            /// This file has handed out all its work: it takes no more threads and no longer counts toward anyone's share. The
            /// threads it holds return as its workers finish their last items.
            /// </summary>
            public void FinishHandingOutWork()
            {
                lock (_budget._lock)
                {
                    HasWork = false;
                    Monitor.PulseAll(_budget._lock);
                }
            }

            /// <summary>
            /// Waits until a thread is returned or a file finishes handing out work, or until <paramref name="timeout"/> passes.
            /// </summary>
            public void WaitForChange(TimeSpan timeout)
            {
                lock (_budget._lock)
                {
                    Monitor.Wait(_budget._lock, timeout);
                }
            }

            /// <summary>
            /// Removes a finished file from the budget once its threads are back.
            /// </summary>
            internal void Leave()
            {
                lock (_budget._lock)
                {
                    HasWork = false;
                    _budget._shares.Remove(this);
                    Monitor.PulseAll(_budget._lock);
                }
            }
        }
    }
}
