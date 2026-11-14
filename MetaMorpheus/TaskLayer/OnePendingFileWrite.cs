using System;
using System.Threading;
using System.Threading.Tasks;

namespace TaskLayer
{
    /// <summary>
    /// Writes a task's exported spectra files on a background thread, one at a time, so the task can go on to its
    /// next spectra file while the last one is still being written. Writing a calibrated or averaged mzML takes 10-30 s
    /// on one core, and calibration and averaging used to wait for each before loading the next file.
    /// </summary>
    /// <remarks>
    /// At most one write is pending: <see cref="Start"/> first waits for the write before it, so one extra file is held
    /// in memory at most. The bookkeeping that reports a file - adding it to the task's results, raising
    /// FinishedWritingFile and FinishedDataFile - is handed in as <c>afterWrite</c> and runs on the task's own thread
    /// only once the file is complete, in the order the files were started. A failed write is surfaced by the next
    /// <see cref="Start"/> or by <see cref="Finish"/>, which rethrow its exception unless the caller handles failures.
    /// The bytes written are the writer's own; only when they are written moves.
    /// </remarks>
    public sealed class OnePendingFileWrite : IDisposable
    {
        private Task _write;
        private Action _afterWrite;
        private Action<Exception> _onWriteFailed;

        /// <summary>
        /// Starts a write on its own thread. A test replaces this to hold a write open.
        /// </summary>
        internal Func<Action, Task> StartWrite = write =>
            Task.Factory.StartNew(write, CancellationToken.None, TaskCreationOptions.LongRunning, TaskScheduler.Default);

        public bool IsWritePending => _write != null;

        /// <summary>
        /// Waits for the pending write and runs its bookkeeping, then starts <paramref name="write"/> in the background.
        /// With no write, <paramref name="afterWrite"/> runs at once, still after the pending file's bookkeeping.
        /// </summary>
        /// <param name="afterWrite"> Runs on the calling thread once the file is complete: at the next Start or at Finish. </param>
        /// <param name="onWriteFailed">
        /// Handles a failed write, after which <paramref name="afterWrite"/> still runs. Without it the next Start or Finish
        /// rethrows the write's exception and <paramref name="afterWrite"/> does not run.
        /// </param>
        public void Start(Action write, Action afterWrite, Action<Exception> onWriteFailed = null)
        {
            Finish();

            if (write == null)
            {
                afterWrite?.Invoke();
                return;
            }

            _afterWrite = afterWrite;
            _onWriteFailed = onWriteFailed;
            _write = StartWrite(write);
        }

        /// <summary>
        /// Waits for the pending write, if there is one, and runs its bookkeeping or surfaces its failure.
        /// </summary>
        public void Finish()
        {
            if (_write == null)
            {
                return;
            }

            Task write = _write;
            Action afterWrite = _afterWrite;
            Action<Exception> onWriteFailed = _onWriteFailed;
            _write = null;
            _afterWrite = null;
            _onWriteFailed = null;

            try
            {
                write.GetAwaiter().GetResult();
            }
            catch (Exception e) when (onWriteFailed != null)
            {
                onWriteFailed(e);
            }

            afterWrite?.Invoke();
        }

        /// <summary>
        /// Waits for a pending write without reporting it or rethrowing its failure: for a task leaving by an exception
        /// of its own, so no write carries on behind it and its own exception is the one that propagates.
        /// </summary>
        public void Dispose()
        {
            Task write = _write;
            _write = null;
            _afterWrite = null;
            _onWriteFailed = null;

            if (write != null)
            {
                try
                {
                    write.Wait();
                }
                catch
                {
                    // observed; the task's own exception is already on its way out
                }
            }
        }
    }
}
