using System.Collections.Generic;

namespace TaskLayer
{
    /// <summary>
    /// In-memory hand-off between tasks in a run list (docs/Truncation-Search.md decision #1). A finishing
    /// task deposits a typed result keyed by its task id; a later task retrieves it by the upstream
    /// task id it was configured with. Used by <see cref="SearchTask"/> to pass its PSM-level
    /// SpectralMatch list — as searched, neither deduped nor reduced to proteoform level — to
    /// <see cref="TruncationSearchTask"/> without a file round-trip; the consumer does the dedup to
    /// proteoform level itself. When the lookup misses, the consumer falls back to parsing from disk.
    ///
    /// The context is handed out by <see cref="EverythingRunnerEngine"/> only when the run list actually
    /// contains a consumer, and <see cref="Clear"/>ed once that consumer has run, so a deposited result
    /// is not kept alive for the remainder of a run that will never read it.
    /// </summary>
    public class TaskChainContext
    {
        private readonly Dictionary<string, object> _results = new();
        private readonly List<string> _depositOrder = new();

        /// <summary>Store <paramref name="result"/> under <paramref name="taskId"/>, overwriting any prior value.</summary>
        public void Deposit<T>(string taskId, T result)
        {
            // Refresh recency on overwrite: a re-deposit under an existing id must count as the most
            // recent for TryGetMostRecent. Remove-then-add keeps _depositOrder in true deposit order.
            // (List.Remove is a no-op when the id isn't present yet.)
            _depositOrder.Remove(taskId);
            _depositOrder.Add(taskId);

            _results[taskId] = result;
        }

        /// <summary>
        /// Retrieve a previously deposited result of type <typeparamref name="T"/> for
        /// <paramref name="taskId"/>. Returns false (and default) when absent or of a different type.
        /// </summary>
        public bool TryGet<T>(string taskId, out T result)
        {
            if (_results.TryGetValue(taskId, out var stored) && stored is T typed)
            {
                result = typed;
                return true;
            }

            result = default;
            return false;
        }

        /// <summary>
        /// Retrieve the most recently deposited result assignable to <typeparamref name="T"/>, regardless
        /// of task id. Lets a consumer (e.g. <see cref="TruncationSearchTask"/>) pick up the preceding
        /// task's output in a run list without being configured with that task's exact id.
        /// </summary>
        public bool TryGetMostRecent<T>(out T result)
        {
            for (int i = _depositOrder.Count - 1; i >= 0; i--)
            {
                if (_results[_depositOrder[i]] is T typed)
                {
                    result = typed;
                    return true;
                }
            }

            result = default;
            return false;
        }

        /// <summary>
        /// Drops every deposited result, releasing the object graphs they pin. The runner calls this once a
        /// consuming task has run, so a search's PSM set is not held for the rest of the run list.
        /// </summary>
        public void Clear()
        {
            _results.Clear();
            _depositOrder.Clear();
        }
    }
}
