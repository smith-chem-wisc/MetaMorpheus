using System.Collections.Generic;

namespace TaskLayer
{
    /// <summary>
    /// In-memory hand-off between tasks in a run list (01_Architecture.md decision #1). A finishing
    /// task deposits a typed result keyed by its task id; a later task retrieves it by the upstream
    /// task id it was configured with. Used by <see cref="SearchTask"/>'s post-search analysis to pass
    /// the deduped proteoform-level SpectralMatch list to <see cref="TruncationSearchTask"/> without a
    /// file round-trip. When the lookup misses, the consumer falls back to parsing from disk.
    /// </summary>
    public class TaskChainContext
    {
        private readonly Dictionary<string, object> _results = new();
        private readonly List<string> _depositOrder = new();

        /// <summary>Store <paramref name="result"/> under <paramref name="taskId"/>, overwriting any prior value.</summary>
        public void Deposit<T>(string taskId, T result)
        {
            if (!_results.ContainsKey(taskId))
            {
                _depositOrder.Add(taskId);
            }

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
    }
}
