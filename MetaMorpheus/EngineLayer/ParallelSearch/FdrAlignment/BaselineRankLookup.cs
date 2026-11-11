#nullable enable

using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.ParallelSearch.FdrAlignment;

/// <summary>
/// Baseline values in the order FDR ranked them, best first. A transient item is given the value of the
/// lowest-ranked baseline entry that still ranks at least as well; one that ranks above every baseline
/// entry gets the first, and one below every baseline entry gets the last.
/// </summary>
public class BaselineRankLookup<TKey, T>
{
    private readonly TKey[] _keys;
    private readonly T[] _values;
    private readonly IComparer<TKey> _betterFirst;

    /// <param name="betterFirst">Orders keys best first: negative when the first key ranks better.</param>
    public BaselineRankLookup(IEnumerable<(TKey Key, T Value)> baseline, IComparer<TKey> betterFirst)
    {
        ArgumentNullException.ThrowIfNull(baseline);
        ArgumentNullException.ThrowIfNull(betterFirst);

        // OrderBy is stable, so equal keys keep the order they were given in
        var sorted = baseline.OrderBy(entry => entry.Key, betterFirst).ToArray();
        _keys = sorted.Select(entry => entry.Key).ToArray();
        _values = sorted.Select(entry => entry.Value).ToArray();
        _betterFirst = betterFirst;
    }

    public int Count => _keys.Length;

    public bool IsEmpty => _keys.Length == 0;

    public T Find(TKey key)
    {
        if (IsEmpty)
            throw new InvalidOperationException("The baseline lookup is empty.");

        // count the baseline entries that rank at least as well as the key
        int low = 0;
        int high = _keys.Length;
        while (low < high)
        {
            int middle = (low + high) / 2;
            if (_betterFirst.Compare(_keys[middle], key) <= 0)
                low = middle + 1;
            else
                high = middle;
        }

        return _values[Math.Max(low - 1, 0)];
    }
}
