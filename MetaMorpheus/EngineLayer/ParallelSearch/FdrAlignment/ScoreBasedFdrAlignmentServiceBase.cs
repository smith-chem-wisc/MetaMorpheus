#nullable enable

using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.ParallelSearch.FdrAlignment;

public abstract class ScoreBasedFdrAlignmentServiceBase<TItem, TBaselineEntry> : IBaselineFdrAlignmentService<TItem>
    where TBaselineEntry : struct, IScoreBaselineEntry
{
    private readonly List<TBaselineEntry> _baselineLookup = [];

    public IReadOnlyList<TBaselineEntry> BaselineLookup => _baselineLookup;

    public bool HasBaselineCache => _baselineLookup.Count > 0;

    public void BuildBaselineCache(IEnumerable<TItem> baselineItems)
    {
        ArgumentNullException.ThrowIfNull(baselineItems);

        _baselineLookup.Clear();

        foreach (var item in baselineItems.Where(p => p is not null).OrderByDescending(GetScore))
        {
            if (TryBuildBaselineEntry(item, out TBaselineEntry entry))
            {
                _baselineLookup.Add(entry);
            }
        }
    }

    public AlignmentStats ApplyBaseline(IList<TItem> items)
    {
        ArgumentNullException.ThrowIfNull(items);

        if (items.Count == 0 || _baselineLookup.Count == 0)
        {
            return AlignmentStats.Empty;
        }

        int baselineIndex = 0;
        int lastBaselineIndex = _baselineLookup.Count - 1;
        double highestBaselineScore = _baselineLookup[0].Score;
        double lowestBaselineScore = _baselineLookup[lastBaselineIndex].Score;
        int clampedHighCount = 0;
        int clampedLowCount = 0;

        foreach (var item in items)
        {
            int selectedIndex;
            double transientScore = GetScore(item);

            if (transientScore > highestBaselineScore)
            {
                selectedIndex = 0;
                clampedHighCount++;
            }
            else
            {
                while (baselineIndex < lastBaselineIndex && _baselineLookup[baselineIndex + 1].Score >= transientScore)
                {
                    baselineIndex++;
                }

                selectedIndex = baselineIndex;
                if (transientScore < lowestBaselineScore)
                {
                    clampedLowCount++;
                }
            }

            ApplyBaselineEntry(item, _baselineLookup[selectedIndex]);
        }

        return new AlignmentStats(items.Count, clampedHighCount, clampedLowCount);
    }

    protected abstract double GetScore(TItem item);

    protected abstract bool TryBuildBaselineEntry(TItem item, out TBaselineEntry entry);

    protected abstract void ApplyBaselineEntry(TItem item, TBaselineEntry entry);
}
