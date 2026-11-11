#nullable enable

using System.Collections.Generic;

namespace EngineLayer.ParallelSearch.FdrAlignment;

public interface IBaselineFdrAlignmentService<TItem>
{
    void BuildBaselineCache(IEnumerable<TItem> baselineItems);

    AlignmentStats ApplyBaseline(IList<TItem> items);
}
