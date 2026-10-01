#nullable enable
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DatabaseLoading;

/// <summary>
/// One encoding of "is there a spectral library in this database list", shared by the tasks, the runner
/// and the loaders. It lives beside <see cref="DbForTask"/> rather than in EngineLayer/Util so that every
/// file able to name the type can already see it, without a new using.
/// </summary>
public static class DbForTaskExtensions
{
    /// <summary>
    /// The precondition LoadSpectralLibraries reads before it opens anything, so a caller can test for a
    /// library without paying for the byte-offset index the SpectralLibrary constructor builds.
    /// </summary>
    public static bool AnySpectralLibrary(this IEnumerable<DbForTask> dbFilenameList)
        => dbFilenameList.Any(p => p.IsSpectralLibrary);

    public static IEnumerable<DbForTask> SpectralLibraries(this IEnumerable<DbForTask> dbFilenameList)
        => dbFilenameList.Where(p => p.IsSpectralLibrary);
}
