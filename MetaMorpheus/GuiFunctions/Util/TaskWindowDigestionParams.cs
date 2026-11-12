using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace GuiFunctions.Util;

/// <summary>
/// The <see cref="DigestionParams"/> a protein task window saves: the values its controls show, and the loaded task's value
/// for every digestion setting it has no control for.
/// </summary>
/// <remarks>
/// <para>Each window builds new digestion parameters on save. Any constructor argument a window left out fell back to its
/// constructor default, so opening a task and pressing Save changed the search: every window turned
/// <see cref="DigestionParams.RespectCleavageBlockingModifications"/> off, the crosslink window reset
/// <see cref="DigestionParams.MaxModsForPeptide"/> to 2, and the calibration, GPTMD, glyco and crosslink windows reset
/// <see cref="DigestionParams.GeneratehUnlabeledProteinsForSilac"/>. A window now passes only what it shows; this class
/// carries the rest. Search mode and terminus are kept through <see cref="TaskWindowSearchMode.Preserve"/>.</para>
/// <para>RNA tasks are not handled here: the windows build <c>RnaDigestionParams</c> separately.</para>
/// </remarks>
public static class TaskWindowDigestionParams
{
    /// <summary>
    /// Builds the digestion parameters a window saves. The five required values are shown by every window. Leave an optional
    /// argument null when the window has no control for that setting; the loaded task's value is kept. A new task, or one
    /// without protein digestion parameters, keeps the constructor defaults.
    /// </summary>
    /// <param name="loaded">The task's digestion parameters when the window opened (TheTask.CommonParameters.DigestionParams).</param>
    /// <param name="protease">The protease the window shows; see <see cref="TaskWindowSearchMode.ProteaseToShow"/>.</param>
    /// <param name="searchMode">For a window with its own search-mode controls (Search, Glyco); null keeps the loaded mode and terminus.</param>
    public static DigestionParams ForSave(IDigestionParams loaded, string protease, int maxMissedCleavages, int minPeptideLength,
        int maxPeptideLength, int maxModificationIsoforms, int? maxModsForPeptides = null,
        InitiatorMethionineBehavior? initiatorMethionineBehavior = null,
        (CleavageSpecificity SearchModeType, FragmentationTerminus Terminus)? searchMode = null,
        bool? generateUnlabeledProteinsForSilac = null)
    {
        // Every setting the window does not show comes from here. Adding a DigestionParams constructor argument means adding
        // it below too; TaskWindowDigestionParamsTests.ASettingNoWindowShows_SurvivesASave fails until it is.
        DigestionParams kept = loaded as DigestionParams ?? new DigestionParams();
        var (searchModeType, terminus) = searchMode ?? TaskWindowSearchMode.Preserve(loaded);

        return new DigestionParams(
            protease: protease,
            maxMissedCleavages: maxMissedCleavages,
            minPeptideLength: minPeptideLength,
            maxPeptideLength: maxPeptideLength,
            maxModificationIsoforms: maxModificationIsoforms,
            initiatorMethionineBehavior: initiatorMethionineBehavior ?? kept.InitiatorMethionineBehavior,
            maxModsForPeptides: maxModsForPeptides ?? kept.MaxMods,
            searchModeType: searchModeType,
            fragmentationTerminus: terminus,
            generateUnlabeledProteinsForSilac: generateUnlabeledProteinsForSilac ?? kept.GeneratehUnlabeledProteinsForSilac,
            keepNGlycopeptide: kept.KeepNGlycopeptide,
            keepOGlycopeptide: kept.KeepOGlycopeptide,
            respectCleavageBlockingModifications: kept.RespectCleavageBlockingModifications,
            respectCleavagePromotingModifications: kept.RespectCleavagePromotingModifications);
    }
}
