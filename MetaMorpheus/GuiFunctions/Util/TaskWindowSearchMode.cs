using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace GuiFunctions.Util;

/// <summary>
/// What a task window loads and saves for <see cref="DigestionParams.SearchModeType"/> and
/// <see cref="DigestionParams.FragmentationTerminus"/>, and which protease it shows.
/// </summary>
/// <remarks>
/// <para>Those two settings together decide whether digestion returns peptides or seeds. Full gives fully specific
/// peptides; Semi with terminus Both gives semi-specific peptides; Semi or None with terminus N or C gives seeds that only
/// the non-specific search engine can use. See <c>DigestionParams.SearchModeType</c> in mzLib for the full table.</para>
/// <para>Windows that build a <see cref="DigestionParams"/> without passing them reset every task to Full and Both on
/// save, which silently turned a semi-specific task into a fully specific one. The Search window has its own semi- and
/// non-specific controls and does not use this class.</para>
/// </remarks>
public static class TaskWindowSearchMode
{
    /// <summary>Whether a window with a semi-specific choice should show it checked for this loaded task.</summary>
    public static bool IsSemiSpecific(IDigestionParams loaded) => loaded?.SearchModeType == CleavageSpecificity.Semi;

    /// <summary>
    /// The search mode and terminus a window with a semi-specific choice (the Glyco window) saves. The terminus is always
    /// Both, which asks for peptides: searches other than the non-specific search cannot use seeds.
    /// </summary>
    public static (CleavageSpecificity SearchModeType, FragmentationTerminus Terminus) ForSemiSpecificChoice(bool semiSpecific) =>
        (semiSpecific ? CleavageSpecificity.Semi : CleavageSpecificity.Full, FragmentationTerminus.Both);

    /// <summary>
    /// The search mode and terminus a window without a semi-specific choice (crosslink, GPTMD, calibration) saves: exactly
    /// what the loaded task had, so that opening and saving a task never changes which peptides it searches. A new task,
    /// or one without protein digestion parameters, gets the defaults. A combination the task cannot use is refused when
    /// the task is run, with a message, rather than silently changed here.
    /// </summary>
    public static (CleavageSpecificity SearchModeType, FragmentationTerminus Terminus) Preserve(IDigestionParams loaded) =>
        loaded is DigestionParams digestionParams
            ? (digestionParams.SearchModeType, digestionParams.FragmentationTerminus)
            : (CleavageSpecificity.Full, FragmentationTerminus.Both);

    /// <summary>
    /// The protease a window shows and saves for a loaded task: the one the user chose. That is
    /// <see cref="DigestionParams.SpecificProtease"/>, because for a non-specific search <see cref="DigestionParams.Protease"/>
    /// is singleN or singleC.
    /// </summary>
    public static Protease ProteaseToShow(DigestionParams loaded) => loaded.SpecificProtease;
}
