using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace TaskLayer
{
    /// <summary>
    /// Refuses a task whose digestion settings ask for seed peptides when the task cannot use seeds.
    /// </summary>
    /// <remarks>
    /// <para><b>Peptides versus seeds.</b> In mzLib, <c>SearchModeType</c> Full gives fully specific peptides, and Semi with
    /// <c>FragmentationTerminus</c> Both gives semi-specific peptides (since mzLib #1303). Semi with N or C, and None with
    /// any terminus, give seeds instead: long stretches fixed at one terminus whose other end the search engine decides
    /// afterwards, from the precursor mass. See <c>DigestionParams.SearchModeType</c> in mzLib for the full table.</para>
    /// <para><b>Who can use seeds.</b> Only the non-specific search engine, which a Search task runs for
    /// <see cref="SearchType.NonSpecific"/> (it makes its own N and C passes). Classic and Modern search, Glyco,
    /// crosslink, GPTMD and calibration score digestion products as they are, so given seeds they finish normally and
    /// report far fewer, wrong identifications. Refusing with a message is better than that silent wrong answer.</para>
    /// <para>Called by <see cref="EverythingRunnerEngine"/> for every task before any runs, and by
    /// <see cref="MetaMorpheusTask.RunTask"/> as a backstop.</para>
    /// </remarks>
    public static class DigestionSearchModeCheck
    {
        /// <summary>
        /// Whether these settings make digestion return seeds rather than peptides: SearchModeType None with any terminus,
        /// or Semi with terminus N or C. The refusal below and the task windows' warnings both use this one rule.
        /// </summary>
        public static bool AsksForSeeds(DigestionParams digestionParams) =>
            digestionParams.SearchModeType == CleavageSpecificity.None
            || (digestionParams.SearchModeType == CleavageSpecificity.Semi && digestionParams.FragmentationTerminus is FragmentationTerminus.N or FragmentationTerminus.C);

        /// <summary>
        /// The message explaining why this task cannot run with its digestion settings, or null when it can.
        /// </summary>
        /// <param name="task">The task to check.</param>
        /// <param name="taskName">The name the user knows the task by, included in the message when given.</param>
        public static string GetRefusal(MetaMorpheusTask task, string taskName = null)
        {
            if (task?.CommonParameters?.DigestionParams is not DigestionParams digestionParams)
            {
                return null; // RNA digestion has no protein seed request
            }
            if (task is SearchTask { SearchParameters.SearchType: SearchType.NonSpecific } || !AsksForSeeds(digestionParams))
            {
                return null; // the non-specific search engine trims seeds; every other task gets peptides
            }

            string which = taskName == null ? "This task" : $"Task \"{taskName}\"";
            if (digestionParams.SearchModeType == CleavageSpecificity.None)
            {
                return $"Cannot proceed. {which} has SearchModeType None (non-specific), which gives seed peptides that only the non-specific search can use. " +
                       "Use a Search task with the non-specific search type, or choose fully or semi-specific digestion.";
            }
            return $"Cannot proceed. {which} has SearchModeType Semi with FragmentationTerminus {digestionParams.FragmentationTerminus}, which gives seed peptides that only the non-specific search can use. " +
                   "For semi-specific peptides, set FragmentationTerminus Both.";
        }
    }
}
