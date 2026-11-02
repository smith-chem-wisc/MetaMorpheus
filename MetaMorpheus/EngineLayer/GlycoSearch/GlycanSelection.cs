using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace EngineLayer.GlycoSearch
{
    /// <summary>
    /// What a saved glycan selection (<c>GlycoSearchParameters.SelectedGlycans</c>) does to one glycan
    /// database: which glycans get searched, and what the user has to be told about it.
    /// </summary>
    /// <remarks>
    /// One place for the rule, so the engine that applies it and the task that reports it cannot drift
    /// apart. The engine is built once per partition per spectra file, so it applies the selection but
    /// leaves the telling to the task, which does it once per run.
    /// </remarks>
    public sealed class GlycanSelection
    {
        /// <summary>
        /// Past this many, the prose gives the count of selected entries without listing them.
        /// </summary>
        public const int MaxEntriesListedInProse = 20;

        private GlycanSelection(string databaseFileName, Glycan[] glycans, int databaseEntries, List<string> found, List<string> missing)
        {
            DatabaseFileName = databaseFileName;
            Glycans = glycans;
            DatabaseEntries = databaseEntries;
            Found = found;
            Missing = missing;
        }

        public string DatabaseFileName { get; }

        /// <summary>
        /// The glycans to search: the selected ones, or the whole database when the selection names none
        /// of it, or names only entries it no longer holds.
        /// </summary>
        public Glycan[] Glycans { get; }

        /// <summary>
        /// Distinct entries (IdWithMotif) in the database, which is what the task window lists and counts.
        /// </summary>
        public int DatabaseEntries { get; }

        /// <summary>
        /// Selected entries the database holds, in the order the selection names them.
        /// </summary>
        public IReadOnlyList<string> Found { get; }

        /// <summary>
        /// Selected entries the database no longer holds, in the order the selection names them.
        /// </summary>
        public IReadOnlyList<string> Missing { get; }

        /// <summary>
        /// The selection named nothing in this database, so it is searched whole and there is nothing to say.
        /// </summary>
        public bool IsWholeDatabaseByDefault => Found.Count == 0 && Missing.Count == 0;

        /// <summary>
        /// The selection named only entries the database no longer holds, so it is searched whole instead.
        /// </summary>
        public bool FellBackToWholeDatabase => Found.Count == 0 && Missing.Count > 0;

        /// <param name="loaded">Every glycan in the database, as loaded.</param>
        /// <param name="databaseFileName">The database's file name, which is how the selection names it.</param>
        /// <param name="selectedGlycans">(database file name, glycan IdWithMotif) pairs; null or empty means none.</param>
        public static GlycanSelection Apply(Glycan[] loaded, string databaseFileName, IEnumerable<(string, string)> selectedGlycans)
        {
            var wanted = (selectedGlycans ?? Enumerable.Empty<(string, string)>())
                .Where(s => s.Item1 == databaseFileName)
                .Select(s => s.Item2)
                .Distinct()
                .ToList();

            var inDatabase = new HashSet<string>(loaded.Select(g => g.IdWithMotif));
            var found = wanted.Where(inDatabase.Contains).ToList();
            var missing = wanted.Where(id => !inDatabase.Contains(id)).ToList();

            // Nothing selected here, or nothing selected that is still here: the whole database. The second
            // beats handing BuildOGlycanBoxes an empty array, whose failure surfaces much later and elsewhere,
            // but it is a different search from the one asked for, so Warning says so.
            var foundSet = new HashSet<string>(found);
            var glycans = found.Count == 0 ? loaded : loaded.Where(g => foundSet.Contains(g.IdWithMotif)).ToArray();

            return new GlycanSelection(databaseFileName, glycans, inDatabase.Count, found, missing);
        }

        /// <summary>
        /// What to add after the database's file name in the methods prose: empty when the whole database
        /// was searched because nothing in it was selected, which leaves the prose as it always was.
        /// </summary>
        public string ProseSuffix()
        {
            if (IsWholeDatabaseByDefault)
            {
                return string.Empty;
            }

            if (FellBackToWholeDatabase)
            {
                return string.Format(CultureInfo.InvariantCulture,
                    " (whole database searched: none of the {0} selected {1} is in it: {2})",
                    Missing.Count, Entries(Missing.Count), string.Join(", ", Missing));
            }

            string suffix = string.Format(CultureInfo.InvariantCulture, " ({0} of {1} entries selected", Found.Count, DatabaseEntries);
            if (Found.Count <= MaxEntriesListedInProse)
            {
                suffix += ": " + string.Join(", ", Found);
            }
            if (Missing.Count > 0)
            {
                suffix += string.Format(CultureInfo.InvariantCulture,
                    "; {0} more selected {1} not in the database and skipped: {2}",
                    Missing.Count, Missing.Count == 1 ? "entry was" : "entries were", string.Join(", ", Missing));
            }
            return suffix + ")";
        }

        /// <summary>
        /// A warning for the user when the selection names entries the database no longer holds; null otherwise.
        /// </summary>
        /// <param name="kind">"O-glycan" or "N-glycan", for the message.</param>
        public string Warning(string kind)
        {
            if (Missing.Count == 0)
            {
                return null;
            }

            if (FellBackToWholeDatabase)
            {
                return string.Format(CultureInfo.InvariantCulture,
                    "None of the {0} {1} selected from the {2} database '{3}' {4} in it any more ({5}), so the whole database is searched instead.",
                    Missing.Count, Entries(Missing.Count), kind, DatabaseFileName, Missing.Count == 1 ? "is" : "are", string.Join(", ", Missing));
            }

            return string.Format(CultureInfo.InvariantCulture,
                "{0} {1} selected from the {2} database '{3}' {4} no longer in it and {5} skipped: {6}. The other {7} selected {8} searched.",
                Missing.Count, Entries(Missing.Count), kind, DatabaseFileName, Missing.Count == 1 ? "is" : "are",
                Missing.Count == 1 ? "was" : "were", string.Join(", ", Missing), Found.Count, Found.Count == 1 ? "entry is" : "entries are");
        }

        private static string Entries(int count) => count == 1 ? "entry" : "entries";
    }
}
