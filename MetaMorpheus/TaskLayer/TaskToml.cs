using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using EngineLayer;
using EngineLayer.SpectrumMatch;
using Nett;

namespace TaskLayer
{
    /// <summary>
    /// What a task TOML's [CommonParameters] table means for the identification filter and its threshold, read once the
    /// whole file is read (<see cref="MetaMorpheusTask.ReadTaskToml{TTask}(string)"/>): the filter an unnamed file gets
    /// depends on its digestion parameters, which may come after the filter's keys, and a file saved before
    /// FilterThreshold existed holds other keys, which CommonParameters has no property for.
    /// </summary>
    /// <remarks>
    /// The mapping, in order:
    /// <list type="number">
    /// <item><description>FilterThreshold present: that threshold, with the filter the file names, or the default for its
    /// analyte when it names none.</description></item>
    /// <item><description>No FilterThreshold, and the file holds an old key (QValueThreshold, PepQValueThreshold, or the
    /// output filters of a file older still, <see cref="LegacyOutputFilterKeys"/>): the file keeps what it meant. A named
    /// filter is held to its old threshold, PepQValueThreshold for the PEP q-value and QValueThreshold otherwise. Naming none,
    /// the PEP q-value at PepQValueThreshold when that is below QValueThreshold (the old GUI's PEP box saved
    /// PepQValueThreshold 0.01 and QValueThreshold 1.0; the search says so, <see cref="LegacySentence"/>), and the q-value
    /// notch at QValueThreshold otherwise. An old key the file lacks means its old default
    /// (<see cref="LegacyQValueThresholdDefault"/>, <see cref="LegacyPepQValueThresholdDefault"/>).</description></item>
    /// <item><description>None of these keys: the filter the file names, or the default for its analyte, at
    /// <see cref="CommonParameters.DefaultFilterThreshold"/>.</description></item>
    /// </list>
    /// Saving the task writes IdentificationFilter and FilterThreshold, so the old keys are read once per old file.
    /// </remarks>
    public static class TaskToml
    {
        /// <summary>The key a task TOML saved before FilterThreshold existed held the q-value threshold under.</summary>
        public const string LegacyQValueThresholdKey = "QValueThreshold";

        /// <summary>The key a task TOML saved before FilterThreshold existed held the PEP q-value threshold under.</summary>
        public const string LegacyPepQValueThresholdKey = "PepQValueThreshold";

        /// <summary>
        /// The keys a task TOML saved before QValueThreshold and PepQValueThreshold existed (MetaMorpheus before June 2023,
        /// #2285) held its output filters under. Their values filtered only what was written, and have been ignored since;
        /// they still are. They mark the file as an old task file, which keeps the q-value notch at 0.01 it has searched
        /// with since.
        /// </summary>
        public static readonly IReadOnlyList<string> LegacyOutputFilterKeys = new[] { "QValueOutputFilter", "PepQValueOutputFilter" };

        /// <summary>What an old task TOML without QValueThreshold meant by it: its old default.</summary>
        public const double LegacyQValueThresholdDefault = 0.01;

        /// <summary>What an old task TOML without PepQValueThreshold meant by it: its old default.</summary>
        public const double LegacyPepQValueThresholdDefault = 1.0;

        /// <summary>
        /// The identification filter and threshold a task TOML's [CommonParameters] table means
        /// (<paramref name="commonParameters"/>, null when the file has none), given the CommonParameters Nett read from it
        /// (<paramref name="read"/>). <c>LegacySentence</c> is the sentence the search says when the old PEP thresholds chose
        /// the filter, and null otherwise. The mapping is this class's.
        /// </summary>
        /// <exception cref="MetaMorpheusException">QValueThreshold or PepQValueThreshold is not a number.</exception>
        public static (FilterType Filter, double Threshold, string LegacySentence) IdentificationFilterOf(TomlTable commonParameters, CommonParameters read)
        {
            bool Holds(string key) => commonParameters?.ContainsKey(key) == true;

            FilterType? named = Holds(nameof(CommonParameters.IdentificationFilter)) ? read.IdentificationFilter : null;
            double? qValueThreshold = Number(commonParameters, LegacyQValueThresholdKey);
            double? pepQValueThreshold = Number(commonParameters, LegacyPepQValueThresholdKey);
            bool oldTaskFile = qValueThreshold.HasValue || pepQValueThreshold.HasValue || LegacyOutputFilterKeys.Any(Holds);

            if (Holds(nameof(CommonParameters.FilterThreshold)) || !oldTaskFile)
            {
                return (named ?? AnalyteDefault(read), read.FilterThreshold, null);
            }

            double q = qValueThreshold ?? LegacyQValueThresholdDefault;
            double pep = pepQValueThreshold ?? LegacyPepQValueThresholdDefault;
            if (named is FilterType filter)
            {
                return (filter, filter == FilterType.PepQValue ? pep : q, null);
            }
            return pep < q
                ? (FilterType.PepQValue, pep, LegacySentence(pep, q))
                : (FilterType.QValueNotch, q, null);
        }

        /// <summary>
        /// The sentence the search warns with, and writes to results.txt, when a task TOML's old PEP thresholds chose its
        /// filter.
        /// </summary>
        public static string LegacySentence(double pepQValueThreshold, double qValueThreshold)
        {
            string pep = pepQValueThreshold.ToString(CultureInfo.InvariantCulture);
            string q = qValueThreshold.ToString(CultureInfo.InvariantCulture);
            return $"IdentificationFilter read from this task's legacy thresholds (PepQValueThreshold {pep} < QValueThreshold {q}): PEP q-value at {pep}";
        }

        /// <summary>
        /// Sets the identification filter and threshold <paramref name="task"/>'s TOML means
        /// (<see cref="IdentificationFilterOf"/>) on its CommonParameters, and records the legacy sentence on the task.
        /// </summary>
        internal static void SetIdentificationFilter(MetaMorpheusTask task, TomlTable commonParameters)
        {
            if (task.CommonParameters == null)
            {
                return;
            }

            var (filter, threshold, legacySentence) = IdentificationFilterOf(commonParameters, task.CommonParameters);
            if (filter != task.CommonParameters.IdentificationFilter || threshold != task.CommonParameters.FilterThreshold)
            {
                task.CommonParameters = task.CommonParameters.CloneWithNewIdentificationFilter(filter, threshold);
            }
            task.SetLegacyFilterSentence(legacySentence);
        }

        private static FilterType AnalyteDefault(CommonParameters read)
        {
            AnalyteType analyte = read.DetermineAnalyteType();
            return CommonParameters.DefaultIdentificationFilter(analyte == AnalyteType.Oligo, analyte == AnalyteType.Proteoform);
        }

        private static double? Number(TomlTable table, string key)
        {
            if (table == null || !table.ContainsKey(key))
            {
                return null;
            }

            TomlObject value = table[key];
            return value switch
            {
                TomlFloat number => number.Value,
                TomlInt number => number.Value,
                _ => throw new MetaMorpheusException($"{key} must be a number, not {value.ReadableTypeName} '{value}'.")
            };
        }
    }
}
