using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;

namespace EngineLayer
{
    /// <summary>
    /// What MetaMorpheus records beside an experimental design it wrote from an SDRF (QuantProject QP-S41, ruled
    /// 2026-11-07): the design file and the SDRF, each by its sha256, and how the condition was built, including any
    /// grouping the user stated (ruled 2026-11-11). It is how a later run tells a design MetaMorpheus wrote from one
    /// the user wrote or edited: while the design's sha256 is the one recorded here, the design is MetaMorpheus's, and
    /// a run writes it again for its own spectra files; once edited, it is the user's and is never touched. The design's
    /// sha256 is of its text with line endings ignored (<see cref="ContentSha256Of"/>, ruled 2026-11-13): the GUI's
    /// design windows save Windows line endings, so saving one unchanged is not an edit.
    ///
    /// The marker also keeps the grouping the user stated, as keys (<see cref="ConditionColumns"/>,
    /// <see cref="SingleCondition"/>, <see cref="UnknownFactorIsOwnCondition"/>), so that writing the design again for a
    /// later run uses it without anyone stating it again: the grouping belongs with the design it made, not with a task.
    ///
    /// A plain text file of tab-separated keys and values beside the design, so a person reading the folder sees why
    /// the design is there.
    /// </summary>
    public sealed class SdrfDesignMarker
    {
        /// <summary>The marker's file name, beside the design.</summary>
        public const string FileName = "DesignFromSdrf.txt";

        private const string Heading = "# MetaMorpheus wrote the experimental design beside this file from an SDRF. While the design is as written, " +
                                       "MetaMorpheus writes it again for each run; once edited, the design is yours and is never touched.";

        private const string GroupingHeading = "# ConditionColumns, SingleCondition and UnknownFactorIsOwnCondition, when present, are the grouping " +
                                               "the user stated; a run that states none uses them again.";

        public SdrfDesignMarker(string designFileName, string designSha256, string sdrfPath, string sdrfSha256, string condition,
            IReadOnlyList<string> conditionColumns = null, string singleCondition = null, bool unknownFactorIsOwnCondition = false)
        {
            DesignFileName = designFileName;
            DesignSha256 = designSha256;
            SdrfPath = sdrfPath;
            SdrfSha256 = sdrfSha256;
            Condition = condition;
            ConditionColumns = conditionColumns?.ToList() ?? new List<string>();
            SingleCondition = singleCondition;
            UnknownFactorIsOwnCondition = unknownFactorIsOwnCondition;
        }

        /// <summary>The design's file name: ExperimentalDesign.tsv or TmtDesign.txt.</summary>
        public string DesignFileName { get; }

        /// <summary>The design's <see cref="ContentSha256Of"/> as MetaMorpheus wrote it, lowercase hex.</summary>
        public string DesignSha256 { get; }

        /// <summary>The SDRF the design was made from, as a full path.</summary>
        public string SdrfPath { get; }

        /// <summary>The SDRF's sha256 when the design was made, lowercase hex.</summary>
        public string SdrfSha256 { get; }

        /// <summary>How the condition was built, in words, with any grouping the user stated.</summary>
        public string Condition { get; }

        /// <summary>The SDRF columns the user named for the condition, in order; empty when none were named.</summary>
        public IReadOnlyList<string> ConditionColumns { get; }

        /// <summary>The one condition the user stated every sample is; null when none was stated.</summary>
        public string SingleCondition { get; }

        /// <summary>True when the user stated that files whose factor nobody stated are a condition of their own.</summary>
        public bool UnknownFactorIsOwnCondition { get; }

        /// <summary>True when the marker keeps a grouping the user stated.</summary>
        public bool StatesGrouping => ConditionColumns.Count > 0 || SingleCondition != null || UnknownFactorIsOwnCondition;

        /// <summary>The marker beside a design in <paramref name="designDirectory"/>.</summary>
        public static string PathIn(string designDirectory) => Path.Combine(designDirectory, FileName);

        /// <summary>The line results.txt carries for a design written from an SDRF.</summary>
        public string Describe() =>
            $"Experimental design ({DesignFileName}) written by MetaMorpheus from the SDRF {SdrfPath} (sha256 {SdrfSha256}); condition: {Condition}.";

        /// <summary>
        /// True when the design this marker names is in <paramref name="designDirectory"/> and its text is what MetaMorpheus
        /// wrote, line endings aside. False when it is missing, unreadable or edited.
        /// </summary>
        public bool DescribesUneditedDesignIn(string designDirectory)
        {
            string design = Path.Combine(designDirectory, DesignFileName);
            try
            {
                return File.Exists(design) && string.Equals(ContentSha256Of(design), DesignSha256, StringComparison.OrdinalIgnoreCase);
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                return false;
            }
        }

        /// <summary>Writes the marker into <paramref name="designDirectory"/>, replacing any there.</summary>
        public void WriteTo(string designDirectory)
        {
            var lines = new List<string>
            {
                Heading,
                "design\t" + DesignFileName,
                "design sha256\t" + DesignSha256,
                "sdrf\t" + SdrfPath,
                "sdrf sha256\t" + SdrfSha256,
                "condition\t" + Condition,
                GroupingHeading,
            };
            if (ConditionColumns.Count > 0)
                lines.Add("ConditionColumns\t" + string.Join("\t", ConditionColumns));
            if (SingleCondition != null)
                lines.Add("SingleCondition\t" + SingleCondition);
            if (UnknownFactorIsOwnCondition)
                lines.Add("UnknownFactorIsOwnCondition\ttrue");
            File.WriteAllText(PathIn(designDirectory), string.Join("\n", lines) + "\n", new UTF8Encoding(false));
        }

        /// <summary>
        /// The marker in <paramref name="designDirectory"/>, or null when there is none or it cannot be read as one: a
        /// marker that is missing a value is no marker, so the design beside it is treated as the user's. A grouping key that
        /// is absent states nothing; ConditionColumns holds the columns tab-separated, and UnknownFactorIsOwnCondition is
        /// stated only by 'true'.
        /// </summary>
        public static SdrfDesignMarker ReadFrom(string designDirectory)
        {
            string path = PathIn(designDirectory);
            Dictionary<string, string> values;
            try
            {
                if (!File.Exists(path))
                    return null;
                values = File.ReadAllLines(path)
                    .Where(l => !l.StartsWith("#") && l.Contains('\t'))
                    .Select(l => l.Split('\t', 2))
                    .GroupBy(kv => kv[0], StringComparer.Ordinal)
                    .ToDictionary(g => g.Key, g => g.First()[1], StringComparer.Ordinal);
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                return null;
            }

            string Value(string key) => values.TryGetValue(key, out var v) && !string.IsNullOrWhiteSpace(v) ? v : null;
            string design = Value("design");
            string designSha = Value("design sha256");
            string sdrf = Value("sdrf");
            string sdrfSha = Value("sdrf sha256");
            string condition = Value("condition");
            if (design == null || designSha == null || sdrf == null || sdrfSha == null || condition == null)
                return null;
            var columns = Value("ConditionColumns")?.Split('\t').Where(c => !string.IsNullOrWhiteSpace(c)).ToList();
            bool unknownIsOwn = string.Equals(Value("UnknownFactorIsOwnCondition"), "true", StringComparison.OrdinalIgnoreCase);
            return new SdrfDesignMarker(design, designSha, sdrf, sdrfSha, condition, columns, Value("SingleCondition"), unknownIsOwn);
        }

        /// <summary>A file's sha256, lowercase hex, read from the bytes on disk.</summary>
        public static string Sha256Of(string path)
        {
            using var stream = File.OpenRead(path);
            return Convert.ToHexString(SHA256.HashData(stream)).ToLowerInvariant();
        }

        /// <summary>
        /// The sha256 of a text file's content, lowercase hex, with its line endings ignored: a byte-order mark dropped, CR LF
        /// and CR read as LF, and trailing line ends dropped. Two files differing only in how lines end give one value.
        /// </summary>
        public static string ContentSha256Of(string path)
        {
            string text = File.ReadAllText(path, Encoding.UTF8).TrimStart('\uFEFF')
                .Replace("\r\n", "\n").Replace('\r', '\n').TrimEnd('\n');
            return Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(text))).ToLowerInvariant();
        }
    }
}
