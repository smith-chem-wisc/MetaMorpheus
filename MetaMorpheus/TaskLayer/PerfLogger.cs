using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    /// <summary>
    /// One row of the rolling, append-only <c>perf_log.tsv</c> performance log (03_Benchmarks.md schema).
    /// A SearchTask run and its paired TruncationSearchTask run each produce one row; truncation-specific
    /// fields are left at their defaults for SearchTask rows.
    /// </summary>
    public class TruncationPerfMetrics
    {
        public string Phase = "Manual";
        public string TaskType = "TruncationSearchTask";
        public string DatasetTag = "";
        public string RunLabel = "";
        public string OutputFolder = "";
        public int NRawFiles;
        public int TotalMs2Scans;
        public double TaskWallSeconds;
        public int NPsmsEmitted;
        public int NProteoformsEmitted;
        public int NPsmsQ01;
        public int NProteoformsQ01;
        public int NPsmsQ05;
        public int NProteoformsQ05;
        // Truncation-specific
        public int NTruncationsTotal;
        public int NTruncationsNterm;
        public int NTruncationsCterm;
        public int NIntactInherited;
        public int NParentsIndexed;
        public int NParentsOversizeExcluded;
        public bool PepQWasUsed;
        public int NDecoysGenerated;
        public double Pass2IndexSeconds;
        public double Pass2ScoringSeconds;
        public double Pass3ChoppingSeconds;
        public double Pass3ScoringSeconds;
        public double FdrPepSeconds;
    }

    /// <summary>
    /// Appends <see cref="TruncationPerfMetrics"/> rows to a tab-separated <c>perf_log.tsv</c>
    /// (03_Benchmarks.md). Writes the header when the file does not yet exist; never overwrites existing
    /// rows (append-only). Best-effort git commit/dirty are captured from the running assembly's directory.
    /// </summary>
    public static class PerfLogger
    {
        private static readonly object FileLock = new object();

        private static readonly string[] Header =
        {
            "timestamp_utc", "git_commit", "git_dirty", "phase", "task_type", "dataset_tag", "run_label",
            "output_folder", "n_raw_files", "total_ms2_scans", "task_wall_seconds", "peak_memory_mb",
            "n_psms_emitted", "n_proteoforms_emitted", "n_psms_q01", "n_proteoforms_q01", "n_psms_q05",
            "n_proteoforms_q05", "n_truncations_total", "n_truncations_nterm", "n_truncations_cterm",
            "n_intact_inherited", "n_parents_indexed", "n_parents_oversize_excluded", "pep_q_was_used",
            "n_decoys_generated", "pass2_index_seconds", "pass2_scoring_seconds", "pass3_chopping_seconds",
            "pass3_scoring_seconds", "fdr_pep_seconds"
        };

        public static void Append(string perfLogPath, TruncationPerfMetrics m)
        {
            (string commit, bool dirty) = TryGetGitState();
            string inv(double d) => d.ToString("0.###", CultureInfo.InvariantCulture);
            double peakMb = Process.GetCurrentProcess().PeakWorkingSet64 / 1024.0 / 1024.0;

            string[] cells =
            {
                DateTime.UtcNow.ToString("yyyy-MM-ddTHH:mm:ssZ", CultureInfo.InvariantCulture),
                commit, dirty ? "true" : "false", m.Phase, m.TaskType, m.DatasetTag, m.RunLabel,
                m.OutputFolder, m.NRawFiles.ToString(CultureInfo.InvariantCulture),
                m.TotalMs2Scans.ToString(CultureInfo.InvariantCulture), inv(m.TaskWallSeconds), inv(peakMb),
                m.NPsmsEmitted.ToString(CultureInfo.InvariantCulture), m.NProteoformsEmitted.ToString(CultureInfo.InvariantCulture),
                m.NPsmsQ01.ToString(CultureInfo.InvariantCulture), m.NProteoformsQ01.ToString(CultureInfo.InvariantCulture),
                m.NPsmsQ05.ToString(CultureInfo.InvariantCulture), m.NProteoformsQ05.ToString(CultureInfo.InvariantCulture),
                m.NTruncationsTotal.ToString(CultureInfo.InvariantCulture), m.NTruncationsNterm.ToString(CultureInfo.InvariantCulture),
                m.NTruncationsCterm.ToString(CultureInfo.InvariantCulture), m.NIntactInherited.ToString(CultureInfo.InvariantCulture),
                m.NParentsIndexed.ToString(CultureInfo.InvariantCulture), m.NParentsOversizeExcluded.ToString(CultureInfo.InvariantCulture),
                m.PepQWasUsed ? "true" : "false", m.NDecoysGenerated.ToString(CultureInfo.InvariantCulture),
                inv(m.Pass2IndexSeconds), inv(m.Pass2ScoringSeconds), inv(m.Pass3ChoppingSeconds),
                inv(m.Pass3ScoringSeconds), inv(m.FdrPepSeconds)
            };

            lock (FileLock)
            {
                bool needHeader = !File.Exists(perfLogPath) || new FileInfo(perfLogPath).Length == 0;
                using var w = new StreamWriter(perfLogPath, append: true);
                if (needHeader)
                {
                    w.WriteLine(string.Join("\t", Header));
                }
                w.WriteLine(string.Join("\t", cells));
            }
        }

        /// <summary>
        /// Parses the 03_Benchmarks run-folder convention &lt;date&gt;_&lt;phase&gt;_&lt;datasetTag&gt;_&lt;runLabel&gt;
        /// into (phase, datasetTag, runLabel). runLabel keeps any trailing underscores. Missing parts stay empty.
        /// </summary>
        public static (string phase, string datasetTag, string runLabel) ParseRunFolderName(string runFolderName)
        {
            if (string.IsNullOrEmpty(runFolderName))
            {
                return ("Manual", "", "");
            }

            string[] parts = runFolderName.Split('_');
            string phase = parts.Length > 1 ? parts[1] : "Manual";
            string datasetTag = parts.Length > 2 ? parts[2] : "";
            string runLabel = parts.Length > 3 ? string.Join("_", parts.Skip(3)) : "";
            return (phase, datasetTag, runLabel);
        }

        private static (string commit, bool dirty) TryGetGitState()
        {
            try
            {
                // The running assembly lives inside the truncation-search worktree (e.g. CMD\bin\...),
                // so git resolves the branch HEAD from there.
                string dir = AppContext.BaseDirectory;
                string commit = RunGit(dir, "rev-parse --short HEAD")?.Trim() ?? "";
                string status = RunGit(dir, "status --porcelain");
                bool dirty = !string.IsNullOrWhiteSpace(status);
                return (commit, dirty);
            }
            catch
            {
                return ("", false);
            }
        }

        private static string RunGit(string workingDir, string args)
        {
            // Use WorkingDirectory rather than `-C "<dir>"`: AppContext.BaseDirectory ends in a separator,
            // and a trailing backslash before a closing quote escapes the quote on Windows.
            var psi = new ProcessStartInfo("git", args)
            {
                WorkingDirectory = workingDir,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                UseShellExecute = false,
                CreateNoWindow = true
            };
            using var p = Process.Start(psi);
            if (p == null)
            {
                return null;
            }
            string output = p.StandardOutput.ReadToEnd();
            p.WaitForExit(5000);
            return p.ExitCode == 0 ? output : null;
        }
    }
}
