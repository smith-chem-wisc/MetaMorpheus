using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// Writes the two truncation-search result files (01_Architecture.md decisions #16, #17) using the
    /// standard MetaMorpheus PSM columns (<see cref="SpectralMatch.GetTabSeparatedHeader"/> /
    /// <see cref="SpectralMatch.ToString(System.Collections.Generic.IReadOnlyDictionary{string,int},bool,bool)"/>),
    /// so the truncation type rides in the standard "Description" column (#13). No q-value/PEP cutoff is
    /// applied at write time — the full pooled set is written, decoys and contaminants included.
    /// </summary>
    public static class TruncationOutput
    {
        /// <summary>AllTruncatedPSMs.psmtsv — one row per pooled PSM (scan × candidate).</summary>
        public static void WritePsms(IEnumerable<SpectralMatch> psms, string filePath, IReadOnlyDictionary<string, int> modsToWrite = null)
        {
            IReadOnlyDictionary<string, int> mods = modsToWrite ?? new Dictionary<string, int>();
            using var output = new StreamWriter(filePath);
            output.WriteLine(SpectralMatch.GetTabSeparatedHeader());
            foreach (SpectralMatch psm in psms.Where(p => p != null).OrderByDescending(p => p))
            {
                output.WriteLine(psm.ToString(mods, writePeptideLevelFdr: false));
            }
        }

        /// <summary>AllTruncatedProteoforms.psmtsv — deduped by truncated FullSequence (best per group), peptide-level FDR columns.</summary>
        public static void WriteProteoforms(IEnumerable<SpectralMatch> psms, string filePath, IReadOnlyDictionary<string, int> modsToWrite = null)
        {
            IReadOnlyDictionary<string, int> mods = modsToWrite ?? new Dictionary<string, int>();
            IEnumerable<SpectralMatch> deduped = psms
                .Where(p => p != null)
                .OrderByDescending(p => p)
                .GroupBy(p => p.FullSequence)
                .Select(g => g.First());

            using var output = new StreamWriter(filePath);
            output.WriteLine(SpectralMatch.GetTabSeparatedHeader());
            foreach (SpectralMatch psm in deduped)
            {
                output.WriteLine(psm.ToString(mods, writePeptideLevelFdr: true));
            }
        }
    }
}
