using System.Collections.Generic;
using System.Linq;
using EngineLayer.FdrAnalysis;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// Pooled FDR + PEP over the union of inherited intact PSMs and Pass 3 truncation PSMs
    /// (01_Architecture.md decisions #15, #18). Mirrors the SearchTask order: resolve ambiguities,
    /// drop duplicate PSMs for the same scan, then run the standard <see cref="FdrAnalysisEngine"/>
    /// (which trains PEP internally when the pool is large enough; top-down PSMs get the top-down
    /// training mode). N-terminal and C-terminal truncations are pooled into one FDR category (#15).
    /// </summary>
    public static class TruncationFdr
    {
        public static List<SpectralMatch> RunPooledFdr(
            List<SpectralMatch> pooledPsms,
            CommonParameters commonParameters,
            int massDiffAcceptorNumNotches,
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            string taskId,
            string outputFolder = null,
            bool doPep = true)
        {
            foreach (SpectralMatch psm in pooledPsms.Where(p => p != null))
            {
                psm.ResolveAllAmbiguities();
            }

            List<SpectralMatch> deduped = pooledPsms
                .Where(p => p != null)
                .OrderByDescending(p => p)
                .GroupBy(p => (p.FullFilePath, p.ScanNumber, p.BioPolymerWithSetModsMonoisotopicMass))
                .Select(g => g.First())
                .ToList();

            new FdrAnalysisEngine(deduped, massDiffAcceptorNumNotches, commonParameters, fileSpecificParameters,
                new List<string> { taskId }, analysisType: "PSM", doPEP: doPep, outputFolder: outputFolder).Run();

            return deduped;
        }
    }
}
