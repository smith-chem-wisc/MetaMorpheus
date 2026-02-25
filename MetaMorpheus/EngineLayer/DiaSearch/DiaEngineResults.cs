// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using System.Collections.Generic;
using System.Linq;
using System.Text;
using MassSpectrometry.Dia;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// Results container for DiaEngine, extending MetaMorpheusEngineResults
    /// for integration with the standard MetaMorpheus reporting pipeline.
    /// </summary>
    public sealed class DiaEngineResults : MetaMorpheusEngineResults
    {
        /// <summary>DIA precursor match results (targets + decoys that passed filters).</summary>
        public List<DiaSearchResult> DiaResults { get; }

        /// <summary>Number of target (non-decoy) library entries provided as input.</summary>
        public int TargetLibraryCount { get; }

        /// <summary>Number of decoy library entries provided as input.</summary>
        public int DecoyLibraryCount { get; }

        public DiaEngineResults(
            DiaEngine engine,
            List<DiaSearchResult> diaResults,
            int targetLibraryCount,
            int decoyLibraryCount)
            : base(engine)
        {
            DiaResults = diaResults ?? new List<DiaSearchResult>();
            TargetLibraryCount = targetLibraryCount;
            DecoyLibraryCount = decoyLibraryCount;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine($"Library: {TargetLibraryCount} targets, {DecoyLibraryCount} decoys");
            sb.AppendLine($"DIA matches: {DiaResults.Count}");

            if (DiaResults.Count > 0)
            {
                int targetHits = 0;
                int decoyHits = 0;
                for (int i = 0; i < DiaResults.Count; i++)
                {
                    if (DiaResults[i].IsDecoy) decoyHits++;
                    else targetHits++;
                }
                sb.AppendLine($"  Target hits: {targetHits}, Decoy hits: {decoyHits}");

                var scored = DiaResults.Where(r => !float.IsNaN(r.DotProductScore)).ToList();
                if (scored.Count > 0)
                {
                    sb.AppendLine($"  Avg dot product: {scored.Average(r => r.DotProductScore):F4}");
                    sb.AppendLine($"  Best dot product: {scored.Max(r => r.DotProductScore):F4}");
                }
            }

            return sb.ToString();
        }
    }
}