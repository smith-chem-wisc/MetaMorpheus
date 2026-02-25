// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using System.Collections.Generic;
using System.Linq;
using System.Text;
using MassSpectrometry.Dia;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// Results from a DIA extraction + scoring engine run.
    /// </summary>
    public class DiaEngineResults : MetaMorpheusEngineResults
    {
        public List<DiaSearchResult> DiaResults { get; }
        public int TargetLibraryCount { get; }
        public int DecoyLibraryCount { get; }
        public int Ms2ScanCount { get; }
        public int WindowCount { get; }

        public int TargetResultCount => DiaResults.Count(r => !r.IsDecoy);
        public int DecoyResultCount => DiaResults.Count(r => r.IsDecoy);

        public DiaEngineResults(
            MetaMorpheusEngine engine,
            List<DiaSearchResult> diaResults,
            int targetLibraryCount,
            int decoyLibraryCount,
            int ms2ScanCount,
            int windowCount)
            : base(engine)
        {
            DiaResults = diaResults;
            TargetLibraryCount = targetLibraryCount;
            DecoyLibraryCount = decoyLibraryCount;
            Ms2ScanCount = ms2ScanCount;
            WindowCount = windowCount;
        }

        public override string ToString()
        {
            var sb = new StringBuilder(base.ToString());
            sb.AppendLine();
            sb.AppendLine($"DIA search: {Ms2ScanCount} MS2 scans, {WindowCount} windows");
            sb.AppendLine($"Library: {TargetLibraryCount} targets + {DecoyLibraryCount} decoys");
            sb.AppendLine($"Results: {TargetResultCount} targets + {DecoyResultCount} decoys");
            return sb.ToString();
        }
    }
}
