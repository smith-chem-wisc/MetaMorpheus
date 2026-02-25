// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry.Dia;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// Contains the results of a DIA engine search for a single raw file.
    /// Wraps a list of <see cref="DiaSearchResult"/> objects from the mzLib DIA engine.
    /// </summary>
    public class DiaEngineResults : MetaMorpheusEngineResults
    {
        /// <summary>
        /// All DIA identifications (targets + decoys) that passed the minimum
        /// fragment and score thresholds during result assembly.
        /// </summary>
        public List<DiaSearchResult> DiaResults { get; }

        /// <summary>Number of target (non-decoy) identifications.</summary>
        public int TargetCount => DiaResults.Count(r => !r.IsDecoy);

        /// <summary>Number of decoy identifications.</summary>
        public int DecoyCount => DiaResults.Count(r => r.IsDecoy);

        /// <summary>Number of target library entries that produced results. Alias for TargetCount.</summary>
        public int TargetLibraryCount => TargetCount;

        /// <summary>Number of decoy library entries that produced results. Alias for DecoyCount.</summary>
        public int DecoyLibraryCount => DecoyCount;

        public DiaEngineResults(MetaMorpheusEngine engine, List<DiaSearchResult> diaResults)
            : base(engine)
        {
            DiaResults = diaResults ?? new List<DiaSearchResult>();
        }

        public override string ToString()
        {
            var sb = new StringBuilder(base.ToString());
            sb.AppendLine();
            sb.Append($"DIA results: {TargetCount:N0} targets, {DecoyCount:N0} decoys " +
                       $"({DiaResults.Count:N0} total)");
            return sb.ToString();
        }
    }
}
