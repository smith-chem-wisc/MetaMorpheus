// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License
using MassSpectrometry.Dia;
using System.Collections.Generic;
using System.Text;
namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// Results returned by DiaEngine after a DIA library search.
    /// Contains the list of DIA search results and the iRT calibration model (if fitted).
    /// </summary>
    public class DiaEngineResults : MetaMorpheusEngineResults
    {
        /// <summary>All DIA search results passing the minimum fragment and score thresholds.</summary>
        public List<DiaSearchResult> DiaResults { get; }
        /// <summary>
        /// The final iRT calibration model from the iterative fitting process.
        /// Null if calibration was disabled, failed, or insufficient anchors were found.
        /// </summary>
        public RtCalibrationModel Calibration { get; }
        /// <summary>Number of target (non-decoy) results.</summary>
        public int TargetCount { get; }
        /// <summary>Number of decoy results.</summary>
        public int DecoyCount { get; }
        /// <summary>Number of target library precursors that were searched.</summary>
        public int TargetLibraryCount { get; }
        /// <summary>Number of decoy library precursors that were searched.</summary>
        public int DecoyLibraryCount { get; }

        public DiaEngineResults(DiaEngine engine, List<DiaSearchResult> results,
            RtCalibrationModel calibration, int targetLibraryCount = 0, int decoyLibraryCount = 0)
            : base(engine)
        {
            DiaResults = results ?? new List<DiaSearchResult>();
            Calibration = calibration;
            TargetLibraryCount = targetLibraryCount;
            DecoyLibraryCount = decoyLibraryCount;
            int targets = 0, decoys = 0;
            for (int i = 0; i < DiaResults.Count; i++)
            {
                if (DiaResults[i].IsDecoy) decoys++;
                else targets++;
            }
            TargetCount = targets;
            DecoyCount = decoys;
        }
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine($"DIA Results: {DiaResults.Count} ({TargetCount} targets, {DecoyCount} decoys)");
            sb.AppendLine($"Library: {TargetLibraryCount} targets, {DecoyLibraryCount} decoys searched");
            if (Calibration != null)
            {
                sb.AppendLine($"iRT Calibration: {Calibration}");
            }
            else
            {
                sb.AppendLine("iRT Calibration: not fitted");
            }
            return sb.ToString();
        }
    }
}