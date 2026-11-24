#nullable enable
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// Precursor q-values across runs searched together (PLAN D4, "precursor, global"): each precursor's best score over the
/// runs, then target-decoy competition with the rule each run's own q-values use. The final list needs both this and the
/// run's q-value at the threshold.
/// </summary>
public static class DiaPrecursorFdr
{
    /// <summary>The global q-value of every target precursor seen in any run, by library precursor index.</summary>
    public static Dictionary<int, double> GlobalQValues(IEnumerable<IReadOnlyList<DiaPrecursorMatch>> runs)
    {
        var best = new BestScores();
        foreach (var run in runs)
            best.Add(run);
        return best.GlobalQValues();
    }

    /// <summary>
    /// Each precursor's best score so far, added a run at a time, so a many-run search need not hold every run's matches.
    /// </summary>
    public sealed class BestScores
    {
        private readonly Dictionary<int, (bool IsDecoy, double Score)> _best = [];

        public void Add(IEnumerable<DiaPrecursorMatch> run)
        {
            foreach (var match in run)
                if (!_best.TryGetValue(match.PrecursorIndex, out var seen) || match.Score > seen.Score)
                    _best[match.PrecursorIndex] = (match.IsDecoy, match.Score);
        }

        /// <summary>The global q-value of every target precursor added, by library precursor index.</summary>
        public Dictionary<int, double> GlobalQValues()
        {
            var targets = _best.Where(p => !p.Value.IsDecoy).ToList();
            var decoyScores = _best.Values.Where(v => v.IsDecoy).Select(v => v.Score).ToList();
            double[] q = DeconvolutionQValueCalculator.AssignQValues(targets.Select(p => p.Value.Score).ToList(), decoyScores);
            // As in each run (DiaLibrarySearchEngine): with no decoy scored, (D+1)/T gives every target 1/T, not mzLib's 0
            if (decoyScores.Count == 0 && targets.Count > 0)
                System.Array.Fill(q, 1.0 / targets.Count);
            return targets.Select((p, i) => (p.Key, Q: q[i])).ToDictionary(t => t.Key, t => t.Q);
        }
    }
}
