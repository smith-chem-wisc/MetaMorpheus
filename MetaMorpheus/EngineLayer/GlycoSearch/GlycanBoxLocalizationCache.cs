using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.GlycoSearch
{
    /// <summary>
    /// The number of glycans needing each motif, for one set of glycan ids. Motifs are listed in order of first appearance.
    /// </summary>
    internal sealed class MotifCount
    {
        public readonly string[] Motifs;
        public readonly int[] Counts;

        public MotifCount(int[] modIds)
        {
            var motifs = new List<string>();
            var counts = new List<int>();
            foreach (var modId in modIds)
            {
                string motif = GlycanBox.MotifOf(modId);
                int i = motifs.IndexOf(motif);
                if (i < 0)
                {
                    motifs.Add(motif);
                    counts.Add(1);
                }
                else
                {
                    counts[i]++;
                }
            }
            Motifs = motifs.ToArray();
            Counts = counts.ToArray();
        }

        /// <summary>
        /// True when sites[0..lastSiteIndex] carry at least as many of each motif as this set needs.
        /// </summary>
        public bool CoveredBy(string[] sites, int lastSiteIndex)
        {
            for (int k = 0; k < Motifs.Length; k++)
            {
                int available = 0;
                for (int i = 0; i <= lastSiteIndex; i++)
                {
                    if (sites[i] == Motifs[k])
                    {
                        available++;
                    }
                }
                if (available < Counts[k])
                {
                    return false;
                }
            }
            return true;
        }
    }

    /// <summary>
    /// What <see cref="LocalizationGraph.LocalizeOGlycan"/> derives from a glycan box's child boxes, which does not depend on
    /// the peptide or the scan. Computing it once per box, rather than once per graph node or edge, is the point.
    /// Each member reproduces exactly what the uncached code computes.
    /// </summary>
    internal sealed class GlycanBoxLocalizationCache
    {
        /// <summary> [y][preY], as <see cref="LocalizationGraph.BuildValidChart"/> returns it. </summary>
        public readonly bool[][] ValidChart;

        /// <summary> Motif counts of each child box, for <see cref="LocalizationGraph.NodeCheck"/>. </summary>
        public readonly MotifCount[] ChildMotifCounts;

        /// <summary>
        /// [y][preY], defined where ValidChart[y][preY]: the motif of the glycan added going from child preY to child y, or null when
        /// nothing is added. <see cref="LocalizationGraph.MotifCheck"/> passes when this is null or equals the site's motif.
        /// </summary>
        public readonly string[][] AddedMotif;

        public GlycanBoxLocalizationCache(GlycanBox box)
        {
            var children = box.ChildGlycanBoxes;
            var chart = LocalizationGraph.BuildValidChart(children);
            ValidChart = new bool[children.Length][];
            AddedMotif = new string[children.Length][];
            ChildMotifCounts = new MotifCount[children.Length];
            for (int y = 0; y < children.Length; y++)
            {
                ValidChart[y] = chart[y];
                ChildMotifCounts[y] = new MotifCount(children[y].ModIds);
                AddedMotif[y] = new string[children.Length];
                for (int preY = 0; preY <= y; preY++)
                {
                    if (!ValidChart[y][preY])
                    {
                        continue;
                    }
                    var diff = LocalizationGraph.GetDiff(children[preY].ModIds, children[y].ModIds);
                    AddedMotif[y][preY] = diff == null || diff.Length == 0 ? null : GlycanBox.MotifOf(diff[0]);
                }
            }
        }
    }
}
