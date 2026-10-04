using System;
using System.Collections.Generic;

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
            ValidChart = new bool[children.Length][];
            AddedMotif = new string[children.Length][];
            ChildMotifCounts = new MotifCount[children.Length];
            for (int y = 0; y < children.Length; y++)
            {
                ValidChart[y] = new bool[children.Length];
                ChildMotifCounts[y] = new MotifCount(children[y].ModIds);
                AddedMotif[y] = new string[children.Length];
                for (int preY = 0; preY <= y; preY++)
                {
                    // As BuildValidChart: the later child holds at most one more mod and contains every mod of the earlier child.
                    if (children[y].NumberOfMods <= children[preY].NumberOfMods + 1
                        && (children[preY].NumberOfMods == 0 || IsGlycanCoveredBy(children[preY].ModIds, children[y].ModIds)))
                    {
                        ValidChart[y][preY] = true;
                        int addedGlycan = AddedGlycan(children[preY].ModIds, children[y].ModIds, out bool anyAdded);
                        AddedMotif[y][preY] = anyAdded ? GlycanBox.MotifOf(addedGlycan) : null;
                    }
                }
            }
        }

        /// <summary>
        /// Whether <paramref name="cover"/> holds every id of <paramref name="glycan"/>, counting repeats; what
        /// LocalizationGraph.TryGetLeft(cover, glycan) returns, without allocating.
        /// </summary>
        internal static bool IsGlycanCoveredBy(int[] glycan, int[] cover)
        {
            for (int i = 0; i < glycan.Length; i++)
            {
                // Check each distinct id once, at its first occurrence in glycan.
                if (Array.IndexOf(glycan, glycan[i], 0, i) >= 0)
                {
                    continue;
                }
                if (Count(cover, glycan[i]) < Count(glycan, glycan[i]))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// The first element of LocalizationGraph.GetDiff(pre, current), without allocating. GetDiff lists the ids of
        /// <paramref name="current"/> left after removing those of <paramref name="pre"/>, grouped by id in the order each id first
        /// appears in current, so its first element is the earliest-appearing id that current holds more times than pre.
        /// <paramref name="pre"/> must be contained in current, as GetDiff also requires.
        /// </summary>
        internal static int AddedGlycan(int[] pre, int[] current, out bool anyAdded)
        {
            for (int i = 0; i < current.Length; i++)
            {
                if (Array.IndexOf(current, current[i], 0, i) >= 0)
                {
                    continue;
                }
                if (Count(current, current[i]) > Count(pre, current[i]))
                {
                    anyAdded = true;
                    return current[i];
                }
            }
            anyAdded = false;
            return 0;
        }

        private static int Count(int[] ids, int id)
        {
            int n = 0;
            foreach (int x in ids)
            {
                if (x == id)
                {
                    n++;
                }
            }
            return n;
        }
    }
}
