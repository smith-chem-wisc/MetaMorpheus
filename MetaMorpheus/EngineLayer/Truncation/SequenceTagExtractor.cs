using System;
using System.Collections.Generic;
using MzLibUtil;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// Extracts short de-novo sequence tags from a deconvoluted MS2 scan: amino-acid strings spelled by runs
    /// of fragment-mass gaps that each match a residue mass. Tags are I/L-normalized (the two residues are
    /// isobaric, both -> 'L'). They are used only as a recall-oriented candidate-protein FILTER for the
    /// database-seeded truncation search (doc §11.2.4 / pTop §3.4.4); Pass 2/3 scoring + FDR police precision.
    ///
    /// Top-down caveat: at large fragment masses the per-peak absolute tolerance is wide relative to
    /// residue-mass differences, so a gap may match several residues (e.g. K vs Q). Such gaps emit one edge
    /// per matching residue, deliberately over-generating tags so the true protein is not missed.
    /// </summary>
    public static class SequenceTagExtractor
    {
        // Standard monoisotopic residue masses (Da). I and L are isobaric -> both reported as 'L'.
        private static readonly (char Aa, double Mass)[] Residues =
        {
            ('G', 57.02146), ('A', 71.03711), ('S', 87.03203), ('P', 97.05276), ('V', 99.06841),
            ('T', 101.04768), ('C', 103.00919), ('L', 113.08406), ('N', 114.04293), ('D', 115.02694),
            ('Q', 128.05858), ('K', 128.09496), ('E', 129.04259), ('M', 131.04049), ('H', 137.05891),
            ('F', 147.06841), ('R', 156.10111), ('Y', 163.06333), ('W', 186.07931)
        };
        private const double MinResidueMass = 57.02146;
        private const double MaxResidueMass = 186.07931;

        /// <summary>I and L are isobaric; normalize both to 'L' so tags and protein k-mers match consistently.</summary>
        public static char Normalize(char aminoAcid) => aminoAcid == 'I' ? 'L' : aminoAcid;

        /// <summary>
        /// Returns the distinct length-<paramref name="tagLength"/> sequence tags found in a scan's
        /// experimental fragment masses (need not be pre-sorted), capped at <paramref name="maxTags"/>.
        /// A gap between two fragment masses that matches a residue mass (within the combined per-peak
        /// tolerance) is an edge; a tag is a run of <paramref name="tagLength"/> consecutive edges.
        /// </summary>
        public static HashSet<string> ExtractTags(IReadOnlyList<double> fragmentMasses, Tolerance tolerance,
            int tagLength = 4, int maxTags = 500)
        {
            var tags = new HashSet<string>();
            if (fragmentMasses == null || fragmentMasses.Count < tagLength + 1)
            {
                return tags;
            }

            double[] masses = new double[fragmentMasses.Count];
            for (int i = 0; i < masses.Length; i++) masses[i] = fragmentMasses[i];
            Array.Sort(masses);
            int n = masses.Length;

            // Edge list per peak: (nextPeakIndex, residueChar) for each gap that matches a residue.
            var edges = new List<(int J, char Aa)>[n];
            for (int i = 0; i < n; i++)
            {
                edges[i] = new List<(int, char)>();
                double mi = masses[i];
                double peakTolI = tolerance.GetMaximumValue(mi) - mi;
                for (int j = i + 1; j < n; j++)
                {
                    double gap = masses[j] - mi;
                    if (gap < MinResidueMass - 0.5) continue;
                    if (gap > MaxResidueMass + 0.5) break; // sorted: no larger residue-sized gap from i
                    // Emit an edge only when the gap matches EXACTLY ONE residue within the window. Gaps that
                    // are ambiguous (two residues in range, e.g. K vs Q at high mass) are dropped rather than
                    // guessed — high precision over recall, so spurious tags do not select the whole database.
                    double window = peakTolI + (tolerance.GetMaximumValue(masses[j]) - masses[j]);
                    char matched = '\0';
                    int matchCount = 0;
                    foreach (var (aa, mass) in Residues)
                    {
                        if (Math.Abs(gap - mass) <= window)
                        {
                            matched = Normalize(aa);
                            if (++matchCount > 1) break;
                        }
                    }
                    if (matchCount == 1)
                    {
                        edges[i].Add((j, matched));
                    }
                }
            }

            // Depth-limited walk emitting tags of exactly tagLength edges.
            var buffer = new char[tagLength];
            void Walk(int node, int depth)
            {
                if (tags.Count >= maxTags) return;
                if (depth == tagLength)
                {
                    tags.Add(new string(buffer));
                    return;
                }
                foreach (var (j, aa) in edges[node])
                {
                    buffer[depth] = aa;
                    Walk(j, depth + 1);
                    if (tags.Count >= maxTags) return;
                }
            }

            for (int i = 0; i < n && tags.Count < maxTags; i++)
            {
                Walk(i, 0);
            }

            return tags;
        }
    }
}
