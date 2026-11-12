using EngineLayer.Gptmd;
using Omics;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Transcriptomics.Digestion;

namespace EngineLayer.Localization
{
    /// <summary>
    /// Which residues a modification may be moved to when its arrangements are enumerated.
    /// </summary>
    public enum LocalizationSiteRule
    {
        /// <summary>
        /// Residues where the modification, or one of the same chemistry, fits by motif and location restriction.
        /// Terminal modifications stay on their terminus.
        /// </summary>
        Motif,

        /// <summary>
        /// Any residue: the modification's chemistry is treated as unknown, so its motif and location restriction are
        /// ignored, and terminal modifications move too. For experiments where a modification may sit anywhere.
        /// </summary>
        AnyResidue
    }

    /// <summary>
    /// A peptide's modifications in one arrangement: the placement of the modifications whose position is in question,
    /// and the full modification dictionary that builds the peptide with them there.
    /// </summary>
    public sealed record PositionalIsoform(LocalizationPlacement Placement, Dictionary<int, Modification> AllModsOneIsNterminus);

    /// <summary>
    /// Enumerates the arrangements of a peptide's modifications over the residues each could occupy.
    /// </summary>
    public static class PositionalIsoforms
    {
        /// <summary>
        /// The most arrangements enumerated for one match. Three phosphates on a peptide with fifteen S, T and Y have 455.
        /// </summary>
        public const int DefaultMaxPlacements = 1024;

        /// <summary>
        /// Enumerates every arrangement of the peptide's movable modifications. Fixed modifications stay where they are
        /// and hold their positions. Under <see cref="LocalizationSiteRule.Motif"/> so do terminal modifications.
        ///
        /// Under <see cref="LocalizationSiteRule.Motif"/>, a movable modification may go to any residue where it, or a
        /// modification of the same chemistry from <paramref name="placeableModifications"/>, fits by motif and location
        /// restriction; at a residue its own motif does not allow, the same-chemistry modification is placed instead
        /// (phosphorylation on S moves to a T as phosphorylation on T). Under <see cref="LocalizationSiteRule.AnyResidue"/>
        /// it may go to any residue as itself.
        ///
        /// Modifications of the same chemistry are interchangeable, so two phosphates on S3 and T5 are one arrangement,
        /// not two. A position holds at most one modification. A modification reported on a peptide terminus stands for
        /// the terminal residue as well (the two fragment identically), so the terminal residue is not offered to it again.
        /// </summary>
        /// <param name="peptide">The peptide in the arrangement the search reported.</param>
        /// <param name="placeableModifications">Modifications a movable one may be swapped for, typically the search's
        /// variable modifications. Unused under <see cref="LocalizationSiteRule.AnyResidue"/>.</param>
        /// <param name="isFixed">Whether a modification is a fixed modification of the search.</param>
        /// <param name="maxPlacements">The most arrangements to enumerate.</param>
        /// <param name="isoforms">Every arrangement, the reported one among them; null when there are more than
        /// <paramref name="maxPlacements"/>.</param>
        /// <param name="siteRule">Which residues a modification may be moved to.</param>
        /// <returns>False when there are more than <paramref name="maxPlacements"/> arrangements.</returns>
        public static bool TryEnumerate(IBioPolymerWithSetMods peptide, IReadOnlyCollection<Modification> placeableModifications,
            Func<Modification, bool> isFixed, int maxPlacements, out List<PositionalIsoform> isoforms,
            LocalizationSiteRule siteRule = LocalizationSiteRule.Motif)
        {
            int length = peptide.Length;
            int nTerminusKey = 1;
            int cTerminusKey = length + 2;
            var stationary = new Dictionary<int, Modification>();
            var groups = new List<ChemistryGroup>();

            foreach (var (key, mod) in peptide.AllModsOneIsNterminus)
            {
                bool onTerminus = key == nTerminusKey || key == cTerminusKey;
                if (isFixed(mod) || (onTerminus && siteRule == LocalizationSiteRule.Motif))
                {
                    stationary.Add(key, mod);
                    continue;
                }

                ChemistryGroup group = groups.FirstOrDefault(g => GptmdEngine.IsSameChemistry(g.Members[0], mod));
                if (group == null)
                {
                    group = new ChemistryGroup();
                    groups.Add(group);
                }
                group.Members.Add(mod);
                group.ReportedKeys.Add(key);
            }

            foreach (ChemistryGroup group in groups)
            {
                // a reported position holds the reported modification, whether or not its motif is the one that fits
                for (int m = 0; m < group.Members.Count; m++)
                {
                    group.Allowed.Add((group.ReportedKeys[m], group.Members[m]));
                }

                for (int r = 0; r < length; r++)
                {
                    int key = r + 2;
                    if (stationary.ContainsKey(key) || group.ReportedKeys.Contains(key))
                    {
                        continue;
                    }

                    // a terminal residue is already offered through its terminus when the group was reported there
                    if ((r == 0 && group.ReportedKeys.Contains(nTerminusKey)) || (r == length - 1 && group.ReportedKeys.Contains(cTerminusKey)))
                    {
                        continue;
                    }

                    Modification placed = siteRule == LocalizationSiteRule.AnyResidue
                        ? group.Members[0]
                        : group.Members.Concat(placeableModifications ?? Array.Empty<Modification>())
                            .FirstOrDefault(m => GptmdEngine.IsSameChemistry(m, group.Members[0]) && Fits(m, peptide, r));

                    if (placed != null)
                    {
                        group.Allowed.Add((key, placed));
                    }
                }

                group.Allowed.Sort((a, b) => ResidueIndexOf(a.Key, length).CompareTo(ResidueIndexOf(b.Key, length)));
            }

            isoforms = new List<PositionalIsoform>();
            var chosen = new List<(int Key, Modification Modification)>();
            var used = new HashSet<int>();
            if (!Enumerate(peptide, groups, 0, 0, chosen, used, stationary, maxPlacements, isoforms))
            {
                isoforms = null;
                return false;
            }
            return true;
        }

        /// <summary>
        /// Builds the peptide or oligo with the given modifications, keeping everything else about it. Null for any
        /// other kind of biopolymer.
        /// </summary>
        public static IBioPolymerWithSetMods WithModifications(IBioPolymerWithSetMods bioPolymer, Dictionary<int, Modification> allModsOneIsNterminus)
        {
            return bioPolymer switch
            {
                PeptideWithSetModifications peptide => new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams,
                    peptide.OneBasedStartResidueInProtein, peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory,
                    peptide.PeptideDescription, peptide.MissedCleavages, allModsOneIsNterminus, peptide.NumFixedMods),
                OligoWithSetMods oligo => new OligoWithSetMods(oligo.NucleicAcid, oligo.DigestionParams as RnaDigestionParams,
                    oligo.OneBasedStartResidue, oligo.OneBasedEndResidue, oligo.MissedCleavages, oligo.CleavageSpecificityForFdrCategory,
                    allModsOneIsNterminus, oligo.NumFixedMods, oligo.FivePrimeTerminus, oligo.ThreePrimeTerminus),
                _ => null
            };
        }

        /// <summary>
        /// The zero-based residue a modification-dictionary key stands for: a terminus stands for its terminal residue.
        /// </summary>
        private static int ResidueIndexOf(int key, int length) => Math.Clamp(key - 2, 0, length - 1);

        /// <summary>
        /// Places each group's modifications on a combination of its allowed positions not already taken, group by group,
        /// in increasing order within a group so that each arrangement is produced once.
        /// </summary>
        private static bool Enumerate(IBioPolymerWithSetMods peptide, List<ChemistryGroup> groups, int groupIndex, int start,
            List<(int Key, Modification Modification)> chosen, HashSet<int> used, Dictionary<int, Modification> stationary,
            int maxPlacements, List<PositionalIsoform> isoforms)
        {
            if (groupIndex == groups.Count)
            {
                if (isoforms.Count == maxPlacements)
                {
                    return false;
                }

                var mods = new Dictionary<int, Modification>(stationary);
                var sites = new List<LocalizationSite>(chosen.Count);
                foreach (var (key, mod) in chosen)
                {
                    mods.Add(key, mod);
                    int r = ResidueIndexOf(key, peptide.Length);
                    sites.Add(new LocalizationSite(r, peptide.BaseSequence[r], mod));
                }
                isoforms.Add(new PositionalIsoform(new LocalizationPlacement(sites.OrderBy(s => s.ResidueIndex).ToArray()), mods));
                return true;
            }

            ChemistryGroup group = groups[groupIndex];
            int placedInGroup = chosen.Count - groups.Take(groupIndex).Sum(g => g.Members.Count);
            if (placedInGroup == group.Members.Count)
            {
                return Enumerate(peptide, groups, groupIndex + 1, 0, chosen, used, stationary, maxPlacements, isoforms);
            }

            for (int a = start; a < group.Allowed.Count; a++)
            {
                var (key, mod) = group.Allowed[a];
                if (!used.Add(key))
                {
                    continue;
                }
                chosen.Add((key, mod));

                bool withinLimit = Enumerate(peptide, groups, groupIndex, a + 1, chosen, used, stationary, maxPlacements, isoforms);

                chosen.RemoveAt(chosen.Count - 1);
                used.Remove(key);
                if (!withinLimit)
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Whether the modification fits the residue at zero-based <paramref name="residueIndex"/> of the peptide, by
        /// motif and location restriction, judged against the parent sequence when there is one.
        /// </summary>
        private static bool Fits(Modification mod, IBioPolymerWithSetMods peptide, int residueIndex)
        {
            if (peptide.Parent?.BaseSequence is { Length: > 0 } parentSequence)
            {
                return ModificationLocalization.ModFits(mod, parentSequence, residueIndex + 1, peptide.Length, peptide.OneBasedStartResidue + residueIndex);
            }
            return ModificationLocalization.ModFits(mod, peptide.BaseSequence, residueIndex + 1, peptide.Length, residueIndex + 1);
        }

        /// <summary>
        /// The movable modifications of one chemistry, the positions (modification-dictionary keys) they were reported
        /// on, and the positions they may occupy.
        /// </summary>
        private sealed class ChemistryGroup
        {
            public List<Modification> Members { get; } = new();
            public List<int> ReportedKeys { get; } = new();
            public List<(int Key, Modification Modification)> Allowed { get; } = new();
        }
    }
}
