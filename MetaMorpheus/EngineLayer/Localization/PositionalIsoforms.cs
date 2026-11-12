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
        /// Enumerates every arrangement of the peptide's movable modifications: those on residues that are not fixed
        /// modifications. Modifications on a peptide terminus and fixed modifications stay where they are, and their
        /// residues are not available to the others.
        ///
        /// A movable modification may go to any residue where it, or a modification of the same chemistry from
        /// <paramref name="placeableModifications"/>, fits by motif and location restriction; at a residue its own
        /// motif does not allow, the same-chemistry modification is placed instead (phosphorylation on S moves to a T as
        /// phosphorylation on T). Modifications of the same chemistry are interchangeable, so two phosphates on S3 and
        /// T5 are one arrangement, not two. A residue holds at most one modification.
        /// </summary>
        /// <param name="peptide">The peptide in the arrangement the search reported.</param>
        /// <param name="placeableModifications">Modifications a movable one may be swapped for, typically the search's
        /// variable modifications.</param>
        /// <param name="isFixed">Whether a modification is a fixed modification of the search.</param>
        /// <param name="maxPlacements">The most arrangements to enumerate.</param>
        /// <param name="isoforms">Every arrangement, the reported one among them; null when there are more than
        /// <paramref name="maxPlacements"/>.</param>
        /// <returns>False when there are more than <paramref name="maxPlacements"/> arrangements.</returns>
        public static bool TryEnumerate(IBioPolymerWithSetMods peptide, IReadOnlyCollection<Modification> placeableModifications,
            Func<Modification, bool> isFixed, int maxPlacements, out List<PositionalIsoform> isoforms)
        {
            int length = peptide.Length;
            var stationary = new Dictionary<int, Modification>();
            var groups = new List<ChemistryGroup>();

            foreach (var (key, mod) in peptide.AllModsOneIsNterminus)
            {
                bool onResidue = key >= 2 && key <= length + 1;
                if (!onResidue || isFixed(mod))
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
                group.ReportedResidues.Add(key - 2);
            }

            foreach (ChemistryGroup group in groups)
            {
                for (int r = 0; r < length; r++)
                {
                    if (stationary.ContainsKey(r + 2))
                    {
                        continue;
                    }

                    // the reported site holds the reported modification, whether or not its motif is the one that fits
                    int reported = group.ReportedResidues.IndexOf(r);
                    Modification fitting = reported >= 0
                        ? group.Members[reported]
                        : group.Members.Concat(placeableModifications ?? Array.Empty<Modification>())
                            .FirstOrDefault(m => GptmdEngine.IsSameChemistry(m, group.Members[0]) && Fits(m, peptide, r));

                    if (fitting != null)
                    {
                        group.Allowed.Add((r, fitting));
                    }
                }
            }

            isoforms = new List<PositionalIsoform>();
            var chosen = new List<LocalizationSite>();
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
        /// Places each group's modifications on a combination of its allowed residues not already taken, group by group,
        /// in increasing residue order within a group so that each arrangement is produced once.
        /// </summary>
        private static bool Enumerate(IBioPolymerWithSetMods peptide, List<ChemistryGroup> groups, int groupIndex, int start,
            List<LocalizationSite> chosen, HashSet<int> used, Dictionary<int, Modification> stationary, int maxPlacements,
            List<PositionalIsoform> isoforms)
        {
            if (groupIndex == groups.Count)
            {
                if (isoforms.Count == maxPlacements)
                {
                    return false;
                }

                var mods = new Dictionary<int, Modification>(stationary);
                foreach (LocalizationSite site in chosen)
                {
                    mods.Add(site.ResidueIndex + 2, site.Modification);
                }
                isoforms.Add(new PositionalIsoform(new LocalizationPlacement(chosen.OrderBy(s => s.ResidueIndex).ToArray()), mods));
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
                var (r, mod) = group.Allowed[a];
                if (!used.Add(r))
                {
                    continue;
                }
                chosen.Add(new LocalizationSite(r, peptide.BaseSequence[r], mod));

                bool withinLimit = Enumerate(peptide, groups, groupIndex, a + 1, chosen, used, stationary, maxPlacements, isoforms);

                chosen.RemoveAt(chosen.Count - 1);
                used.Remove(r);
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
        /// The movable modifications of one chemistry, the residues they were reported on, and the residues they may occupy.
        /// </summary>
        private sealed class ChemistryGroup
        {
            public List<Modification> Members { get; } = new();
            public List<int> ReportedResidues { get; } = new();
            public List<(int Residue, Modification Modification)> Allowed { get; } = new();
        }
    }
}
