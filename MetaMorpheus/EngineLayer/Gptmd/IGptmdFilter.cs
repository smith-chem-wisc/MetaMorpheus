#nullable enable
using System;
using Omics.Fragmentation;
using System.Collections.Generic;
using Omics;
using System.Linq;
using Omics.Digestion;
using Omics.Modifications;

namespace EngineLayer;

public interface IGptmdFilter
{
    public static string GetFilterTypeName(IGptmdFilter filter) => filter.GetType().Name;

    bool Passes(
        IBioPolymerWithSetMods candidatePeptide,
        SpectralMatch psm,
        double newScore,
        double originalScore,
        List<MatchedFragmentIon> matchedIons,
        int peptideOneBasedModSite,
        int peptideLength, 
        Modification modAttemptingToAdd);
}

/// <summary>
/// Requires that the new score is greater than the original score.
/// </summary>
public sealed class ImprovedScoreFilter : IGptmdFilter
{
    public bool Passes(
        IBioPolymerWithSetMods candidatePeptide,
        SpectralMatch psm,
        double newScore,
        double originalScore,
        List<MatchedFragmentIon> matchedIons,
        int peptideOneBasedModSite,
        int peptideLength,
        Modification modAttemptingToAdd)
    {
        return newScore > originalScore;
    }
}

/// <summary>
/// Requires the mod site to be covered by at least one N-terminal and one C-terminal ion.
/// That is, ions from both directions must include the mod, even if not flanking it.
/// </summary>
public sealed class DualDirectionalIonCoverageFilter : IGptmdFilter
{
    public bool Passes(
        IBioPolymerWithSetMods candidatePeptide,
        SpectralMatch psm,
        double newScore,
        double originalScore,
        List<MatchedFragmentIon> matchedIons,
        int peptideOneBasedModSite,
        int peptideLength,
        Modification modAttemptingToAdd)
    {
        if (matchedIons == null || matchedIons.Count == 0)
            return false;

        int site = peptideOneBasedModSite;

        bool coveredFromNTerm = matchedIons.Any(m =>
           m.NeutralTheoreticalProduct.ProductType == ProductType.M ||
           (m.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.N or FragmentationTerminus.FivePrime &&
            m.NeutralTheoreticalProduct.ResiduePosition >= site)
        );

        bool coveredFromCTerm = matchedIons.Any(m =>
           m.NeutralTheoreticalProduct.ProductType == ProductType.M ||
           (m.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.C or FragmentationTerminus.ThreePrime &&
            m.NeutralTheoreticalProduct.ResiduePosition <= site)
        );

        if (modAttemptingToAdd.LocationRestriction.Contains("terminal", StringComparison.InvariantCultureIgnoreCase))
            return coveredFromCTerm || coveredFromNTerm;

        return coveredFromNTerm && coveredFromCTerm;
    }
}

/// <summary>
/// Requires the mod site to be covered by at least one N-terminal or one C-terminal ion.
/// That is, ions from one direction must include the mod, even if not flanking it.
/// </summary>
public sealed class UniDirectionalIonCoverageFilter : IGptmdFilter
{
    public bool Passes(
        IBioPolymerWithSetMods candidatePeptide,
        SpectralMatch psm,
        double newScore,
        double originalScore,
        List<MatchedFragmentIon> matchedIons,
        int peptideOneBasedModSite,
        int peptideLength,
        Modification modAttemptingToAdd)
    {
        if (matchedIons == null || matchedIons.Count == 0)
            return false;

        int site = peptideOneBasedModSite;

        bool coveredFromNTerm = matchedIons.Any(m =>
            m.NeutralTheoreticalProduct.ProductType == ProductType.M ||
            (m.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.N or FragmentationTerminus.FivePrime &&
             m.NeutralTheoreticalProduct.ResiduePosition >= site)
        );

        bool coveredFromCTerm = matchedIons.Any(m =>
            m.NeutralTheoreticalProduct.ProductType == ProductType.M ||
            (m.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.C or FragmentationTerminus.ThreePrime &&
             m.NeutralTheoreticalProduct.ResiduePosition <= site)
        );

        return coveredFromNTerm || coveredFromCTerm;
    }
}

/// <summary>
/// Requires flanking ions — a fragment from *before* and one from *after* the mod site,
/// regardless of fragmentation direction.
/// </summary>
public sealed class FlankingIonCoverageFilter : IGptmdFilter
{
    public bool Passes(
        IBioPolymerWithSetMods candidatePeptide,
        SpectralMatch psm,
        double newScore,
        double originalScore,
        List<MatchedFragmentIon> matchedIons,
        int peptideOneBasedModSite,
        int peptideLength,
        Modification modAttemptingToAdd)
    {
        if (matchedIons == null || matchedIons.Count == 0)
            return false;

        int site = peptideOneBasedModSite;

        bool leftFlank = matchedIons.Any(m => CleavesAt(m, site - 1));

        bool rightFlank = matchedIons.Any(m => CleavesAt(m, site));

        if (modAttemptingToAdd.LocationRestriction.Contains("terminal", StringComparison.InvariantCultureIgnoreCase))
            return leftFlank || rightFlank;

        return leftFlank && rightFlank;
    }

    private static bool CleavesAt(MatchedFragmentIon ion, int position)
    {
        Product product = ion.NeutralTheoreticalProduct;
        return product.Terminus switch
        {
            FragmentationTerminus.N or FragmentationTerminus.FivePrime => product.ResiduePosition == position,
            FragmentationTerminus.C or FragmentationTerminus.ThreePrime => product.ResiduePosition - 1 == position,
            _ => false
        };
    }
}

/// <summary>
/// Rejects an amino acid substitution that removes a cut site defining one of the peptide's termini
/// or creates one inside it. Either way the peptide being scored is not produced, or not produced
/// cleanly, so its mass is not a reliable unmodified-mass-plus-substitution-delta hypothesis.
/// Modifications that are not substitutions always pass.
/// </summary>
public sealed class CleavageSiteFilter : IGptmdFilter
{
    public bool Passes(
        IBioPolymerWithSetMods candidatePeptide,
        SpectralMatch psm,
        double newScore,
        double originalScore,
        List<MatchedFragmentIon> matchedIons,
        int peptideOneBasedModSite,
        int peptideLength,
        Modification modAttemptingToAdd)
    {
        if (!TryGetSubstitution(modAttemptingToAdd, out char original, out char substituted))
            return true;

        IDigestionParams? digestionParams = candidatePeptide?.DigestionParams;
        DigestionAgent? agent = digestionParams?.DigestionAgent;
        if (agent?.DigestionMotifs is not { Count: > 0 })
            return true;

        // Only a fully specific search requires the peptide's termini to be cleavage sites. In a
        // semi- or non-specific search the peptide can exist wherever the substitution puts the
        // sites, so moving one says nothing about whether the candidate is real.
        if (digestionParams!.SearchModeType != CleavageSpecificity.Full)
            return true;

        int indexInPeptide = peptideOneBasedModSite - 1;
        if (indexInPeptide < 0 || indexInPeptide >= candidatePeptide!.BaseSequence.Length)
            return true;

        string parentSequence = candidatePeptide.Parent?.BaseSequence ?? candidatePeptide.BaseSequence;
        int peptideStart = candidatePeptide.Parent is null ? 0 : candidatePeptide.OneBasedStartResidue - 1;
        int indexInParent = peptideStart + indexInPeptide;
        if (indexInParent < 0 || indexInParent >= parentSequence.Length
            || parentSequence[indexInParent] != original)
            return true;

        // A cut decision reads at most InducingCleavage + PreventingCleavage residues around its own
        // location, so a window that wide past each end of the peptide contains every site the
        // substitution could move. Doubling it keeps the window's own edges out of the comparison.
        int radius = 2 * agent.DigestionMotifs.Max(m =>
            m.InducingCleavage.Length + (m.PreventingCleavage?.Length ?? 0));
        int windowStart = Math.Max(0, peptideStart - radius);
        int windowEnd = Math.Min(parentSequence.Length,
            peptideStart + candidatePeptide.BaseSequence.Length + radius);

        string window = parentSequence.Substring(windowStart, windowEnd - windowStart);
        char[] substitutedWindow = window.ToCharArray();
        substitutedWindow[indexInParent - windowStart] = substituted;

        var sitesBefore = new HashSet<int>(agent.GetDigestionSiteIndices(window));
        var sitesAfter = new HashSet<int>(agent.GetDigestionSiteIndices(new string(substitutedWindow)));
        int start = peptideStart - windowStart;
        int end = start + candidatePeptide.BaseSequence.Length;

        // Reject losing a cut that defines a terminus, or gaining one inside the peptide. Losing an
        // internal missed cleavage still yields this peptide, so it passes.
        bool terminusLost = (sitesBefore.Contains(start) && !sitesAfter.Contains(start))
            || (sitesBefore.Contains(end) && !sitesAfter.Contains(end));
        bool internalSiteGained = sitesAfter.Any(s => s > start && s < end && !sitesBefore.Contains(s));
        return !terminusLost && !internalSiteGained;
    }

    /// <summary>
    /// Identifies substitutions the way mzLib does when it converts them to sequence variants:
    /// a "nucleotide substitution" modification type with an "X-&gt;Y" identifier.
    /// </summary>
    private static bool TryGetSubstitution(Modification mod, out char original, out char substituted)
    {
        original = default;
        substituted = default;

        if (mod?.ModificationType is null || mod.OriginalId is null
            || !mod.ModificationType.Contains("nucleotide substitution"))
            return false;

        string[] residues = mod.OriginalId.Split("->", StringSplitOptions.RemoveEmptyEntries);
        if (residues.Length != 2 || residues[0].Length != 1 || residues[1].Length != 1)
            return false;

        original = residues[0][0];
        substituted = residues[1][0];
        return true;
    }
}
