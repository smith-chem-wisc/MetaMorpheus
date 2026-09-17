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
            m.NeutralTheoreticalProduct.ResiduePosition < site)
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
             m.NeutralTheoreticalProduct.ResiduePosition < site)
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

        bool leftFlank = matchedIons.Any(m =>
            m.NeutralTheoreticalProduct.ResiduePosition == site - 1);

        bool rightFlank = matchedIons.Any(m =>
            m.NeutralTheoreticalProduct.ResiduePosition == site);

        if (modAttemptingToAdd.LocationRestriction.Contains("terminal", StringComparison.InvariantCultureIgnoreCase))
            return leftFlank || rightFlank;

        return leftFlank && rightFlank;
    }
}


/// <summary>
/// Rejects an amino acid substitution that moves the digestion agent's cut sites. Such a
/// substitution changes where the protein is cleaved, so the peptide actually produced is not the
/// one being scored and its mass is not the unmodified peptide's mass plus the substitution delta.
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

        DigestionAgent? agent = candidatePeptide?.DigestionParams?.DigestionAgent;
        if (agent?.DigestionMotifs is not { Count: > 0 })
            return true;

        // GptmdEngine passes the site in "one is N-terminus" space, so residue one arrives as two.
        int indexInPeptide = peptideOneBasedModSite - 2;
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

        // GetDigestionSiteIndices builds its result from a hash set, so compare as sets.
        return new HashSet<int>(agent.GetDigestionSiteIndices(window))
            .SetEquals(agent.GetDigestionSiteIndices(new string(substitutedWindow)));
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
