#nullable enable
using System;
using Omics.Fragmentation;
using System.Collections.Generic;
using Omics;
using System.Linq;
using Omics.Modifications;

namespace EngineLayer;

public interface IGptmdFilter : IEquatable<IGptmdFilter>
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

    bool IEquatable<IGptmdFilter>.Equals(IGptmdFilter? other)
    {
        return other != null && GetType() == other.GetType();
    }
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
    // TODO: Consider N-terminal and C-terminal mods as special cases, where the flanking ion not possible except for the M ion. 
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

