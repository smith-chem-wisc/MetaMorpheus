#nullable enable
using Omics;
using System.Collections.Generic;
using Omics.Fragmentation;
using System;

namespace EngineLayer.SpectrumMatch;

/// <summary>
/// Class to represent a possible spectral match to be analyzed at a later point
/// </summary>
/// <param name="notch"></param>
/// <param name="pwsm"></param>
/// <param name="matchedIons"></param>
public class SpectralMatchHypothesis(int notch, IBioPolymerWithSetMods pwsm, List<MatchedFragmentIon> matchedIons, double score) 
    : IEquatable<SpectralMatchHypothesis>, ISearchAttempt
{
    public double Score { get; } = score;
    public int Notch { get; } = notch;
    public readonly IBioPolymerWithSetMods SpecificBioPolymer = pwsm;
    public readonly List<MatchedFragmentIon> MatchedIons = matchedIons;

    public bool IsDecoy => SpecificBioPolymer.Parent.IsDecoy;
    public string FullSequence => SpecificBioPolymer.FullSequence;

    #region Notch Ambiguous Q-Value Properties

    // If we have an ambiguous spectral match that contains multiple notches, we store the q-values and cumulative counts for each notch here.
    // We can then use this information to disambiguate at the end. 
    // If we end up with a notch ambiguous PSM and these fields are not set after FDR analysis engine, something went wrong. 

    public double? QValueNotch { get; set; } = null;
    public double? CumulativeTargetNotch { get; set; } = null;
    public double? CumulativeDecoyNotch { get; set; } = null;

    public double? PeptideQValueNotch { get; set; } = null;
    public double? PeptideCumulativeTargetNotch { get; set; } = null;
    public double? PeptideCumulativeDecoyNotch { get; set; } = null;

    #endregion

    public bool Equals(ISearchAttempt? other)
    {
        if (other is SpectralMatchHypothesis hypothesis) return Equals(hypothesis);
        if (other is null) return false;
        return IsDecoy == other.IsDecoy
               && Notch == other.Notch
               && Math.Abs(Score - other.Score) < SpectralMatch.ToleranceForScoreDifferentiation;
    }

    public bool Equals(SpectralMatchHypothesis? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        return Notch == other.Notch
            && IsDecoy == other.IsDecoy
            && Equals(MatchedIons.Count, other.MatchedIons.Count)
            && Math.Abs(Score - other.Score) < SpectralMatch.ToleranceForScoreDifferentiation
            && SpecificBioPolymer.Equals(other.SpecificBioPolymer);
    }

    public override bool Equals(object? obj)
    {
        if (obj is null) return false;
        if (ReferenceEquals(this, obj)) return true;
        if (obj.GetType() != GetType()) return false;
        return Equals((SpectralMatchHypothesis)obj);
    }

    public override int GetHashCode()
    {
        return HashCode.Combine(SpecificBioPolymer, MatchedIons, Score, Notch);
    }
}