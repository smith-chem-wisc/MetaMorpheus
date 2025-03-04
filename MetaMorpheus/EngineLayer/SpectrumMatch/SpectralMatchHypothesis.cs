#nullable enable
using Omics;
using System.Collections.Generic;
using Omics.Fragmentation;
using System;

namespace EngineLayer.SpectrumMatch;

/// <summary>
/// Class to represent a possible spectral match
/// </summary>
/// <param name="notch"></param>
/// <param name="pwsm"></param>
/// <param name="matchedIons"></param>
public class SpectralMatchHypothesis(int notch, IBioPolymerWithSetMods pwsm, List<MatchedFragmentIon> matchedIons) : IEquatable<SpectralMatchHypothesis>
{
    public readonly int Notch = notch;
    public readonly IBioPolymerWithSetMods WithSetMods = pwsm;
    public readonly List<MatchedFragmentIon> MatchedIons = matchedIons;

    public bool IsDecoy => WithSetMods.Parent.IsDecoy;
    public string FullSequence => WithSetMods.FullSequence;

    public bool Equals(SpectralMatchHypothesis? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        return Notch == other.Notch && WithSetMods.Equals(other.WithSetMods) && Equals(MatchedIons.Count, other.MatchedIons.Count);
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
        return HashCode.Combine(Notch, WithSetMods, MatchedIons);
    }
}