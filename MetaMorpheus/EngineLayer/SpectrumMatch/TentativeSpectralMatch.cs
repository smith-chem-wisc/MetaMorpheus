using Omics;
using System.Collections.Generic;
using Omics.Fragmentation;

namespace EngineLayer.SpectrumMatch;

public class TentativeSpectralMatch(int notch, IBioPolymerWithSetMods pwsm, List<MatchedFragmentIon> matchedIons)
{
    public readonly int Notch = notch;
    public readonly IBioPolymerWithSetMods WithSetMods = pwsm;
    public readonly List<MatchedFragmentIon> MatchedIons = matchedIons;

    public bool IsDecoy => WithSetMods.Parent.IsDecoy;
    public string FullSequence => WithSetMods.FullSequence;
}