using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using EngineLayer.SpectrumMatch;
using Omics;
using Omics.Fragmentation;

namespace EngineLayer;

/// <summary>
/// Placeholder class for when RNA code gets fully incorporated
/// </summary>
[ExcludeFromCodeCoverage]
public class OligoSpectralMatch : SpectralMatch
{
    public OligoSpectralMatch(IBioPolymerWithSetMods peptide, int notch, double score, int scanIndex,
        Ms2ScanWithSpecificMass scan, CommonParameters commonParameters, List<MatchedFragmentIon> matchedFragmentIons) 
        : base(peptide, notch, score, scanIndex, scan, commonParameters, matchedFragmentIons)
    {

    }

    protected OligoSpectralMatch(SpectralMatch psm, List<SpectralMatchHypothesis> bestMatchingPeptides)
        : base(psm, bestMatchingPeptides)
    {
    }
}