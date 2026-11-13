#nullable enable
using System;
using System.Collections.Generic;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;

namespace TaskLayer.ParallelSearch.Analysis;

/// <summary>
/// The value and threshold a match is held to when the collectors count it as confident. The choice is the one
/// <see cref="FilteredPsms"/> makes for a normal search: the PEP q-value at the PEP q-value threshold when PEP filtering
/// is asked for and can be done, and the q-value otherwise (at the threshold FilteredPsms falls back to when PEP
/// filtering is asked for but cannot be done).
/// </summary>
public sealed class ConfidenceCutoff
{
    public bool UsesPepQValue { get; }
    public double Threshold { get; }
    public bool PeptideLevel { get; }

    private ConfidenceCutoff(bool usesPepQValue, double threshold, bool peptideLevel)
    {
        UsesPepQValue = usesPepQValue;
        Threshold = threshold;
        PeptideLevel = peptideLevel;
    }

    /// <summary>
    /// Asks FilteredPsms which value and threshold it would filter these matches on.
    /// </summary>
    public static ConfidenceCutoff For(IEnumerable<SpectralMatch> matches, CommonParameters commonParameters, bool peptideLevel)
    {
        // without PEP filtering FilteredPsms always filters on the q-value at the smaller threshold; skip the work
        if (!PepFilteringRequested(commonParameters))
        {
            return new ConfidenceCutoff(false,
                Math.Min(commonParameters.QValueThreshold, commonParameters.PepQValueThreshold), peptideLevel);
        }

        var filtered = FilteredPsms.Filter(matches, commonParameters,
            includeDecoys: true,
            includeContaminants: true,
            includeAmbiguous: true,
            includeHighQValuePsms: false,
            filterAtPeptideLevel: peptideLevel);
        return new ConfidenceCutoff(filtered.FilterType == FilterType.PepQValue, filtered.FilterThreshold, peptideLevel);
    }

    /// <summary>
    /// The protein group q-value threshold a normal search writes protein groups at, given how its PSMs were filtered:
    /// the PEP q-value threshold under PEP filtering, the fallback threshold when PEP filtering was asked for but not
    /// done, and the q-value threshold otherwise.
    /// </summary>
    public static double ProteinGroupThreshold(IEnumerable<SpectralMatch> psms, CommonParameters commonParameters)
    {
        if (!PepFilteringRequested(commonParameters))
        {
            return commonParameters.QValueThreshold;
        }

        var filtered = FilteredPsms.Filter(psms, commonParameters,
            includeDecoys: true,
            includeContaminants: true,
            includeAmbiguous: true,
            includeHighQValuePsms: false);
        return filtered.FilterType switch
        {
            FilterType.PepQValue => commonParameters.PepQValueThreshold,
            _ when filtered.PepFilteringNotDone != null && !filtered.FilteringNotPerformed => filtered.FilterThreshold,
            _ => commonParameters.QValueThreshold
        };
    }

    private static bool PepFilteringRequested(CommonParameters commonParameters) =>
        commonParameters.PepQValueThreshold < commonParameters.QValueThreshold;

    public bool Passes(SpectralMatch match) => Passes(match.GetFdrInfo(PeptideLevel));

    // on the q-value path FilteredPsms also holds the notch q-value to the threshold
    public bool Passes(FdrInfo? fdrInfo) =>
        fdrInfo != null && (UsesPepQValue
            ? fdrInfo.PEP_QValue <= Threshold
            : fdrInfo.QValue <= Threshold && fdrInfo.QValueNotch <= Threshold);
}
