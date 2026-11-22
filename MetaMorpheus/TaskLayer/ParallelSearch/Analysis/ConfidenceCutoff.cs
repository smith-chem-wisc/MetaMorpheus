#nullable enable
using System;
using System.Collections.Generic;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;

namespace TaskLayer.ParallelSearch.Analysis;

/// <summary>
/// The value and threshold a match is held to when the collectors count it as confident: the search's own
/// <see cref="CommonParameters.IdentificationFilter"/> at <see cref="CommonParameters.FilterThreshold"/>, judged as
/// <see cref="FilteredPsms.Passes"/> judges it for a normal search (under the PEP q-value filter, a match whose PEP was
/// not trained is judged by its q-value notch).
/// </summary>
public sealed class ConfidenceCutoff
{
    public FilterType FilterType { get; }
    public bool UsesPepQValue => FilterType == FilterType.PepQValue;
    public double Threshold { get; }
    public bool PeptideLevel { get; }

    private ConfidenceCutoff(FilterType filterType, double threshold, bool peptideLevel)
    {
        FilterType = filterType;
        Threshold = threshold;
        PeptideLevel = peptideLevel;
    }

    /// <summary>The search's identification filter and threshold, at the PSM or the peptide level.</summary>
    public static ConfidenceCutoff For(IEnumerable<SpectralMatch> matches, CommonParameters commonParameters, bool peptideLevel) =>
        new(commonParameters.IdentificationFilter, commonParameters.FilterThreshold, peptideLevel);

    /// <summary>
    /// The protein group q-value threshold a normal search writes protein groups at: the search's
    /// <see cref="CommonParameters.FilterThreshold"/>, as ProteinScoringAndFdrEngine holds them to it.
    /// </summary>
    public static double ProteinGroupThreshold(IEnumerable<SpectralMatch> psms, CommonParameters commonParameters) =>
        commonParameters.FilterThreshold;

    /// <summary>What the matches were held to, worded as a normal search's results.txt words it, e.g. "q-value notch &lt;= 0.01".</summary>
    public string Label => FilteredPsms.Describe(FilterType, Threshold);

    public bool Passes(SpectralMatch match) => FilteredPsms.Passes(match, PeptideLevel, FilterType, Threshold);

    public bool Passes(FdrInfo? fdrInfo) => fdrInfo != null && Value(fdrInfo) <= Threshold;

    // The same value FilteredPsms judges a match by
    private double Value(FdrInfo fdrInfo) => FilterType switch
    {
        FilterType.PepQValue => fdrInfo.PepTrained ? fdrInfo.PEP_QValue : fdrInfo.QValueNotch,
        FilterType.QValueNotch => fdrInfo.QValueNotch,
        FilterType.QValue => fdrInfo.QValue,
        _ => throw new ArgumentOutOfRangeException(nameof(FilterType), FilterType, null)
    };
}
