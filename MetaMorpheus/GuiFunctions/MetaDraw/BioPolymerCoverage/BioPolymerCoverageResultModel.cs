using Readers;

namespace GuiFunctions;

public class BioPolymerCoverageResultModel(SpectrumMatchFromTsv match, string baseSequence, int start, int end, BioPolymerCoverageType coverageType)
{
    public readonly SpectrumMatchFromTsv Match = match;
    public string BaseSequence { get; } = baseSequence;
    public int Start { get; } = start;
    public int End { get; } = end;
    public BioPolymerCoverageType CoverageType { get; } = coverageType;
    public double QValue => Match.QValue;
    public string Ambiguity => Match.AmbiguityLevel;
    public string FileName => Match.FileName;
}