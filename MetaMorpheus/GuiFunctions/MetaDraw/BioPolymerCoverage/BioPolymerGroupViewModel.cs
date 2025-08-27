using System.Collections.Generic;
using System.Collections.ObjectModel;
using Readers;

namespace GuiFunctions;

public enum BioPolymerCoverageType
{
    Unique,
    UniqueMissedCleavage, 
    TandemRepeat,               // Shared within the same biopolymer
    TandemRepeatMissedCleavage,
    Shared,                     // Shared between different biopolymers
    SharedMissedCleavage,
}


public class BioPolymerGroupViewModel(string accession, string proteinName, string sequence, IEnumerable<BioPolymerCoverageResultModel> results) : BaseViewModel
{
    public string Accession { get; } = accession;
    public string ProteinName { get; } = proteinName;
    public string Sequence { get; } = sequence;
    public int Length => Sequence.Length;
    public int GroupCount => CoverageResults.Count;
    public double SequenceCoverage { get; } = 0; // TODO: Calculate this

    public ObservableCollection<BioPolymerCoverageResultModel> CoverageResults { get; } = new(results);
}

public class BioPolymerCoverageResultModel(SpectrumMatchFromTsv match, string baseSequence, int start, int end, BioPolymerCoverageType coverageType)
{
    public string BaseSequence { get; } = baseSequence;
    public int Start { get; } = start;
    public int End { get; } = end;
    public BioPolymerCoverageType CoverageType { get; } = coverageType;
    public double QValue => match.QValue;
    public string Ambiguity => match.AmbiguityLevel;
    public string FileName => match.FileName;
}

