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

    public ObservableCollection<BioPolymerCoverageResultModel> CoverageResults { get; } = new(results);
}

public class BioPolymerCoverageResultModel : BaseViewModel
{
    private SpectrumMatchFromTsv _spectrumMatch;

    public BioPolymerCoverageResultModel(SpectrumMatchFromTsv match, string baseSequence, int start, int end, BioPolymerCoverageType coverageType)
    {
        _spectrumMatch = match;
        BaseSequence = baseSequence;
        Start = start;
        End = end;
        CoverageType = coverageType;
    }

    public string BaseSequence { get; }
    public int Start { get; }
    public int End { get; }
    public BioPolymerCoverageType CoverageType { get; }
    public double QValue => _spectrumMatch.QValue;
    public string Ambiguity => _spectrumMatch.AmbiguityLevel;
    public string FileName => _spectrumMatch.FileName;
}

