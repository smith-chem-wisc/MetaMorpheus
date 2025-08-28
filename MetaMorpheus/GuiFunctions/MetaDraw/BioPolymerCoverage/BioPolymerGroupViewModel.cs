using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
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


public class BioPolymerGroupViewModel : BaseViewModel
{
    public string Accession { get; }
    public string ProteinName { get; }
    public string Sequence { get; }
    public int Length => Sequence.Length;
    public int GroupCount { get; private set; }
    public double UniqueSequenceCoverage { get; private set; }
    public double MaximumSequenceCoverage { get; private set; }
    public ObservableCollection<BioPolymerCoverageResultModel> CoverageResults { get; }

    public BioPolymerGroupViewModel(string accession, string proteinName, string sequence, IEnumerable<BioPolymerCoverageResultModel> results)
    {
        Accession = accession;
        ProteinName = proteinName;
        Sequence = sequence;
        CoverageResults = new(results);
        UpdatePropertiesAfterFilter();
    }

    public void UpdatePropertiesAfterFilter()
    {
        int count = 0;
        bool[] uniqueCovered = new bool[Sequence.Length];
        bool[] maxCovered = new bool[Sequence.Length];

        foreach (var res in CoverageResults.Where(p => MetaDrawSettings.FilterAcceptsPsm(p.Match)))
        {
            count++;
            int start = Math.Max(0, res.Start - 1);
            int end = Math.Min(Sequence.Length - 1, res.End - 1);

            // Maximum coverage: all results
            for (int i = start; i <= end; i++)
                maxCovered[i] = true;

            // Unique coverage: only Unique or UniqueMissedCleavage
            if (res.CoverageType == BioPolymerCoverageType.Unique || res.CoverageType == BioPolymerCoverageType.UniqueMissedCleavage)
            {
                for (int i = start; i <= end; i++)
                    uniqueCovered[i] = true;
            }
        }

        GroupCount = count;
        UniqueSequenceCoverage = uniqueCovered.Count(b => b) / (double)Sequence.Length;
        MaximumSequenceCoverage = maxCovered.Count(b => b) / (double)Sequence.Length;
        OnPropertyChanged(nameof(GroupCount));
        OnPropertyChanged(nameof(MaximumSequenceCoverage));
        OnPropertyChanged(nameof(MaximumSequenceCoverage));
    }
}

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

