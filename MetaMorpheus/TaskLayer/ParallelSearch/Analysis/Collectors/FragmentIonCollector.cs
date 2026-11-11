using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using MathNet.Numerics.Statistics;

namespace TaskLayer.ParallelSearch.Analysis.Collectors;
public class FragmentIonCollector : IMetricCollector
{
    // PSM Column Names
    public const string PSM_LongestIonSeriesBidirectionalTargets = "PSM_LongestIonSeriesBidirectionalTargets";
    public const string PSM_ComplementaryIonCountTargets = "PSM_ComplementaryIonCountTargets";
    public const string PSM_SequenceCoverageFractionTargets = "PSM_SequenceCoverageFractionTargets";
    public const string PSM_LongestIonSeriesBidirectionalDecoys = "PSM_LongestIonSeriesBidirectionalDecoys";
    public const string PSM_ComplementaryIonCountDecoys = "PSM_ComplementaryIonCountDecoys";
    public const string PSM_SequenceCoverageFractionDecoys = "PSM_SequenceCoverageFractionDecoys";
    public const string PSM_LongestIonSeriesBidirectional_AllTargets = "PSM_LongestIonSeriesBidirectional_AllTargets";
    public const string PSM_ComplementaryIonCount_AllTargets = "PSM_ComplementaryIonCount_AllTargets";
    public const string PSM_SequenceCoverageFraction_AllTargets = "PSM_SequenceCoverageFraction_AllTargets";
    public const string PSM_LongestIonSeriesBidirectional_AllDecoys = "PSM_LongestIonSeriesBidirectional_AllDecoys";
    public const string PSM_ComplementaryIonCount_AllDecoys = "PSM_ComplementaryIonCount_AllDecoys";
    public const string PSM_SequenceCoverageFraction_AllDecoys = "PSM_SequenceCoverageFraction_AllDecoys";
    public const string PSM_FragmentPPMErrors = "PSM_PPMErrors";

    // Peptide Column Names
    public const string Peptide_LongestIonSeriesBidirectionalTargets = "Peptide_LongestIonSeriesBidirectionalTargets";
    public const string Peptide_ComplementaryIonCountTargets = "Peptide_ComplementaryIonCountTargets";
    public const string Peptide_SequenceCoverageFractionTargets = "Peptide_SequenceCoverageFractionTargets";
    public const string Peptide_LongestIonSeriesBidirectionalDecoys = "Peptide_LongestIonSeriesBidirectionalDecoys";
    public const string Peptide_ComplementaryIonCountDecoys = "Peptide_ComplementaryIonCountDecoys";
    public const string Peptide_SequenceCoverageFractionDecoys = "Peptide_SequenceCoverageFractionDecoys";
    public const string Peptide_LongestIonSeriesBidirectional_AllTargets = "Peptide_LongestIonSeriesBidirectional_AllTargets";
    public const string Peptide_ComplementaryIonCount_AllTargets = "Peptide_ComplementaryIonCount_AllTargets";
    public const string Peptide_SequenceCoverageFraction_AllTargets = "Peptide_SequenceCoverageFraction_AllTargets";
    public const string Peptide_LongestIonSeriesBidirectional_AllDecoys = "Peptide_LongestIonSeriesBidirectional_AllDecoys";
    public const string Peptide_ComplementaryIonCount_AllDecoys = "Peptide_ComplementaryIonCount_AllDecoys";
    public const string Peptide_SequenceCoverageFraction_AllDecoys = "Peptide_SequenceCoverageFraction_AllDecoys";
    public const string Peptide_FragmentPPMErrors = "Peptide_PPMErrors";

    public string CollectorName => "FragmentIons";
    
    public IEnumerable<string> GetOutputColumns()
    {
        yield return PSM_LongestIonSeriesBidirectionalTargets;
        yield return PSM_ComplementaryIonCountTargets;
        yield return PSM_SequenceCoverageFractionTargets;
        yield return PSM_LongestIonSeriesBidirectionalDecoys;
        yield return PSM_ComplementaryIonCountDecoys;
        yield return PSM_SequenceCoverageFractionDecoys;
        yield return PSM_LongestIonSeriesBidirectional_AllTargets;
        yield return PSM_ComplementaryIonCount_AllTargets;
        yield return PSM_SequenceCoverageFraction_AllTargets;
        yield return PSM_LongestIonSeriesBidirectional_AllDecoys;
        yield return PSM_ComplementaryIonCount_AllDecoys;
        yield return PSM_SequenceCoverageFraction_AllDecoys;

        yield return Peptide_LongestIonSeriesBidirectionalTargets;
        yield return Peptide_ComplementaryIonCountTargets;
        yield return Peptide_SequenceCoverageFractionTargets;
        yield return Peptide_LongestIonSeriesBidirectionalDecoys;
        yield return Peptide_ComplementaryIonCountDecoys;
        yield return Peptide_SequenceCoverageFractionDecoys;
        yield return Peptide_LongestIonSeriesBidirectional_AllTargets;
        yield return Peptide_ComplementaryIonCount_AllTargets;
        yield return Peptide_SequenceCoverageFraction_AllTargets;
        yield return Peptide_LongestIonSeriesBidirectional_AllDecoys;
        yield return Peptide_ComplementaryIonCount_AllDecoys;
        yield return Peptide_SequenceCoverageFraction_AllDecoys;
    }

    public Dictionary<string, object> CollectData(TransientDatabaseContext context)
    {
        double qValueThreshold = Math.Min(context.CommonParameters.QValueThreshold, context.CommonParameters.PepQValueThreshold);

        // Psms
        var confidentPsms = context.TransientPsms
            .Where(p => p.FdrInfo.QValue <= qValueThreshold)
            .ToList();

        var psmPpmErrors = confidentPsms.SelectMany(p => p.MatchedFragmentIons.Select(m => m.MassErrorPpm)).ToArray();

        foreach (var psm in confidentPsms)
            if (psm.FragmentCoveragePositionInPeptide is null)
                psm.GetAminoAcidCoverage();

        var psmBidirectional = confidentPsms
            .Select(p => (p.IsDecoy,
                (double)SpectralMatch.GetLongestIonSeriesBidirectional(p.BestMatchingBioPolymersWithSetMods.First()) / p.BestMatchingBioPolymersWithSetMods.Average(pep => pep.SpecificBioPolymer.Length)))
            .ToList();

        var psmComplementaryCounts = confidentPsms
            .Select(p => (p.IsDecoy,
                (double)SpectralMatch.GetCountComplementaryIons(p.BestMatchingBioPolymersWithSetMods.First()) / p.BestMatchingBioPolymersWithSetMods.Average(pep => pep.SpecificBioPolymer.Length)))
            .ToList();

        var psmSequenceCoverage = confidentPsms
            .Select(p => (p.IsDecoy, p.FragmentCoveragePositionInPeptide is null ? 0 :
                p.FragmentCoveragePositionInPeptide.Count(p => p != -1) / (double)p.BaseSequence.Length))
            .ToList();

        // peptides
        var confidentPeptides = context.TransientPeptides
            .Where(p => p.FdrInfo.QValue <= qValueThreshold)
            .ToList();

        var peptidePpmErrors = confidentPeptides.SelectMany(p => p.MatchedFragmentIons.Select(m => m.MassErrorPpm)).ToArray();

        foreach (var psm in confidentPeptides)
            if (psm.FragmentCoveragePositionInPeptide is null)
                psm.GetAminoAcidCoverage();

        var peptidesBidirectional = confidentPeptides
            .Select(p => (p.IsDecoy,
                (double)SpectralMatch.GetLongestIonSeriesBidirectional(p.BestMatchingBioPolymersWithSetMods.First()) / p.BestMatchingBioPolymersWithSetMods.Average(pep => pep.SpecificBioPolymer.Length)))
            .ToList();

        var peptidesComplementaryCounts = confidentPeptides
            .Select(p => (p.IsDecoy,
                (double)SpectralMatch.GetCountComplementaryIons(p.BestMatchingBioPolymersWithSetMods.First()) / p.BestMatchingBioPolymersWithSetMods.Average(pep => pep.SpecificBioPolymer.Length)))
            .ToList();

        var peptidesSequenceCoverage = confidentPeptides
            .Select(p => (p.IsDecoy, p.FragmentCoveragePositionInPeptide is null ? 0 :
                p.FragmentCoveragePositionInPeptide.Count(p => p != -1) / (double)p.BaseSequence.Length))
            .ToList();

        return new Dictionary<string, object>
        {
            [PSM_LongestIonSeriesBidirectionalTargets] = psmBidirectional.Where(p => !p.Item1).Select(p => p.Item2).Median(),
            [PSM_ComplementaryIonCountTargets] = psmComplementaryCounts.Where(p => !p.Item1).Select(p => p.Item2).Median(),
            [PSM_SequenceCoverageFractionTargets] = psmSequenceCoverage.Where(p => !p.Item1).Select(p => p.Item2).Median(),
            [PSM_LongestIonSeriesBidirectionalDecoys] = psmBidirectional.Where(p => p.Item1).Select(p => p.Item2).Median(),
            [PSM_ComplementaryIonCountDecoys] = psmComplementaryCounts.Where(p => p.Item1).Select(p => p.Item2).Median(),
            [PSM_SequenceCoverageFractionDecoys] = psmSequenceCoverage.Where(p => p.Item1).Select(p => p.Item2).Median(),
            [PSM_LongestIonSeriesBidirectional_AllTargets] = psmBidirectional.Where(p => !p.Item1).Select(p => p.Item2).ToArray(),
            [PSM_ComplementaryIonCount_AllTargets] = psmComplementaryCounts.Where(p => !p.Item1).Select(p => p.Item2).ToArray(),
            [PSM_SequenceCoverageFraction_AllTargets] = psmSequenceCoverage.Where(p => !p.Item1).Select(p => p.Item2).ToArray(),
            [PSM_LongestIonSeriesBidirectional_AllDecoys] = psmBidirectional.Where(p => p.Item1).Select(p => p.Item2).ToArray(),
            [PSM_ComplementaryIonCount_AllDecoys] = psmComplementaryCounts.Where(p => p.Item1).Select(p => p.Item2).ToArray(),
            [PSM_SequenceCoverageFraction_AllDecoys] = psmSequenceCoverage.Where(p => p.Item1).Select(p => p.Item2).ToArray(),
            [PSM_FragmentPPMErrors] = psmPpmErrors,

            [Peptide_LongestIonSeriesBidirectionalTargets] = peptidesBidirectional.Where(p => !p.Item1).Select(p => p.Item2).Median(),
            [Peptide_ComplementaryIonCountTargets] = peptidesComplementaryCounts.Where(p => !p.Item1).Select(p => p.Item2).Median(),
            [Peptide_SequenceCoverageFractionTargets] = peptidesSequenceCoverage.Where(p => !p.Item1).Select(p => p.Item2).Median(),
            [Peptide_LongestIonSeriesBidirectionalDecoys] = peptidesBidirectional.Where(p => p.Item1).Select(p => p.Item2).Median(),
            [Peptide_ComplementaryIonCountDecoys] = peptidesComplementaryCounts.Where(p => p.Item1).Select(p => p.Item2).Median(),
            [Peptide_SequenceCoverageFractionDecoys] = peptidesSequenceCoverage.Where(p => p.Item1).Select(p => p.Item2).Median(),
            [Peptide_LongestIonSeriesBidirectional_AllTargets] = peptidesBidirectional.Where(p => !p.Item1).Select(p => p.Item2).ToArray(),
            [Peptide_ComplementaryIonCount_AllTargets] = peptidesComplementaryCounts.Where(p => !p.Item1).Select(p => p.Item2).ToArray(),
            [Peptide_SequenceCoverageFraction_AllTargets] = peptidesSequenceCoverage.Where(p => !p.Item1).Select(p => p.Item2).ToArray(),
            [Peptide_LongestIonSeriesBidirectional_AllDecoys] = peptidesBidirectional.Where(p => p.Item1).Select(p => p.Item2).ToArray(),
            [Peptide_ComplementaryIonCount_AllDecoys] = peptidesComplementaryCounts.Where(p => p.Item1).Select(p => p.Item2).ToArray(),
            [Peptide_SequenceCoverageFraction_AllDecoys] = peptidesSequenceCoverage.Where(p => p.Item1).Select(p => p.Item2).ToArray(),
            [Peptide_FragmentPPMErrors] = peptidePpmErrors
        };
    }

    public bool CanCollectData(TransientDatabaseContext context)
    {
        return context.TransientPeptides != null && context.TransientPsms != null;
    }
}
