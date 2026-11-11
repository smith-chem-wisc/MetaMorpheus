#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using MathNet.Numerics.Statistics;

namespace TaskLayer.ParallelSearch.Analysis.Collectors;

/// <summary>
/// Analyzer for organism-specific metrics (bacterial vs human ambiguity)
/// </summary>
public class PsmPeptideSearchCollector : IMetricCollector
{
    // PSM Column Names
    public const string PsmTargets = "PsmTargets";
    public const string PsmDecoys = "PsmDecoys";
    public const string PsmBacterialTargets = "PsmBacterialTargets";
    public const string PsmBacterialDecoys = "PsmBacterialDecoys";
    public const string PsmBacterialAmbiguous = "PsmBacterialAmbiguous";
    public const string PsmBacterialUnambiguousTargets = "PsmBacterialUnambiguousTargets";
    public const string PsmBacterialUnambiguousDecoys = "PsmBacterialUnambiguousDecoys";
    public const string PsmBacterialUnambiguousTargetScores = "PsmBacterialUnambiguousTargetScores";
    public const string PsmBacterialUnambiguousDecoyScores = "PsmBacterialUnambiguousDecoyScores";
    public const string PsmBacteriaTargetlDeltaScores = "PsmBacterialTargetDeltaScores";
    public const string PsmPrecursorMassErrors = "PsmPrecursorMassErrors";
    public const string PsmPrecursorEnvelopePeakCount = "PsmPrecursorEnvelopePeakCount";
    public const string PsmPrecusorDeconScores = "PsmPrecursorDeconScores";
    public const string PsmPrecusorFractionalIntensity = "PsmPrecursorFractionalIntensity";

    // Peptide Column Names
    public const string PeptideTargets = "PeptideTargets";
    public const string PeptideDecoys = "PeptideDecoys";
    public const string PeptideBacterialTargets = "PeptideBacterialTargets";
    public const string PeptideBacterialDecoys = "PeptideBacterialDecoys";
    public const string PeptideBacterialAmbiguous = "PeptideBacterialAmbiguous";
    public const string PeptideBacterialUnambiguousTargets = "PeptideBacterialUnambiguousTargets";
    public const string PeptideBacterialUnambiguousDecoys = "PeptideBacterialUnambiguousDecoys";
    public const string PeptideBacterialUnambiguousTargetScores = "PeptideBacterialUnambiguousTargetScores";
    public const string PeptideBacterialUnambiguousDecoyScores = "PeptideBacterialUnambiguousDecoyScores";
    public const string PeptideBacterialTargetDeltaScores = "PeptideBacterialTargetDeltaScores";
    public const string PeptidePrecursorMassErrors = "PeptidePrecursorMassErrors";
    public const string PeptidePrecursorEnvelopePeakCount = "PeptidePrecursorEnvelopePeakCount";
    public const string PeptidePrecusorDeconScores = "PeptidePrecursorDeconScores";
    public const string PeptidePrecusorFractionalIntensity = "PeptidePrecursorFractionalIntensity";

    private readonly string _targetOrganism;

    public PsmPeptideSearchCollector(string targetOrganism = "Homo sapiens")
    {
        _targetOrganism = targetOrganism;
    }

    public string CollectorName => "OrganismSpecificity";

    public IEnumerable<string> GetOutputColumns()
    {
        // PSM metrics
        yield return PsmTargets;
        yield return PsmDecoys;
        yield return PsmBacterialTargets;
        yield return PsmBacterialDecoys;
        yield return PsmBacterialAmbiguous;
        yield return PsmBacterialUnambiguousTargets;
        yield return PsmBacterialUnambiguousDecoys;
        yield return PsmBacterialUnambiguousTargetScores;
        yield return PsmBacterialUnambiguousDecoyScores;
        yield return PsmBacteriaTargetlDeltaScores;
        yield return PsmPrecursorMassErrors;
        yield return PsmPrecursorEnvelopePeakCount;
        yield return PsmPrecusorDeconScores;
        yield return PsmPrecusorFractionalIntensity;

        // Peptide metrics
        yield return PeptideTargets;
        yield return PeptideDecoys;
        yield return PeptideBacterialTargets;
        yield return PeptideBacterialDecoys;
        yield return PeptideBacterialAmbiguous;
        yield return PeptideBacterialUnambiguousTargets;
        yield return PeptideBacterialUnambiguousDecoys;
        yield return PeptideBacterialUnambiguousTargetScores;
        yield return PeptideBacterialUnambiguousDecoyScores;
        yield return PeptideBacterialTargetDeltaScores;
        yield return PeptidePrecursorMassErrors;
        yield return PeptidePrecursorEnvelopePeakCount;
        yield return PeptidePrecusorDeconScores;
        yield return PeptidePrecusorFractionalIntensity;
    }

    public bool CanCollectData(TransientDatabaseContext context)
    {
        return context.AllPsms != null
               && context.TransientPsms != null
               && context.AllPeptides != null
               && context.TransientPeptides != null;
    }

    public Dictionary<string, object> CollectData(TransientDatabaseContext context)
    {
        double qValueThreshold = Math.Min(context.CommonParameters.QValueThreshold, context.CommonParameters.PepQValueThreshold);
        var globalPsmMetrics = AnalyzeSpectralMatches(context.AllPsms, qValueThreshold);
        var globalPeptideMetrics = AnalyzeSpectralMatches(context.AllPeptides, qValueThreshold, isPeptideLevel: true);
        var psmMetrics = AnalyzeSpectralMatches(context.TransientPsms, qValueThreshold);
        var peptideMetrics = AnalyzeSpectralMatches(context.TransientPeptides, qValueThreshold, isPeptideLevel: true);

        return new Dictionary<string, object>
        {
            [PsmTargets] = globalPsmMetrics.Targets,
            [PsmDecoys] = globalPsmMetrics.Decoys,
            [PsmBacterialTargets] = psmMetrics.Targets,
            [PsmBacterialDecoys] = psmMetrics.Decoys,
            [PsmBacterialAmbiguous] = psmMetrics.Targets - psmMetrics.UnambiguousTargets,
            [PsmBacterialUnambiguousTargets] = psmMetrics.UnambiguousTargets,
            [PsmBacterialUnambiguousDecoys] = psmMetrics.UnambiguousDecoys,
            [PsmBacterialUnambiguousTargetScores] = psmMetrics.TargetScores.ToArray(),
            [PsmBacterialUnambiguousDecoyScores] = psmMetrics.DecoyScores.ToArray(),
            [PsmBacteriaTargetlDeltaScores] = psmMetrics.DeltaScores.ToArray(),
            [PsmPrecursorMassErrors] = psmMetrics.PrecursorMassErrors.ToArray(),
            [PsmPrecursorEnvelopePeakCount] = psmMetrics.PrecursorPeakCounts.Select(v => (double)v).ToArray(),
            [PsmPrecusorDeconScores] = psmMetrics.PrecursorDeconScores.ToArray(),
            [PsmPrecusorFractionalIntensity] = psmMetrics.PrecursorFractionalIntensity.ToArray(),

            [PeptideTargets] = globalPeptideMetrics.Targets,
            [PeptideDecoys] = globalPeptideMetrics.Decoys,
            [PeptideBacterialTargets] = peptideMetrics.Targets,
            [PeptideBacterialDecoys] = peptideMetrics.Decoys,
            [PeptideBacterialAmbiguous] = peptideMetrics.Targets - peptideMetrics.UnambiguousTargets,
            [PeptideBacterialUnambiguousTargets] = peptideMetrics.UnambiguousTargets,
            [PeptideBacterialUnambiguousDecoys] = peptideMetrics.UnambiguousDecoys,
            [PeptideBacterialUnambiguousTargetScores] = peptideMetrics.TargetScores.ToArray(),
            [PeptideBacterialUnambiguousDecoyScores] = peptideMetrics.DecoyScores.ToArray(),
            [PeptideBacterialTargetDeltaScores] = peptideMetrics.DeltaScores.ToArray(),
            [PeptidePrecursorMassErrors] = peptideMetrics.PrecursorMassErrors.ToArray(),
            [PeptidePrecursorEnvelopePeakCount] = peptideMetrics.PrecursorPeakCounts.Select(v => (double)v).ToArray(),
            [PeptidePrecusorDeconScores] = peptideMetrics.PrecursorDeconScores.ToArray(),
            [PeptidePrecusorFractionalIntensity] = peptideMetrics.PrecursorFractionalIntensity.ToArray()
        };
    }

    private (int Targets, int Decoys, int UnambiguousTargets, int UnambiguousDecoys,
        List<double> TargetScores, List<double> DecoyScores, List<double> DeltaScores, List<double> PrecursorMassErrors, List<int> PrecursorPeakCounts, List<double> PrecursorDeconScores, List<double> PrecursorFractionalIntensity)
        AnalyzeSpectralMatches(List<SpectralMatch> spectralMatches, double qValueThreshold, bool isPeptideLevel = false)
    {
        int targets = 0, decoys = 0, unambiguousTargets = 0, unambiguousDecoys = 0;
        List<double> targetScores = [];
        List<double> decoyScores = [];
        List<double> deltaScores = [];
        List<double> precursorMassErrors = [];
        List<int> precursorPeakCounts = [];
        List<double> precursorDeconScores = [];
        List<double> precursorFractionalIntensity = [];

        foreach (var psm in spectralMatches)
        {
            double qValue = isPeptideLevel
                ? psm.PeptideFdrInfo?.QValue ?? 1.0
                : psm.GetFdrInfo(false)?.QValue ?? 1.0;

            if (qValue > qValueThreshold) continue;

            bool isOrganismAmbiguous = false;
            foreach (var match in psm.BestMatchingBioPolymersWithSetMods)
            {
                string? organism = match.SpecificBioPolymer.Parent.Organism;
                if (organism != null && organism.Contains(_targetOrganism, StringComparison.Ordinal))
                {
                    isOrganismAmbiguous = true;
                    break;
                }
            }

            if (psm.IsDecoy)
            {
                decoys++;
                if (!isOrganismAmbiguous)
                {
                    unambiguousDecoys++;
                    decoyScores.Add(psm.Score);
                }
            }
            else
            {
                targets++;
                if (!isOrganismAmbiguous)
                {
                    unambiguousTargets++;
                    targetScores.Add(psm.Score);
                }
                deltaScores.Add(psm.DeltaScore);
                precursorMassErrors.Add(psm.PrecursorMassErrorPpm.Mean());
                precursorPeakCounts.Add(psm.PrecursorScanEnvelopePeakCount);
                precursorDeconScores.Add(psm.PrecursorScanDeconvolutionScore);
                precursorFractionalIntensity.Add(psm.PrecursorFractionalIntensity);
            }
        }

        return (targets, decoys, unambiguousTargets, unambiguousDecoys, targetScores, decoyScores, deltaScores, precursorMassErrors, precursorPeakCounts, precursorDeconScores, precursorFractionalIntensity);
    }
}
