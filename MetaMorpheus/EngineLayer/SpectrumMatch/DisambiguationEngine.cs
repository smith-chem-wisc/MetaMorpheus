using EngineLayer.FdrAnalysis;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.SpectrumMatch;

/// <summary>
/// Engine designed to disambiguate spectral matches.
/// <remarks>
/// This is currently done in several locations and should be consolidated to this engine.
/// Places in which disambiguation occurs:
/// SearchTask -> Internal Ions
/// ProteinParsimonyEngine -> remove non-parsimonious peptides.
/// </remarks>
/// </summary>
public class DisambiguationEngine : MetaMorpheusEngine
{
    private readonly double _qvalueNotchDisambiguationThreshold = 0.05;
    private readonly List<SpectralMatch> _allSpectralMatches;
    private readonly IPepDisambiguationRule _pepRule;

    /// <param name="pepRule">Removes hypotheses by their PEP. Null removes none by PEP.</param>
    public DisambiguationEngine(List<SpectralMatch> allPsms, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds,
        IPepDisambiguationRule pepRule = null)
        : base(commonParameters, fileSpecificParameters, nestedIds)
    {
        _allSpectralMatches = allPsms;
        _pepRule = pepRule;
    }

    protected override MetaMorpheusEngineResults RunSpecific()
    {
        Status("Running Disambiguation Engine...");

        // By PEP first, so the notch q-values judged next are computed without the hypotheses it removed.
        int? removedByPep = null;
        if (_pepRule != null)
        {
            removedByPep = DisambiguateByPep();
            if (removedByPep > 0)
            {
                // Each removal already resolved its own match's ambiguities. Recalculate Q-Values.
                FdrAnalysisEngine.DoFalseDiscoveryRateAnalysis(_allSpectralMatches, false, FileSpecificParameters, null, null, CommonParameters);
            }
        }

        // Remove ambiguous PSMs by various methods.
        int removedQValueNotch = DisambiguateByQValueNotch();

        if (removedQValueNotch > 0)
        {
            // Resolve all remaining ambiguities
            foreach (var psm in _allSpectralMatches)
                psm.ResolveAllAmbiguities();

            // Recalculate Q-Values
            FdrAnalysisEngine.DoFalseDiscoveryRateAnalysis(_allSpectralMatches, false, FileSpecificParameters, null, null, CommonParameters);
        }

        Status("Done.");
        return new DisambiguationEngineResults(this)
        {
            RemovedByPEP = removedByPep,
            RemovedByQValueNotch = removedQValueNotch,
        };
    }

    /// <summary>
    /// Removes the hypotheses <see cref="_pepRule"/> selects from each ambiguous match. A match is skipped unless
    /// every hypothesis has a PEP, i.e. unless the PEP engine scored it.
    /// </summary>
    private int DisambiguateByPep()
    {
        int removed = 0;
        foreach (var psm in _allSpectralMatches.Where(p => p != null && p.BestMatchingBioPolymersWithSetMods.Count() > 1))
        {
            var hypotheses = psm.BestMatchingBioPolymersWithSetMods.ToList();
            if (hypotheses.Any(h => !h.PEP.HasValue))
                continue; // can't disambiguate by PEP if we don't have PEPs for all of them.

            // Materialized before the first removal, which changes the match's hypotheses.
            foreach (var remove in _pepRule.HypothesesToRemove(hypotheses).ToList())
            {
                psm.RemoveThisAmbiguousPeptide(remove);
                removed++;
            }
        }
        return removed;
    }

    private int DisambiguateByQValueNotch()
    {
        int removed = 0;
        foreach (var psm in _allSpectralMatches.Where(p => p.Notch == null && p.BestMatchingBioPolymersWithSetMods.Count() > 1))
        {
            if (psm.BestMatchingBioPolymersWithSetMods.Any(b => !b.QValueNotch.HasValue))
                continue; // can't disambiguate by q-value if we don't have q-values for all of them.
            var bestQValue = psm.BestMatchingBioPolymersWithSetMods.Min(b => b.QValueNotch!.Value);
            var toRemove = psm.BestMatchingBioPolymersWithSetMods.Where(b => Math.Abs(b.QValueNotch!.Value - bestQValue) > _qvalueNotchDisambiguationThreshold);

            foreach (var remove in toRemove)
            {
                psm.RemoveThisAmbiguousPeptide(remove);
                removed++;
            }
        }
        return removed;
    }
}

public class DisambiguationEngineResults : MetaMorpheusEngineResults
{
    /// <summary>
    /// Null when the engine was given no PEP rule.
    /// </summary>
    public int? RemovedByPEP { get; set; }
    public int RemovedByQValueNotch { get; set; }
    //public int RemovedByInternalIonCount { get; set; }

    public DisambiguationEngineResults(DisambiguationEngine s) : base(s)
    {
    }

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine(base.ToString());
        if (RemovedByPEP.HasValue)
            sb.AppendLine($"Ambiguous {GlobalVariables.AnalyteType.GetUniqueFormLabel()}s removed by PEP: {RemovedByPEP}");
        sb.AppendLine($"Ambiguous {GlobalVariables.AnalyteType.GetUniqueFormLabel()}s removed QValueNotch: {RemovedByQValueNotch}");
        //sb.AppendLine($"Ambiguous {GlobalVariables.AnalyteType.GetUniqueFormLabel()}s removed Internal Ion Count: {RemovedByInternalIonCount}");
        return sb.ToString();
    }
}
