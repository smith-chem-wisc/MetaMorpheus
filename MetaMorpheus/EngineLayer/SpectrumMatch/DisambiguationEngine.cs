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
/// PEPAnalysisEngine -> By PEP
/// ProteinParsimonyEngine -> remove non-parsimonious peptides. 
/// </remarks>
/// </summary>
public class DisambiguationEngine : MetaMorpheusEngine
{
    private readonly double _qvalueNotchDisambiguationThreshold = 0.05;
    private readonly List<SpectralMatch> _allSpectralMatches;

    public DisambiguationEngine(List<SpectralMatch> allPsms, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds)
        : base(commonParameters, fileSpecificParameters, nestedIds)
    {
        _allSpectralMatches = allPsms;
    }

    protected override MetaMorpheusEngineResults RunSpecific()
    {
        Status("Running Disambiguation Engine...");

        // Remove ambiguous PSMs by various methods.
        int removedQValueNotch = DisambiguateByQValueNotch();

        if (removedQValueNotch > 0)
        {        
            // Resolve all remaining ambiguities
            foreach (var psm in _allSpectralMatches)
                psm.ResolveAllAmbiguities();

            // Recalculate Q-Values 
            FdrAnalysisEngine.DoFalseDiscoveryRateAnalysis(_allSpectralMatches, false, FileSpecificParameters, null, null);
        }

        Status("Done.");
        return new DisambiguationEngineResults(this)
        {
            RemovedByQValueNotch = removedQValueNotch,
        };
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
    //public int RemovedByPEP { get; set; }
    public int RemovedByQValueNotch { get; set; }
    //public int RemovedByInternalIonCount { get; set; }

    public DisambiguationEngineResults(DisambiguationEngine s) : base(s)
    {
    }

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine(base.ToString());
        //sb.AppendLine($"Ambiguous {GlobalVariables.AnalyteType.GetUniqueFormLabel()}s removed by PEP: {RemovedByPEP}");
        sb.AppendLine($"Ambiguous {GlobalVariables.AnalyteType.GetUniqueFormLabel()}s removed QValueNotch: {RemovedByQValueNotch}");
        //sb.AppendLine($"Ambiguous {GlobalVariables.AnalyteType.GetUniqueFormLabel()}s removed Internal Ion Count: {RemovedByInternalIonCount}");
        return sb.ToString();
    }
}