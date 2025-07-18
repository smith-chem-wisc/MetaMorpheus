using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using MassSpectrometry;
using Omics;
using Omics.Modifications;
using OxyPlot;
using Readers;

namespace GuiFunctions.MetaDraw;

/// <summary>
/// View model for a single chimeric spectral match
/// </summary>
public class ChimericSpectralMatchModel : BaseViewModel
{
    public SpectrumMatchFromTsv Psm { get; set; }
    public IsotopicEnvelope PrecursorEnvelope { get; set; }
    public OxyColor Color { get; set; }
    public OxyColor ProteinColor { get; set; }
    public string Letter { get; set; }
    public string ModString { get; set; }
    public Dictionary<int, Modification> AllModsOneIsNterminus { get; set; }

    public ChimericSpectralMatchModel(SpectrumMatchFromTsv psm, IsotopicEnvelope precursorEnvelope, OxyColor color, OxyColor proteinColor)
    {
        Psm = psm;
        PrecursorEnvelope = precursorEnvelope;
        Color = color;
        ProteinColor = proteinColor;

        var modDict = GlobalVariables.AnalyteType switch {
            AnalyteType.Oligo => GlobalVariables.AllRnaModsKnownDictionary,
            _ => GlobalVariables.AllModsKnownDictionary
        };

        // If ambiguous, process all sequences and concatenate their mod strings with '|'
        var sequences = psm.FullSequence.Split('|');
        var modStrings = new List<string>();
        foreach (var seq in sequences)
        {
            var mods = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(seq, modDict);
            modStrings.Add(string.Join(", ", mods.Select(m => $"{m.Key} - {m.Value.IdWithMotif}")));
        }
        AllModsOneIsNterminus = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(sequences[0], modDict);
        ModString = string.Join(" | ", modStrings);
    }
}