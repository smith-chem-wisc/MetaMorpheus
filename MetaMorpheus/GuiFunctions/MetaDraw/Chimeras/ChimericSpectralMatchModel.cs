using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using MassSpectrometry;
using Omics;
using Omics.Modifications;
using OxyPlot;
using Readers;

namespace GuiFunctions.MetaDraw.Chimeras;

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

        // If ambiguous, takes the first sequence. 
        AllModsOneIsNterminus = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(psm.FullSequence.Split('|')[0], modDict);
        ModString = string.Join(", ", AllModsOneIsNterminus.Select(m => $"{m.Key} - {m.Value.IdWithMotif}"));
    }
}