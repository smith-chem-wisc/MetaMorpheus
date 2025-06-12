using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using OxyPlot;
using Readers;

namespace GuiFunctions;

/// <summary>
/// All data and user triggered manipulations for the Chimera Analysis tab
/// </summary>
public class ChimeraAnalysisTabViewModel : BaseViewModel
{
    public ChimeraAnalysisTabViewModel(List<SpectrumMatchFromTsv> allPsms, Dictionary<string, MsDataFile> dataFiles)
    {

    }
}

/// <summary>
/// View model for a group of chimeric IDs from a single MS2 scan
/// </summary>
public class ChimeraGroupViewModel : BaseViewModel
{
    public string FileNameWithoutExtension { get; set; }
    public int OneBasedPrecursorScanNumber { get; set; }
    public int Ms2ScanNumber { get; set; }
    public int Count => ChimericPsms.Count;
    public int ProteinCount { get; }
    public ObservableCollection<ChimericSpectralMatchModel> ChimericPsms { get; set; }
}

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
        var mods = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(psm.FullSequence, modDict);
        ModString = string.Join(", ", mods.Select(m => $"{m.Key} - {m.Value.IdWithMotif}"));
    }
}
