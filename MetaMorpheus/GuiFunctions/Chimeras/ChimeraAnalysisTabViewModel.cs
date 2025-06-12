using System.Collections.Generic;
using System.Linq;
using System.Windows;
using Easy.Common.Extensions;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Modifications;
using OxyPlot;
using Readers;

namespace GuiFunctions;

/// <summary>
/// All data and user triggered manipulations for the Chimera Analysis tab
/// </summary>
public class ChimeraAnalysisTabViewModel : BaseViewModel
{
    #region Displayed in GUI
    public ChimeraSpectrumMatchPlot ChimeraSpectrumMatchPlot { get; set; }
    public Ms1ChimeraPlot Ms1ChimeraPlot { get; set; }
    public List<ChimeraGroupViewModel> ChimeraGroupViewModels { get; set; }

    private ChimeraGroupViewModel _selectedChimeraGroup;
    public ChimeraGroupViewModel SelectedChimeraGroup
    {
        get => _selectedChimeraGroup;
        set
        {
            _selectedChimeraGroup = value;
            ChimeraLegendViewModel.ChimeraLegendItems = value.LegendItems;
            OnPropertyChanged(nameof(SelectedChimeraGroup));
        }
    }

    private ChimeraLegendViewModel _chimeraLegendViewModel;
    public ChimeraLegendViewModel ChimeraLegendViewModel
    {
        get => _chimeraLegendViewModel;
        set
        {
            _chimeraLegendViewModel = value;
            OnPropertyChanged(nameof(ChimeraLegendViewModel));
        }
    }

    private ChimeraDrawnSequence _chimeraDrawnSequence;
    public ChimeraDrawnSequence ChimeraDrawnSequence
    {
        get => _chimeraDrawnSequence;
        set
        {
            _chimeraDrawnSequence = value;
            OnPropertyChanged(nameof(ChimeraDrawnSequence));
        }
    }

    #endregion

    #region Settings that change behavior

    private bool _groupProteinsInSequenceAnnotation;
    public bool GroupProteinsInSequenceAnnotation
    {
        get => _groupProteinsInSequenceAnnotation;
        set
        {
            _groupProteinsInSequenceAnnotation = value;
            OnPropertyChanged(nameof(GroupProteinsInSequenceAnnotation));
        }
    }

    private bool useLetterOnly;
    public bool UseLetterOnly
    {
        get => useLetterOnly;
        set
        {
            if (useLetterOnly == value)
                return;

            useLetterOnly = value;
            ChimeraGroupViewModels.ForEach(p => p.AssignIonColors(useLetterOnly));
            OnPropertyChanged(nameof(UseLetterOnly));
        }
    }

    #endregion

    public ChimeraAnalysisTabViewModel(List<SpectrumMatchFromTsv> allPsms, Dictionary<string, MsDataFile> dataFiles)
    {
        ChimeraLegendViewModel = new ChimeraLegendViewModel();
        ChimeraGroupViewModels = ConstructChimericPsms(allPsms, dataFiles)
            .OrderByDescending(p => p.Count)
            .ToList();
    }

    private static IEnumerable<ChimeraGroupViewModel> ConstructChimericPsms(List<SpectrumMatchFromTsv> psms, Dictionary<string, MsDataFile> dataFiles)
    {
        return psms
            .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T")
            .GroupBy(p => (p.FileNameWithoutExtension, p.Ms2ScanNumber))
            .Where(p => p.Count() > 1)
            .Select(group =>
            {
                if (!dataFiles.TryGetValue(group.First().FileNameWithoutExtension, out MsDataFile spectraFile))
                    return null;

                var ms1Scan = spectraFile.GetOneBasedScanFromDynamicConnection(group.First().PrecursorScanNum);
                var ms2Scan = spectraFile.GetOneBasedScanFromDynamicConnection(group.First().Ms2ScanNumber);

                if (ms1Scan == null || ms2Scan == null)
                    return null;

                var groupVm = new ChimeraGroupViewModel(group.OrderBy(p => p.PrecursorMz), ms1Scan, ms2Scan);
                return groupVm.ChimericPsms.Count > 0 ? groupVm : null;
            })
            .Where(groupVm => groupVm != null);
    }
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
        AllModsOneIsNterminus = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(psm.FullSequence, modDict);
        ModString = string.Join(", ", AllModsOneIsNterminus.Select(m => $"{m.Key} - {m.Value.IdWithMotif}"));
    }
}


internal static class Extensions
{
    /// <summary>
    /// Determines if a majority of values are within a range
    /// </summary>
    /// <param name="range"></param>
    /// <param name="values"></param>
    /// <returns></returns>
    internal static bool MajorityWithin(this MzRange range, IEnumerable<double> values)
    {
        int within = values.Count(p => p >= range.Minimum && p <= range.Maximum);
        return within > values.Count() / 2;
    }

    // method to add a value to a list in a dictionary if the key is present, and craete a new list if the key is not present
    public static void AddOrReplace<TKey, TValue, TValue2>(this Dictionary<TKey, List<(TValue, TValue2)>> dictionary, TKey key,
        TValue value, TValue2 value2)
    {
        if (dictionary.ContainsKey(key))
        {
            var previousVersion = dictionary[key].FirstOrDefault(p => p.Item1.Equals(value));
            if (!previousVersion.GetType().IsDefault())
            {
                dictionary[key].Remove(previousVersion);
            }
            dictionary[key].Add((value, value2));

        }
        else
            dictionary.Add(key, new List<(TValue, TValue2)> { (value, value2) });

    }
}