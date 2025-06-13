using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using MassSpectrometry;
using Readers;
using EngineLayer;

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

    public ChimeraAnalysisTabViewModel(List<SpectrumMatchFromTsv> allPsms, Dictionary<string, MsDataFile> dataFiles, string? exportDirectory = null)
    {
        ChimeraLegendViewModel = new ChimeraLegendViewModel();
        ChimeraGroupViewModels = ConstructChimericPsms(allPsms, dataFiles)
            .OrderByDescending(p => p.Count)
            .ToList();
        ExportDirectory = exportDirectory ?? Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
        SelectedExportType = "Png";
        ExportTypes = [ "Pdf", "Png", "Svg"];
        OnPropertyChanged(nameof(ExportTypes));
    }

    private static List<ChimeraGroupViewModel> ConstructChimericPsms(List<SpectrumMatchFromTsv> psms, Dictionary<string, MsDataFile> dataFiles)
    {
        // Groups are made from psms that pass the MetaDraw quality filter and only include decoys if we are showing decoys.
        var filteredPsms = new List<SpectrumMatchFromTsv>(psms.Count);
        foreach (var p in psms)
        {
            if (p.QValue <= MetaDrawSettings.QValueFilter && (MetaDrawSettings.ShowDecoys || !p.DecoyContamTarget.Contains('D')))
                filteredPsms.Add(p);
        }

        var groupDict = new Dictionary<(string, int), IList<SpectrumMatchFromTsv>>(filteredPsms.Count);
        foreach (var p in filteredPsms)
        {
            var key = (p.FileNameWithoutExtension, p.Ms2ScanNumber);
            groupDict.AddOrCreate(key, p);
        }

        List<ChimeraGroupViewModel> toReturn = new(groupDict.Count);
        foreach (var group in groupDict.Values)
        {
            if (group.Count <= 1)
                continue;

            var first = group[0];
            if (!dataFiles.TryGetValue(first.FileNameWithoutExtension, out MsDataFile spectraFile))
                continue;

            var ms1Scan = spectraFile.GetOneBasedScanFromDynamicConnection(first.PrecursorScanNum);
            var ms2Scan = spectraFile.GetOneBasedScanFromDynamicConnection(first.Ms2ScanNumber);

            if (ms1Scan == null || ms2Scan == null)
                continue;

            var orderedGroup = group.OrderBy(p => p.PrecursorMz);
            var groupVm = new ChimeraGroupViewModel(orderedGroup, ms1Scan, ms2Scan);
            if (groupVm.ChimericPsms.Count > 0)
                toReturn.Add(groupVm);
        }
        return toReturn;
    }

    #region IO

    private string _exportDirectory;
    public string ExportDirectory
    {
        get
        {
            if (!Directory.Exists(_exportDirectory))
                Directory.CreateDirectory(_exportDirectory);
            return _exportDirectory;
        }
        set
        {
            _exportDirectory = value;
            OnPropertyChanged(nameof(ExportDirectory));
        }
    }
    public ObservableCollection<string> ExportTypes { get; set; }

    private string _selectedExportType;
    public string SelectedExportType
    {
        get => _selectedExportType;
        set
        {
            _selectedExportType = value;
            OnPropertyChanged(nameof(SelectedExportType));
        }
    }

    #endregion
}