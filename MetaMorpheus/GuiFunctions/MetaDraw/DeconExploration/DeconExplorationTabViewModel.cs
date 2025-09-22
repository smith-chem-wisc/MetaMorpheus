#nullable enable
using EngineLayer;
using GuiFunctions.MetaDraw;
using MassSpectrometry;
using MathNet.Numerics;
using MzLibUtil;
using OxyPlot.Wpf;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows.Data;
using System.Windows.Input;

namespace GuiFunctions;

public enum DeconvolutionMode
{
    FullSpectrum,
    IsolationRegion
}

public class DeconExplorationTabViewModel : MetaDrawTabViewModel
{
    private readonly MetaDrawLogic _metaDrawLogic;
    public override string TabHeader { get; init; } = "Deconvolution Exploration";
    public ObservableCollection<MsDataFile> MsDataFiles { get; set; } = new();
    public ObservableCollection<MsDataScan> Scans { get; set; } = new();
    public ObservableCollection<DeconvolutedSpeciesViewModel> DeconvolutedSpecies { get; set; } = new();
    public DeconvolutionPlot? Plot { get; set; }
    public List<DeconvolutionMode> DeconvolutionModes { get; } = System.Enum.GetValues<DeconvolutionMode>().ToList();
    public DeconHostViewModel DeconHostViewModel { get; set; } = new();

    private MsDataFile? _selectedMsDataFile;
    public MsDataFile? SelectedMsDataFile
    {
        get => _selectedMsDataFile;
        set
        {
            if (_selectedMsDataFile == value) return;
            _selectedMsDataFile = value;

            PopulateScansCollection();
            OnPropertyChanged(nameof(SelectedMsDataFile));
        }
    }

    private MsDataScan? _selectedMsDataScan;
    public MsDataScan? SelectedMsDataScan
    {
        get => _selectedMsDataScan;
        set
        {
            if (_selectedMsDataScan == value) return;
            _selectedMsDataScan = value;
            OnPropertyChanged(nameof(SelectedMsDataScan));
        }
    }

    private DeconvolutionMode _mode;
    public DeconvolutionMode Mode
    {
        get => _mode;
        set
        {
            if (_mode == value) return;
            _mode = value;

            PopulateScansCollection();
            OnPropertyChanged(nameof(Mode));
        }
    }

    private bool _onlyIdentifiedScans;
    public bool OnlyIdentifiedScans
    {
        get => _onlyIdentifiedScans;
        set
        {
            if (_onlyIdentifiedScans == value) return;
            _onlyIdentifiedScans = value;

            PopulateScansCollection();
            OnPropertyChanged(nameof(OnlyIdentifiedScans));
        }
    }

    public ICommand RunDeconvolutionCommand { get; }

    public DeconExplorationTabViewModel(MetaDrawLogic metaDrawLogic)
    {
        Mode = DeconvolutionMode.IsolationRegion;
        RunDeconvolutionCommand = new DelegateCommand(pv => RunDeconvolution((pv as PlotView)!));

        BindingOperations.EnableCollectionSynchronization(MsDataFiles, ThreadLocker);
        _metaDrawLogic = metaDrawLogic;
    }

    private void RunDeconvolution(PlotView plotView)
    {
        DeconvolutedSpecies.Clear();
        if (SelectedMsDataScan == null) return;

        IEnumerable<DeconvolutedSpeciesViewModel> results;
        MzRange? isolationRange;
        MsDataScan scanToPlot;

        if (Mode == DeconvolutionMode.FullSpectrum)
        {
            isolationRange = null;
            scanToPlot = SelectedMsDataScan;
            results = DeconvoluteFullSpectrum(scanToPlot);
        }
        else
        {
            isolationRange = SelectedMsDataScan.IsolationRange;
            scanToPlot = SelectedMsDataFile!.GetOneBasedScan(SelectedMsDataScan.OneBasedPrecursorScanNumber!.Value);
            results = DeconvoluteIsolationRegion(SelectedMsDataScan, scanToPlot);
        }

        var sortedSpecies = results
            .OrderByDescending(p => p.MonoisotopicMass.Round(2))
            .ThenByDescending(p => p.Charge)
            .ToList();

        foreach (var species in sortedSpecies)
            DeconvolutedSpecies.Add(species);

        Plot = new DeconvolutionPlot(plotView, scanToPlot, sortedSpecies, isolationRange);
    }

    private IEnumerable<DeconvolutedSpeciesViewModel> DeconvoluteFullSpectrum(MsDataScan scan) 
        => Deconvoluter.Deconvolute(scan, scan.MsnOrder == 1 ? DeconHostViewModel.PrecursorDeconvolutionParameters.Parameters : DeconHostViewModel.ProductDeconvolutionParameters.Parameters)
        .Select(envelope => new DeconvolutedSpeciesViewModel(envelope));

    private IEnumerable<DeconvolutedSpeciesViewModel> DeconvoluteIsolationRegion(MsDataScan scan, MsDataScan precursorScan) 
        => scan.GetIsolatedMassesAndCharges(precursorScan, DeconHostViewModel.PrecursorDeconvolutionParameters.Parameters)
        .Select(envelope => new DeconvolutedSpeciesViewModel(envelope));

    private void PopulateScansCollection()
    {
        Scans.Clear();

        if (SelectedMsDataFile == null)
            return;

        IEnumerable<MsDataScan> scansToDisplay = _mode switch
        {
            DeconvolutionMode.IsolationRegion => SelectedMsDataFile.GetMsDataScans().Where(scan => scan.MsnOrder == 2),
            DeconvolutionMode.FullSpectrum => SelectedMsDataFile.GetMsDataScans(),
            _ => []
        };

        if (OnlyIdentifiedScans)
        {
            var key = Path.GetFileName(SelectedMsDataFile.FilePath.Replace(GlobalVariables.GetFileExtension(SelectedMsDataFile.FilePath), string.Empty));
            if (_metaDrawLogic.SpectralMatchesGroupedByFile.TryGetValue(key, out var matches))
            {
                HashSet<int> scanNumbers;
                if (_mode == DeconvolutionMode.FullSpectrum)
                {
                    scanNumbers = matches.Select(m => m.Ms2ScanNumber)
                        .Concat(matches.Select(m => m.PrecursorScanNum))
                        .ToHashSet();
                }
                else
                {
                    scanNumbers = matches.Select(m => m.Ms2ScanNumber)
                        .ToHashSet();
                }
                scansToDisplay = scansToDisplay.Where(scan => scanNumbers.Contains(scan.OneBasedScanNumber));
            }
        }

        foreach (var scan in scansToDisplay)
            Scans.Add(scan);
    }
}