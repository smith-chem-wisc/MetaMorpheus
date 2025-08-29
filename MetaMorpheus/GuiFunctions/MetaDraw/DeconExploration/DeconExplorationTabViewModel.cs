#nullable enable
using MassSpectrometry;
using MathNet.Numerics;
using MzLibUtil;
using OxyPlot.Wpf;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows.Input;

namespace GuiFunctions;

public enum DeconvolutionMode
{
    FullSpectrum,
    IsolationRegion
}

public class DeconExplorationTabViewModel : BaseViewModel
{
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

    public ICommand RunDeconvolutionCommand { get; }

    public DeconExplorationTabViewModel()
    {
        Mode = DeconvolutionMode.IsolationRegion;
        RunDeconvolutionCommand = new DelegateCommand(pv => RunDeconvolution((pv as PlotView)!));
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

        switch (Mode)
        {
            // Display only MS2
            case DeconvolutionMode.IsolationRegion:
            {
                foreach (var scan in SelectedMsDataFile.GetMsDataScans())
                    if (scan.MsnOrder == 2)
                        Scans.Add(scan);
                break;
            }
            case DeconvolutionMode.FullSpectrum:
            {
                foreach (var scan in SelectedMsDataFile.GetMsDataScans())
                    Scans.Add(scan);
                break;
            }
        }
    }
}