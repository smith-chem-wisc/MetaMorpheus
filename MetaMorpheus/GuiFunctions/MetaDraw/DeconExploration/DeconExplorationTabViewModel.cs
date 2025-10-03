#nullable enable
using EngineLayer;
using GuiFunctions.MetaDraw;
using MassSpectrometry;
using MathNet.Numerics;
using MzLibUtil;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Data;
using System.Windows.Input;
using ThermoFisher.CommonCore.Data.Business;

namespace GuiFunctions;

public enum DeconvolutionMode
{
    FullSpectrum,
    IsolationRegion
}

public class DeconExplorationTabViewModel : MetaDrawTabViewModel
{
    private readonly MetaDrawLogic _metaDrawLogic;
    private readonly SemaphoreSlim _populateScansSemaphore = new(1, 1);

    public override string TabHeader { get; init; } = "Deconvolution Exploration";
    public ObservableCollection<MsDataFile> MsDataFiles { get; set; } = new();
    public ObservableCollection<MsDataScan> Scans { get; set; } = new();
    public ObservableCollection<DeconvolutedSpeciesViewModel> DeconvolutedSpecies { get; set; } = new();
    public DeconvolutionPlot? Plot { get; set; }
    public List<DeconvolutionMode> DeconvolutionModes { get; } = System.Enum.GetValues<DeconvolutionMode>().ToList();

    private DeconHostViewModel _deconHostViewModel = new();
    public DeconHostViewModel DeconHostViewModel
    {
        get => _deconHostViewModel;
        set
        {
            _deconHostViewModel = value;
            OnPropertyChanged(nameof(DeconHostViewModel));
        }
    }

    private MsDataFile? _selectedMsDataFile;
    public MsDataFile? SelectedMsDataFile
    {
        get => _selectedMsDataFile;
        set
        {
            if (_selectedMsDataFile == value) return;
            _selectedMsDataFile = value;

            _ = PopulateScansCollection();
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

    private DeconvolutionMode _mode = DeconvolutionMode.IsolationRegion;
    public DeconvolutionMode Mode
    {
        get => _mode;
        set
        {
            if (_mode == value) return;
            _mode = value;

            _ = PopulateScansCollection();
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

            _ = PopulateScansCollection();
            OnPropertyChanged(nameof(OnlyIdentifiedScans));
        }
    }

    private double minMzToPlot;
    public double MinMzToPlot
    {
        get => minMzToPlot;
        set
        {
            minMzToPlot = value;
            OnPropertyChanged(nameof(MinMzToPlot));
        }
    }
    private double maxMzToPlot;
    public double MaxMzToPlot
    {
        get => maxMzToPlot;
        set
        {
            maxMzToPlot = value;
            OnPropertyChanged(nameof(MaxMzToPlot));
        }
    }

    public ICommand RunDeconvolutionCommand { get; }

    public DeconExplorationTabViewModel(MetaDrawLogic? metaDrawLogic)
    {
        RunDeconvolutionCommand = new DelegateCommand(pv => RunDeconvolution((pv as PlotView)!));

        BindingOperations.EnableCollectionSynchronization(MsDataFiles, ThreadLocker);
        BindingOperations.EnableCollectionSynchronization(Scans, ThreadLocker);
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
            var min = MinMzToPlot > 0 ? MinMzToPlot : SelectedMsDataScan.MassSpectrum.Range.Minimum;
            var max = MaxMzToPlot > 0 ? MaxMzToPlot : SelectedMsDataScan.MassSpectrum.Range.Maximum;
            isolationRange = new MzRange(min, max);

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
            .Where(p => isolationRange is null || p.Envelope.Peaks.Any(peak => isolationRange.Contains(peak.mz)))
            .OrderByDescending(p => p.MonoisotopicMass.Round(2))
            .ThenByDescending(p => p.Charge)
            .ToList();

        foreach (var species in sortedSpecies)
            DeconvolutedSpecies.Add(species);

        Plot = new DeconvolutionPlot(plotView, scanToPlot, sortedSpecies, Mode, isolationRange);
    }

    private IEnumerable<DeconvolutedSpeciesViewModel> DeconvoluteFullSpectrum(MsDataScan scan) 
        => Deconvoluter.Deconvolute(scan, scan.MsnOrder == 1 ? DeconHostViewModel.PrecursorDeconvolutionParameters.Parameters : DeconHostViewModel.ProductDeconvolutionParameters.Parameters)
        .Select(envelope => new DeconvolutedSpeciesViewModel(envelope));

    private IEnumerable<DeconvolutedSpeciesViewModel> DeconvoluteIsolationRegion(MsDataScan scan, MsDataScan precursorScan) 
        => scan.GetIsolatedMassesAndCharges(precursorScan, DeconHostViewModel.PrecursorDeconvolutionParameters.Parameters)
        .Select(envelope => new DeconvolutedSpeciesViewModel(envelope));

    private async Task PopulateScansCollection()
    {
        await _populateScansSemaphore.WaitAsync();
        try
        {
            Scans.Clear();

            if (SelectedMsDataFile == null)
                return;

            IsLoading = true;
            try
            {
                var selectedMsDataFile = SelectedMsDataFile;
                var mode = _mode;
                var onlyIdentifiedScans = OnlyIdentifiedScans;
                var metaDrawLogic = _metaDrawLogic;

                // Run scan filtering on a background thread
                var scansToDisplay = await Task.Run(() =>
                {
                    IEnumerable<MsDataScan> scans = mode switch
                    {
                        DeconvolutionMode.IsolationRegion => selectedMsDataFile.GetMsDataScans().Where(scan => scan.MsnOrder == 2),
                        DeconvolutionMode.FullSpectrum => selectedMsDataFile.GetMsDataScans(),
                        _ => Enumerable.Empty<MsDataScan>()
                    };

                    if (onlyIdentifiedScans)
                    {
                        var key = Path.GetFileName(selectedMsDataFile.FilePath.Replace(GlobalVariables.GetFileExtension(selectedMsDataFile.FilePath), string.Empty));
                        if (metaDrawLogic.SpectralMatchesGroupedByFile.TryGetValue(key, out var matches))
                        {
                            HashSet<int> scanNumbers;
                            if (mode == DeconvolutionMode.FullSpectrum)
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
                            scans = scans.Where(scan => scanNumbers.Contains(scan.OneBasedScanNumber));
                        }
                    }

                    return scans;
                });

                foreach (var scan in scansToDisplay)
                    Scans.Add(scan);
            }
            catch (Exception e) { MessageBoxHelper.Warn($"An error occurred while populating the scan collection:\n{e}"); }
            finally
            {
                IsLoading = false;
            }
        }
        finally
        {
            _populateScansSemaphore.Release();
        }
    }
}


[ExcludeFromCodeCoverage] // for design time display
public class DeconTabModel() : DeconExplorationTabViewModel(new MetaDrawLogic())
{
    public static DeconTabModel Instance => new();
}