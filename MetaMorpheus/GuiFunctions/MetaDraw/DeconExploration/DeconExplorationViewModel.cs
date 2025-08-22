using MassSpectrometry;
using MathNet.Numerics;
using MzLibUtil;
using OxyPlot.Wpf;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows.Input;
using OxyPlot;

namespace GuiFunctions;

public enum DeconvolutionMode
{
    FullSpectrum,
    IsolationRegion
}

public class DeconExplorationViewModel : BaseViewModel
{
    public ObservableCollection<MsDataFile> MsDataFiles { get; set; } = new();
    public ObservableCollection<MsDataScan> Scans { get; set; } = new();
    public ObservableCollection<DeconvolutedSpeciesViewModel> DeconvolutedSpecies { get; set; } = new();
    public DeconvolutionPlot Plot { get; set; }
    public ICommand RunDeconvolutionCommand { get; }

    public DeconExplorationViewModel()
    {
        Mode = DeconvolutionMode.FullSpectrum;
        RunDeconvolutionCommand = new DelegateCommand((pv) => RunDeconvolution((PlotView)pv));
    }

    private MsDataFile _selectedMsDataFile;
    public MsDataFile SelectedMsDataFile
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

    private MsDataScan _selectedMsDataScan;
    public MsDataScan SelectedMsDataScan
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

    public List<DeconvolutionMode> DeconvolutionModes { get; } = System.Enum.GetValues<DeconvolutionMode>().ToList();
    public DeconHostViewModel DeconHostViewModel { get; set; } = new();

    private void RunDeconvolution(PlotView plotView)
    {
        DeconvolutedSpecies.Clear();
        if (SelectedMsDataScan == null) return;

        IEnumerable<DeconvolutedSpeciesViewModel> results;
        MzRange? range = null;
        MsDataScan scan;

        if (Mode == DeconvolutionMode.FullSpectrum)
        {
            scan = SelectedMsDataScan;
            results = DeconvoluteFullSpectrum(scan);
        }
        else
        {
            range = SelectedMsDataScan.IsolationRange;
            var precursorScan = SelectedMsDataFile.GetOneBasedScan(SelectedMsDataScan.OneBasedPrecursorScanNumber!.Value);
            scan = precursorScan;
            results = DeconvoluteIsolationRegion(SelectedMsDataScan, precursorScan);
        }

        foreach (var species in results.OrderByDescending(p => p.MonoisotopicMass.Round(2)).ThenByDescending(p => p.Charge))
            DeconvolutedSpecies.Add(species);

        Plot = new DeconvolutionPlot(plotView, scan, DeconvolutedSpecies.ToList(), range);
    }

    private IEnumerable<DeconvolutedSpeciesViewModel> DeconvoluteFullSpectrum(MsDataScan scan)
    {
        IEnumerable<IsotopicEnvelope> deconResults;
        if (scan.MsnOrder == 1)
            deconResults = Deconvoluter.Deconvolute(scan, DeconHostViewModel.PrecursorDeconvolutionParameters.Parameters);
        else
            deconResults = Deconvoluter.Deconvolute(scan, DeconHostViewModel.ProductDeconvolutionParameters.Parameters);

        foreach (var envelope in deconResults)
        {
            yield return new DeconvolutedSpeciesViewModel(envelope);
        }
    }

    private IEnumerable<DeconvolutedSpeciesViewModel> DeconvoluteIsolationRegion(MsDataScan scan, MsDataScan precursorScan)
    {
        return scan.GetIsolatedMassesAndCharges(precursorScan, DeconHostViewModel.PrecursorDeconvolutionParameters.Parameters).Select(envelope => new DeconvolutedSpeciesViewModel(envelope));
    }

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
                foreach (var scan in _selectedMsDataFile.GetMsDataScans())
                    Scans.Add(scan);
                break;
            }
        }
    }
}

public class DeconvolutedSpeciesViewModel(IsotopicEnvelope envelope) : BaseViewModel
{
    private string? _annotation;
    private double? _totalIntensity;
    private double? _mostAbundantMz;
    private string? _peakMzs;

    public IsotopicEnvelope Envelope { get; } = envelope;
    public OxyColor Color { get; set; } = OxyColors.Automatic;
    public double MonoisotopicMass => Envelope.MonoisotopicMass;
    public int Charge => Envelope.Charge;
    public int PeakCount => Envelope.Peaks.Count;
    public double Intensity => _totalIntensity ??= Envelope.Peaks.Sum(p => p.intensity);
    public string Annotation => _annotation ??= $"M={MonoisotopicMass.Round(2)}\nz={Charge}";
    public double MostAbundantMz => _mostAbundantMz ?? Envelope.Peaks.MaxBy(p => p.intensity).mz;
    public string PeakMzs => _peakMzs ?? string.Join(',', Envelope.Peaks.Select(p => p.mz.Round(2)).OrderBy(p => p));
}