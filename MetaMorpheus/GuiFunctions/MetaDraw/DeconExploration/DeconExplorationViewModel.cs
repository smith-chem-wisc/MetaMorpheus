using MassSpectrometry;
using MzLibUtil;
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

public class DeconExplorationViewModel : BaseViewModel
{
    public ObservableCollection<MsDataFile> MsDataFiles { get; set; } = new();
    public ObservableCollection<MsDataScan> Scans { get; set; } = new();
    public ObservableCollection<DeconvolutedSpeciesViewModel> DeconvolutedSpecies { get; set; } = new();
    public ICommand RunDeconvolutionCommand { get; }

    public DeconExplorationViewModel()
    {
        Mode = DeconvolutionMode.FullSpectrum;
        RunDeconvolutionCommand = new RelayCommand(RunDeconvolution);
    }

    private MsDataFile _selectedMsDataFile;
    public MsDataFile SelectedMsDataFile
    {
        get => _selectedMsDataFile;
        set
        {
            if (_selectedMsDataFile == value) return;
            _selectedMsDataFile = value;
            Scans.Clear();
            if (_selectedMsDataFile != null)
            {
                // Load scans from the selected file
                foreach (var scan in _selectedMsDataFile.GetAllScansList())
                {
                    Scans.Add(scan);
                }
            }
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
            OnPropertyChanged(nameof(Mode));
        }
    }

    public List<DeconvolutionMode> DeconvolutionModes { get; } = System.Enum.GetValues<DeconvolutionMode>().ToList();
    public DeconHostViewModel DeconHostViewModel { get; set; } = new();

    private void RunDeconvolution()
    {
        DeconvolutedSpecies.Clear();
        if (SelectedMsDataScan == null) return;

        IEnumerable<DeconvolutedSpeciesViewModel> results;

        if (Mode == DeconvolutionMode.FullSpectrum)
        {
            results = DeconvoluteFullSpectrum(SelectedMsDataScan);
        }
        else
        {
            var precursorScan = Scans.First(p => p.OneBasedScanNumber == SelectedMsDataScan.OneBasedPrecursorScanNumber);
            results = DeconvoluteIsolationRegion(SelectedMsDataScan, precursorScan);
        }

        foreach (var species in results)
            DeconvolutedSpecies.Add(species);
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
}

public class DeconvolutedSpeciesViewModel(IsotopicEnvelope envelope) : BaseViewModel
{
    private string? _annotation;
    private double? _totalIntensity;
    private double? _mostAbundantMz;

    public IsotopicEnvelope Envelope { get; } = envelope;
    public double MonoisotopicMass => Envelope.MonoisotopicMass;
    public int Charge => Envelope.Charge;
    public int PeakCount => Envelope.Peaks.Count;
    public double Intensity => _totalIntensity ??= Envelope.Peaks.Sum(p => p.intensity);
    public string Annotation => _annotation ??= $"M={MonoisotopicMass}\nz={Charge}";
    public double MostAbundantMz => _mostAbundantMz ?? Envelope.Peaks.MaxBy(p => p.intensity).mz;
}