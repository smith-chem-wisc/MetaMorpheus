using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using MassSpectrometry;
using Readers;
using EngineLayer;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media.Imaging;
using System.Windows.Media;
using System.Drawing;
using System.Drawing.Imaging;
using Brushes = System.Windows.Media.Brushes;
using Point = System.Windows.Point;

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

        ExportMs1Command = new RelayCommand(ExportMs1);
        ExportMs2Command = new RelayCommand(ExportMs2);
        ExportSequenceCoverageCommand = new RelayCommand(ExportSequenceCoverage);
        ExportLegendCommand = new DelegateCommand(ExportLegend);
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

    public ICommand ExportMs1Command { get; set; }
    public ICommand ExportMs2Command { get; set; }
    public ICommand ExportSequenceCoverageCommand { get; set; }
    public ICommand ExportLegendCommand { get; set; }

    private void ExportMs1()
    {
        if (SelectedChimeraGroup == null)
        {
            MessageBox.Show("No chimera group selected for export.");
            return;
        }

        string path = Path.Combine(ExportDirectory,
            $"{SelectedChimeraGroup.FileNameWithoutExtension}_{SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_MS1.{SelectedExportType.ToLower()}");

        switch (SelectedExportType)
        {
            case "Pdf":
                Ms1ChimeraPlot.ExportToPdf(path, (int)Ms1ChimeraPlot.Model.Width, (int)Ms1ChimeraPlot.Model.Height);
                break;
            case "Png":
                Ms1ChimeraPlot.ExportToPng(path, (int)Ms1ChimeraPlot.Model.Width, (int)Ms1ChimeraPlot.Model.Height);
                break;
            case "Svg":
                Ms1ChimeraPlot.ExportToSvg(path, (int)Ms1ChimeraPlot.Model.Width, (int)Ms1ChimeraPlot.Model.Height);
                break;
            default:
                throw new ArgumentOutOfRangeException();
        }

        MessageBox.Show(MetaDrawSettings.ExportType + "(s) exported to: " + path);
    }

    private void ExportMs2()
    {
        if (SelectedChimeraGroup == null)
        {
            MessageBox.Show("No chimera group selected for export.");
            return;
        }

        string path = Path.Combine(ExportDirectory,
            $"{SelectedChimeraGroup.FileNameWithoutExtension}_{SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_MS2.{SelectedExportType.ToLower()}");

        switch (SelectedExportType)
        {
            case "Pdf":
                ChimeraSpectrumMatchPlot.ExportToPdf(path, (int)ChimeraSpectrumMatchPlot.Model.Width, (int)ChimeraSpectrumMatchPlot.Model.Height);
                break;
            case "Png":
                ChimeraSpectrumMatchPlot.ExportToPng(path, (int)ChimeraSpectrumMatchPlot.Model.Width, (int)ChimeraSpectrumMatchPlot.Model.Height);
                break;
            case "Svg":
                ChimeraSpectrumMatchPlot.ExportToSvg(path, (int)ChimeraSpectrumMatchPlot.Model.Width, (int)ChimeraSpectrumMatchPlot.Model.Height);
                break;
            default:
                throw new ArgumentOutOfRangeException();
        }
        MessageBox.Show(MetaDrawSettings.ExportType + "(s) exported to: " + path);
    }

    private void ExportSequenceCoverage()
    {
        string path = Path.Combine(ExportDirectory,
            $"{SelectedChimeraGroup.FileNameWithoutExtension}_{SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_SequenceCoverage.{SelectedExportType.ToLower()}");
        // change path to .png
        path = Path.ChangeExtension(path, "png");

        // convert canvas to bitmap
        Rect bounds = VisualTreeHelper.GetDescendantBounds(ChimeraDrawnSequence.SequenceDrawingCanvas);
        double dpi = 96d;

        // Defaults in case the canvas has not been expanded. 
        if (double.IsNegativeInfinity(bounds.Width))
        {
            bounds.Width = 800;
            bounds.Height = 80 * ChimeraDrawnSequence.ChimeraGroupViewModel.Count;
        }

        RenderTargetBitmap rtb = new(
            (int)bounds.Width, //width
            (int)bounds.Height, //height
            dpi, //dpi x
            dpi, //dpi y
            System.Windows.Media.PixelFormats.Default // pixelformat
        );

        DrawingVisual dv = new();
        using (DrawingContext dc = dv.RenderOpen())
        {
            VisualBrush vb = new(ChimeraDrawnSequence.SequenceDrawingCanvas);
            dc.DrawRectangle(vb, null, new Rect(new System.Windows.Point(), bounds.Size));
        }

        rtb.Render(dv);

        // export
        using (FileStream stream = new(path, FileMode.Create))
        {
            PngBitmapEncoder encoder = new();
            encoder.Frames.Add(BitmapFrame.Create(rtb));
            encoder.Save(stream);
        }
        MessageBox.Show(MetaDrawSettings.ExportType + "(s) exported to: " + path);
    }

    private void ExportLegend(object frameworkElement)
    {
        var element = frameworkElement as FrameworkElement;
        if (element == null)
        {
            MessageBox.Show("No legend available for export.");
            return;
        }

        string path = Path.Combine(ExportDirectory,
            $"{SelectedChimeraGroup.FileNameWithoutExtension}_{SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_Legend.{SelectedExportType.ToLower()}");


        var bounds = VisualTreeHelper.GetDescendantBounds(element);
        double dpi = 96d;

        RenderTargetBitmap rtb = new RenderTargetBitmap((int)(bounds.Width),
            (int)(bounds.Height),
            dpi,
            dpi,
            PixelFormats.Pbgra32);

        DrawingVisual dv = new DrawingVisual();
        using (DrawingContext ctx = dv.RenderOpen())
        {
            VisualBrush vb = new VisualBrush(element);
            ctx.DrawRectangle(vb, null, new Rect(new Point(), bounds.Size));
        }

        rtb.Render(dv);

        // export
        using (FileStream stream = new(path, FileMode.Create))
        {
            PngBitmapEncoder encoder = new();
            encoder.Frames.Add(BitmapFrame.Create(rtb));
            encoder.Save(stream);
        }
        MessageBox.Show(MetaDrawSettings.ExportType + "(s) exported to: " + path);
    }

    #endregion
}