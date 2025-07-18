using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using MassSpectrometry;
using Readers;
using EngineLayer;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media.Imaging;
using System.Windows.Media;
using Easy.Common.Extensions;
using Point = System.Windows.Point;
using System.Reflection.Emit;
using OxyPlot.Wpf;
using OxyPlot;

namespace GuiFunctions.MetaDraw;

/// <summary>
/// All data and user triggered manipulations for the Chimera Analysis tab
/// </summary>
public class ChimeraAnalysisTabViewModel : BaseViewModel
{
    #region Displayed in GUI
    public ChimeraSpectrumMatchPlot ChimeraSpectrumMatchPlot { get; set; }
    public Ms1ChimeraPlot Ms1ChimeraPlot { get; set; }
    public ObservableCollection<ChimeraGroupViewModel> ChimeraGroupViewModels { get; set; }

    private ChimeraGroupViewModel _selectedChimeraGroup;
    public ChimeraGroupViewModel SelectedChimeraGroup
    {
        get => _selectedChimeraGroup;
        set
        {
            _selectedChimeraGroup = value;
            if (value != null && MetaDrawSettings.DisplayChimeraLegend)
                LegendCanvas = new ChimeraLegendCanvas(value);
            OnPropertyChanged(nameof(SelectedChimeraGroup));
        }
    }

    private ChimeraLegendCanvas _legendCanvas;
    public ChimeraLegendCanvas LegendCanvas
    {
        get => _legendCanvas;
        set
        {
            _legendCanvas = value;
            OnPropertyChanged(nameof(LegendCanvas));
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

    public ChimeraAnalysisTabViewModel(List<SpectrumMatchFromTsv> allPsms, Dictionary<string, MsDataFile> dataFiles, string exportDirectory = null)
    {
        ChimeraGroupViewModels = [..ConstructChimericPsms(allPsms, dataFiles)
            .OrderByDescending(p => p.Count)
            .ThenByDescending(p => p.UniqueFragments)
            .ThenByDescending(p => p.TotalFragments)];
        ExportDirectory = exportDirectory ?? Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

        ExportMs1Command = new RelayCommand(ExportMs1);
        ExportMs2Command = new RelayCommand(ExportMs2);
        ExportSequenceCoverageCommand = new RelayCommand(ExportSequenceCoverage);
        ExportLegendCommand = new DelegateCommand(ExportLegend);
        ExportAllCommand = new DelegateCommand(ExportAll);
    }

    public static List<ChimeraGroupViewModel> ConstructChimericPsms(List<SpectrumMatchFromTsv> psms, Dictionary<string, MsDataFile> dataFiles)
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

            var orderedGroup = group.OrderByDescending(p => p.Score);
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

    public ICommand ExportMs1Command { get; set; }
    public ICommand ExportMs2Command { get; set; }
    public ICommand ExportSequenceCoverageCommand { get; set; }
    public ICommand ExportLegendCommand { get; set; }
    public ICommand ExportAllCommand { get; set; }
    private void ExportSequenceCoverage()
    {
        string path = Path.Combine(ExportDirectory,
            $"{SelectedChimeraGroup.FileNameWithoutExtension}_{SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_SequenceCoverage.{MetaDrawSettings.ExportType.ToLower()}");
        // change path to .png
        path = Path.ChangeExtension(path, "png");

        // convert canvas to bitmap
        Rect bounds = new Rect(new Point(0, 0), new Point(ChimeraDrawnSequence.SequenceDrawingCanvas.Width, ChimeraDrawnSequence.SequenceDrawingCanvas.Height));
        double dpi = 96d;

        var width = double.IsNegativeInfinity(bounds.Width)
            ? MetaDrawSettings.AnnotatedSequenceTextSpacing * ChimeraDrawnSequence.ChimeraGroupViewModel.ChimericPsms.Max(p => p.Psm.BaseSeq.Length)
            : bounds.Width;
        var height = double.IsNegativeInfinity(bounds.Height)
            ? 80 * ChimeraDrawnSequence.ChimeraGroupViewModel.Count
            : bounds.Height;

        RenderTargetBitmap rtb = new(
            (int)width, //width
            (int)height, //height
            dpi, //dpi x
            dpi, //dpi y
            PixelFormats.Default // pixelformat
        );
        var size = new Size(width, height);

        DrawingVisual dv = new();
        using (DrawingContext dc = dv.RenderOpen())
        {
            VisualBrush vb = new(ChimeraDrawnSequence.SequenceDrawingCanvas);
            dc.DrawRectangle(vb, null, new Rect(new Point(0, 0), size));
        }

        rtb.Render(dv);

        // export
        using (FileStream stream = new(path, FileMode.Create))
        {
            PngBitmapEncoder encoder = new();
            encoder.Frames.Add(BitmapFrame.Create(rtb));
            encoder.Save(stream);
        }
        MessageBoxHelper.Show(MetaDrawSettings.ExportType + "(s) exported to: " + path);
    }

    private void ExportLegend(object frameworkElement)
    {
        var element = frameworkElement as FrameworkElement;
        if (element == null)
        {
            MessageBoxHelper.Warn("No legend available for export.");
            return;
        }

        string path = Path.Combine(ExportDirectory,
            $"{SelectedChimeraGroup.FileNameWithoutExtension}_{SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_Legend.{MetaDrawSettings.ExportType.ToLower()}");


        var bounds = VisualTreeHelper.GetDescendantBounds(element);
        double dpi = 96d;

        RenderTargetBitmap rtb = new RenderTargetBitmap((int)bounds.Width,
            (int)bounds.Height,
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
        MessageBoxHelper.Show(MetaDrawSettings.ExportType + "(s) exported to: " + path);
    }

    private void ExportMs1()
    {
        if (SelectedChimeraGroup == null)
        {
            MessageBoxHelper.Warn("No chimera group selected for export.");
            return;
        }

        string path = Path.Combine(ExportDirectory,
            $"{SelectedChimeraGroup.FileNameWithoutExtension}_{SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_MS1.{MetaDrawSettings.ExportType.ToLower()}");

        var bitmap = CombinePlotAndLegend(Ms1ChimeraPlot.PlotView,
            MetaDrawSettings.DisplayChimeraLegend ? LegendCanvas : null);

        int width = Ms1ChimeraPlot.Model.Width > 0 ? (int)Ms1ChimeraPlot.Model.Width : 700;
        int height = Ms1ChimeraPlot.Model.Height > 0 ? (int)Ms1ChimeraPlot.Model.Height : 300;

        SpectrumMatchPlot.ExportPlot(path, bitmap, width, height);
        MessageBoxHelper.Show(MetaDrawSettings.ExportType + "(s) exported to: " + path);
    }

    private void ExportMs2()
    {
        if (SelectedChimeraGroup == null)
        {
            MessageBoxHelper.Warn("No chimera group selected for export.");
            return;
        }

        string path = Path.Combine(ExportDirectory,
            $"{SelectedChimeraGroup.FileNameWithoutExtension}_{SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_MS2.{MetaDrawSettings.ExportType.ToLower()}");

        var bitmap = CombinePlotAndLegend(ChimeraSpectrumMatchPlot.PlotView,
            MetaDrawSettings.DisplayChimeraLegend ? LegendCanvas : null);

        int width = ChimeraSpectrumMatchPlot.Model.Width > 0 ? (int)ChimeraSpectrumMatchPlot.Model.Width : 700;
        int height = ChimeraSpectrumMatchPlot.Model.Height > 0 ? (int)ChimeraSpectrumMatchPlot.Model.Height : 300;

        SpectrumMatchPlot.ExportPlot(path, bitmap, width, height);
        MessageBoxHelper.Show(MetaDrawSettings.ExportType + "(s) exported to: " + path);
    }

    private void ExportAll(object frameworkElement) { }

    private System.Drawing.Bitmap CombinePlotAndLegend(FrameworkElement plotView, FrameworkElement legend = null)
    {
        if (plotView == null)
            throw new ArgumentNullException(nameof(plotView));

        // DPI scale
        var dpiScale = MetaDrawSettings.CanvasPdfExportDpi / 96.0;

        // Render the base plot at its actual size
        var (plotBitmap, _) = RenderFrameworkElementToBitmap(plotView);
        if (plotBitmap == null)
            throw new InvalidOperationException("Failed to render plot view.");

        // No legend? Just return the plot bitmap
        if (legend == null)
            return plotBitmap;

        // Render legend + capture visual offset
        var (legendBitmap, legendVisualOffset) = RenderFrameworkElementToBitmap(legend);
        if (legendBitmap == null)
            return plotBitmap;

        try
        {
            // Find a common ancestor to compute relative layout offset
            var commonAncestor = FindCommonAncestor(plotView, legend) as Visual;
            if (commonAncestor == null)
                return plotBitmap;

            // Get layout positions in screen space relative to the ancestor
            Point plotPosition = plotView.TransformToAncestor(commonAncestor).Transform(new Point(0, 0));
            Point legendPosition = legend.TransformToAncestor(commonAncestor).Transform(new Point(0, 0));

            // Offset of legend relative to plot (assumed same coordinate space)
            Vector layoutOffset = legendPosition - plotPosition;

            // Adjust for the fact that the visual starts rendering at legendVisualOffset
            Point legendOffsetInBitmap = new(
                (layoutOffset.X - legendVisualOffset.X + 20) * dpiScale,
                (layoutOffset.Y - legendVisualOffset.Y + 20) * dpiScale);

            return MetaDrawLogic.CombineBitmap(
                new List<System.Drawing.Bitmap> { plotBitmap, legendBitmap },
                new List<Point> { new(0, 0), legendOffsetInBitmap },
                overlap: true);
        }
        catch
        {
            return plotBitmap;
        }
    }

    private static (System.Drawing.Bitmap bitmap, Point visualOffset) RenderFrameworkElementToBitmap(FrameworkElement element)
    {
        if (element == null) return (null, new Point(0, 0));

        var dpi = MetaDrawSettings.CanvasPdfExportDpi;
        var dpiScale = dpi / 96.0;

        // Get the true visual bounds of the element and all children
        Rect bounds = VisualTreeHelper.GetDescendantBounds(element);

        int pixelWidth = (int)Math.Ceiling(bounds.Width * dpiScale);
        int pixelHeight = (int)Math.Ceiling(bounds.Height * dpiScale);

        if (pixelWidth <= 0 || pixelHeight <= 0)
            return (null, new Point(0, 0));

        var rtb = new RenderTargetBitmap(pixelWidth, pixelHeight, dpi, dpi, PixelFormats.Pbgra32);

        var dv = new DrawingVisual();
        using (DrawingContext dc = dv.RenderOpen())
        {
            var vb = new VisualBrush(element);
            dc.DrawRectangle(vb, null, new Rect(new Point(), bounds.Size));
        }

        rtb.Render(dv);

        var encoder = new PngBitmapEncoder();
        encoder.Frames.Add(BitmapFrame.Create(rtb));

        using var ms = new MemoryStream();
        encoder.Save(ms);

        using var bmp = new System.Drawing.Bitmap(ms);
        return (new System.Drawing.Bitmap(bmp), bounds.TopLeft); // deep copy of map 
    }

    private static DependencyObject FindCommonAncestor(DependencyObject a, DependencyObject b)
    {
        var ancestorsA = new HashSet<DependencyObject>();
        for (var current = a; current != null; current = VisualTreeHelper.GetParent(current))
            ancestorsA.Add(current);

        for (var current = b; current != null; current = VisualTreeHelper.GetParent(current))
            if (ancestorsA.Contains(current))
                return current;

        return null;
    }


    #endregion
}