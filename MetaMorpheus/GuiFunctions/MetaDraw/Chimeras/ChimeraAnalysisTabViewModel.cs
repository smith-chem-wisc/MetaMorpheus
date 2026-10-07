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
using Easy.Common.Extensions;
using Point = System.Windows.Point;
using System.Windows.Data;
using System.Threading;

namespace GuiFunctions.MetaDraw;

/// <summary>
/// All data and user triggered manipulations for the Chimera Analysis tab
/// </summary>
public class ChimeraAnalysisTabViewModel : MetaDrawTabViewModel
{
    #region Displayed in GUI
    public override string TabHeader { get; init; } = "Chimera Annotation";
    public ChimeraSpectrumMatchPlot ChimeraSpectrumMatchPlot { get; set; }
    public Ms1ChimeraPlot Ms1ChimeraPlot { get; set; }
    public ObservableCollection<ChimeraGroupViewModel> ChimeraGroupViewModels { get; set; }
    public Dictionary<string, MsDataFile> MsDataFiles { get; private set; }

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
    public ChimeraAnalysisTabViewModel(string exportDirectory = null)
    {
        ChimeraGroupViewModels = new();
        ExportDirectory = exportDirectory ?? Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
        ExportMs1Command = new RelayCommand(ExportMs1);
        ExportMs2Command = new RelayCommand(ExportMs2);
        ExportSequenceCoverageCommand = new RelayCommand(ExportSequenceCoverage);
        ExportLegendCommand = new DelegateCommand(ExportLegend);
        ExportAllCommand = new DelegateCommand(ExportAll);

        BindingOperations.EnableCollectionSynchronization(ChimeraGroupViewModels, ThreadLocker);
    }

    public void ProcessChimeraData(List<SpectrumMatchFromTsv> allPsms, Dictionary<string, MsDataFile> dataFiles, CancellationToken cancellationToken = default)
    {
        MsDataFiles = dataFiles;
        ChimeraGroupViewModels.Clear();

        foreach (var chimeraGroup in ConstructChimericPsms(allPsms, dataFiles, cancellationToken)
                     .OrderByDescending(p => p.Count)
                     .ThenByDescending(p => p.UniqueFragments)
                     .ThenByDescending(p => p.TotalFragments))
        {
            ChimeraGroupViewModels.Add(chimeraGroup);
        }
    }

    public static List<ChimeraGroupViewModel> ConstructChimericPsms(List<SpectrumMatchFromTsv> psms, Dictionary<string, MsDataFile> dataFiles, CancellationToken cancellationToken = default)
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
            if (cancellationToken.IsCancellationRequested)
                throw new OperationCanceledException();

            if (group.Count <= 1)
                continue;

            var first = group[0];
            if (!dataFiles.TryGetValue(first.FileNameWithoutExtension, out MsDataFile spectraFile))
                continue;

            var ms1Scan = spectraFile.GetOneBasedScanFromDynamicConnection(first.PrecursorScanNum);
            var ms2Scan = spectraFile.GetOneBasedScanFromDynamicConnection(first.Ms2ScanNumber);

            if (ms1Scan == null || ms2Scan == null)
                continue;

            // Attempt to construct ChimeraGroupViewModel, but catch any exceptions that may occur
            // (e.g., due to deconvolution failures or mismatched isotopic envelopes. This could be caused by different decon params in the visualization and the search).
            // This prevents the program from crashing and skips problematic groups.
            try
            {
                var orderedGroup = group.OrderByDescending(p => p.Score); 
                var groupVm = new ChimeraGroupViewModel(orderedGroup, ms1Scan, ms2Scan);
                if (groupVm.ChimericPsms.Count > 0)
                    toReturn.Add(groupVm);
            }
            catch (Exception)
            {
                // Log or handle the exception as needed (e.g., for debugging)
                // For now, just skip this group and continue processing others.
            }
        }
        return toReturn;
    }

    #region IO

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

    private void ExportAll(object frameworkElement) {
        if (SelectedChimeraGroup == null)
        {
            MessageBoxHelper.Warn("No chimera group selected for export.");
            return;
        }

        string path = Path.Combine(ExportDirectory,
            $"{SelectedChimeraGroup.FileNameWithoutExtension}_{SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_All.{MetaDrawSettings.ExportType.ToLower()}");

        var dpiScale = MetaDrawSettings.CanvasPdfExportDpi / 96.0;

        // Render individual components
        var (ms1Bitmap, _) = RenderFrameworkElementToBitmap(Ms1ChimeraPlot.PlotView);
        var (ms2Bitmap, _) = RenderFrameworkElementToBitmap(ChimeraSpectrumMatchPlot.PlotView);
        var (sequenceBitmap, _) = RenderFrameworkElementToBitmap(ChimeraDrawnSequence.SequenceDrawingCanvas);

        var bitmaps = new List<System.Drawing.Bitmap> { ms1Bitmap, ms2Bitmap, sequenceBitmap };

        // Stack vertically
        var stackedImage = MetaDrawLogic.CombineBitmap(bitmaps, new List<Point>
        {
            new(0, 0),
            new(0, ms1Bitmap.Height),
            new(0, ms1Bitmap.Height + ms2Bitmap.Height)
        }, overlap: false);

        // Legend overlay
        if (MetaDrawSettings.DisplayChimeraLegend && LegendCanvas != null)
        {
            var (legendBitmap, legendVisualOffset) = RenderFrameworkElementToBitmap(LegendCanvas);
            var commonAncestor = FindCommonAncestor(Ms1ChimeraPlot.PlotView, LegendCanvas) as Visual;

            if (commonAncestor != null)
            {
                Point ms1Position = Ms1ChimeraPlot.PlotView.TransformToAncestor(commonAncestor).Transform(new Point(0, 0));
                Point legendPosition = LegendCanvas.TransformToAncestor(commonAncestor).Transform(new Point(0, 0));

                Vector offsetInLayout = legendPosition - ms1Position;

                Point legendOffsetInStackedImage = new(
                    (offsetInLayout.X - legendVisualOffset.X + 20) * dpiScale,
                    (offsetInLayout.Y - legendVisualOffset.Y + 20) * dpiScale); // already in ms1 space

                // Overlay the legend on the stacked image
                stackedImage = MetaDrawLogic.CombineBitmap(
                    new List<System.Drawing.Bitmap> { stackedImage, legendBitmap },
                    new List<Point> { new(0, 0), legendOffsetInStackedImage },
                    overlap: true);
            }
        }

        // Export
        SpectrumMatchPlot.ExportPlot(path, stackedImage, stackedImage.Width, stackedImage.Height);
        MessageBoxHelper.Show("Exported all to: " + path);
    }

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

        int pixelWidth, pixelHeight;
        if (double.IsNegativeInfinity(bounds.Width) || double.IsNaN(bounds.Width) || bounds.Width <= 0)
            pixelWidth = (int)Math.Ceiling(700 * dpiScale);
        else
            pixelWidth = (int)Math.Ceiling(bounds.Width * dpiScale);

        if (double.IsNegativeInfinity(bounds.Height) || double.IsNaN(bounds.Height) || bounds.Height <= 0)
            pixelHeight = (int)Math.Ceiling(400 * dpiScale);
        else
            pixelHeight = (int)Math.Ceiling(bounds.Height * dpiScale);

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