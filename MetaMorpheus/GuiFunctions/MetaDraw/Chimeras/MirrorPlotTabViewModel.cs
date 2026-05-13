using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using GuiFunctions;
using MassSpectrometry;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers;

namespace GuiFunctions.MetaDraw;

public class MirrorPlotTabViewModel : MetaDrawTabViewModel
{
    public override string TabHeader { get; init; } = "Mirror Plot";

    public ObservableCollection<SpectrumMatchFromTsv> AllFilteredPsms { get; set; }

    private SpectrumMatchFromTsv _selectedLeftPsm;
    public SpectrumMatchFromTsv SelectedLeftPsm
    {
        get => _selectedLeftPsm;
        set
        {
            _selectedLeftPsm = value;
            OnPropertyChanged(nameof(SelectedLeftPsm));
            RebuildMirrorPlot();
        }
    }

    private SpectrumMatchFromTsv _selectedRightPsm;
    public SpectrumMatchFromTsv SelectedRightPsm
    {
        get => _selectedRightPsm;
        set
        {
            _selectedRightPsm = value;
            OnPropertyChanged(nameof(SelectedRightPsm));
            RebuildMirrorPlot();
        }
    }

    private PeptideSpectrumMatchPlot _mirrorPlot;
    public PeptideSpectrumMatchPlot MirrorPlot
    {
        get => _mirrorPlot;
        set
        {
            _mirrorPlot = value;
            OnPropertyChanged(nameof(MirrorPlot));
        }
    }

    public Dictionary<string, MsDataFile> MsDataFiles { get; private set; }

    public ICommand ExportCommand { get; set; }

    public MirrorPlotTabViewModel(string exportDirectory = null)
    {
        AllFilteredPsms = new ObservableCollection<SpectrumMatchFromTsv>();
        ExportDirectory = exportDirectory ?? Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
        ExportCommand = new DelegateCommand(obj => ExportPlot());
        BindingOperations.EnableCollectionSynchronization(AllFilteredPsms, ThreadLocker);
    }

    public void ProcessMirrorData(List<SpectrumMatchFromTsv> allPsms, Dictionary<string, MsDataFile> dataFiles)
    {
        MsDataFiles = dataFiles;
        AllFilteredPsms.Clear();

        foreach (var psm in allPsms)
        {
            AllFilteredPsms.Add(psm);
        }
    }

    public void RebuildMirrorPlot()
    {
        if (SelectedLeftPsm == null || SelectedRightPsm == null)
            return;

        if (!MsDataFiles.TryGetValue(SelectedLeftPsm.FileNameWithoutExtension, out var dataFileA))
            return;
        if (!MsDataFiles.TryGetValue(SelectedRightPsm.FileNameWithoutExtension, out var dataFileB))
            return;

        var scanA = dataFileA.GetOneBasedScan(SelectedLeftPsm.Ms2ScanNumber);
        var scanB = dataFileB.GetOneBasedScan(SelectedRightPsm.Ms2ScanNumber);

        if (scanA == null || scanB == null)
            return;

        var plotViewRef = MirrorPlotView;
        if (plotViewRef == null)
            return;

        MirrorPlot = new PeptideSpectrumMatchPlot(
            plotViewRef,
            SelectedLeftPsm,
            scanA,
            SelectedLeftPsm.MatchedIons,
            annotateProperties: false,
            mirrorIons: SelectedRightPsm.MatchedIons);

        MirrorPlot.AnnotatePropertiesBottom(SelectedRightPsm);
    }

    public void ExportPlot()
    {
        if (MirrorPlot == null || TopCanvasExport == null || BottomCanvasExport == null)
            return;

        string path = Path.Combine(ExportDirectory,
            $"MirrorPlot_{SelectedLeftPsm?.Ms2ScanNumber}_{SelectedRightPsm?.Ms2ScanNumber}.{MetaDrawSettings.ExportType.ToLower()}");

        double width = 700;
        double height = 370;
        double dpiScale = MetaDrawSettings.CanvasPdfExportDpi / 96.0;

        string tempModelPath = Path.Combine(Path.GetDirectoryName(path), Path.GetRandomFileName() + "." + MetaDrawSettings.ExportType);
        string tempTopSeqPath = Path.Combine(Path.GetDirectoryName(path), Path.GetRandomFileName() + "_topSeq.png");
        string tempBottomSeqPath = Path.Combine(Path.GetDirectoryName(path), Path.GetRandomFileName() + "_bottomSeq.png");

        List<System.Drawing.Bitmap> bitmaps = new();
        List<Point> points = new();

        MirrorPlot.ExportToPng(tempModelPath, (int)width, (int)height);
        bitmaps.Add(new System.Drawing.Bitmap(tempModelPath));
        points.Add(new Point(0, 0));

        RenderCanvasToFile(TopCanvasExport, tempTopSeqPath, dpiScale);
        System.Drawing.Bitmap topBitmap = new(tempTopSeqPath);
        bitmaps.Add(topBitmap);
        var topOffset = (Vector)TopCanvasExport.GetType()
            .GetProperty("VisualOffset", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            .GetValue(TopCanvasExport);
        points.Add(new Point(topOffset.X, topOffset.Y));

        RenderCanvasToFile(BottomCanvasExport, tempBottomSeqPath, dpiScale);
        System.Drawing.Bitmap bottomBitmap = new(tempBottomSeqPath);
        bitmaps.Add(bottomBitmap);
        var bottomOffset = (Vector)BottomCanvasExport.GetType()
            .GetProperty("VisualOffset", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
            .GetValue(BottomCanvasExport);
        points.Add(new Point(bottomOffset.X, bottomOffset.Y));

        System.Drawing.Bitmap combined = MetaDrawLogic.CombineBitmap(bitmaps, points);

        topBitmap.Dispose();
        bottomBitmap.Dispose();
        File.Delete(tempModelPath);
        File.Delete(tempTopSeqPath);
        File.Delete(tempBottomSeqPath);

        SpectrumMatchPlot.ExportPlot(path, combined, width, height);
    }

    private static void RenderCanvasToFile(Canvas canvas, string path, double dpiScale)
    {
        canvas.Height += 30;
        canvas.Width += 30;
        Size size = new((int)canvas.Width, (int)canvas.Height);
        canvas.Measure(size);
        canvas.Arrange(new Rect(size));

        RenderTargetBitmap rtb = new(
            (int)(dpiScale * canvas.Width),
            (int)(dpiScale * canvas.Height),
            MetaDrawSettings.CanvasPdfExportDpi,
            MetaDrawSettings.CanvasPdfExportDpi,
            PixelFormats.Pbgra32);
        rtb.Render(canvas);

        PngBitmapEncoder encoder = new();
        encoder.Frames.Add(BitmapFrame.Create(rtb));
        using FileStream file = File.Create(path);
        encoder.Save(file);
    }

    public OxyPlot.Wpf.PlotView MirrorPlotView { get; set; }
    public Canvas TopCanvasExport { get; set; }
    public Canvas BottomCanvasExport { get; set; }
}
