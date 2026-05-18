using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Input;
using System.Windows;
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

    private MirrorSpectrumMatchPlot _mirrorPlot;
    public MirrorSpectrumMatchPlot MirrorPlot
    {
        get => _mirrorPlot;
        set
        {
            _mirrorPlot = value;
            OnPropertyChanged(nameof(MirrorPlot));
        }
    }

    public Dictionary<string, MsDataFile> MsDataFiles { get; private set; }

    private bool _useRelativeIntensity;
    public bool UseRelativeIntensity
    {
        get => _useRelativeIntensity;
        set
        {
            _useRelativeIntensity = value;
            OnPropertyChanged(nameof(UseRelativeIntensity));
            RefreshPlot();
        }
    }

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

        var ionsA = SelectedLeftPsm.MatchedIons;
        var ionsB = SelectedRightPsm.MatchedIons;

        if (UseRelativeIntensity)
        {
            var yArrayA = scanA.MassSpectrum.YArray;
            var yArrayB = scanB.MassSpectrum.YArray;
            double maxA = yArrayA.Max();
            double maxB = yArrayB.Max();
            if (maxA > 0 && maxB > 0)
            {
                scanA = NormalizeScan(scanA, maxA);
                scanB = NormalizeScan(scanB, maxB);
                ionsA = NormalizeIons(ionsA, maxA);
                ionsB = NormalizeIons(ionsB, maxB);
            }
        }

        var plotViewRef = MirrorPlotView;
        if (plotViewRef == null)
            return;

        MirrorPlot = new MirrorSpectrumMatchPlot(
            plotViewRef,
            SelectedLeftPsm, scanA, ionsA,
            SelectedRightPsm, scanB, ionsB,
            annotateProperties: true);
    }

    public void RefreshPlot()
    {
        RebuildMirrorPlot();
    }

    public void ExportPlot()
    {
        if (MirrorPlot == null || MirrorPlotExportElement == null)
            return;

        string directoryPath = Path.Combine(ExportDirectory, "MetaDrawExport",
            DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));
        Directory.CreateDirectory(directoryPath);

        string path = Path.Combine(directoryPath,
            $"MirrorPlot_Scan{SelectedLeftPsm?.Ms2ScanNumber}_vs_Scan{SelectedRightPsm?.Ms2ScanNumber}.{MetaDrawSettings.ExportType.ToLower()}");

        using var bitmap = RenderFrameworkElementToBitmap(MirrorPlotExportElement);
        SpectrumMatchPlot.ExportPlot(path, bitmap, bitmap.Width, bitmap.Height);
        MessageBoxHelper.Show(MetaDrawSettings.ExportType + "(s) exported to: " + path);
    }

    public OxyPlot.Wpf.PlotView MirrorPlotView { get; set; }
    public FrameworkElement MirrorPlotExportElement { get; set; }
    public Canvas TopCanvasExport { get; set; }
    public Canvas BottomCanvasExport { get; set; }

    private static MsDataScan NormalizeScan(MsDataScan original, double maxIntensity)
    {
        var xArray = original.MassSpectrum.XArray;
        var yArray = original.MassSpectrum.YArray;
        var normalizedY = new double[xArray.Length];
        for (int i = 0; i < xArray.Length; i++)
            normalizedY[i] = yArray[i] / maxIntensity;

        var spectrum = new MzSpectrum(xArray, normalizedY, false);
        return new MsDataScan(spectrum, original.OneBasedScanNumber, original.MsnOrder,
            original.IsCentroid, original.Polarity, original.RetentionTime,
            original.ScanWindowRange, original.ScanFilter, original.MzAnalyzer,
            original.TotalIonCurrent, original.InjectionTime, original.NoiseData,
            original.NativeId, original.SelectedIonMZ, original.SelectedIonChargeStateGuess,
            original.SelectedIonIntensity, original.IsolationMz, original.IsolationWidth,
            original.DissociationType, original.OneBasedPrecursorScanNumber,
            original.SelectedIonMonoisotopicGuessMz, original.HcdEnergy,
            original.ScanDescription, original.CompensationVoltage);
    }

    private static List<MatchedFragmentIon> NormalizeIons(List<MatchedFragmentIon> ions, double maxIntensity)
    {
        var normalized = new List<MatchedFragmentIon>(ions.Count);
        foreach (var ion in ions)
        {
            var product = new Product(
                ion.NeutralTheoreticalProduct.ProductType,
                ion.NeutralTheoreticalProduct.Terminus,
                ion.NeutralTheoreticalProduct.NeutralMass,
                ion.NeutralTheoreticalProduct.FragmentNumber,
                ion.NeutralTheoreticalProduct.ResiduePosition,
                ion.NeutralTheoreticalProduct.NeutralLoss,
                ion.NeutralTheoreticalProduct.SecondaryProductType,
                ion.NeutralTheoreticalProduct.SecondaryFragmentNumber);

            double normalizedIntensity = ion.Intensity / maxIntensity;

            if (ion is MatchedFragmentIonWithEnvelope envIon)
            {
                var normalizedEnvelope = envIon.Envelope is null
                    ? null
                    : new IsotopicEnvelope(
                        envIon.Envelope.Peaks.Select(p => (p.mz, p.intensity / maxIntensity)).ToList(),
                        envIon.Envelope.MonoisotopicMass,
                        envIon.Envelope.Charge,
                        envIon.Envelope.TotalIntensity / maxIntensity,
                        envIon.Envelope.Peaks.Count);

                normalized.Add(new MatchedFragmentIonWithEnvelope(product, ion.Mz,
                    normalizedIntensity, ion.Charge, normalizedEnvelope));
            }
            else
            {
                normalized.Add(new MatchedFragmentIon(product, ion.Mz,
                    normalizedIntensity, ion.Charge));
            }
        }
        return normalized;
    }

    private static System.Drawing.Bitmap RenderFrameworkElementToBitmap(FrameworkElement element)
    {
        var dpi = MetaDrawSettings.CanvasPdfExportDpi;
        var dpiScale = dpi / 96.0;

        int width = (int)Math.Ceiling(element.ActualWidth > 0 ? element.ActualWidth : element.Width);
        int height = (int)Math.Ceiling(element.ActualHeight > 0 ? element.ActualHeight : element.Height);

        if (width <= 0)
            width = 700;
        if (height <= 0)
            height = 490;

        var size = new Size(width, height);
        element.Measure(size);
        element.Arrange(new Rect(size));
        element.UpdateLayout();

        var renderTargetBitmap = new RenderTargetBitmap(
            (int)Math.Ceiling(width * dpiScale),
            (int)Math.Ceiling(height * dpiScale),
            dpi,
            dpi,
            PixelFormats.Pbgra32);

        var drawingVisual = new DrawingVisual();
        using (var drawingContext = drawingVisual.RenderOpen())
        {
            var visualBrush = new VisualBrush(element);
            drawingContext.DrawRectangle(visualBrush, null, new Rect(new Point(0, 0), size));
        }

        renderTargetBitmap.Render(drawingVisual);

        var encoder = new PngBitmapEncoder();
        encoder.Frames.Add(BitmapFrame.Create(renderTargetBitmap));

        using var memoryStream = new System.IO.MemoryStream();
        encoder.Save(memoryStream);
        using var temporaryBitmap = new System.Drawing.Bitmap(memoryStream);
        return new System.Drawing.Bitmap(temporaryBitmap);
    }
}
