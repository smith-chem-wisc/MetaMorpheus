using GuiFunctions;
using System;
using System.IO;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media.Imaging;
using System.Windows.Media;
using System.Collections.Generic;
using System.Drawing;
using Brushes = System.Windows.Media.Brushes;
using ImageFormat = System.Drawing.Imaging.ImageFormat;
using Point = System.Windows.Point;
using Size = System.Windows.Size;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ChimeraAnalysisTabView.xaml
    /// </summary>
    public partial class ChimeraAnalysisTabView : UserControl
    {
        public ChimeraAnalysisTabView()
        {
            InitializeComponent();
        }

        private void ChermicDataGrid_OnSelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            if (chermicDataGrid.SelectedItem == null || sender == null || dataContext == null)
            {
                return;
            }

            ChimeraGroupViewModel chimeraGroup = dataContext.SelectedChimeraGroup;
            if (chimeraGroup == null)
            {
                MessageBox.Show("No chimera group found for this PSM");
                return;
            }
            dataContext.Ms1ChimeraPlot = new Ms1ChimeraPlot(ms1ChimeraOverlaPlot, chimeraGroup);
            dataContext.ChimeraSpectrumMatchPlot = new ChimeraSpectrumMatchPlot(ms2ChimeraPlot, chimeraGroup);
            dataContext.ChimeraDrawnSequence = new ChimeraDrawnSequence(chimeraSequenceCanvas, chimeraGroup, dataContext);
        }

        private void ExportMs1(object sender, RoutedEventArgs e)
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            if (dataContext == null || dataContext.SelectedChimeraGroup == null)
            {
                MessageBox.Show("No chimera group selected for export.");
                return;
            }

            string path = Path.Combine(dataContext.ExportDirectory,
                $"{dataContext.SelectedChimeraGroup.FileNameWithoutExtension}_{dataContext.SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{dataContext.SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_MS1.{dataContext.SelectedExportType.ToLower()}");

            switch (dataContext!.SelectedExportType)
            {
                case "Pdf":
                    dataContext.Ms1ChimeraPlot.ExportToPdf(path, (int)dataContext.Ms1ChimeraPlot.Model.Width, (int)dataContext.Ms1ChimeraPlot.Model.Height);
                    break;
                case "Png":
                    dataContext.Ms1ChimeraPlot.ExportToPng(path, (int)dataContext.Ms1ChimeraPlot.Model.Width, (int)dataContext.Ms1ChimeraPlot.Model.Height);
                    break;
                case "Svg":
                    dataContext.Ms1ChimeraPlot.ExportToSvg(path, (int)dataContext.Ms1ChimeraPlot.Model.Width, (int)dataContext.Ms1ChimeraPlot.Model.Height);
                    break;
                default:
                    throw new ArgumentOutOfRangeException();
            }
        }

        private void ExportMs2(object sender, RoutedEventArgs e)
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            if (dataContext == null || dataContext.SelectedChimeraGroup == null)
            {
                MessageBox.Show("No chimera group selected for export.");
                return;
            }

            string path = Path.Combine(dataContext.ExportDirectory,
                $"{dataContext.SelectedChimeraGroup.FileNameWithoutExtension}_{dataContext.SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{dataContext.SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_MS2.{dataContext.SelectedExportType.ToLower()}");

            switch (dataContext!.SelectedExportType)
            {
                case "Pdf":
                    dataContext.ChimeraSpectrumMatchPlot.ExportToPdf(path, (int)dataContext.ChimeraSpectrumMatchPlot.Model.Width, (int)dataContext.ChimeraSpectrumMatchPlot.Model.Height);
                    break;
                case "Png":
                    dataContext.ChimeraSpectrumMatchPlot.ExportToPng(path, (int)dataContext.ChimeraSpectrumMatchPlot.Model.Width, (int)dataContext.ChimeraSpectrumMatchPlot.Model.Height);
                    break;
                case "Svg":
                    dataContext.ChimeraSpectrumMatchPlot.ExportToSvg(path, (int)dataContext.ChimeraSpectrumMatchPlot.Model.Width, (int)dataContext.ChimeraSpectrumMatchPlot.Model.Height);
                    break;
                default:
                    throw new ArgumentOutOfRangeException();
            }
        }

        private void ExportSequenceCoverage(object sender, RoutedEventArgs e)
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            string path = Path.Combine(dataContext.ExportDirectory,
                $"{dataContext.SelectedChimeraGroup.FileNameWithoutExtension}_{dataContext.SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{dataContext.SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_SequenceCoverage.{dataContext.SelectedExportType.ToLower()}");
            // change path to .png
            path = Path.ChangeExtension(path, "png");

            // convert canvas to bitmap
            Rect bounds = VisualTreeHelper.GetDescendantBounds(dataContext.ChimeraDrawnSequence.SequenceDrawingCanvas);
            double dpi = 96d;

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
                VisualBrush vb = new(dataContext.ChimeraDrawnSequence.SequenceDrawingCanvas);
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
        }

        private void ExportLegend(object sender, RoutedEventArgs e)
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            var element = (sender as Button)!.CommandParameter as FrameworkElement;
            if (dataContext == null || element == null)
            {
                MessageBox.Show("No legend available for export.");
                return;
            }
            string path = Path.Combine(dataContext.ExportDirectory,
                $"{dataContext.SelectedChimeraGroup.FileNameWithoutExtension}_{dataContext.SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{dataContext.SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_Legend.{dataContext.SelectedExportType.ToLower()}");


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
        }

        private void ExportAll(object sender, RoutedEventArgs e)
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            var element = (sender as Button)!.CommandParameter as FrameworkElement;

            // initialize all 
            string path = Path.Combine(dataContext!.ExportDirectory,
                $"{dataContext!.SelectedChimeraGroup.FileNameWithoutExtension}_{dataContext.SelectedChimeraGroup.Ms1Scan.OneBasedScanNumber}_{dataContext.SelectedChimeraGroup.Ms2Scan.OneBasedScanNumber}_Combined.{dataContext.SelectedExportType.ToLower()}");

            List<System.Drawing.Bitmap> bitmaps = new();
            List<Point> points = new();
            double dpi = MetaDrawSettings.CanvasPdfExportDpi;
            double outterBuffer = 10;


            // scale and rotate sequence annotation
            Rect annotationBounds = VisualTreeHelper.GetDescendantBounds(dataContext.ChimeraDrawnSequence.SequenceDrawingCanvas);
            var rtb = new RenderTargetBitmap((int)(annotationBounds.Width * dpi / 96), (int)(annotationBounds.Height *
                dpi / 96), dpi, dpi, PixelFormats.Default);

            var dv = new DrawingVisual();
            using (DrawingContext ctx = dv.RenderOpen())
            {
                ctx.DrawRectangle(Brushes.White, null, new Rect(new Point(), annotationBounds.Size));

                VisualBrush vb = new VisualBrush(dataContext.ChimeraDrawnSequence.SequenceDrawingCanvas);
                ctx.DrawRectangle(vb, null, new Rect(new Point(), annotationBounds.Size));
            }
            rtb.Render(dv);

            TransformedBitmap rotatedBitmap = new TransformedBitmap(rtb, new RotateTransform(270));
            TransformedBitmap scaledBitmap = new TransformedBitmap(rotatedBitmap, new ScaleTransform(1, 0.6));
            using (var memory = new MemoryStream())
            {
                BitmapEncoder encoder = new BmpBitmapEncoder();
                encoder.Frames.Add(BitmapFrame.Create(scaledBitmap));
                encoder.Save(memory);
                bitmaps.Add(new System.Drawing.Bitmap(memory));
            }
            var annotationWidth = bitmaps[0].Width;
            var annotationHeight = bitmaps[0].Height;
            points.Add(new Point(outterBuffer, outterBuffer));

            // legend
            Rect legendBounds = VisualTreeHelper.GetDescendantBounds(element);
            rtb = new RenderTargetBitmap((int)(legendBounds.Width * dpi / 96), (int)(legendBounds.Height * dpi / 96), dpi, dpi, PixelFormats.Default);
            dv = new DrawingVisual();
            using (DrawingContext ctx = dv.RenderOpen())
            {
                VisualBrush vb = new VisualBrush(element);
                ctx.DrawRectangle(vb, null, new Rect(new Point(), legendBounds.Size));
            }
            rtb.Render(dv);
            using (var memory = new MemoryStream())
            {
                BitmapEncoder encoder = new BmpBitmapEncoder();
                encoder.Frames.Add(BitmapFrame.Create(rtb));
                encoder.Save(memory);
                bitmaps.Add(new System.Drawing.Bitmap(memory));
            }
            var legendHeight = bitmaps[1].Height;
            var legendWidth = bitmaps[1].Width;
            points.Add(new Point(bitmaps[0].Width + 4 * outterBuffer, annotationHeight - legendHeight + outterBuffer));

            // create temporary exports for spectraS
            // give spectra the remaining space
            string tempDir = Path.Combine(dataContext.ExportDirectory, "temp");
            if (!Directory.Exists(tempDir))
                Directory.CreateDirectory(tempDir);
            int remainingX = (int)(legendWidth);
            double remainingY = annotationHeight - legendHeight;
            int specHeight = (int)(remainingY / 2);

            var ms1TempPath = Path.Combine(tempDir, "ms1.png");
            dataContext.Ms1ChimeraPlot.ExportToPng(ms1TempPath, remainingX, specHeight);
            bitmaps.Add(new Bitmap(ms1TempPath));
            points.Add(new Point(annotationWidth + outterBuffer, outterBuffer));


            var ms2TempPath = Path.Combine(tempDir, "ms2.png");
            dataContext.ChimeraSpectrumMatchPlot.ExportToPng(ms2TempPath, remainingX, specHeight);
            bitmaps.Add(new Bitmap(ms2TempPath));
            points.Add(new Point(annotationWidth + outterBuffer, specHeight));

            var combinedBitmap = new System.Drawing.Bitmap((int)(legendWidth + annotationWidth + 3 * outterBuffer),
                (int)(annotationHeight + 2 * outterBuffer));
            using (var g = System.Drawing.Graphics.FromImage(combinedBitmap))
            {
                //g.ScaleTransform(scaleX, scaleY);
                g.Clear(System.Drawing.Color.White);
                for (int i = 0; i < bitmaps.Count; i++)
                {
                    g.DrawImage(bitmaps[i], (float)points[i].X, (float)points[i].Y);
                }
            }
            combinedBitmap.Save(path, ImageFormat.Png);

            // clean up
            bitmaps.ForEach(b => b.Dispose());
            Directory.Delete(tempDir, true);
        }
    }
}
