using System;
using System.Drawing;
using System.IO;
using System.Linq;
using iText.IO.Image;
using iText.Kernel.Pdf;
using iText.Layout;
using MassSpectrometry;
using MzLibUtil;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Wpf;
using LinearAxis = OxyPlot.Axes.LinearAxis;
using LineSeries = OxyPlot.Series.LineSeries;
using Plot = mzPlot.Plot;
using TextAnnotation = OxyPlot.Annotations.TextAnnotation;

namespace GuiFunctions;

public class MassSpectrumPlot : Plot
{
    public PlotView PlotView { get; protected set; }
    public MsDataScan Scan { get; private set; }
    public MassSpectrumPlot(PlotView plotView, MsDataScan scan) : base(plotView)
    {
        PlotView = plotView;
        Model.Title = string.Empty;
        Model.Subtitle = string.Empty;
        Scan = scan;
        DrawSpectrum();
        RefreshChart();
    }

    /// <summary>
    /// Adds the spectrum from the MSDataScan to the Model
    /// </summary>
    protected void DrawSpectrum()
    {
        var yArray = Scan.MassSpectrum.YArray;
        var xArray = Scan.MassSpectrum.XArray;
        double yMax = 0;
        for (int i = 0; i < yArray.Length; i++)
        {
            if (yArray[i] > yMax) yMax = yArray[i];
        }
        // set up axes
        Model.Axes.Add(new LinearAxis
        {
            Position = AxisPosition.Bottom,
            Title = "m/z",
            Minimum = Scan.ScanWindowRange.Minimum,
            Maximum = Scan.ScanWindowRange.Maximum,
            AbsoluteMinimum = Math.Max(0, Scan.ScanWindowRange.Minimum - 100),
            AbsoluteMaximum = Scan.ScanWindowRange.Maximum + 100,
            MajorTickSize = 2,
            TitleFontWeight = FontWeights.Bold,
            TitleFontSize = MetaDrawSettings.AxisTitleTextSize,
            FontSize = MetaDrawSettings.AxisLabelTextSize,
        });

        Model.Axes.Add(new LinearAxis
        {
            Position = AxisPosition.Left,
            Title = "Intensity",
            Minimum = 0,
            Maximum = yMax,
            AbsoluteMinimum = 0,
            AbsoluteMaximum = yMax * 2,
            StringFormat = "0e-0",
            MajorTickSize = 2,
            TitleFontWeight = FontWeights.Bold,
            TitleFontSize = MetaDrawSettings.AxisTitleTextSize,
            FontSize = MetaDrawSettings.AxisLabelTextSize,
            AxisTitleDistance = 10,
            ExtraGridlines = new double[] { 0 },
            ExtraGridlineColor = OxyColors.Black,
            ExtraGridlineThickness = 1
        });

        // Batch peaks into a single LineSeries for unannotated peaks
        var unannotatedSeries = new LineSeries
        {
            Color = MetaDrawSettings.UnannotatedPeakColor,
            StrokeThickness = MetaDrawSettings.StrokeThicknessUnannotated,
            LineStyle = LineStyle.Solid
        };

        for (int i = 0; i < xArray.Length; i++)
        {
            double mz = xArray[i];
            double intensity = yArray[i];
            // Add as vertical line (two points per peak)
            unannotatedSeries.Points.Add(new DataPoint(mz - 0.0001, 0));
            unannotatedSeries.Points.Add(new DataPoint(mz, intensity));
            unannotatedSeries.Points.Add(new DataPoint(mz + 0.0001, 0));
        }
        Model.Series.Add(unannotatedSeries);
    }

    /// <summary>
    /// Draws a peak on the spectrum
    /// </summary>
    /// <param name="mz">x value of peak to draw</param>
    /// <param name="intensity">y max of peak to draw</param>
    /// <param name="strokeWidth"></param>
    /// <param name="color">Color to draw peak</param>
    /// <param name="annotation">text to display above the peak</param>
    protected void DrawPeak(double mz, double intensity, double strokeWidth, OxyColor color,
        TextAnnotation annotation)
    {
        // peak line
        var line = new LineSeries
        {
            Color = color,
            StrokeThickness = strokeWidth
        };
        line.Points.Add(new DataPoint(mz, 0));
        line.Points.Add(new DataPoint(mz, intensity));

        // Miso is a tag used in chimeric ms1 plotting to indicate that we should not annotate this peak with a label. 
        if (annotation != null && !annotation.Text.Contains("Miso"))
        {
            var x = annotation.TextPosition.X;
            var y = annotation.TextPosition.Y + 20;
            var splits = annotation.Text.Split('\n');

            // Calculate y step for annotation lines based on Y-axis range
            var yAxis = Model.Axes.FirstOrDefault(a => a.Position == AxisPosition.Left);
            double yRange = yAxis != null ? Math.Abs(yAxis.ActualMaximum - yAxis.ActualMinimum) : 1.0;
            double yStep = yRange * 0.05; // 5% of the y-range per line 

            for (int j = splits.Length - 1; j >= 0; j--)
            {
                var split = splits[j];
                var annotationLine = new TextAnnotation
                {
                    Font = annotation.Font,
                    FontSize = annotation.FontSize,
                    FontWeight = annotation.FontWeight,
                    TextColor = annotation.TextColor,
                    StrokeThickness = annotation.StrokeThickness,
                    Text = split,
                    TextPosition = new DataPoint(x, y),
                    TextVerticalAlignment = annotation.TextVerticalAlignment,
                    TextHorizontalAlignment = HorizontalAlignment.Center
                };
                Model.Annotations.Add(annotationLine);
                y += yStep;
            }
        }

        Model.Series.Add(line);
    }

    protected void AnnotateIsolationWindow(DoubleRange isolationRange)
    {
        var yMax = Model.Axes[1].AbsoluteMaximum;
        var points = new DataPoint[]
        {
            new(isolationRange.Minimum, yMax),
            new(isolationRange.Minimum, 0),
            new(isolationRange.Maximum, 0),
            new(isolationRange.Maximum, yMax),
        };
        var lineSeries = new LineSeries
        {
            LineStyle = LineStyle.Dash,
            Color = OxyColors.Red
        };
        lineSeries.Points.AddRange(points);
        Model.Series.Add(lineSeries);
    }

    /// <summary>
    /// Exports plot from a combined bitmap created in the children classes
    /// </summary>
    /// <param name="path"></param>
    /// <param name="combinedBitmaps"></param>
    /// <param name="width"></param>
    /// <param name="height"></param>
    public static void ExportPlot(string path, Bitmap combinedBitmaps, double width = 700, double height = 370)
    {
        width = width > 0 ? width : 700;
        height = height > 0 ? height : 300;
        switch (MetaDrawSettings.ExportType)
        {
            case "Pdf":
                string tempCombinedPath =
                    System.IO.Path.Combine(System.IO.Path.GetDirectoryName(path), "tempCombined.png");
                combinedBitmaps.Save(tempCombinedPath, System.Drawing.Imaging.ImageFormat.Png);

                PdfDocument pdfDoc = new(new PdfWriter(path));
                Document document = new(pdfDoc,
                    new iText.Kernel.Geom.PageSize((float)width - 30, (float)height - 30));

                ImageData sequenceAndLegendImageData = ImageDataFactory.Create(tempCombinedPath);
                iText.Layout.Element.Image sequenceAndPtmLegendImage = new(sequenceAndLegendImageData);
                sequenceAndPtmLegendImage.SetMarginLeft(-30);
                sequenceAndPtmLegendImage.SetMarginTop(-30);
                sequenceAndPtmLegendImage.ScaleToFit((float)width, (float)height);
                document.Add(sequenceAndPtmLegendImage);

                pdfDoc.Close();
                document.Close();
                File.Delete(tempCombinedPath);
                break;

            case "Png":
                combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Png);
                break;

            case "Jpeg":
                combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Jpeg);
                break;

            case "Tiff":
                combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Tiff);
                break;

            case "Wmf":
                combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Wmf);
                break;

            case "Bmp":
                combinedBitmaps.Save(path, System.Drawing.Imaging.ImageFormat.Bmp);
                break;
        }
    }
}