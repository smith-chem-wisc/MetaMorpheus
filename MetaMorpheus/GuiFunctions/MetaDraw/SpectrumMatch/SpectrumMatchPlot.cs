using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Chemistry;
using iText.IO.Image;
using iText.Kernel.Pdf;
using iText.Layout;
using MassSpectrometry;
using mzPlot;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using Readers;
using FontWeights = OxyPlot.FontWeights;
using HorizontalAlignment = OxyPlot.HorizontalAlignment;
using VerticalAlignment = OxyPlot.VerticalAlignment;

namespace GuiFunctions
{
    public class SpectrumMatchPlot : Plot
    {
        public static int MaxCharactersPerDescriptionLine = 32;
        protected List<MatchedFragmentIon> matchedFragmentIons;
        public MsDataScan Scan { get; protected set; }
        public SpectrumMatchFromTsv SpectrumMatch { get; set; }

        /// <summary>
        /// Base Spectrum match constructor
        /// Matched fragment ions should only be passed in for glyco parent/child scan
        /// </summary>
        /// <param name="plotView">Where the plot is going</param>
        /// <param name="sm">sm to plot</param>
        /// <param name="scan">spectrum to plot</param>
        /// <param name="matchedIons">glyco ONLY child matched ions</param>
        public SpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, SpectrumMatchFromTsv sm,
            MsDataScan scan, List<MatchedFragmentIon> matchedIons = null) : base(plotView)
        {
            Model.Title = string.Empty;
            Model.Subtitle = string.Empty;
            Scan = scan;
            matchedFragmentIons = new();

            DrawSpectrum();
            if (matchedIons is null && sm is not null)
            {
                SpectrumMatch = sm;
                matchedFragmentIons = SpectrumMatch.MatchedIons;
                AnnotateMatchedIons(isBetaPeptide: false, matchedFragmentIons);
            }
            else if (matchedIons is not null && sm is not null)
            {
                SpectrumMatch = sm;
                matchedFragmentIons = matchedIons;
                AnnotateMatchedIons(false, matchedIons);
            }

            ZoomAxes();
            RefreshChart();
        }

        /// <summary>
        /// Adds the spectrum from the MSDataScan to the Model
        /// </summary>
        protected void DrawSpectrum()
        {
            // set up axes
            Model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "m/z",
                Minimum = Scan.ScanWindowRange.Minimum,
                Maximum = Scan.ScanWindowRange.Maximum,
                AbsoluteMinimum = Math.Max(0, Scan.ScanWindowRange.Minimum - 100),
                AbsoluteMaximum = Scan.ScanWindowRange.Maximum + 100,
                MajorStep = 200,
                MinorStep = 200,
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
                Maximum = Scan.MassSpectrum.YArray.Max(),
                AbsoluteMinimum = 0,
                AbsoluteMaximum = Scan.MassSpectrum.YArray.Max() * 2,
                MajorStep = Scan.MassSpectrum.YArray.Max() / 10,
                MinorStep = Scan.MassSpectrum.YArray.Max() / 10,
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

            // draw all peaks in the scan
            for (int i = 0; i < Scan.MassSpectrum.XArray.Length; i++)
            {
                double mz = Scan.MassSpectrum.XArray[i];
                double intensity = Scan.MassSpectrum.YArray[i];

                DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessUnannotated,
                    MetaDrawSettings.UnannotatedPeakColor, null);
            }
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
            var line = new LineSeries();
            line.Color = color;
            line.StrokeThickness = strokeWidth;
            line.Points.Add(new DataPoint(mz, 0));
            line.Points.Add(new DataPoint(mz, intensity));

            if (annotation != null)
            {
                Model.Annotations.Add(annotation);
            }

            Model.Series.Add(line);
        }

        /// <summary>
        /// Annotates all matched ion peaks
        /// </summary>
        /// <param name="isBetaPeptide"></param>
        /// <param name="matchedFragmentIons"></param>
        /// <param name="useLiteralPassedValues"></param>
        protected void AnnotateMatchedIons(bool isBetaPeptide, List<MatchedFragmentIon> matchedFragmentIons,
            bool useLiteralPassedValues = false)
        {
            List<MatchedFragmentIon> ionsToDisplay = !MetaDrawSettings.DisplayInternalIons
                ? matchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.SecondaryProductType == null).ToList()
                : matchedFragmentIons;

            foreach (MatchedFragmentIon matchedIon in ionsToDisplay)
            {
                AnnotatePeak(matchedIon, isBetaPeptide, useLiteralPassedValues);
            }
        }

        /// <summary>
        /// Annotates a single matched ion peak
        /// </summary>
        /// <param name="matchedIon">matched ion to annotate</param>
        /// <param name="isBetaPeptide">is a beta x-linked peptide</param>
        /// <param name="useLiteralPassedValues"></param>
        protected void AnnotatePeak(MatchedFragmentIon matchedIon, bool isBetaPeptide,
            bool useLiteralPassedValues = false, OxyColor? ionColorNullable = null)
        {
            OxyColor ionColor;
            if (ionColorNullable == null)
            {
                if (SpectrumMatch.VariantCrossingIons.Contains(matchedIon))
                {
                    ionColor = MetaDrawSettings.VariantCrossColor;
                }
                else if (matchedIon.NeutralTheoreticalProduct.SecondaryProductType != null) //if internal fragment
                {
                    ionColor = MetaDrawSettings.InternalIonColor;
                }
                else if (isBetaPeptide)
                {
                    ionColor = MetaDrawSettings.BetaProductTypeToColor[
                        matchedIon.NeutralTheoreticalProduct.ProductType];
                }
                else
                {
                    ionColor = MetaDrawSettings.ProductTypeToColor[matchedIon.NeutralTheoreticalProduct.ProductType];
                }
            }
            else
            {
                ionColor = (OxyColor)ionColorNullable;
            }

            int i = Scan.MassSpectrum.GetClosestPeakIndex(
                matchedIon.NeutralTheoreticalProduct.NeutralMass.ToMz(matchedIon.Charge));
            double mz = Scan.MassSpectrum.XArray[i];
            double intensity = Scan.MassSpectrum.YArray[i];

            if (useLiteralPassedValues)
            {
                mz = matchedIon.Mz;
                intensity = matchedIon.Intensity;
            }

            // peak annotation
            string prefix = "";
            if (SpectrumMatch.IsCrossLinkedPeptide())
            {
                if (isBetaPeptide)
                {
                    //prefix = "β";
                    prefix = "B-";
                }
                else
                {
                    //prefix = "α";
                    prefix = "A-";
                }
            }

            var peakAnnotation = new TextAnnotation();
            if (MetaDrawSettings.DisplayIonAnnotations)
            {
                string peakAnnotationText = prefix;

                if (MetaDrawSettings.SubAndSuperScriptIons)
                    foreach (var character in matchedIon.NeutralTheoreticalProduct.Annotation)
                    {
                        if (char.IsDigit(character))
                            peakAnnotationText += MetaDrawSettings.SubScriptNumbers[character - '0'];
                        else switch (character)
                        {
                            case '-':
                                peakAnnotationText += "\u208B"; // sub scripted Hyphen
                                break;
                            case '[':
                            case ']':
                                continue;
                            default:
                                peakAnnotationText += character;
                                break;
                        }
                    }
                else
                    peakAnnotationText += matchedIon.NeutralTheoreticalProduct.Annotation;
                
                if (matchedIon.NeutralTheoreticalProduct.NeutralLoss != 0 &&
                    !peakAnnotationText.Contains("-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss.ToString("F2")))
                {
                    peakAnnotationText += "-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss.ToString("F2");
                }

                if (MetaDrawSettings.AnnotateCharges)
                {
                    char chargeAnnotation = matchedIon.Charge > 0 ? '+' : '-';
                    if (MetaDrawSettings.SubAndSuperScriptIons)
                    {
                        var superScript = new string(Math.Abs(matchedIon.Charge).ToString()
                            .Select(digit => MetaDrawSettings.SuperScriptNumbers[digit - '0'])
                            .ToArray());

                        peakAnnotationText += superScript;
                        if (chargeAnnotation == '+')
                            peakAnnotationText += (char)(chargeAnnotation + 8271);
                        else
                            peakAnnotationText += (char)(chargeAnnotation + 8270);
                    }
                    else
                        peakAnnotationText += chargeAnnotation.ToString() + matchedIon.Charge;
                }

                if (MetaDrawSettings.AnnotateMzValues)
                {
                    peakAnnotationText += " (" + matchedIon.Mz.ToString("F3") + ")";
                }


                peakAnnotation.Font = "Arial";
                peakAnnotation.FontSize = MetaDrawSettings.AnnotatedFontSize;
                peakAnnotation.FontWeight = MetaDrawSettings.AnnotationBold ? FontWeights.Bold : 2.0;
                peakAnnotation.TextColor = ionColor;
                peakAnnotation.StrokeThickness = 0;
                peakAnnotation.Text = peakAnnotationText;
                peakAnnotation.TextPosition = new DataPoint(mz, intensity);
                peakAnnotation.TextVerticalAlignment = intensity < 0 ? VerticalAlignment.Top : VerticalAlignment.Bottom;
                peakAnnotation.TextHorizontalAlignment = HorizontalAlignment.Center;
            }
            else
            {
                peakAnnotation.Text = string.Empty;
            }

            if (matchedIon.NeutralTheoreticalProduct.SecondaryProductType != null &&
                !MetaDrawSettings.DisplayInternalIonAnnotations) //if internal fragment
            {
                peakAnnotation.Text = string.Empty;
            }

            DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessAnnotated, ionColor, peakAnnotation);
        }

        /// <summary>
        /// Zooms the axis of the graph to the matched ions
        /// </summary>
        /// <param name="yZoom"></param>
        /// <param name="matchedFramgentIons">ions to zoom to. if null, it will used the stored protected matchedFragmentIons</param>
        protected void ZoomAxes(IEnumerable<MatchedFragmentIon> matchedFramgentIons = null, double yZoom = 1.2)
        {
            matchedFramgentIons ??= matchedFragmentIons;
            double highestAnnotatedIntensity = 0;
            double highestAnnotatedMz = double.MinValue;
            double lowestAnnotatedMz = double.MaxValue;

            foreach (var ion in matchedFragmentIons)
            {
                double mz = ion.NeutralTheoreticalProduct.NeutralMass.ToMz(ion.Charge);
                int i = Scan.MassSpectrum.GetClosestPeakIndex(mz);
                double intensity = Scan.MassSpectrum.YArray[i];

                if (intensity > highestAnnotatedIntensity)
                {
                    highestAnnotatedIntensity = intensity;
                }

                if (highestAnnotatedMz < mz)
                {
                    highestAnnotatedMz = mz;
                }

                if (mz < lowestAnnotatedMz)
                {
                    lowestAnnotatedMz = mz;
                }
            }

            if (highestAnnotatedIntensity > 0)
            {
                Model.Axes[1].Zoom(0, highestAnnotatedIntensity * yZoom);
            }

            if (highestAnnotatedMz > double.MinValue && lowestAnnotatedMz < double.MaxValue)
            {
                Model.Axes[0].Zoom(lowestAnnotatedMz - 100, highestAnnotatedMz + 100);
            }
        }

        /// <summary>
        /// Exports plot from a combined bitmap created in the children classes
        /// </summary>
        /// <param name="path"></param>
        /// <param name="combinedBitmaps"></param>
        /// <param name="width"></param>
        /// <param name="height"></param>
        public void ExportPlot(string path, Bitmap combinedBitmaps, double width = 700, double height = 370)
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

        protected void AnnotateLibraryIons(bool isBetaPeptide, List<MatchedFragmentIon> libraryIons)
        {
            // figure out the sum of the intensities of the matched fragment ions
            double sumOfMatchedIonIntensities = 0;
            double sumOfLibraryIntensities = 0;
            foreach (var libraryIon in libraryIons)
            {
                var matchedIon = SpectrumMatch.MatchedIons.FirstOrDefault(p =>
                    p.NeutralTheoreticalProduct.ProductType == libraryIon.NeutralTheoreticalProduct.ProductType
                    && p.NeutralTheoreticalProduct.FragmentNumber ==
                    libraryIon.NeutralTheoreticalProduct.FragmentNumber);

                if (matchedIon == null)
                {
                    continue;
                }

                int i = Scan.MassSpectrum.GetClosestPeakIndex(libraryIon.Mz);
                double intensity = Scan.MassSpectrum.YArray[i];
                sumOfMatchedIonIntensities += intensity;
                sumOfLibraryIntensities += libraryIon.Intensity;
            }

            double multiplier = -1 * sumOfMatchedIonIntensities / sumOfLibraryIntensities;

            List<MatchedFragmentIon> mirroredLibraryIons = new List<MatchedFragmentIon>();

            foreach (MatchedFragmentIon libraryIon in libraryIons)
            {
                var neutralProduct = new Product(libraryIon.NeutralTheoreticalProduct.ProductType,
                    libraryIon.NeutralTheoreticalProduct.Terminus,
                    libraryIon.NeutralTheoreticalProduct.NeutralMass,
                    libraryIon.NeutralTheoreticalProduct.FragmentNumber,
                    libraryIon.NeutralTheoreticalProduct.ResiduePosition,
                    libraryIon.NeutralTheoreticalProduct.NeutralLoss);

                mirroredLibraryIons.Add(new MatchedFragmentIon(neutralProduct, libraryIon.Mz,
                    multiplier * libraryIon.Intensity, libraryIon.Charge));
            }

            AnnotateMatchedIons(isBetaPeptide, mirroredLibraryIons, useLiteralPassedValues: true);

            // zoom to accomodate the mirror plot
            double min = mirroredLibraryIons.Min(p => p.Intensity) * 1.2;
            Model.Axes[1].AbsoluteMinimum = min * 2;
            Model.Axes[1].AbsoluteMaximum = -min * 2;
            Model.Axes[1].Zoom(min, -min);
            Model.Axes[1].LabelFormatter = DrawnSequence.YAxisLabelFormatter;
        }

        protected void AnnotateProperties(LibrarySpectrum librarySpectrum = null)
        {
            StringBuilder text = new StringBuilder();
            if (MetaDrawSettings.SpectrumDescription["Precursor Charge: "])
            {
                text.Append("Precursor Charge: ");
                text.Append(SpectrumMatch.PrecursorCharge);
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["Precursor Mass: "])
            {
                text.Append("Precursor Mass: ");
                text.Append(SpectrumMatch.PrecursorMass.ToString("F3"));
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["Theoretical Mass: "])
            {
                text.Append("Theoretical Mass: ");
                text.Append(
                    double.TryParse(SpectrumMatch.MonoisotopicMassString, NumberStyles.Any, CultureInfo.InvariantCulture,
                        out var monoMass)
                        ? monoMass.ToString("F3")
                        : SpectrumMatch.MonoisotopicMassString);
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["Protein Accession: "])
            {
                text.Append("Protein Accession: ");
                if (SpectrumMatch.Accession.Length > 10)
                {
                    text.Append("\r\n   " + SpectrumMatch.Accession);
                }
                else
                    text.Append(SpectrumMatch.Accession);

                text.Append("\r\n");
            }

            if (SpectrumMatch.Name != null && MetaDrawSettings.SpectrumDescription["Protein: "])
            {
                text.Append("Protein: ");
                text.Append(SpectrumMatch.Name);
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["Retention Time: "])
            {
                text.Append("Retention Time: ");
                text.Append(SpectrumMatch.RetentionTime);
                text.Append("\r\n");
            }

            if (SpectrumMatch.OneOverK0 != null && MetaDrawSettings.SpectrumDescription["1/K\u2080: "])
            {
                text.Append("1/K\u2080: ");
                text.Append(SpectrumMatch.OneOverK0);
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["Decoy/Contaminant/Target: "])
            {
                text.Append("Decoy/Contaminant/Target: ");
                text.Append(SpectrumMatch.DecoyContamTarget);
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["Sequence Length: "])
            {
                text.Append("Sequence Length: ");
                text.Append(SpectrumMatch.BaseSeq.Length.ToString("F3").Split('.')[0]);
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["Ambiguity Level: "])
            {
                text.Append("Ambiguity Level: ");
                text.Append(SpectrumMatch.AmbiguityLevel);
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["Spectral Angle: "])
            {
                if (SpectrumMatch.SpectralAngle != null)
                {
                    text.Append("Original Spectral Angle: ");
                    text.Append(SpectrumMatch.SpectralAngle.ToString() + "\r\n");
                }

                if (librarySpectrum != null)
                {
                    text.Append("Displayed Spectral Angle: ");
                    text.Append(librarySpectrum.CalculateSpectralAngleOnTheFly(this.matchedFragmentIons) + "\r\n");
                }
            }

            if (MetaDrawSettings.SpectrumDescription["Score: "])
            {
                text.Append("Score: ");
                text.Append(SpectrumMatch.Score.ToString("F3"));
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["Q-Value: "])
            {
                text.Append("Q-Value: ");
                text.Append(SpectrumMatch.QValue.ToString("F3"));
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["PEP: "])
            {
                text.Append("PEP: ");
                text.Append(SpectrumMatch.PEP.ToString("F3"));
                text.Append("\r\n");
            }

            if (MetaDrawSettings.SpectrumDescription["PEP Q-Value: "])
            {
                text.Append("PEP Q-Value: ");
                text.Append(SpectrumMatch.PEP_QValue.ToString("F3"));
                text.Append("\r\n");
            }

            double fontSize = MetaDrawSettings.SpectrumDescriptionFontSize;
            string annotationText = text.ToString();
            double averageCharWidth = 0.48; // Magic number determined by trial and error
            int maxLineLength = Math.Min(annotationText.Split('\n').Max(line => line.Length), MaxCharactersPerDescriptionLine);
            double estimatedWidth = fontSize * averageCharWidth * maxLineLength;

            // Set X offset to negative estimated width minus a small margin
            double xOffset = -estimatedWidth - 10 - (fontSize / 3);

            var annotation = new PlotTextAnnotation()
            {
                Text = text.ToString(),
                XPosition = PlotTextAnnotation.RelativeToX.Right,
                YPosition = PlotTextAnnotation.RelativeToY.Top,
                FontWeight = 450,
                X = xOffset,
                Font = "Arial",
                FontSize = MetaDrawSettings.SpectrumDescriptionFontSize,
                TextColor = OxyColors.Black, 
            };

            Model.Annotations.Add(annotation);
        }
    }

    //TODO: move this to mzLib (https://github.com/smith-chem-wisc/mzLib/blob/master/mzPlot/Annotations/PlotTextAnnotation.cs)
    public class PlotTextAnnotation : Annotation
    {
        public enum RelativeToX { Left, Center, Right }
        public enum RelativeToY { Top, Center, Bottom }

        /// <summary>
        /// Creates a "floating" text annotation, i.e., the annotation's x,y position is coupled to the chart area and not a data point.
        /// </summary>
        public PlotTextAnnotation()
        {

        }

        public string Text { get; set; }

        /// <summary>
        /// The x-position of the annotation, in pixels. Relative to the left of the chart (higher X is more to the right). Can be negative.
        /// </summary>
        public double X { get; set; }

        /// <summary>
        /// The y-position of the annotation, in pixels. Relative to the top of the chart (higher Y is lower on the chart). Can be negative.
        /// </summary>
        public double Y { get; set; }

        public RelativeToX XPosition { get; set; } = RelativeToX.Left;
        public RelativeToY YPosition { get; set; } = RelativeToY.Top;

        public override void Render(IRenderContext rc)
        {
            base.Render(rc);
            double pX = 0;
            double pY = 0;

            switch (XPosition)
            {
                case RelativeToX.Left: pX = PlotModel.PlotArea.Left; break;
                case RelativeToX.Center: pX = PlotModel.PlotArea.Center.X; break;
                case RelativeToX.Right: pX = PlotModel.PlotArea.Right; break;
            }

            switch (YPosition)
            {
                case RelativeToY.Top: pY = PlotModel.PlotArea.Top; break;
                case RelativeToY.Center: pY = PlotModel.PlotArea.Center.Y; break;
                case RelativeToY.Bottom: pY = PlotModel.PlotArea.Bottom; break;
            }

            pX += X;
            pY += Y;

            // Ensure text does not overlay in Y dimension. 
            double lineSpacing = 1.2; // 20% extra space
            double lineHeight = FontSize * lineSpacing;

            string processedText = WrapLongLines(Text, SpectrumMatchPlot.MaxCharactersPerDescriptionLine, SpectrumMatchPlot.MaxCharactersPerDescriptionLine);
            var lines = (processedText ?? string.Empty).Split(new[] { "\r\n", "\n" }, StringSplitOptions.None);

            for (int i = 0; i < lines.Length; i++)
            {
                double yOffset = i * lineHeight;
                rc.DrawText(
                    new ScreenPoint(pX, pY + yOffset),
                    lines[i],
                    TextColor,
                    Font,
                    FontSize,
                    FontWeight);
            }
        }

        private string WrapLongLines(string text, int firstLineLimit = 32, int wrapLineLimit = 29, string indent = "   ")
        {
            var lines = text.Split(new[] { "\r\n", "\n" }, StringSplitOptions.None);
            var result = new StringBuilder();

            foreach (var line in lines)
            {
                if (line.Length > firstLineLimit)
                {
                    int written = 0;
                    int remaining = line.Length;
                    // First line
                    result.Append(line.Substring(0, firstLineLimit));
                    written += firstLineLimit;
                    remaining -= firstLineLimit;

                    // Subsequent lines
                    while (remaining > 0)
                    {
                        int take = Math.Min(wrapLineLimit - indent.Length, remaining);
                        result.Append("\r\n" + indent + line.Substring(written, take));
                        written += take;
                        remaining -= take;
                    }
                    result.Append("\r\n");
                }
                else
                {
                    result.Append(line + "\r\n");
                }
            }
            return result.ToString().TrimEnd('\r', '\n');
        }

    }
}
