using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Chemistry;
using EngineLayer;
using MassSpectrometry;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using OxyPlot;
using OxyPlot.Wpf;
using Readers;
using Annotation = OxyPlot.Annotations.Annotation;
using FontWeights = OxyPlot.FontWeights;
using HorizontalAlignment = OxyPlot.HorizontalAlignment;
using TextAnnotation = OxyPlot.Annotations.TextAnnotation;
using VerticalAlignment = OxyPlot.VerticalAlignment;

namespace GuiFunctions
{
    public class SpectrumMatchPlot : MassSpectrumPlot
    {
        public static int MaxCharactersPerDescriptionLine = 32;
        public List<MatchedFragmentIon> MatchedFragmentIons { get; protected set; }
        public SpectrumMatchFromTsv SpectrumMatch { get; set; }

        /// <summary>
        /// Base Spectrum match constructor
        /// Matched fragment ions should only be passed in for glyco parent/child scan
        /// </summary>
        /// <param name="plotView">Where the plot is going</param>
        /// <param name="sm">sm to plot</param>
        /// <param name="scan">spectrum to plot</param>
        /// <param name="matchedIons">glyco ONLY child matched ions</param>
        public SpectrumMatchPlot(PlotView plotView, SpectrumMatchFromTsv sm,
            MsDataScan scan, List<MatchedFragmentIon> matchedIons = null) : base(plotView, scan)
        {
            Model.Title = string.Empty;
            Model.Subtitle = string.Empty;

            if (matchedIons is null && sm is not null)
            {
                SpectrumMatch = sm;
                MatchedFragmentIons = SpectrumMatch.MatchedIons;
                AnnotateMatchedIons(isBetaPeptide: false, MatchedFragmentIons);
            }
            else if (matchedIons is not null && sm is not null)
            {
                SpectrumMatch = sm;
                MatchedFragmentIons = matchedIons;
                AnnotateMatchedIons(false, matchedIons);
            }
            else
                MatchedFragmentIons = new();

            ZoomAxes();
            RefreshChart();
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
            if (MetaDrawSettings.AnnotateIsotopicEnvelopes && Scan != null)
            {
                PopulateEnvelopesForScanIfNeeded(Scan, matchedFragmentIons);
            }

            List<MatchedFragmentIon> ionsToDisplay;
            if (!MetaDrawSettings.DisplayInternalIons)
            {
                ionsToDisplay = new List<MatchedFragmentIon>(matchedFragmentIons.Count);
                foreach (var p in matchedFragmentIons)
                {
                    if (p.NeutralTheoreticalProduct.SecondaryProductType == null)
                        ionsToDisplay.Add(p);
                }
            }
            else
            {
                ionsToDisplay = matchedFragmentIons;
            }

            foreach (MatchedFragmentIon matchedIon in ionsToDisplay)
            {
                AnnotatePeak(matchedIon, isBetaPeptide, useLiteralPassedValues);
            }
        }

        /// <summary>
        /// Lazily populates envelopes for matched ions if needed
        /// </summary>
        private void PopulateEnvelopesForScanIfNeeded(MsDataScan scan, List<MatchedFragmentIon> ions)
        {
            var ionsNeedingEnvelopes = ions
                .OfType<MatchedFragmentIonWithEnvelope>()
                .Where(p => p.Envelope == null)
                .ToList();

            if (ionsNeedingEnvelopes.Count == 0)
                return;

            var products = ionsNeedingEnvelopes
                .Select(p => p.NeutralTheoreticalProduct)
                .Distinct()
                .ToList();

            var commonParams = new CommonParameters(
                precursorDeconParams: MetaDrawSettingsViewModel.Instance?.DeconHostViewModel?.PrecursorDeconvolutionParameters?.Parameters,
                productDeconParams: MetaDrawSettingsViewModel.Instance?.DeconHostViewModel?.ProductDeconvolutionParameters?.Parameters);

            var specificMass = new Ms2ScanWithSpecificMass(Scan, SpectrumMatch.PrecursorMz,
                SpectrumMatch.PrecursorCharge, SpectrumMatch.FileNameWithoutExtension, commonParams);

            var allMatches = MetaMorpheusEngine.MatchFragmentIons(specificMass, products, commonParams, matchAllCharges: false, includeExperimentalEnvelope: true);

            var productToMatch = allMatches
                .OfType<MatchedFragmentIonWithEnvelope>()
                .Where(p => p.Envelope != null)
                .ToDictionary(p => p.NeutralTheoreticalProduct);

            foreach (var ion in ionsNeedingEnvelopes)
            {
                if (productToMatch.TryGetValue(ion.NeutralTheoreticalProduct, out var match))
                {
                    ion.Envelope = match.Envelope;
                }
            }
        }

        /// <summary>
        /// Annotates a single matched ion peak
        /// </summary>
        /// <param name="matchedIon">matched ion to annotate</param>
        /// <param name="isBetaPeptide">is a beta x-linked peptide</param>
        /// <param name="useLiteralPassedValues"></param>
        protected void AnnotatePeak(MatchedFragmentIon matchedIon, bool isBetaPeptide,
            bool useLiteralPassedValues = false, OxyColor? ionColorNullable = null, string annotation = null)
        {
            OxyColor ionColor;
            if (ionColorNullable == null)
            {
                if (SpectrumMatch.VariantCrossingIons is not null && SpectrumMatch.VariantCrossingIons.Contains(matchedIon))
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

            // peak annotation
            string prefix = string.Empty;
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

            // Check for isotopic envelope
            if (MetaDrawSettings.AnnotateIsotopicEnvelopes && 
                matchedIon is MatchedFragmentIonWithEnvelope envIon && 
                envIon.Envelope != null)
            {
                AnnotateEnvelopePeaks(envIon, ionColor, isBetaPeptide, useLiteralPassedValues, prefix, annotation);
                return;
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

            var peakAnnotation = new TextAnnotation();
            string peakAnnotationText;

            // Fast path: direct annotation
            if (annotation != null)
            {
                peakAnnotationText = prefix + annotation;
                intensity += intensity * 0.05;
                peakAnnotation.TextColor = ionColor;
            }
            else
            {
                peakAnnotationText = BuildAnnotationText(matchedIon, isBetaPeptide, prefix, null);
                peakAnnotation.TextColor = ionColor;
            }

            // Hide internal fragment annotation if not displaying
            if (matchedIon.NeutralTheoreticalProduct.SecondaryProductType != null &&
                !MetaDrawSettings.DisplayInternalIonAnnotations)
            {
                peakAnnotationText = string.Empty;
            }

            // Set annotation properties
            peakAnnotation.Text = peakAnnotationText;
            peakAnnotation.Font = "Arial";
            peakAnnotation.FontSize = MetaDrawSettings.AnnotatedFontSize;
            peakAnnotation.FontWeight = MetaDrawSettings.AnnotationBold ? FontWeights.Bold : 2.0;
            peakAnnotation.StrokeThickness = 0;
            peakAnnotation.TextPosition = new DataPoint(mz, intensity);
            peakAnnotation.TextVerticalAlignment = intensity < 0 ? VerticalAlignment.Top : VerticalAlignment.Bottom;
            peakAnnotation.TextHorizontalAlignment = HorizontalAlignment.Center;

            DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessAnnotated, ionColor, peakAnnotation);
        }

        /// <summary>
        /// Annotates all peaks in an isotopic envelope, with text annotation only on the tallest peak
        /// </summary>
        private void AnnotateEnvelopePeaks(MatchedFragmentIonWithEnvelope envIon, OxyColor ionColor, 
            bool isBetaPeptide, bool useLiteralPassedValues, string prefix, string directAnnotation)
        {
            var envelope = envIon.Envelope!;
            var product = envIon.NeutralTheoreticalProduct;

            // Find the tallest peak in the envelope
            var tallestPeak = envelope.Peaks.MaxBy(p => p.intensity);
            if (tallestPeak.intensity == 0)
                return;

            // Annotate all peaks in the envelope
            foreach (var (envMz, envIntensity) in envelope.Peaks)
            {
                bool isTallest = Math.Abs(envMz - tallestPeak.mz) < 0.001 && Math.Abs(envIntensity - tallestPeak.intensity) < 0.001;

                // Draw peak marker (colored)
                DrawPeak(envMz, envIntensity, MetaDrawSettings.StrokeThicknessAnnotated, ionColor, null);

                // Draw TEXT annotation only on tallest peak
                if (isTallest)
                {                // Hide internal fragment annotation if not displaying
                    string peakAnnotationText;
                    if (envIon.NeutralTheoreticalProduct.SecondaryProductType != null &&
                        !MetaDrawSettings.DisplayInternalIonAnnotations)
                    {
                        peakAnnotationText = string.Empty;
                    }
                    else
                    {
                        peakAnnotationText = BuildAnnotationText(envIon, isBetaPeptide, prefix, directAnnotation);
                    }
                    
                    var peakAnnotation = new TextAnnotation
                    {
                        Text = peakAnnotationText,
                        Font = "Arial",
                        FontSize = MetaDrawSettings.AnnotatedFontSize,
                        FontWeight = MetaDrawSettings.AnnotationBold ? FontWeights.Bold : 2.0,
                        StrokeThickness = 0,
                        TextPosition = new DataPoint(envMz, envIntensity),
                        TextVerticalAlignment = envIntensity < 0 ? VerticalAlignment.Top : VerticalAlignment.Bottom,
                        TextHorizontalAlignment = HorizontalAlignment.Center,
                        TextColor = ionColor
                    };

                    DrawPeak(envMz, envIntensity, MetaDrawSettings.StrokeThicknessAnnotated, ionColor, peakAnnotation);
                }
            }
        }

        /// <summary>
        /// Builds the annotation text for a matched fragment ion
        /// </summary>
        private string BuildAnnotationText(MatchedFragmentIon matchedIon, bool isBetaPeptide, string prefix, string directAnnotation)
        {
            if (directAnnotation != null)
                return prefix + directAnnotation;

            if (!MetaDrawSettings.DisplayIonAnnotations)
                return string.Empty;

            string peakAnnotationText = prefix;

            // Fragment Number annotation
            if (MetaDrawSettings.SubAndSuperScriptIons)
            {
                foreach (var character in matchedIon.NeutralTheoreticalProduct.Annotation)
                {
                    if (char.IsDigit(character))
                        peakAnnotationText += MetaDrawSettings.SubScriptNumbers[character - '0'];
                    else switch (character)
                    {
                        case '-':
                            peakAnnotationText += "\u208B";
                            break;
                        case '[':
                        case ']':
                            continue;
                        default:
                            peakAnnotationText += character;
                            break;
                    }
                }
            }
            else
            {
                peakAnnotationText += matchedIon.NeutralTheoreticalProduct.Annotation;
            }

            // Charge annotation
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

            // m/z annotation
            if (MetaDrawSettings.AnnotateMzValues)
            {
                peakAnnotationText += " (" + matchedIon.Mz.ToString("F3") + ")";
            }

            return peakAnnotationText;
        }

        /// <summary>
        /// Zooms the axis of the graph to the matched ions
        /// </summary>
        /// <param name="yZoom"></param>
        /// <param name="matchedFramgentIons">ions to zoom to. if null, it will used the stored protected MatchedFragmentIons</param>
        protected void ZoomAxes(IEnumerable<MatchedFragmentIon> matchedFramgentIons = null, double yZoom = 1.2)
        {
            matchedFramgentIons ??= MatchedFragmentIons;
            double highestAnnotatedIntensity = 0;
            double highestAnnotatedMz = double.MinValue;
            double lowestAnnotatedMz = double.MaxValue;
            var yArray = Scan.MassSpectrum.YArray;

            foreach (var ion in matchedFramgentIons)
            {
                double mz = ion.NeutralTheoreticalProduct.NeutralMass.ToMz(ion.Charge);
                int i = Scan.MassSpectrum.GetClosestPeakIndex(mz);
                double intensity = yArray[i];

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
                    text.Append(librarySpectrum.CalculateSpectralAngleOnTheFly(this.MatchedFragmentIons) + "\r\n");
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
