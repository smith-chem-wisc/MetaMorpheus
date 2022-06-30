using EngineLayer;
using OxyPlot;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace GuiFunctions
{

    /// <summary>
    /// Class for drawing the sequences within metadraw, both scrollable and stationary
    /// </summary>
    public class DrawnSequence
    {
        public Canvas SequenceDrawingCanvas;
        public bool Stationary;
        public bool Annotation;
        public PsmFromTsv SpectrumMatch;
        public DrawnSequence(Canvas sequenceDrawingCanvas, PsmFromTsv psm, bool stationary, bool annotation = false)
        {
            SequenceDrawingCanvas = sequenceDrawingCanvas;
            SpectrumMatch = psm;
            Stationary = stationary;
            Annotation = annotation;
            SequenceDrawingCanvas.Width = 600;
            SequenceDrawingCanvas.Height = 60;

            if (Annotation)
            {
                DrawSequenceAnnotation(psm, this);
            }
            else if (Stationary)
            {
                DrawStationarySequence(psm, this, 10);
            }
            else
            {
                AnnotateBaseSequence(psm.BaseSeq, psm.FullSequence, 10, psm.MatchedIons, SpectrumMatch);
            }
        }

        /// <summary>
        /// Draws the Letters and matched ions for each sequence
        /// </summary>
        /// <param name="baseSequence"></param>
        /// <param name="fullSequence"></param>
        /// <param name="yLoc"></param>
        /// <param name="matchedFragmentIons"></param>
        /// <param name="canvas"></param>
        /// <param name="psm"></param>
        public void AnnotateBaseSequence(string baseSequence, string fullSequence, int yLoc, List<MatchedFragmentIon> matchedFragmentIons, PsmFromTsv psm, 
            bool stationary = false, int annotationRow = 0, int chunkPositionInRow = 0)
        {
            if (!Annotation && (psm.BetaPeptideBaseSequence == null || !psm.BetaPeptideBaseSequence.Equals(baseSequence)))
            {
                ClearCanvas(SequenceDrawingCanvas);
            }
            double canvasWidth = SequenceDrawingCanvas.Width;
            int spacing = 12;

            // draw initial amino acid number
            if (stationary && MetaDrawSettings.DrawNumbersUnderStationary)
            {
                int psmStartResidue = int.Parse(psm.StartAndEndResiduesInProtein.Split("to")[0].Replace("[", ""));
                var startAA = (MetaDrawSettings.FirstAAonScreenIndex + psmStartResidue).ToString().ToCharArray().Reverse().ToArray();
                double x = 22;

                Polyline line = new Polyline();
                line.Points = new PointCollection() { new Point(x, yLoc + 25), new Point(x, yLoc + 35) };
                line.Stroke = new SolidColorBrush(Colors.Black);
                line.StrokeThickness = 3;
                SequenceDrawingCanvas.Children.Add(line);
                for (int i = 0; i < startAA.Length; i++)
                {
                    x = startAA.Length < 2 ? (i + 1) * -spacing + 22 : i  * -spacing + 14;
                    if (startAA.Length > 2)
                        x += 14;
                    DrawText(SequenceDrawingCanvas, new Point(x, yLoc + 35), startAA[i].ToString(), Brushes.Black);
                }
            }

            // draw base sequence
            for (int r = 0; r < baseSequence.Length; r++)
            {
                double x = r * MetaDrawSettings.AnnotatedSequenceTextSpacing + 10;

                // adjust for spacing in sequence annotation
                if (Annotation)
                {
                    x = ((r + (MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment + 1) * chunkPositionInRow) * MetaDrawSettings.AnnotatedSequenceTextSpacing) + 10 ;
                }

                DrawText(SequenceDrawingCanvas, new Point(x, yLoc), baseSequence[r].ToString(), Brushes.Black);
                
                canvasWidth = x;
            }
            SequenceDrawingCanvas.Width = Math.Max(SequenceDrawingCanvas.Width, canvasWidth) + 185; // this number is the width of the grayed out box

            // draw final amino acid number
            if (stationary && MetaDrawSettings.DrawNumbersUnderStationary)
            {
                int psmStartResidue = int.Parse(psm.StartAndEndResiduesInProtein.Split("to")[0].Replace("[", ""));
                var endAA = (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen + psmStartResidue - 1).ToString();
                canvasWidth += spacing;
                double x = canvasWidth;

                Polyline line = new Polyline();
                line.Points = new PointCollection() { new Point(x, yLoc + 25), new Point(x, yLoc + 35) };
                line.Stroke = new SolidColorBrush(Colors.Black);
                line.StrokeThickness = 3;
                SequenceDrawingCanvas.Children.Add(line);
                for (int i = 0; i < endAA.Length; i++)
                {
                    x = spacing * (i) + canvasWidth - 22;
                    if (endAA.Length > 2)
                        x -= 14;
                    DrawText(SequenceDrawingCanvas, new Point(x, yLoc + 35), endAA[i].ToString(), Brushes.Black);
                }
            }

            if (MetaDrawSettings.DrawMatchedIons)
            {
                // draw the fragment ion annotations on the base sequence
                foreach (var ion in matchedFragmentIons)
                {
                    //if it's not an internal fragment
                    if (ion.NeutralTheoreticalProduct.SecondaryProductType == null)
                    {
                        string annotation = ion.NeutralTheoreticalProduct.ProductType + "" + ion.NeutralTheoreticalProduct.FragmentNumber;
                        OxyColor oxycolor = psm.VariantCrossingIons.Contains(ion) ?
                            MetaDrawSettings.VariantCrossColor : MetaDrawSettings.ProductTypeToColor[ion.NeutralTheoreticalProduct.ProductType];
                        Color color = Color.FromArgb(oxycolor.A, oxycolor.R, oxycolor.G, oxycolor.B);

                        if (ion.NeutralTheoreticalProduct.NeutralLoss != 0)
                        {
                            annotation += "-" + ion.NeutralTheoreticalProduct.NeutralLoss;
                        }

                        int residue;
                        if (Stationary)
                        {
                            residue = ion.NeutralTheoreticalProduct.AminoAcidPosition - MetaDrawSettings.FirstAAonScreenIndex;
                        }
                        else if (Annotation)
                        {
                            residue = ion.NeutralTheoreticalProduct.AminoAcidPosition + (chunkPositionInRow) - (MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment * MetaDrawSettings.SequenceAnnotationSegmentPerRow * annotationRow);
                        }
                        else
                        {
                            residue = ion.NeutralTheoreticalProduct.AminoAcidPosition;
                        }
                        
                        double x = residue * MetaDrawSettings.AnnotatedSequenceTextSpacing + 11;
                        double y = yLoc + MetaDrawSettings.ProductTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType];

                        if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.C)
                        {
                            DrawCTermIon(SequenceDrawingCanvas, new Point(x, y), color, annotation);
                        }
                        else if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N)
                        {
                            DrawNTermIon(SequenceDrawingCanvas, new Point(x, y), color, annotation);
                        }
                        // don't draw diagnostic ions, precursor ions, etc
                    }
                }
            }
            AnnotateModifications(psm, SequenceDrawingCanvas, fullSequence, yLoc);
        }

        /// <summary>
        /// Draws the Modifications for each sequence
        /// </summary>
        /// <param name="spectrumMatch"></param>
        /// <param name="sequenceDrawingCanvas"></param>
        /// <param name="fullSequence"></param>
        /// <param name="yLoc"></param>
        /// <param name="spacer"></param>
        /// <param name="xShift"></param>
        public static void AnnotateModifications(PsmFromTsv spectrumMatch, Canvas sequenceDrawingCanvas, string fullSequence, int yLoc, double? spacer = null, int xShift = 12)
        {
            var peptide = new PeptideWithSetModifications(fullSequence, GlobalVariables.AllModsKnownDictionary);

            // read glycans if applicable
            List<Tuple<int, string, double>> localGlycans = null;
            if (spectrumMatch.GlycanLocalizationLevel != null)
            {
                localGlycans = PsmFromTsv.ReadLocalizedGlycan(spectrumMatch.LocalizedGlycan);
            }

            // annotate mods
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                double xLocation = (mod.Key - 1) * (spacer ?? MetaDrawSettings.AnnotatedSequenceTextSpacing) - xShift;
                double yLocation = yLoc + 2;

                if (mod.Value.ModificationType == "O-Glycosylation")
                {
                    if (localGlycans.Where(p => p.Item1 + 1 == mod.Key).Count() > 0)
                    {
                        DrawCircle(sequenceDrawingCanvas, new Point(xLocation, yLocation), ParseColorBrushFromOxyColor(MetaDrawSettings.GetValueOrDefault(MetaDrawSettings.ModificationTypeToColor, mod.Value.IdWithMotif, OxyColors.Blue)));
                    }
                    else
                    {
                        DrawCircle(sequenceDrawingCanvas, new Point(xLocation, yLocation), Brushes.Gray);
                    }
                }
                else
                {
                    DrawCircle(sequenceDrawingCanvas, new Point(xLocation, yLocation), ParseColorBrushFromOxyColor(MetaDrawSettings.ModificationTypeToColor[mod.Value.IdWithMotif]));
                }
            }
        }


        /// <summary>
        /// Redraws the Stationary Sequence on the spectrum in reference to the position of the scrollable sequence
        /// </summary>
        /// <param name="lettersOnScreen"></param>
        /// <param name="firstLetterOnScreen"></param>
        /// <param name="psm"></param>
        /// <param name="canvas"></param>
        public static void DrawStationarySequence(PsmFromTsv psm, DrawnSequence stationarySequence, int yLoc)
        {
            ClearCanvas(stationarySequence.SequenceDrawingCanvas); 
            string baseSequence = psm.BaseSeq.Substring(MetaDrawSettings.FirstAAonScreenIndex, MetaDrawSettings.NumberOfAAOnScreen);
            string fullSequence = baseSequence;

            // Trim full sequences selectively based upon what is show in scrollable sequence
            Dictionary<int, List<string>> modDictionary = PsmFromTsv.ParseModifications(psm.FullSequence);
            foreach (var mod in modDictionary.OrderByDescending(p => p.Key))
            {
                // if modification is within the visible region
                if (mod.Key >= MetaDrawSettings.FirstAAonScreenIndex && mod.Key < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen))
                {
                    // account for multiple modifications on the same amino acid
                    for (int i = mod.Value.Count - 1; i > -1; i--)
                    {
                        fullSequence = fullSequence.Insert(mod.Key - MetaDrawSettings.FirstAAonScreenIndex, "[" + mod.Value[i] + "]");
                        if (i >= 1)
                        {
                            fullSequence = fullSequence.Insert(mod.Key, "|");
                        }
                    }
                }
            }

            List<MatchedFragmentIon> matchedIons = psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition > MetaDrawSettings.FirstAAonScreenIndex &&
                                                   p.NeutralTheoreticalProduct.AminoAcidPosition < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen)).ToList();
            stationarySequence.AnnotateBaseSequence(baseSequence, fullSequence, yLoc, matchedIons, psm, true);
        }

        /// <summary>
        /// Draws the annotated sequence located below sequence coverage view taking into account the mutable display settings
        /// </summary>
        /// <param name="psm"></param>
        /// <param name="sequence"></param>
        public static void DrawSequenceAnnotation(PsmFromTsv psm, DrawnSequence sequence)
        {
            ClearCanvas(sequence.SequenceDrawingCanvas); 
            int segmentsPerRow = MetaDrawSettings.SequenceAnnotationSegmentPerRow;
            int residuesPerSegment = MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment;
            Dictionary<int, List<string>> modDictionary = PsmFromTsv.ParseModifications(psm.FullSequence);
            int numberOfRows = (int)Math.Ceiling(((double)psm.BaseSeq.Length / residuesPerSegment) / segmentsPerRow);
            int remaining = psm.BaseSeq.Length;

            sequence.SequenceDrawingCanvas.Height = 42 * numberOfRows + 10;
            // create an individual psm for each chunk to be drawn
            List<PsmFromTsv> segments = new();
            List<List<MatchedFragmentIon>> matchedIonSegments = new();
            for (int i = 0; i < psm.BaseSeq.Length; i += residuesPerSegment)
            {
                // split base seq
                string baseSequence;
                List<MatchedFragmentIon> ions = new();
                if (i + residuesPerSegment < psm.BaseSeq.Length)
                {
                    baseSequence = psm.BaseSeq.Substring(i, residuesPerSegment);
                    ions = psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition > i && p.NeutralTheoreticalProduct.AminoAcidPosition < (i + residuesPerSegment)).ToList();
                    ions.AddRange(psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition == i && p.NeutralTheoreticalProduct.Annotation.Contains('y')));
                    ions.AddRange(psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition == (i + residuesPerSegment) && p.NeutralTheoreticalProduct.Annotation.Contains('b')));
                    remaining -= residuesPerSegment;
                }
                else
                {
                    baseSequence = psm.BaseSeq.Substring(i, remaining);
                    ions = psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition > i && p.NeutralTheoreticalProduct.AminoAcidPosition < (i + remaining)).ToList();
                    ions.AddRange(psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition == i && p.NeutralTheoreticalProduct.Annotation.Contains('y')));
                    ions.AddRange(psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition == (i + residuesPerSegment) && p.NeutralTheoreticalProduct.Annotation.Contains('b')));
                    remaining -= remaining;
                }

                // add the mods onto the trimmed full sequence
                string fullSequence = baseSequence;
                foreach (var mod in modDictionary.OrderByDescending(p => p.Key))
                {
                    // account for multiple modifications on the same amino acid
                    for (int k = mod.Value.Count - 1; k > -1; k--)
                    {
                        // if modification is within the visible region
                        if (mod.Key >= i && mod.Key <= i + residuesPerSegment)
                        {
                            fullSequence = fullSequence.Insert(mod.Key - i, "[" + mod.Value[k] + "]");
                            if (k >= 1)
                            {
                                fullSequence = fullSequence.Insert(mod.Key - i, "|");
                            }
                        }
                    }
                }
                PsmFromTsv tempPsm = new(psm, fullSequence, baseSequence: baseSequence);
                segments.Add(tempPsm);
                matchedIonSegments.Add(ions);
            }

            // draw each resulting psm
            for (int i = 0; i < segments.Count; i++)
            {
                int startPosition = i * residuesPerSegment;
                int endPosition = (i + 1) * residuesPerSegment;
                List<MatchedFragmentIon> matchedIons = psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition > startPosition &&
                                                       p.NeutralTheoreticalProduct.AminoAcidPosition < endPosition).ToList();
                int currentRowZeroIndexed = i / segmentsPerRow;
                int yLoc = 10 + (currentRowZeroIndexed * 42);
                int chunkPositionInRow = (i % segmentsPerRow);

                sequence.AnnotateBaseSequence(segments[i].BaseSeq, segments[i].FullSequence, yLoc, matchedIonSegments[i], segments[i], false, currentRowZeroIndexed, chunkPositionInRow);
            }
        }


        public void DrawCrossLinkSequence()
        {
            this.AnnotateBaseSequence(SpectrumMatch.BetaPeptideBaseSequence, SpectrumMatch.BetaPeptideFullSequence, 100, SpectrumMatch.BetaPeptideMatchedIons, SpectrumMatch);
            // annotate crosslinker
            int alphaSite = int.Parse(Regex.Match(SpectrumMatch.FullSequence, @"\d+").Value);
            int betaSite = int.Parse(Regex.Match(SpectrumMatch.BetaPeptideFullSequence, @"\d+").Value);

            AnnotateCrosslinker(SequenceDrawingCanvas,
                new Point(alphaSite * MetaDrawSettings.AnnotatedSequenceTextSpacing, 50),
                new Point(betaSite * MetaDrawSettings.AnnotatedSequenceTextSpacing, 90),
                Colors.Black);
        }


        /// <summary>
        /// This method exists because of the mirror plotting of spectral libraries. Library spectral ions are displayed
        /// as having negative intensities for easy visualization, but obviously the ions do not actually have negative
        /// intensities. This formatter is used on the Y-axis (intensity) to turn negative values into positive ones
        /// so the Y-axis doesn't display negative intensities.
        /// </summary>
        public static string YAxisLabelFormatter(double d)
        {
            if (d < 0)
            {
                return (-d).ToString("0e-0", CultureInfo.InvariantCulture);
            }

            return d.ToString("0e-0", CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Clear canvas board
        /// </summary>
        public static void ClearCanvas(Canvas cav)
        {
            cav.Children.Clear();
        }

        public static SolidColorBrush ParseColorBrushFromOxyColor(OxyColor color)
        {
            var colorVal = color.ToByteString().Split(',');
            return new SolidColorBrush(System.Windows.Media.Color.FromArgb(Byte.Parse(colorVal[0]), Byte.Parse(colorVal[1]), Byte.Parse(colorVal[2]), Byte.Parse(colorVal[3])));
        }

        public static SolidColorBrush ParseColorBrushFromName(string name)
        {
            OxyColor color = MetaDrawSettings.PossibleColors.Keys.Where(p => p.GetColorName().Equals(name.Replace(" ", ""))).First();
            return ParseColorBrushFromOxyColor(color);
        }

        public static OxyColor ParseOxyColorFromName(string name)
        {
            return MetaDrawSettings.PossibleColors.Keys.Where(p => p.GetColorName().Equals(name.Replace(" ", ""))).First();
        }

        public static Color ParseColorFromOxyColor(OxyColor color)
        {
            var colorVal = color.ToByteString().Split(',');
            return System.Windows.Media.Color.FromArgb(Byte.Parse(colorVal[0]), Byte.Parse(colorVal[1]), Byte.Parse(colorVal[2]), Byte.Parse(colorVal[3]));
        }

        /// <summary>
        /// Draw the line seperator @ top
        /// </summary>
        private static void DrawCTermIon(Canvas cav, Point topLoc, Color clr, string footnote)
        {
            double x = topLoc.X, y = topLoc.Y;
            Polyline bot = new Polyline();
            bot.Points = new PointCollection() { new Point(x + 10, y + 10), new Point(x, y + 10), new Point(x, y + 24) };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 1;
            cav.Children.Add(bot);
            Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
        }

        /// <summary>
        /// Draw the line seperator @ bottom
        /// </summary>
        private static void DrawNTermIon(Canvas cav, Point botLoc, Color clr, string footnote)
        {
            double x = botLoc.X, y = botLoc.Y;
            Polyline bot = new Polyline();
            bot.Points = new PointCollection() { new Point(x - 10, y - 10), new Point(x, y - 10), new Point(x, y - 24) };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 1;
            Canvas.SetZIndex(bot, 1); //on top of any other things in canvas
            cav.Children.Add(bot);
        }

        /// <summary>
        /// Create text blocks on canvas
        /// </summary>
        private static void DrawText(Canvas cav, Point loc, string txt, Brush clr)
        {
            TextBlock tb = new TextBlock();
            tb.Foreground = clr;
            tb.Text = txt;
            tb.Height = 30;
            tb.FontSize = 25;
            tb.FontWeight = System.Windows.FontWeights.Bold;
            tb.FontFamily = new FontFamily("Arial");
            tb.TextAlignment = TextAlignment.Center;
            tb.HorizontalAlignment = System.Windows.HorizontalAlignment.Center;
            tb.Width = 24; // W (tryptophan) seems to be widest letter, make sure it fits if you're editing this

            Canvas.SetTop(tb, loc.Y);
            Canvas.SetLeft(tb, loc.X);
            Panel.SetZIndex(tb, 2); //lower priority
            cav.Children.Add(tb);
            cav.UpdateLayout();
        }

        /// <summary>
        /// Draws a circle
        /// </summary>
        private static void DrawCircle(Canvas cav, Point loc, SolidColorBrush clr)
        {
            Ellipse circle = new Ellipse()
            {
                Width = 24,
                Height = 24,
                Stroke = clr,
                StrokeThickness = 1,
                Fill = clr,
                Opacity = 0.7
            };
            Canvas.SetLeft(circle, loc.X);
            Canvas.SetTop(circle, loc.Y);
            Panel.SetZIndex(circle, 1);
            cav.Children.Add(circle);
        }

        private static void AnnotateCrosslinker(Canvas cav, Point alphaBotLoc, Point betaBotLoc, Color clr)
        {
            Polyline bot = new Polyline();
            double distance = (betaBotLoc.Y - alphaBotLoc.Y) / 2;
            bot.Points = new PointCollection() { alphaBotLoc, new Point(alphaBotLoc.X, alphaBotLoc.Y + distance), new Point(betaBotLoc.X, alphaBotLoc.Y + distance), betaBotLoc };
            bot.Stroke = new SolidColorBrush(clr);
            bot.StrokeThickness = 2;
            cav.Children.Add(bot);
        }
    }
}
