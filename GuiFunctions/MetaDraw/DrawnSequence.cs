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

namespace GuiFunctions.MetaDraw
{

    /// <summary>
    /// Class for drawing the sequences within metadraw, both scrollable and stationary
    /// </summary>
    public class DrawnSequence
    {
        public Canvas SequenceDrawingCanvas;
        public bool Stationary;
        public PsmFromTsv SpectrumMatch;
        public DrawnSequence(Canvas sequenceDrawingCanvas, PsmFromTsv psm, bool stationary)
        {
            SequenceDrawingCanvas = sequenceDrawingCanvas;
            SpectrumMatch = psm;
            Stationary = stationary;
            SequenceDrawingCanvas.Width = 600;
            SequenceDrawingCanvas.Height = 60;


            if (stationary)
            {
                DrawStationarySequence(psm, this);
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
        public void AnnotateBaseSequence(string baseSequence, string fullSequence, int yLoc, List<MatchedFragmentIon> matchedFragmentIons, PsmFromTsv psm, bool stationary = false)
        {
            if (psm.BetaPeptideBaseSequence == null || !psm.BetaPeptideBaseSequence.Equals(baseSequence))
            {
                ClearCanvas(SequenceDrawingCanvas);
            }
            double canvasWidth = SequenceDrawingCanvas.Width;
            int spacing = 12;

            // draw initial amino acid number
            if (stationary)
            {
                var startAA = (MetaDrawSettings.FirstAAonScreenIndex + 1).ToString().ToCharArray().Reverse().ToArray();
                double x = 22;

                Polyline line = new Polyline();
                line.Points = new PointCollection() { new Point(x, yLoc + 25), new Point(x, yLoc + 35) };
                line.Stroke = new SolidColorBrush(Colors.Black);
                line.StrokeThickness = 3;
                SequenceDrawingCanvas.Children.Add(line);
                for (int i = 0; i < startAA.Length; i++)
                {
                    x = MetaDrawSettings.FirstAAonScreenIndex < 10 ? (i + 1) * -spacing + 22 : i  * -spacing + 14;
                    DrawText(SequenceDrawingCanvas, new Point(x, yLoc + 35), startAA[i].ToString(), Brushes.Black);
                }
            }

            // draw base sequence
            for (int r = 0; r < baseSequence.Length; r++)
            {
                double x = r * MetaDrawSettings.AnnotatedSequenceTextSpacing + 10;
                DrawText(SequenceDrawingCanvas, new Point(x, yLoc), baseSequence[r].ToString(), Brushes.Black);
                
                canvasWidth = x;
            }
            SequenceDrawingCanvas.Width = Math.Max(SequenceDrawingCanvas.Width, canvasWidth) + 185; // this number is the width of the grayed out box

            // draw final amino acid number
            if (stationary)
            {
                var endAA = (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen).ToString();
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
                        int residue;
                        if (Stationary)
                        {
                            residue = ion.NeutralTheoreticalProduct.AminoAcidPosition - MetaDrawSettings.FirstAAonScreenIndex;
                        }
                        else
                        {
                            residue = ion.NeutralTheoreticalProduct.AminoAcidPosition;
                        }

                        string annotation = ion.NeutralTheoreticalProduct.ProductType + "" + ion.NeutralTheoreticalProduct.FragmentNumber;
                        OxyColor oxycolor = psm.VariantCrossingIons.Contains(ion) ?
                            MetaDrawSettings.VariantCrossColor : MetaDrawSettings.ProductTypeToColor[ion.NeutralTheoreticalProduct.ProductType];
                        Color color = Color.FromArgb(oxycolor.A, oxycolor.R, oxycolor.G, oxycolor.B);

                        if (ion.NeutralTheoreticalProduct.NeutralLoss != 0)
                        {
                            annotation += "-" + ion.NeutralTheoreticalProduct.NeutralLoss;
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
                        DrawCircle(sequenceDrawingCanvas, new Point(xLocation, yLocation), MetaDrawSettings.ModificationAnnotationColor);
                    }
                    else
                    {
                        DrawCircle(sequenceDrawingCanvas, new Point(xLocation, yLocation), Brushes.Gray);
                    }
                }
                else
                {
                    DrawCircle(sequenceDrawingCanvas, new Point(xLocation, yLocation), MetaDrawSettings.ModificationAnnotationColor);
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
        public static void DrawStationarySequence(PsmFromTsv psm, DrawnSequence stationarySequence)
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
            stationarySequence.AnnotateBaseSequence(baseSequence, fullSequence, 10, matchedIons, psm, true);
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
