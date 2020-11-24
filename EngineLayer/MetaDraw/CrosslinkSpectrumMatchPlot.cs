using MassSpectrometry;
using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace EngineLayer
{
    public class CrosslinkSpectrumMatchPlot : PeptideSpectrumMatchPlot
    {
        public CrosslinkSpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, Canvas sequenceDrawingCanvas, PsmFromTsv csm, MsDataScan scan) 
            : base(plotView, sequenceDrawingCanvas, csm, scan, csm.MatchedIons)
        {
            SequenceDrawingCanvas.Height = 150;

            // annotate beta peptide base sequence
            AnnotateBaseSequence(csm.BetaPeptideBaseSequence, csm.BetaPeptideFullSequence, 100, csm.BetaPeptideMatchedIons);

            // annotate beta peptide matched ions
            AnnotateMatchedIons(isBetaPeptide: true, csm.BetaPeptideMatchedIons);

            // annotate crosslinker
            int alphaSite = int.Parse(Regex.Match(SpectrumMatch.FullSequence, @"\d+").Value);
            int betaSite = int.Parse(Regex.Match(SpectrumMatch.BetaPeptideFullSequence, @"\d+").Value);

            AnnotateCrosslinker(SequenceDrawingCanvas, 
                new Point(alphaSite * MetaDrawSettings.AnnotatedSequenceTextSpacing, 50), 
                new Point(betaSite * MetaDrawSettings.AnnotatedSequenceTextSpacing, 90), 
                Colors.Black);

            ZoomAxes(csm.MatchedIons.Concat(csm.BetaPeptideMatchedIons), yZoom: 1.5);
            RefreshChart();
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