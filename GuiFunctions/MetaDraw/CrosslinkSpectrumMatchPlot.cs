using EngineLayer;
using GuiFunctions.MetaDraw;
using MassSpectrometry;
using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;

namespace GuiFunctions
{
    /// <summary>
    /// Class for displaying CrossLinked spectra in MetaDraw
    /// </summary>
    public class CrosslinkSpectrumMatchPlot : PeptideSpectrumMatchPlot
    {
        public CrosslinkSpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, PsmFromTsv csm, MsDataScan scan, Canvas stationaryCanvas)
            : base(plotView, csm, scan, csm.MatchedIons)
        {
            // annotate beta peptide matched ions
            AnnotateMatchedIons(isBetaPeptide: true, csm.BetaPeptideMatchedIons);

            ZoomAxes(csm.MatchedIons.Concat(csm.BetaPeptideMatchedIons), yZoom: 1.5);
            RefreshChart();
        }
    }
}