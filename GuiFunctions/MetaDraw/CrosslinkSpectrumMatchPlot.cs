using EngineLayer;
using MassSpectrometry;
using System.Linq;
using System.Windows.Controls;
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