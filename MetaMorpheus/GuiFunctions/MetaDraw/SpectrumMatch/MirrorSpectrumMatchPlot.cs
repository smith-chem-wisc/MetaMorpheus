using System.Collections.Generic;
using MassSpectrometry;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using OxyPlot.Wpf;
using Readers;

namespace GuiFunctions
{
    /// <summary>
    /// Creates a mirror plot comparing two different spectrum matches (PSMs) from two different scans.
    /// The primary PSM's spectrum and matched fragment ions are displayed above the axis (positive Y),
    /// while the secondary PSM's spectrum and matched fragment ions are mirrored below the axis (negative Y).
    /// Useful for visually comparing RNA or peptide sequences to highlight differentiating fragment ions.
    /// </summary>
    public class MirrorSpectrumMatchPlot : SpectrumMatchPlot
    {
        public MirrorSpectrumMatchPlot(
            PlotView plotView,
            SpectrumMatchFromTsv psmA,
            MsDataScan scanA,
            List<MatchedFragmentIon> ionsA,
            SpectrumMatchFromTsv psmB,
            MsDataScan scanB,
            List<MatchedFragmentIon> ionsB,
            bool annotateProperties = true)
            : base(plotView, psmA, scanA, ionsA)
        {
            DrawInvertedSpectrum(scanB);

            AnnotateMirrorIons(isBetaPeptide: false, ionsB);

            AdjustMirrorAxes(scanA, scanB);

            if (annotateProperties)
            {
                AnnotateProperties();
            }

            RefreshChart();
        }
    }
}
