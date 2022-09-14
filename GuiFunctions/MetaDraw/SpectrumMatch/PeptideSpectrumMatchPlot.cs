using Chemistry;
using EngineLayer;
using iText.IO.Image;
using iText.Kernel.Pdf;
using iText.Layout;
using MassSpectrometry;
using mzPlot;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using iText.Kernel.Geom;
using Point = System.Windows.Point;
using Vector = System.Windows.Vector;
using Canvas = System.Windows.Controls.Canvas;

namespace GuiFunctions
{
    /// <summary>
    /// Class for the peptide spectrum match plot within the metadraw window
    /// </summary>
    public class PeptideSpectrumMatchPlot : SpectrumMatchPlot
    {
        public PeptideSpectrumMatchPlot(OxyPlot.Wpf.PlotView plotView, PsmFromTsv psm, MsDataScan scan,
            List<MatchedFragmentIon> matchedFragmentIons, bool annotateProperties = true, LibrarySpectrum librarySpectrum = null, bool stationarySequence = false)
            : base(plotView, psm, scan)
        {
            if (annotateProperties)
            {
                AnnotateProperties();
            }

            ZoomAxes(matchedFragmentIons);

            if (librarySpectrum != null)
            {
                AnnotateLibraryIons(isBetaPeptide: false, librarySpectrum.MatchedFragmentIons);
            }

            RefreshChart();
        }
    }

}