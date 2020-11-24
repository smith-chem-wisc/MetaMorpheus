using System.ComponentModel;
using System.Collections.ObjectModel;
using System.Windows.Controls;
using mzPlot;
using EngineLayer;

namespace EngineLayer
{
    public class ParentChildScanPlotsView
    {
        public ObservableCollection<ParentChildScanPlotTemplate> Plots { get; set; }
        
        public ParentChildScanPlotsView()
        {
            Plots = new ObservableCollection<ParentChildScanPlotTemplate>();
        }
    }

    public class ParentChildScanPlotTemplate
    {
        public PeptideSpectrumMatchPlot Plot { get; set; }

        public string SpectrumLabel { get; set; }

        public Canvas TheCanvas { get; set; }
    }
}