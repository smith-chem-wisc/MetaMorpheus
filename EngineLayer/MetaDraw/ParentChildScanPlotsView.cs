using System.ComponentModel;
using System.Collections.ObjectModel;
using System.Windows.Controls;
using mzPlot;
using EngineLayer;

namespace EngineLayer
{
    public class ParentChildScanPlotsView : INotifyPropertyChanged
    {
        public ObservableCollection<ParentChildScanPlotTemplate> Plots { get; set; }
        private int _myColumnCount = 1;
        public event PropertyChangedEventHandler PropertyChanged = (s, e) => { };
        public int MyColumnCount
        {
            get { return _myColumnCount; }
            set
            {
                _myColumnCount = value;
                this.PropertyChanged(this, new PropertyChangedEventArgs(nameof(MyColumnCount)));
            }
        }

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