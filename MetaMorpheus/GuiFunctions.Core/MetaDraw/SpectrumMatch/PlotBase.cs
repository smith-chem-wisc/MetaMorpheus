using System.ComponentModel;
using OxyPlot;

namespace GuiFunctions
{
    /// <summary>
    /// The base every MetaDraw plot derives from: an OxyPlot model plus change notification.
    ///
    /// Replaces mzPlot.Plot, which is the same idea but lives in an assembly referencing
    /// PresentationFramework and OxyPlot.Wpf - so inheriting from it made the whole
    /// SpectrumMatchPlot hierarchy Windows-bound even though the hierarchy's own source is
    /// WPF-free. Only Model, RefreshChart and ExportToPng were ever used from it; ExportToPng
    /// is WPF rendering and stays with the WPF layer.
    /// </summary>
    public abstract class PlotBase : INotifyPropertyChanged
    {
        private PlotModel _model = new();

        /// <summary>The OxyPlot model for the chart.</summary>
        public PlotModel Model
        {
            get => _model;
            set
            {
                _model = value;
                RefreshChart();
            }
        }

        public void RefreshChart()
        {
            Model.InvalidatePlot(true);
            NotifyPropertyChanged(nameof(Model));
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected void NotifyPropertyChanged(string propertyName)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}
