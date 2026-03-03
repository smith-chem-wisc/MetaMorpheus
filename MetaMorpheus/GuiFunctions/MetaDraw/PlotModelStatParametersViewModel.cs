using System;
using System.Collections.ObjectModel;

namespace GuiFunctions.MetaDraw
{
    /// <summary>
    /// Singleton ViewModel for Data Visualization plot parameters
    /// </summary>
    public class PlotModelStatParametersViewModel : BaseViewModel
    {
        private static readonly Lazy<PlotModelStatParametersViewModel> _instance =
            new(() => new PlotModelStatParametersViewModel());

        public static PlotModelStatParametersViewModel Instance => _instance.Value;

        private PlotModelStatParameters _current;

        private PlotModelStatParametersViewModel()
        {
            _current = new PlotModelStatParameters();
        }

        public ObservableCollection<string> GroupingProperties { get; } = new(MetaDrawSettings.GroupingProperties);

        public bool UseLogScaleYAxis
        {
            get => _current.UseLogScaleYAxis;
            set
            {
                _current.UseLogScaleYAxis = value;
                OnPropertyChanged(nameof(UseLogScaleYAxis));
            }
        }

        public string GroupingProperty
        {
            get => _current.GroupingProperty;
            set
            {
                _current.GroupingProperty = value ?? "None";
                OnPropertyChanged(nameof(GroupingProperty));
            }
        }

        public double MinRelativeCutoff
        {
            get => _current.MinRelativeCutoff;
            set
            {
                if (value < 0) value = 0;
                if (value > 100) value = 100;
                if (value > MaxRelativeCutoff) value = MaxRelativeCutoff;
                _current.MinRelativeCutoff = value;
                OnPropertyChanged(nameof(MinRelativeCutoff));
            }
        }

        public double MaxRelativeCutoff
        {
            get => _current.MaxRelativeCutoff;
            set
            {
                if (value < 0) value = 0;
                if (value > 100) value = 100;
                if (value < MinRelativeCutoff) value = MinRelativeCutoff;
                _current.MaxRelativeCutoff = value;
                OnPropertyChanged(nameof(MaxRelativeCutoff));
            }
        }

        public bool NormalizeHistogramToFile
        {
            get => _current.NormalizeHistogramToFile;
            set
            {
                _current.NormalizeHistogramToFile = value;
                OnPropertyChanged(nameof(NormalizeHistogramToFile));
            }
        }

        public bool DisplayFilteredOnly
        {
            get => _current.DisplayFilteredOnly;
            set
            {
                _current.DisplayFilteredOnly = value;
                OnPropertyChanged(nameof(DisplayFilteredOnly));
            }
        }

        /// <summary>
        /// Gets a clone of the current parameters for passing to PlotModelStat
        /// </summary>
        public PlotModelStatParameters GetParameters()
        {
            return _current.Clone();
        }

        /// <summary>
        /// Loads parameters from a snapshot
        /// </summary>
        public void LoadFromSnapshot(PlotModelStatParameters parameters)
        {
            _current = parameters.Clone();
            OnPropertyChanged(nameof(UseLogScaleYAxis));
            OnPropertyChanged(nameof(GroupingProperty));
            OnPropertyChanged(nameof(MinRelativeCutoff));
            OnPropertyChanged(nameof(MaxRelativeCutoff));
            OnPropertyChanged(nameof(NormalizeHistogramToFile));
            OnPropertyChanged(nameof(DisplayFilteredOnly));
        }
    }
}
