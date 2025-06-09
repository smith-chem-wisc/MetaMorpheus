using System.Collections.ObjectModel;

namespace GuiFunctions
{
    public class LegendViewModel : BaseViewModel
    {
        private bool visibility;
        private double topOffset;
        private ObservableCollection<LegendItemViewModel> legendItemViewModels = new();
        public string Header { get; set; } = "Legend";
        public int HeaderSize { get; set; } = 12;


        public bool Visibility
        {
            get { return visibility; }
            set
            {
                visibility = value;
                OnPropertyChanged(nameof(Visibility));
            }
        }

        public double TopOffset
        {
            get => topOffset;
            set { topOffset = value; OnPropertyChanged(nameof(TopOffset)); }
        }

        public ObservableCollection<LegendItemViewModel> LegendItemViewModels
        {
            get => legendItemViewModels;
            set => value = LegendItemViewModels;
        }

        public LegendViewModel()
        {
            LegendItemViewModels = new();
            Visibility = MetaDrawSettings.ShowLegend;
        }
    }
}
