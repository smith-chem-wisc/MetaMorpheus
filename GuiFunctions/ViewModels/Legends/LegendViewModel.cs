using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace GuiFunctions
{
    public class LegendViewModel : BaseViewModel
    {
        
        private Visibility visibility;
        private double topOffset;
        private ObservableCollection<LegendItemViewModel> legendItemViewModels = new();
        public string Header { get; set; } = "Legend";
        public int HeaderSize { get; set; } = 12;


        public Visibility Visibility
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
        }

        public LegendViewModel()
        {

        }
    }
}
