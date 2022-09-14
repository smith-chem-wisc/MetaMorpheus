using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using GuiFunctions.ViewModels.Legends;

namespace GuiFunctions
{
    public class ChimeraLegendViewModel : LegendViewModel
    {
        public ChimeraLegendItemViewModel SharedIons { get; set; }
        public ChimeraLegendViewModel(List<ChimeraLegendItemViewModel> legendItems, string ascession, double offset = 0) : base()
        {
            Header = ascession;
            var sharedIons = legendItems.First(p => p.Name == "Shared Ions");
            if (legendItems.Count > 2)
            {
                SharedIons = sharedIons;
                SharedIonStackPanelVisibility = Visibility.Visible;
            }
            else
            {
                SharedIonStackPanelVisibility = Visibility.Collapsed;
            }
            legendItems.Remove(sharedIons);
            legendItems.ForEach(p => LegendItemViewModels.Add(p));
            TopOffset = offset;
            Visibility = Visibility.Visible;
        }
    }
}
