using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// View Model class for the ptm color legend
    /// </summary>
    public class PtmLegendViewModel : BaseViewModel
    {
        public string Header { get; set; } = "Legend";
        public int HeaderSize { get; set; } = 12;
        public ObservableCollection<PtmLegendItemViewModel> LegendItems { get; }

        public PtmLegendViewModel(List<Modification> mods)
        {
            LegendItems = new ObservableCollection<PtmLegendItemViewModel>();
            foreach (var mod in mods)
            {
                var modItem = new PtmLegendItemViewModel(mod.IdWithMotif);
                LegendItems.Add(modItem);
            }
        }
    }
}
