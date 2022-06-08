using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace GuiFunctions
{
    /// <summary>
    /// View Model class for the ptm color legend
    /// </summary>
    public class PtmLegendView : BaseView
    {
        public string Header { get; set; } = "Legend";
        public int HeaderSize { get; set; } = 12;
        public ObservableCollection<PtmLegendItemView> LegendItems { get; }
        public Visibility Visibility { get; set; }

        public PtmLegendView(List<Modification> mods)
        {
            LegendItems = new ObservableCollection<PtmLegendItemView>();
            foreach (var mod in mods.Distinct())
            {
                var modItem = new PtmLegendItemView(mod.IdWithMotif);
                LegendItems.Add(modItem);
            }

            Visibility = mods.Count > 0 ? Visibility.Visible : Visibility.Hidden;
        }
    }
}
