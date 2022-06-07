using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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

        public PtmLegendView(List<Modification> mods)
        {
            LegendItems = new ObservableCollection<PtmLegendItemView>();
            foreach (var mod in mods)
            {
                var modItem = new PtmLegendItemView(mod.IdWithMotif);
                LegendItems.Add(modItem);
            }
        }
    }
}
