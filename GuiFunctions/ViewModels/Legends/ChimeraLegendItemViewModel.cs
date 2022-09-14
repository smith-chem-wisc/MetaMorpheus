using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;
using OxyPlot;

namespace GuiFunctions.ViewModels.Legends
{
    public class ChimeraLegendItemViewModel : LegendItemViewModel
    {
        #region Constructor

        public ChimeraLegendItemViewModel(string mods, OxyColor color)
        {
            if (mods == null || mods == "")
                Name = "No Modifications";
            else
                Name = mods;
            ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);
        }

        #endregion

    }
}
