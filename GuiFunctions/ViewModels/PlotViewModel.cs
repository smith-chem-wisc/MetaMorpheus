using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using mzPlot;

namespace GuiFunctions
{
    public class PlotViewModel : BaseViewModel
    {
        private Plot plot;

        public Plot Plot
        {
            get { return plot; }
            set { plot = value; OnPropertyChanged(nameof(plot)); plot.RefreshChart();  }
        }
    }
}
