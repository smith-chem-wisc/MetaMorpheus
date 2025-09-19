using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI.ForDisplayingInDataGrids
{
    public class SelectedPrecursorForDataGrid 
    {
        public SelectedPrecursorForDataGrid(double mz, int charge, double rtStartInMinutes, double rtEndInMinutes)
        {
            Mz = mz;
            Charge = charge;
            RtStartInMinutes = rtStartInMinutes;
            RtEndInMinutes = rtEndInMinutes;
        }

        public SelectedPrecursorForDataGrid()
        {
            Mz = 500;
            Charge = 1;
            RtStartInMinutes = 0;
            RtEndInMinutes = 120;
        }

        public double Mz { get; set; }
        public int Charge { get; set; }
        public double RtStartInMinutes { get; set; }
        public double RtEndInMinutes { get; set; }
        public override string ToString()
        {
            return $"{Mz}\t{Charge}\t{RtStartInMinutes}\t{RtEndInMinutes}";
        }
    }
}
