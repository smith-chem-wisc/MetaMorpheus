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
        public double Mz { get; private set; }
        public int Charge { get; private set; }
        public double RtStartInMinutes { get; private set; }
        public double RtEndInMinutes { get; private set; }
        public override string ToString()
        {
            return $"{Mz}\t{Charge}\t{RtStartInMinutes}\t{RtEndInMinutes}";
        }
    }
}
