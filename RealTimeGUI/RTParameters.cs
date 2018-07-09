using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RealTimeGUI
{
    public class RTParameters
    {
        public RTParameters()
        {
            TimeScale = 10000;
        }
        public int TimeScale { get; set; }
    }
}
