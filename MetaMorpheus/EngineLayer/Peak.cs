using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class Peak
    {
        public Peak(double mz, double rt, double intensity, int scanNumber=0,int index=0)
        {
            Mz = mz;
            Intensity = intensity;
            RT = rt;
            SanNumber = scanNumber;
            Index = index;
        }

        public double Mz { get; set; }
        public double Intensity { get; set; }
        public double RT { get; set; }
        public int SanNumber { get; set; }
        public int Index { get; set; }
    }
}
