using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;

namespace EngineLayer.DIA
{
    public class PeakCurve
    {
        public List<(double,double,double)> PeakList { get; set; }
        public double RetentionTime { get; set; }
        public double PrecursorMz { get; set; }
        public int ChargeState { get; set; }
        public virtual List<MatchedFragmentIon> MatchedFragmentIons { get; set; }
        public bool IsDecoy { get; set; }

        public double[] XArray { get; private set; }
        public double[] YArray { get; private set; }


        public PeakCurve()
        {
            PeakList = new List<(double, double, double)>();
            //PrecursorMz = precursorMz;
            //MatchedFragmentIons = peaks;
            //ChargeState = chargeState;
            //IsDecoy = isDecoy;
            //RetentionTime = rt;
            //XArray = peaks.Select(p => p.Mz).ToArray();
            //YArray = peaks.Select(p => p.Intensity).ToArray();
            //Array.Sort(XArray, YArray);
        }
         
        public class XYZDataList
        {
            public XYZDataList() { }
            List<(double, double, double)> xyzdata = new List<(double, double, double)>();
            public void addXYZ((double, double, double) xyz)
            {
                xyzdata.Add(xyz);
            }

        }

        public void AddPeak(double x, double y, double z)
        {
            (double, double, double) newPeak = (x, y, z);
            PeakList.Add(newPeak);
            //TotalIntMzF += y * z * z;
            //TotalIntF += z * z;
            //if (z > ApexInt)
            //{
            //    ApexInt = z;
            //    ApexRT = x;
            //}
            //if (z < minIntF)
            //{
            //    minIntF = z;
            //}
            //TargetMz = TotalIntMzF / TotalIntF;
        }
    }
}
