using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using MzLibUtil;
using Newtonsoft.Json.Linq;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;


namespace EngineLayer
{
    public class PeakCurve
    {
        public PeakCurve(double average_mz, List<Peak> peaks,int msLevel, MzRange isolationWindow =null)
        {
            Average_Mz = average_mz;
            Peaks = peaks;
            MsLevel = msLevel;
            IsolationWindow = isolationWindow;
        }

        public double Average_Mz { get; set; }
        public List<Peak> Peaks { get; set; }
        public int MsLevel { get; set; }
        public MzRange IsolationWindow { get; set; }


        public static double CalPeakCorr(PeakCurve peakCurve1, PeakCurve peakCurve2)
        {
            var matchedPeaks = RTPairs(peakCurve1,peakCurve2);
            var dataSet1 = matchedPeaks.Select(p => p.Item1).ToArray();
            var dataSet2 = matchedPeaks.Select(p => p.Item2).ToArray();
            if (matchedPeaks.Count>=5)
            {
                // Calculate Pearson correlation
                double correlation = CalculatePearsonCorrelation(dataSet1, dataSet2);
                return correlation;
            }
            else
            {
                return double.NaN;
            }
        }
        public static double CalculatePearsonCorrelation(double[] x, double[] y)
        {
            int n = x.Length;
            double sumX = 0, sumY = 0, sumXy = 0, sumX2 = 0, sumY2 = 0;
            for (int i = 0; i < n; i++)
            {
                sumX += x[i];
                sumY += y[i];
                sumXy += x[i] * y[i];
                sumX2 += x[i] * x[i];
                sumY2 += y[i] * y[i];
            }
            double numerator = n * sumXy - sumX * sumY;
            double denominator = Math.Sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));
            
            return numerator / denominator;
        }

        private static List<(double, double)> RTPairs(PeakCurve peakCurve1, PeakCurve peakCurve2)
        {
            double[] RT_1 = peakCurve1.Peaks.Select(p=>p.RT).ToArray();
            double[] Intensity_1 = peakCurve1.Peaks.Select(p => p.Intensity).ToArray();
            double[] RT_2 = peakCurve2.Peaks.Select(p => p.RT).ToArray();
            double[] Intensity_2 = peakCurve2.Peaks.Select(p => p.Intensity).ToArray();     
            
            if (RT_1 == null || RT_2 == null)
            {
                return new List<(double, double)> { (-1.0, -1.0) };
            }

            List<(double, double)> list = new List<(double, double)>();
            List<(double, double)> list2 = new List<(double, double)>();
            List<(double, double)> list3 = new List<(double, double)>();

            for (int j = 0; j < RT_1.Length; j++)
            {
                list2.Add((RT_1[j], Intensity_1[j]));
            }

            for (int k = 0; k < RT_2.Length; k++)
            {
                list3.Add((RT_2[k], Intensity_2[k]));
            }

            list2 = list2.OrderByDescending(((double, double) i) => i.Item2).ToList();
            list3 = list3.OrderByDescending(((double, double) i) => i.Item2).ToList();
            foreach (var item in list3)
            {
                int num = 0;
                while (list2.Count > 0 && num < list2.Count)
                {
                    if (Math.Abs(list2[num].Item1-item.Item1)/( list2[num].Item1+ item.Item1 )< 0.01)
                    {
                        list.Add((list2[num].Item2, item.Item2));
                        list2.RemoveAt(num);
                        num = -1;
                        break;
                    }

                    num++;
                }

                if (list2.Count == 0)
                {
                    num++;
                }

                if (num > 0)
                {
                    list.Add((0.0, item.Item2));
                }
            }

            return list;
        }
    }
}
