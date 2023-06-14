using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Proteomics;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;

namespace EngineLayer
{
    public enum MultiplexLabel
    {
        TMT6,
        TMT10,
        TMT11,
        TMT18,
        iTRAQ4,
        iTRAQ8,
        DiLeu4,
        DiLeu12
    }
    
    public class MultiplexMod
    {
        public Modification Label { get; }
        public MultiplexLabel LabelType { get; }
        public MultiplexMod(Modification label, MultiplexLabel labelType)
        {
            LabelType = labelType;
            Label = label;
        }

        public List<string> MultiplexLabels()
        {
            List<string> ionLabels = new();
            var labelGroups = Label.DiagnosticIons.First().Value
                .Select(x => x.ToMz(1))
                .OrderBy(x => x)
                .GroupBy(x => (int)Math.Floor(x));

            if (Label.IdWithMotif.Contains("TMT"))
            {
                // TMT 126 contains no heavy isotopes. TMT 127N has one N15, TMT127C has one C15.
                // The "N" labels are slightly lighter than the "C" labels. 
                // Labels for the diagnostic ions are created accordingly
                foreach (var group in labelGroups)
                {
                    if (group.Count() == 1)
                    {
                        ionLabels.Add(group.Key.ToString());
                    }
                    else if (group.Count() == 2)
                    {
                        ionLabels.Add(group.Key + "N");
                        ionLabels.Add(group.Key + "C");
                    }
                }
            }
            else
            {
                foreach (var group in labelGroups)
                {
                    if (group.Count() == 1)
                    {
                        ionLabels.Add(group.Key.ToString());
                    }
                    else
                    {
                        ionLabels.AddRange(group.Select(mz => Math.Round(mz, 3).ToString(CultureInfo.CurrentCulture)));
                    }
                }
            }
            return ionLabels;
        }

        public double[] GetMultiplexIonIntensities(MzSpectrum scan, double ppmTolerance = 10)
        {
            Tolerance tolerance = new PpmTolerance(ppmTolerance);
            double[] theoreticalIonMzs = Label.DiagnosticIons.First().Value
                .Select(x => x.ToMz(1))
                .OrderBy(x => x)
                .ToArray();
            int peakIndex = scan.GetClosestPeakIndex(theoreticalIonMzs[0]);
            int lastPeakIndex = Math.Min(scan.GetClosestPeakIndex(theoreticalIonMzs.Last()) + 1, scan.XArray.Length - 1);
            double[] ionIntensities = new double[theoreticalIonMzs.Length];

            for (int ionIndex = 0; ionIndex < ionIntensities.Length; ionIndex++)
            {
                while (peakIndex <= lastPeakIndex &&
                       scan.XArray[peakIndex] < tolerance.GetMinimumValue(theoreticalIonMzs[ionIndex]))
                {
                    peakIndex++;
                }
                if (peakIndex > lastPeakIndex)
                {
                    break;
                }
                if (tolerance.Within(scan.XArray[peakIndex], theoreticalIonMzs[ionIndex]))
                {
                    ionIntensities[ionIndex] = scan.YArray[peakIndex];
                    peakIndex++;
                }
            }

            return ionIntensities;
        }
    }
}
