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
    public class MultiplexMod
    {
        public Modification Label { get; }
        public List<string> DiagnosticIonLabels { get; }
        public double[] DiagnosticIonMzs { get; }

        /// <summary>
        /// Functions for finding and writing the diagnostic ions associated with a specific multiplex label.
        /// </summary>
        /// <param name="label">Modification of type "Multiplex Label"</param>
        public MultiplexMod(Modification label)
        {
            if(label == null || !label.ModificationType.Equals("Multiplex Label"))
            {
                throw new ArgumentException(
                    "MultiplexMod can only be initialized with a Multiplex Label type modification");
            }
            Label = label;
            DiagnosticIonLabels = GetMultiplexLabels(label);
            DiagnosticIonMzs = Label.DiagnosticIons.First().Value
                .Select(x => x.ToMz(1))
                .OrderBy(x => x)
                .ToArray();
        }
        
        public static List<string> GetMultiplexLabels(Modification label)
        {
            List<string> ionLabels = new();
            var labelGroups = label.DiagnosticIons.First().Value
                .Select(x => x.ToMz(1))
                .OrderBy(x => x)
                .GroupBy(x => (int)Math.Floor(x));

            if (label.IdWithMotif.Contains("TMT"))
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
            int peakIndex = scan.GetClosestPeakIndex(DiagnosticIonMzs[0]);
            int lastPeakIndex = Math.Min(scan.GetClosestPeakIndex(DiagnosticIonMzs.Last()) + 1, scan.XArray.Length - 1);
            double[] ionIntensities = new double[DiagnosticIonMzs.Length];

            for (int ionIndex = 0; ionIndex < ionIntensities.Length; ionIndex++)
            {
                while (peakIndex <= lastPeakIndex &&
                       scan.XArray[peakIndex] < tolerance.GetMinimumValue(DiagnosticIonMzs[ionIndex]))
                {
                    peakIndex++;
                }
                if (peakIndex > lastPeakIndex)
                {
                    break;
                }
                if (tolerance.Within(scan.XArray[peakIndex], DiagnosticIonMzs[ionIndex]))
                {
                    ionIntensities[ionIndex] = scan.YArray[peakIndex];
                    peakIndex++;
                }
            }

            return ionIntensities;
        }
    }
}
