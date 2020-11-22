using MassSpectrometry;
using MzLibUtil;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public static class SpectralLibrarySearchFunction
    {
        public static double CalculateNormalizedSpectralAngle(List<MatchedFragmentIon> theoreticalLibraryIons, MsDataScan scan, CommonParameters commonParameters)
        {
            double mzCutoff = 300;
            int fragmentNumberCutoff = 3;

            // if the spectrum has no peaks
            if (scan.MassSpectrum.XArray.Length == 0)
            {
                return 0;
            }

            // L2 norm
            double normalizer = Math.Sqrt(scan.MassSpectrum.YArray.Sum(p => Math.Pow(p, 2)));

            var normalizedIntensities = scan.MassSpectrum.YArray.Select(p => Math.Sqrt(p) / normalizer).ToArray();

            Dictionary<int, MatchedFragmentIon> matchedIndices = new Dictionary<int, MatchedFragmentIon>();

            // search for each theoretical ion
            for (int i = 0; i < theoreticalLibraryIons.Count; i++)
            {
                var libraryIon = theoreticalLibraryIons[i];

                // see https://www.nature.com/articles/s41592-019-0426-7
                // "All non-zero fragment ions (m/z > 300, ion >3, no neutral loss fragment ions) were considered for spectral angle calculation"
                if (libraryIon.Mz <= mzCutoff || libraryIon.NeutralTheoreticalProduct.FragmentNumber <= fragmentNumberCutoff)
                {
                    continue;
                }

                // get the closest peak in the spectrum to the library peak
                var closestPeakIndex = scan.MassSpectrum.GetClosestPeakIndex(libraryIon.Mz);
                double mz = scan.MassSpectrum.XArray[closestPeakIndex];

                // is the mass error acceptable?
                if (commonParameters.ProductMassTolerance.Within(mz, libraryIon.Mz))
                {
                    matchedIndices.TryAdd(i, libraryIon);
                }
            }

            double numerator = 0;
            double denominatorA = 0;
            double denominatorB = 0;

            for (int i = 0; i < normalizedIntensities.Length; i++)
            {
                double mz = scan.MassSpectrum.XArray[i];
                if (mz <= mzCutoff)
                {
                    continue;
                }

                matchedIndices.TryGetValue(i, out var matchedIon);

                double libraryIntensity = matchedIon == null ? 0 : matchedIon.Intensity;
                double normalizedExperimentalIntensity = normalizedIntensities[i];

                numerator += libraryIntensity * normalizedExperimentalIntensity;

                denominatorA += Math.Pow(libraryIntensity, 2);
                denominatorB += Math.Pow(normalizedExperimentalIntensity, 2);
            }

            double denominator = Math.Sqrt(denominatorA) * Math.Sqrt(denominatorB);
            double dotProduct = denominator > 0 ? numerator / denominator : 0;

            double normalizedSpectralAngle = 1 - (2 * Math.Acos(dotProduct) / Math.PI);

            return normalizedSpectralAngle;
        }
    }
}
