using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public static class SpectralLibrarySearchFunction
    {
        public static void CalculateSpectralAngles(SpectralLibrary spectralLibrary, PeptideSpectralMatch[] peptideSpectralMatches,
            Ms2ScanWithSpecificMass[] arrayOfSortedMs2Scans, CommonParameters commonParameters)
        {
            foreach (PeptideSpectralMatch psm in peptideSpectralMatches.Where(p => p != null))
            {
                Ms2ScanWithSpecificMass scan = arrayOfSortedMs2Scans[psm.ScanIndex];

                //TODO: spectral angle could be used to disambiguate PSMs. right now for ambiguous PSMs, the spectral angle for only one peptide option is saved
                foreach (var peptide in psm.PeptidesToMatchingFragments)
                {
                    if (spectralLibrary == null || !spectralLibrary.TryGetSpectrum(peptide.Key.FullSequence, scan.PrecursorCharge, out var librarySpectrum))
                    {
                        continue;
                    }

                    double spectralAngle = CalculateNormalizedSpectralAngle(librarySpectrum.MatchedFragmentIons, scan.TheScan, commonParameters);

                    psm.SpectralAngle = spectralAngle;
                }
            }
        }

        /// <summary>
        /// Calculates the spectral angle, as described by Prosit ( https://www.nature.com/articles/s41592-019-0426-7 ).
        /// </summary>
        public static double CalculateNormalizedSpectralAngle(List<MatchedFragmentIon> theoreticalLibraryIons, MsDataScan scan, CommonParameters commonParameters)
        {
            double mzCutoff = 300;
            int fragmentNumberCutoff = 3;

            // if the spectrum has no peaks
            if (scan.MassSpectrum.XArray.Length == 0)
            {
                return 0;
            }

            Dictionary<MatchedFragmentIon, MatchedFragmentIon> matchedIons = new Dictionary<MatchedFragmentIon, MatchedFragmentIon>();

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
                double experimentalIntensity = scan.MassSpectrum.YArray[closestPeakIndex];

                // is the mass error acceptable?
                if (commonParameters.ProductMassTolerance.Within(mz.ToMass(libraryIon.Charge), libraryIon.Mz.ToMass(libraryIon.Charge)))
                {
                    var test = new Product(libraryIon.NeutralTheoreticalProduct.ProductType, libraryIon.NeutralTheoreticalProduct.Terminus,
                        libraryIon.NeutralTheoreticalProduct.NeutralMass, libraryIon.NeutralTheoreticalProduct.FragmentNumber,
                        libraryIon.NeutralTheoreticalProduct.AminoAcidPosition, libraryIon.NeutralTheoreticalProduct.NeutralLoss);

                    matchedIons.Add(libraryIon, new MatchedFragmentIon(ref test, mz, experimentalIntensity, libraryIon.Charge));
                }
            }

            // L2 norm
            double expNormalizer = Math.Sqrt(matchedIons.Sum(p => Math.Pow(p.Value.Intensity, 2)));
            double theorNormalizer = Math.Sqrt(theoreticalLibraryIons.Sum(p => Math.Pow(p.Intensity, 2)));

            double dotProduct = 0;

            foreach (var libraryIon in theoreticalLibraryIons)
            {
                if (matchedIons.TryGetValue(libraryIon, out var experIon))
                {
                    dotProduct += (libraryIon.Intensity / theorNormalizer) * (experIon.Intensity / expNormalizer);
                }
            }

            double normalizedSpectralAngle = 1 - (2 * Math.Acos(dotProduct) / Math.PI);

            return normalizedSpectralAngle;
        }
    }
}
