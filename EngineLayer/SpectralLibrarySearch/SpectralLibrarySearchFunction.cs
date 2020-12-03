using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
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
                    if (spectralLibrary != null && peptide.Key.Protein.IsDecoy && !spectralLibrary.ContainsSpectrum(peptide.Key.FullSequence, scan.PrecursorCharge))
                    {
                        double spectralAngle = CalculateDecoyNormalizedSpectralAngle(peptide.Key, scan, commonParameters);
                        psm.SpectralAngle = spectralAngle;
                    }
                    else if (spectralLibrary != null && spectralLibrary.TryGetSpectrum(peptide.Key.FullSequence, scan.PrecursorCharge, out var librarySpectrum))
                    {
                        double spectralAngle = CalculateNormalizedSpectralAngle(librarySpectrum.MatchedFragmentIons, scan.TheScan, commonParameters);
                        psm.SpectralAngle = spectralAngle;
                    }
  
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

        //1024TestDoneCompareFunction
        public static double MatchedSpectraCompare(List<MatchedFragmentIon> standardSpectra, List<MatchedFragmentIon> spectraToCompare)
        {

            double[] mz1 = standardSpectra.Select(b => b.Mz).ToArray();
            double intensitySum1 = standardSpectra.Select(b => b.Intensity).Sum();
            double[] intensity1 = standardSpectra.Select(b => Math.Sqrt(b.Intensity / intensitySum1)).ToArray();
            //Console.WriteLine(mz1.Length + "  " + intensity1.Length);
            Array.Sort(mz1, intensity1);

            double[] mz2 = spectraToCompare.Select(b => b.Mz).ToArray();
            double intensitySum2 = spectraToCompare.Select(b => b.Intensity).Sum();
            double[] intensity2 = spectraToCompare.Select(b => Math.Sqrt(b.Intensity / intensitySum2)).ToArray();
            Array.Sort(mz2, intensity2);
            //Console.WriteLine(mz2.Length + "  " + intensity2.Length);

            var commonNumbers = mz1.Union(mz2).ToArray();
            double min = commonNumbers.Min();
            double max = commonNumbers.Max();
            int roundMin = (int)min;
            int roundMax = (int)max + 1;
            //Console.WriteLine(roundMin + "  " + roundMax);

            //convert spectra to vectors
            List<double> vector1 = new List<double>();
            List<double> vector2 = new List<double>();

            int i = 0; //iterate through mz1
            int k = 0; //iterate through bin
            double oneMz = mz1[0];
            double oneIntensity = intensity1[0];
            //find where peaks match
            while (roundMin + k * 0.5 < roundMax)
            {
                List<double> x1 = new List<double>();
                while (i < mz1.Length && roundMin + k * 0.5 <= oneMz && oneMz < roundMin + k * 0.5 + 0.5)
                {
                    x1.Add(oneIntensity);
                    i++;
                    if (i != mz1.Length)
                    {
                        oneMz = mz1[i];
                        oneIntensity = intensity1[i];
                    }
                }
                vector1.Add(x1.Sum());
                k++;
            }

            int j = 0; //iterate through mz2
            int n = 0; //iterate through bin
            double twoMz = mz2[0];
            double twoIntensity = intensity2[0];
            while (roundMin + n * 0.5 < roundMax)
            {
                List<double> x2 = new List<double>();
                while (j < mz2.Length && roundMin + n * 0.5 <= twoMz && twoMz < roundMin + n * 0.5 + 0.5)
                {
                    x2.Add(twoIntensity);
                    j++;
                    if (j != mz2.Length)
                    {
                        twoMz = mz2[j];
                        twoIntensity = intensity2[j];
                    }
                }
                vector2.Add(x2.Sum());
                n++;
            }

            //numerator of dot product
            double numerator = 0;
            for (i = 0; i < vector1.Count; i++)
            {
                numerator += vector1[i] * vector2[i];
            }

            //denominator of dot product
            double denominator = Math.Sqrt(vector1.Sum(x => x * x)) * Math.Sqrt(vector2.Sum(x => x * x));

            var score = numerator / denominator;
            return score;
        }

        public static List<MatchedFragmentIon> AverageTwoSpectra(List<MatchedFragmentIon> spectraOne, List<MatchedFragmentIon> spectraTwo)
        {
            Dictionary<String, MatchedFragmentIon> averagedPeaksDictionary = new Dictionary<String, MatchedFragmentIon>();
            var averageTwoSpectraResult = new List<MatchedFragmentIon>();
            Dictionary<String, MatchedFragmentIon> spectraOneDictionary = new Dictionary<String, MatchedFragmentIon>();
            Dictionary<String, MatchedFragmentIon> spectraTwoDictionary = new Dictionary<String, MatchedFragmentIon>();
            double intensitySum1 = spectraOne.Select(b => b.Intensity).Sum();
            double intensitySum2 = spectraTwo.Select(b => b.Intensity).Sum();
            foreach (var ionOne in spectraOne)
            {
                var productWithChargeInOne = ionOne.NeutralTheoreticalProduct.ProductType.ToString() + ionOne.NeutralTheoreticalProduct.FragmentNumber + "^" + ionOne.Charge;
                spectraOneDictionary.Add(productWithChargeInOne, ionOne);
            }
            foreach (var ionTwo in spectraTwo)
            {
                var productWithChargeInTwo = ionTwo.NeutralTheoreticalProduct.ProductType.ToString() + ionTwo.NeutralTheoreticalProduct.FragmentNumber + "^" + ionTwo.Charge;
                spectraTwoDictionary.Add(productWithChargeInTwo, ionTwo);
            }
            var oneKeyList = spectraOneDictionary.Keys.ToList();
            var twoKeyList = spectraTwoDictionary.Keys.ToList();
            List<String> AllproductWithCharge = new List<string>();
            AllproductWithCharge.AddRange(oneKeyList);
            foreach (var x in twoKeyList)
            {
                if (!AllproductWithCharge.Contains(x))
                {
                    AllproductWithCharge.Add(x);
                }

            }

            //Tolerance tolerance = new PpmTolerance(40);
            foreach (var productWithCharge in AllproductWithCharge)
            {
                if (spectraOneDictionary.ContainsKey(productWithCharge) && spectraTwoDictionary.ContainsKey(productWithCharge))
                {
                    var newMz = (spectraOneDictionary[productWithCharge].Mz + spectraTwoDictionary[productWithCharge].Mz) / 2;
                    var newNorIntensity = (spectraOneDictionary[productWithCharge].Intensity / intensitySum1 + spectraTwoDictionary[productWithCharge].Intensity / intensitySum2) / 2;
                    Product product = spectraOneDictionary[productWithCharge].NeutralTheoreticalProduct;
                    var newIon = new MatchedFragmentIon(ref product, newMz, newNorIntensity, spectraOneDictionary[productWithCharge].Charge);
                    averageTwoSpectraResult.Add(newIon);
                }
                else if (spectraOneDictionary.ContainsKey(productWithCharge))
                {
                    Product productOne = spectraOneDictionary[productWithCharge].NeutralTheoreticalProduct;
                    var oneIon = new MatchedFragmentIon(ref productOne, spectraOneDictionary[productWithCharge].Mz, spectraOneDictionary[productWithCharge].Intensity / intensitySum1, spectraOneDictionary[productWithCharge].Charge);
                    averageTwoSpectraResult.Add(oneIon);
                }
                else if (spectraTwoDictionary.ContainsKey(productWithCharge))
                {
                    Product productTwo = spectraTwoDictionary[productWithCharge].NeutralTheoreticalProduct;
                    var twoIon = new MatchedFragmentIon(ref productTwo, spectraTwoDictionary[productWithCharge].Mz, spectraTwoDictionary[productWithCharge].Intensity / intensitySum2, spectraTwoDictionary[productWithCharge].Charge);
                    averageTwoSpectraResult.Add(twoIon);
                }
            }
            var testmz1 = averageTwoSpectraResult.Select(b => b.Mz).ToArray();
            var testin1 = averageTwoSpectraResult.Select(b => b.Intensity).ToArray();
            return averageTwoSpectraResult;
        }

        public static double CalculateDecoyNormalizedSpectralAngle(PeptideWithSetModifications peptide, Ms2ScanWithSpecificMass scan, CommonParameters commonParameters)

        {
            var theoreticalLibraryIons = new List<Product>();
            peptide.Fragment(commonParameters.DissociationType, commonParameters.DigestionParams.FragmentationTerminus, theoreticalLibraryIons);
            double mzCutoff = 300;
            int fragmentNumberCutoff = 3;

            // if the spectrum has no peaks
            if (scan.TheScan.MassSpectrum.XArray.Length == 0)
            {
                return 0;
            }
            List<MatchedFragmentIon> libraryDecoyIons = new List<MatchedFragmentIon>();
            Dictionary<MatchedFragmentIon, MatchedFragmentIon> matchedIons = new Dictionary<MatchedFragmentIon, MatchedFragmentIon>();

            // search for each theoretical ion
            for (int i = 0; i < theoreticalLibraryIons.Count; i++)
            {
                var libraryIon = theoreticalLibraryIons[i];
                double theoreticalFragmentMz = Math.Round(libraryIon.NeutralMass.ToMz(1) / 1.0005079, 0) * 1.0005079;

                // see https://www.nature.com/articles/s41592-019-0426-7
                // "All non-zero fragment ions (m/z > 300, ion >3, no neutral loss fragment ions) were considered for spectral angle calculation"
                if (libraryIon.NeutralMass <= mzCutoff || libraryIon.FragmentNumber <= fragmentNumberCutoff)
                {
                    continue;
                }



                // get the closest peak in the spectrum to the library peak
                var closestExperimentalMass = scan.GetClosestExperimentalIsotopicEnvelope(libraryIon.NeutralMass);
                if (closestExperimentalMass != null && commonParameters.ProductMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, libraryIon.NeutralMass) && closestExperimentalMass.Charge <= scan.PrecursorCharge)//TODO apply this filter before picking the envelope
                {
                    var test = new Product(libraryIon.ProductType, libraryIon.Terminus,
                        libraryIon.NeutralMass, libraryIon.FragmentNumber,
                        libraryIon.AminoAcidPosition, libraryIon.NeutralLoss);
                    libraryDecoyIons.Add(new MatchedFragmentIon(ref libraryIon, theoreticalFragmentMz, 1, 1));
                    matchedIons.Add(new MatchedFragmentIon(ref libraryIon, theoreticalFragmentMz, 1, 1), new MatchedFragmentIon(ref test, closestExperimentalMass.MonoisotopicMass.ToMz(closestExperimentalMass.Charge), closestExperimentalMass.Peaks.First().intensity, closestExperimentalMass.Charge));
                }

            }

            // L2 norm
            double expNormalizer = Math.Sqrt(matchedIons.Sum(p => Math.Pow(p.Value.Intensity, 2)));
            double theorNormalizer = Math.Sqrt(theoreticalLibraryIons.Sum(p => Math.Pow(1, 2)));

            double dotProduct = 0;

            foreach (var libraryIon in libraryDecoyIons)
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
