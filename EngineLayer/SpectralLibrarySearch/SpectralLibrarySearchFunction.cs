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
                    if (spectralLibrary != null && !peptide.Key.Protein.IsDecoy && spectralLibrary.TryGetSpectrum(peptide.Key.FullSequence, scan.PrecursorCharge, out var librarySpectrum))
                    {
                        double spectralAngle = CalculateSquareIntensitySpectralAngle(librarySpectrum.MatchedFragmentIons, scan.TheScan, psm, commonParameters);
                        psm.SpectralAngle = spectralAngle;
                    }
                    else if (spectralLibrary != null && peptide.Key.Protein.IsDecoy && spectralLibrary.TryGetSpectrum(peptide.Key.PeptideDescription, scan.PrecursorCharge, out var targetlibrarySpectrum))
                    {
                        var decoyPeptideTheorProducts = new List<Product>();
                        peptide.Key.Fragment(commonParameters.DissociationType, commonParameters.DigestionParams.FragmentationTerminus, decoyPeptideTheorProducts);
                        var decoylibrarySpectrum = GetDecoyLibrarySpectrumFromTargetByReverse(targetlibrarySpectrum, decoyPeptideTheorProducts);
                        double spectralAngle = CalculateSquareIntensitySpectralAngle(decoylibrarySpectrum, scan.TheScan, psm, commonParameters);
                        psm.SpectralAngle = spectralAngle;
                    }
                }
            }
        }

        /// <summary>
        /// Calculates the spectral angle, as described by Prosit ( https://www.nature.com/articles/s41592-019-0426-7 ).
        /// </summary>
        public static double CalculateSquareIntensitySpectralAngle(List<MatchedFragmentIon> theoreticalLibraryIons, MsDataScan scan, PeptideSpectralMatch psm, CommonParameters commonParameters)
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
                    if (!matchedIons.ContainsKey(libraryIon))
                    {
                        matchedIons.Add(libraryIon, new MatchedFragmentIon(ref test, mz, experimentalIntensity, libraryIon.Charge));
                    }
                }
            }
            psm.LibraryMatchedFragments = matchedIons.Values.ToList();
            // L2 norm
            double expNormalizer = Math.Sqrt(matchedIons.Sum(p => Math.Pow(Math.Sqrt(p.Value.Intensity), 2)));
            double theorNormalizer = Math.Sqrt(theoreticalLibraryIons.Sum(p => Math.Pow(Math.Sqrt(p.Intensity), 2)));

            double dotProduct = 0;

            foreach (var libraryIon in theoreticalLibraryIons)
            {
                if (matchedIons.TryGetValue(libraryIon, out var experIon))
                {
                    dotProduct += (Math.Sqrt(libraryIon.Intensity) / theorNormalizer) * (Math.Sqrt(experIon.Intensity) / expNormalizer);
                }
            }

            double normalizedSpectralAngle = 1 - (2 * Math.Acos(dotProduct) / Math.PI);

            return normalizedSpectralAngle;
        }


        public static List<MatchedFragmentIon> GetDecoyLibrarySpectrumFromTargetByReverse(LibrarySpectrum targetSpectrum, List<Product> decoyPeptideTheorProducts)
        {
            var decoyFragmentIons = new List<MatchedFragmentIon>();
            foreach (var targetIon in targetSpectrum.MatchedFragmentIons)
            {
                foreach (var decoyPeptideTheorIon in decoyPeptideTheorProducts)
                {
                    if (targetIon.NeutralTheoreticalProduct.ProductType == decoyPeptideTheorIon.ProductType && targetIon.NeutralTheoreticalProduct.FragmentNumber == decoyPeptideTheorIon.FragmentNumber)
                    {
                        double decoyFragmentMz = decoyPeptideTheorIon.NeutralMass.ToMz(targetIon.Charge);
                        Product temProduct = decoyPeptideTheorIon;
                        decoyFragmentIons.Add(new MatchedFragmentIon(ref temProduct, decoyFragmentMz, targetIon.Intensity, targetIon.Charge));
                    }
                }
            }
            return decoyFragmentIons;
        }

        // This function is used to converting several spectra matching the same sequence with same charge to a consensus library spectrum. In this function, ll the mz, intensities
        //are averaged with weight. The MetaMopeheus score which comes from database search are treated as weight.If a peak exists in less than 50% of all spectra, it will be discarded
        public static LibrarySpectrum ConvertingPsmsToConcensusSpectrumWithWeight(List<PeptideSpectralMatch> mutiplePsms)
        {
            Dictionary<Tuple<Product, int>, double[]> allIntensityBeforeNormalize = new Dictionary<Tuple<Product, int>, double[]>();
            Dictionary<Tuple<Product, int>, double> allIntensityAfterNormalize = new Dictionary<Tuple<Product, int>, double>();
            Dictionary<Tuple<Product, int>, double[]> allMz = new Dictionary<Tuple<Product, int>, double[]>();
            Dictionary<Tuple<Product, int>, double> allAverageMz = new Dictionary<Tuple<Product, int>, double>();
            var allIntensitySumForNormalize = new double[mutiplePsms.Count];
            Dictionary<Tuple<Product, int>, int> numberOfEachIons = new Dictionary<Tuple<Product, int>, int>();
            List<MatchedFragmentIon> libraryIons = new List<MatchedFragmentIon>();

            Dictionary<Tuple<Product, int>, double> weightScoreForEachIon = new Dictionary<Tuple<Product, int>, double>();

            for (int i = 0; i < mutiplePsms.Count; i++)
            {
                var matchedIons = mutiplePsms[i].LibraryMatchedFragments.ToList();
                for (int j = 0; j < matchedIons.Count; j++)
                {
                    Tuple<Product, int> newKey = new Tuple<Product, int>(matchedIons[j].NeutralTheoreticalProduct, matchedIons[j].Charge);
                    if (!allIntensityBeforeNormalize.ContainsKey(newKey))
                    {
                        var listOfIntensity = new double[mutiplePsms.Count];
                        var listOfMz = new double[mutiplePsms.Count];
                        listOfIntensity[i] = matchedIons[j].Intensity;
                        listOfMz[i] = matchedIons[j].Mz;
                        allIntensityBeforeNormalize.Add(newKey, listOfIntensity);
                        allMz.Add(newKey, listOfMz);
                        numberOfEachIons.Add(newKey, 1);
                        weightScoreForEachIon.Add(newKey, mutiplePsms[i].Score);
                    }
                    else
                    {
                        allIntensityBeforeNormalize[newKey][i] = matchedIons[j].Intensity;
                        allMz[newKey][i] = matchedIons[j].Mz;
                        int num = numberOfEachIons[newKey];
                        numberOfEachIons[newKey] = num + 1;
                        var tem = weightScoreForEachIon[newKey];
                        weightScoreForEachIon[newKey] = tem + mutiplePsms[i].Score;
                    }
                }
            }

            List<Tuple<Product, int>> commonIons = numberOfEachIons.Where(p => p.Value == mutiplePsms.Count).Select(q => q.Key).ToList();

            //If a peak exists in less than 50% of all spectra, it will be discarded
            List<Tuple<Product, int>> unqualifiedIons = numberOfEachIons.Where(p => p.Value < ((double)mutiplePsms.Count / (double)2)).Select(q => q.Key).ToList();

            // store all intensity values before normalize
            for (int k = 0; k < mutiplePsms.Count; k++)
            {
                double intensitySum = 0;
                foreach (var x in commonIons)
                {
                    intensitySum = intensitySum + allIntensityBeforeNormalize[x][k];
                }
                allIntensitySumForNormalize[k] = intensitySum;
            }

            // store all intensity values after normalize with weight. The MetaMopeheus score which comes from database search are treated as weight.
            foreach (var y in allIntensityBeforeNormalize.Where(p => !unqualifiedIons.Contains(p.Key)))
            {
                for (int l = 0; l < y.Value.Length; l++)
                {
                    double intensityBeforeNormalize = y.Value[l];
                    double intensityAfterNormalize = intensityBeforeNormalize / allIntensitySumForNormalize[l];
                    double intensityAfterWeight = (intensityAfterNormalize * mutiplePsms[l].Score) / weightScoreForEachIon[y.Key];
                    y.Value[l] = intensityAfterWeight;
                }
                allIntensityAfterNormalize.Add(y.Key, y.Value.Where(p => p > 0).ToList().Sum());
            }

            // second normalize by max of the intensity values
            var max = allIntensityAfterNormalize.Select(p => p.Value).Max();
            foreach (var eachMz in allMz.Where(p => !unqualifiedIons.Contains(p.Key)))
            {
                Product temProduct = eachMz.Key.Item1;
                libraryIons.Add(new MatchedFragmentIon(ref temProduct, eachMz.Value.Average(), allIntensityAfterNormalize[eachMz.Key] / max, eachMz.Key.Item2));
            }
            return new LibrarySpectrum(mutiplePsms[0].FullSequence, mutiplePsms.Select(p => p.ScanPrecursorMonoisotopicPeakMz).Average(), mutiplePsms[0].ScanPrecursorCharge, libraryIons, mutiplePsms.Select(p => p.ScanRetentionTime).Average(), false);
        }
    }
}
