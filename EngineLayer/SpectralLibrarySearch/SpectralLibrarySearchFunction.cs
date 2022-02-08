using Chemistry;
using MassSpectrometry;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using Proteomics.ProteolyticDigestion;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using MassSpectrometry.MzSpectra;

namespace EngineLayer
{
    public static class SpectralLibrarySearchFunction
    {
        public static void CalculateSpectralAngles(SpectralLibrary spectralLibrary, PeptideSpectralMatch[] psms,
            Ms2ScanWithSpecificMass[] arrayOfSortedMs2Scans, CommonParameters commonParameters)
        {
            if (spectralLibrary != null)
            {
                // one lock for each MS2 scan; a scan can only be accessed by one thread at a time
                var myLocks = new object[psms.Length];
                for (int i = 0; i < myLocks.Length; i++)
                {
                    myLocks[i] = new object();
                }

                int maxThreadsPerFile = commonParameters.MaxThreadsToUsePerFile;
                int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
                Parallel.ForEach(threads, (i) =>
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }
                    for (; i < psms.Length; i += maxThreadsPerFile)
                    {

                        lock (myLocks[i])
                        {
                            if (psms[i] != null)
                            {
                                Ms2ScanWithSpecificMass scan = arrayOfSortedMs2Scans[psms[i].ScanIndex];
                                List<(int, PeptideWithSetModifications)> pwsms = new();
                                List<double> pwsmSpectralAngles = new();
                                foreach (var (Notch, Peptide) in psms[i].BestMatchingPeptides)
                                {
                                    //if peptide is target, directly look for the target's spectrum in the spectral library
                                    if (!Peptide.Protein.IsDecoy && spectralLibrary.TryGetSpectrum(Peptide.FullSequence, scan.PrecursorCharge, out var librarySpectrum))
                                    {
                                        SpectralSimilarity s = new SpectralSimilarity(scan.TheScan.MassSpectrum, librarySpectrum.XArray, librarySpectrum.YArray, SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, commonParameters.ProductMassTolerance.Value, false);
                                        if (s.SpectralContrastAngle().HasValue)
                                        {
                                            pwsms.Add((Notch, Peptide));
                                            pwsmSpectralAngles.Add((double)s.SpectralContrastAngle());
                                        }

                                    }

                                    //if peptide is decoy, look for the decoy's corresponding target's spectrum in the spectral library and generate decoy spectrum by function GetDecoyLibrarySpectrumFromTargetByRevers
                                    else if (Peptide.Protein.IsDecoy && spectralLibrary.TryGetSpectrum(Peptide.PeptideDescription, scan.PrecursorCharge, out var targetlibrarySpectrum))
                                    {
                                        var decoyPeptideTheorProducts = new List<Product>();
                                        Peptide.Fragment(commonParameters.DissociationType, commonParameters.DigestionParams.FragmentationTerminus, decoyPeptideTheorProducts);
                                        var decoylibrarySpectrum = GetDecoyLibrarySpectrumFromTargetByReverse(targetlibrarySpectrum, decoyPeptideTheorProducts);
                                        SpectralSimilarity s = new SpectralSimilarity(scan.TheScan.MassSpectrum, decoylibrarySpectrum.Select(x => x.Mz).ToArray(),decoylibrarySpectrum.Select(x => x.Intensity).ToArray(), SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, commonParameters.ProductMassTolerance.Value, false);
                                        if (s.SpectralContrastAngle().HasValue)
                                        {
                                            pwsms.Add((Notch, Peptide));
                                            pwsmSpectralAngles.Add((double)s.SpectralContrastAngle());
                                        }
                                    }
                                }
                                if (pwsmSpectralAngles.Count > 0 && !pwsmSpectralAngles.Max().Equals(null))
                                {
                                    psms[i].SpectralAngle = pwsmSpectralAngles.Max();
                                }
                                else
                                {
                                    psms[i].SpectralAngle = -1;
                                }
                            }
                        }

                    }
                });
            }
                
        }


        /// <summary>
        /// Calculates the spectral angle, as described by Prosit ( https://www.nature.com/articles/s41592-019-0426-7 ).
        /// </summary>
        public static double CalculateSquareIntensitySpectralAngle(List<MatchedFragmentIon> spectrumLibraryIons, MsDataScan scan,  CommonParameters commonParameters)
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
            for (int i = 0; i < spectrumLibraryIons.Count; i++)
            {
                var libraryIon = spectrumLibraryIons[i];

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
            //psm.LibraryMatchedFragments = matchedIons.Values.ToList();
            // L2 norm
            double expNormalizer = Math.Sqrt(matchedIons.Sum(p => Math.Pow(Math.Sqrt(p.Value.Intensity), 2)));
            double theorNormalizer = Math.Sqrt(spectrumLibraryIons.Sum(p => Math.Pow(Math.Sqrt(p.Intensity), 2)));

            double dotProduct = 0;

            foreach (var libraryIon in spectrumLibraryIons)
            {
                if (matchedIons.TryGetValue(libraryIon, out var experIon))
                {
                    dotProduct += (Math.Sqrt(libraryIon.Intensity) / theorNormalizer) * (Math.Sqrt(experIon.Intensity) / expNormalizer);
                }
            }

            double normalizedSpectralAngle = 1 - (2 * Math.Acos(dotProduct) / Math.PI);

            return normalizedSpectralAngle;
        }

        // For decoy library spectrum generation, we use the predicted m/z valuse of the decoy sequence and we use the decoy's corresponding target's library spectrum's intensity values as decoy's intensities
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
                var matchedIons = mutiplePsms[i].MatchedFragmentIons.ToList();
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
                libraryIons.Add(new MatchedFragmentIon(ref temProduct, eachMz.Value.Where(p=>p!=0).Average(), allIntensityAfterNormalize[eachMz.Key] / max, eachMz.Key.Item2));
            }
            return new LibrarySpectrum(mutiplePsms[0].FullSequence, mutiplePsms.Select(p => p.ScanPrecursorMonoisotopicPeakMz).Average(), mutiplePsms[0].ScanPrecursorCharge, libraryIons, mutiplePsms.Select(p => p.ScanRetentionTime).Average(), false);
        }
    }
}
