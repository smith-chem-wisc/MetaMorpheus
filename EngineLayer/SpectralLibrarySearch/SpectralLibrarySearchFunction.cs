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
    }
}
