using Chemistry;
using Easy.Common.Extensions;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
namespace EngineLayer
{
    public static class SpectralLibrarySearchFunction
    {
        public static void CalculateSpectralAngles(SpectralLibrary spectralLibrary, SpectralMatch[] psms,
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
                // Spectrum similarity calculations often use ppm tolerance for mass spectra
                // If tolerance is Absolute, then we need to convert it to PpmTolerance.
                var isPpm = commonParameters.ProductMassTolerance is PpmTolerance;

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
                                // If ppm, keep as-is; if Absolute (Da), convert to ppm using the scan's precursor m/z
                                double tolerance = isPpm
                                    ? commonParameters.ProductMassTolerance.Value
                                    : commonParameters.ProductMassTolerance.Value / scan.PrecursorMonoisotopicPeakMz * 1e6;
                                List<double> pwsmSpectralAngles = new();
                                foreach (var bestMatch in psms[i].BestMatchingBioPolymersWithSetMods)
                                {
                                    if (spectralLibrary.TryGetSpectrum(bestMatch.FullSequence, scan.PrecursorCharge, out var librarySpectrum))
                                    {
                                        SpectralSimilarity s = new SpectralSimilarity(scan.TheScan.MassSpectrum, librarySpectrum.XArray, librarySpectrum.YArray,
                                            SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, tolerance, false);
                                            SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, commonParameters.ProductMassTolerance.Value, false);
>>>>>>>>> Temporary merge branch 2
                                        if (s.SpectralContrastAngle().HasValue)
                                        {
                                            pwsmSpectralAngles.Add((double)s.SpectralContrastAngle());
                                        }
                                    }
                                    // Fallback: if peptide is decoy without library spectrum, look for the target's spectrum and generate decoy spectrum by reversing
                                    else if (bestMatch.IsDecoy && spectralLibrary.TryGetSpectrum(bestMatch.SpecificBioPolymer.Description, scan.PrecursorCharge, out var targetlibrarySpectrum))
                                    {
                                        var decoyPeptideTheorProducts = new List<Product>();
                                        bestMatch.SpecificBioPolymer.Fragment(commonParameters.DissociationType, commonParameters.DigestionParams.FragmentationTerminus, decoyPeptideTheorProducts);
                                        var decoylibrarySpectrum = LibrarySpectrum.GetDecoyLibrarySpectrumFromTargetByReverse(targetlibrarySpectrum, decoyPeptideTheorProducts);
                                        SpectralSimilarity s = new SpectralSimilarity(scan.TheScan.MassSpectrum, decoylibrarySpectrum.Select(x => x.Mz).ToArray(),
                                            decoylibrarySpectrum.Select(x => x.Intensity).ToArray(), SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
                                            tolerance, false);
                                        if (s.SpectralContrastAngle().HasValue)
                                        {
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

    }
}
