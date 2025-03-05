using EngineLayer.Util;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.ClassicSearch
{
    public class MiniClassicSearchEngine
    {
        private readonly SpectralLibrary SpectralLibrary;
        private readonly MassDiffAcceptor SearchMode;
        private readonly Ms2ScanWithSpecificMass[] MS2ByRetentionTime;
        private readonly double[] ArrayOfRTs;
        private CommonParameters myCommonParameters;
        private CommonParameters FileSpecificParameters;

        public MiniClassicSearchEngine(
            Ms2ScanWithSpecificMass[] arrayOfRTSortedMS2Scans,
            MassDiffAcceptor searchMode,
            CommonParameters commonParameters,
            SpectralLibrary spectralLibrary,
            CommonParameters fileSpecificParameters)
        {
            MS2ByRetentionTime = arrayOfRTSortedMS2Scans;
            ArrayOfRTs = MS2ByRetentionTime.Select(p => p.RetentionTime).ToArray();
            SearchMode = searchMode;
            SpectralLibrary = spectralLibrary;
            myCommonParameters = commonParameters;
            FileSpecificParameters = fileSpecificParameters ?? commonParameters;
            // Each instance of MCSE will be specific to one file.
        }

        /// <summary>
        /// Searches all Ms2 scans in a 2 minute window around a peak apex for possible PSMs.
        /// Psms are scored against a peptideWithSetModifications that acted as a donor in MBR.
        /// Calculates traditional and spectral contrast angle scores
        /// </summary>
        /// <param name="donorPwsm"> Ms2 scans in window are searched for matches to this donor peptide</param>
        /// <param name="peakApexRT"> The center of the 2 minute window where the search occurs</param>
        /// <returns></returns>
        public IEnumerable<SpectralMatch> SearchAroundPeak(IBioPolymerWithSetMods donorPwsm, double peakApexRT)
        {
            var targetFragmentsForEachDissociationType = new Dictionary<DissociationType, List<Product>>();

            // check if we're supposed to autodetect dissociation type from the scan header or not
            if (FileSpecificParameters.DissociationType == DissociationType.Autodetect)
            {
                foreach (var item in GlobalVariables.AllSupportedDissociationTypes.Where(p => p.Value != DissociationType.Autodetect))
                {
                    targetFragmentsForEachDissociationType.Add(item.Value, new List<Product>());
                }
            }
            else
            {
                targetFragmentsForEachDissociationType.Add(FileSpecificParameters.DissociationType, new List<Product>());
            }

            // score each scan that has an acceptable precursor mass
            IEnumerable<ExtendedScanWithIndexAndNotchInfo> acceptableScans = GetAcceptableScans(donorPwsm.MonoisotopicMass, peakApexRT, SearchMode); // GetAcceptableScans is asynchronous, in case you care
            if (!acceptableScans.Any())
            {
                return null;
            }

            List<SpectralMatch> acceptablePsms = new();
            foreach (ExtendedScanWithIndexAndNotchInfo scan in acceptableScans)
            {
                var dissociationType = FileSpecificParameters.DissociationType == DissociationType.Autodetect ?
                    scan.TheScan.TheScan.DissociationType.Value : FileSpecificParameters.DissociationType;

                if (!targetFragmentsForEachDissociationType.TryGetValue(dissociationType, out var peptideTheorProducts))
                {
                    //TODO: print some kind of warning here. the scan header dissociation type was unknown
                    continue;
                }

                // check if we've already generated theoretical fragments for this peptide+dissociation type
                if (peptideTheorProducts.Count == 0)
                {
                    donorPwsm.Fragment(dissociationType, FileSpecificParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);
                    targetFragmentsForEachDissociationType[dissociationType] = peptideTheorProducts;
                }

                // match theoretical target ions to spectrum
                List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan.TheScan, peptideTheorProducts, FileSpecificParameters);

                // calculate the peptide's score
                double thisScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan.TheScan, matchedIons);

                // Add psm to list
                acceptablePsms.Add(new PeptideSpectralMatch(donorPwsm, scan.Notch, thisScore, scan.ScanIndex, scan.TheScan, FileSpecificParameters, matchedIons, 0));
            }

            IEnumerable<SpectralMatch> matchedSpectra = acceptablePsms.Where(p => p != null);
            foreach (SpectralMatch psm in matchedSpectra)
            {
                psm.ResolveAllAmbiguities();
            }

            CalculateSpectralAngles(SpectralLibrary, acceptablePsms.ToArray(), GetScansInWindow(peakApexRT), myCommonParameters, FileSpecificParameters);

            return matchedSpectra;
        }

        public static void CalculateSpectralAngles(SpectralLibrary spectralLibrary, SpectralMatch[] psms,
            Ms2ScanWithSpecificMass[] arrayOfSortedMs2Scans, CommonParameters commonParameters, CommonParameters fileSpecificParameters)
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
                for (int i = 0; i < psms.Length; i++)
                {

                    lock (myLocks[i])
                    {
                        if (psms[i] != null)
                        {
                            Ms2ScanWithSpecificMass scan = arrayOfSortedMs2Scans[psms[i].ScanIndex];
                            List<double> pwsmSpectralAngles = new();
                            foreach (var bestMatch in psms[i].BestMatchingBioPolymersWithSetMods)
                            {
                                //if peptide is target, directly look for the target's spectrum in the spectral library
                                if (!bestMatch.IsDecoy && spectralLibrary.TryGetSpectrum(bestMatch.FullSequence, scan.PrecursorCharge, out var librarySpectrum))
                                {
                                    SpectralSimilarity s = new SpectralSimilarity(scan.TheScan.MassSpectrum, librarySpectrum.XArray, librarySpectrum.YArray,
                                        SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, fileSpecificParameters.ProductMassTolerance.Value, false);
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
            }
        }

        private IEnumerable<ExtendedScanWithIndexAndNotchInfo> GetAcceptableScans(double peptideMonoisotopicMass, double apexRT, MassDiffAcceptor searchMode)
        {
            Ms2ScanWithSpecificMass[] arrayOfSortedMs2Scans = GetScansInWindow(apexRT);
            double[] myScanPrecursorMasses = arrayOfSortedMs2Scans.Select(p => p.PrecursorMass).ToArray();
            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMonoisotopicMass).ToList())
            {
                int scanIndex = GetFirstScanWithMassOverOrEqual(allowedIntervalWithNotch.Minimum, myScanPrecursorMasses);
                if (scanIndex < arrayOfSortedMs2Scans.Length)
                {
                    var scanMass = myScanPrecursorMasses[scanIndex];
                    while (scanMass <= allowedIntervalWithNotch.Maximum)
                    {
                        var scan = arrayOfSortedMs2Scans[scanIndex];
                        yield return new ExtendedScanWithIndexAndNotchInfo(scan, allowedIntervalWithNotch.Notch, scanIndex);
                        scanIndex++;
                        if (scanIndex == arrayOfSortedMs2Scans.Length)
                        {
                            break;
                        }

                        scanMass = myScanPrecursorMasses[scanIndex];
                    }
                }
            }
        }

        /// <summary>
        /// Finds MS2 scans falling within the relevant time window
        /// </summary>
        /// <param name="mbrPeak"> Acceptor peak </param>
        /// <param name="arrayOfRTs"> </param>
        /// <param name="arrayOfMs2ScansSortedByRT"> </param>
        /// <returns> An array of MS2 scans falling within a 2 minute retention time window of the Acceptor peak apex.
        ///           This array is sorted by precursor mass. </returns>
        private Ms2ScanWithSpecificMass[] GetScansInWindow(double apexRT)
        {
            double peakHalfWidth = 1.0; //Placeholder value to determine retention time window
            int startIndex = Array.BinarySearch(ArrayOfRTs, apexRT - peakHalfWidth);
            if (startIndex < 0)
                startIndex = ~startIndex;
            int endIndex = Array.BinarySearch(ArrayOfRTs, apexRT + peakHalfWidth);
            if (endIndex < 0)
                endIndex = ~endIndex;

            return MS2ByRetentionTime[startIndex..endIndex].OrderBy(b => b.PrecursorMass).ToArray();
        }

        private int GetFirstScanWithMassOverOrEqual(double minimum, double[] precursorMasses)
        {
            int index = Array.BinarySearch(precursorMasses, minimum);
            if (index < 0)
            {
                index = ~index;
            }

            // index of the first element that is larger than value
            return index;
        }
    }
}