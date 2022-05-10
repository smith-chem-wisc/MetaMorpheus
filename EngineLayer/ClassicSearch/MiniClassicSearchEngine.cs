using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.ClassicSearch
{
    public class MiniClassicSearchEngine : MetaMorpheusEngine
    {
        private readonly SpectralLibrary SpectralLibrary;
        private readonly MassDiffAcceptor SearchMode;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        private readonly PeptideWithSetModifications PeptideWithSetMods;
        private readonly PeptideSpectralMatch[] PeptideSpectralMatches;
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private readonly double[] MyScanPrecursorMasses;

        public MiniClassicSearchEngine(PeptideWithSetModifications pwsm, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans,
            List<Modification> variableModifications, List<Modification> fixedModifications,
            MassDiffAcceptor searchMode, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            SpectralLibrary spectralLibrary, List<string> nestedIds)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            PeptideWithSetMods = pwsm;
            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            PeptideSpectralMatches = new PeptideSpectralMatch[arrayOfSortedMS2Scans.Length];
            MyScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;


            SearchMode = searchMode;
            SpectralLibrary = spectralLibrary;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Getting ms2 scans...");

            // one lock for each MS2 scan; a scan can only be accessed by one thread at a time
            var myLocks = new object[ArrayOfSortedMS2Scans.Length]; // There aren't multiple PSMs, just the one we're looking at. 
            for (int i = 0; i < myLocks.Length; i++)
            {
                myLocks[i] = new object();
            }

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();

            Status("Performing mini classic search...");

            var targetFragmentsForEachDissociationType = new Dictionary<DissociationType, List<Product>>();

            // check if we're supposed to autodetect dissociation type from the scan header or not
            if (CommonParameters.DissociationType == DissociationType.Autodetect)
            {
                foreach (var item in GlobalVariables.AllSupportedDissociationTypes.Where(p => p.Value != DissociationType.Autodetect))
                {
                    targetFragmentsForEachDissociationType.Add(item.Value, new List<Product>());
                }
            }
            else
            {
                targetFragmentsForEachDissociationType.Add(CommonParameters.DissociationType, new List<Product>());
            }

            // clear fragments from the last peptide
            foreach (var fragmentSet in targetFragmentsForEachDissociationType)
            {
                fragmentSet.Value.Clear();
            }

            // score each scan that has an acceptable precursor mass
            // TODO: Parallelize this section of the code
            IEnumerable<ScanWithIndexAndNotchInfo> acceptableScans = GetAcceptableScans(PeptideWithSetMods.MonoisotopicMass, SearchMode); // GetAcceptableScans is asynchronous, in case you care
            foreach (ScanWithIndexAndNotchInfo scan in acceptableScans)
            {
                var dissociationType = CommonParameters.DissociationType == DissociationType.Autodetect ?
                    scan.TheScan.TheScan.DissociationType.Value : CommonParameters.DissociationType;

                if (!targetFragmentsForEachDissociationType.TryGetValue(dissociationType, out var peptideTheorProducts))
                {
                    //TODO: print some kind of warning here. the scan header dissociation type was unknown
                    continue;
                }

                // check if we've already generated theoretical fragments for this peptide+dissociation type
                if (peptideTheorProducts.Count == 0)
                {
                    PeptideWithSetMods.Fragment(dissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);
                }

                // match theoretical target ions to spectrum
                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan, peptideTheorProducts, CommonParameters);

                // calculate the peptide's score
                double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons);

                AddPeptideCandidateToPsm(scan, myLocks, thisScore, PeptideWithSetMods, matchedIons);
            }


            //foreach (PeptideSpectralMatch psm in PeptideSpectralMatches.Where(p => p != null))
            IEnumerable<PeptideSpectralMatch> matchedSpectra = PeptideSpectralMatches.Where(p => p != null);
            int numMatches = matchedSpectra.Count();
            int matchScanNumber = 0;
            foreach (PeptideSpectralMatch psm in matchedSpectra)
            {
                psm.ResolveAllAmbiguities();
                matchScanNumber = psm.ScanNumber;
            }

            CalculateSpectralAngles(SpectralLibrary, PeptideSpectralMatches, ArrayOfSortedMS2Scans, CommonParameters);

            return new MetaMorpheusEngineResults(this);
        }

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

        private void AddPeptideCandidateToPsm(ScanWithIndexAndNotchInfo scan, object[] myLocks, double thisScore, PeptideWithSetModifications peptide, List<MatchedFragmentIon> matchedIons)
        {
            // this is thread-safe because even if the score improves from another thread writing to this PSM,
            // the lock combined with AddOrReplace method will ensure thread safety

            // valid hit (met the cutoff score); lock the scan to prevent other threads from accessing it
            lock (myLocks[scan.ScanIndex])
            {
                bool scoreImprovement = PeptideSpectralMatches[scan.ScanIndex] == null || (thisScore - PeptideSpectralMatches[scan.ScanIndex].RunnerUpScore) > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

                if (scoreImprovement)
                {
                    if (PeptideSpectralMatches[scan.ScanIndex] == null)
                    {
                        PeptideSpectralMatches[scan.ScanIndex] = new PeptideSpectralMatch(peptide, scan.Notch, thisScore, scan.ScanIndex, scan.TheScan, CommonParameters, matchedIons, 0);
                    }
                    else
                    {
                        PeptideSpectralMatches[scan.ScanIndex].AddOrReplace(peptide, thisScore, scan.Notch, CommonParameters.ReportAllAmbiguity, matchedIons, 0);
                    }
                }
            }
        }

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(double peptideMonoisotopicMass, MassDiffAcceptor searchMode)
        {
            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMonoisotopicMass).ToList())
            {
                DoubleRange allowedInterval = allowedIntervalWithNotch.AllowedInterval;
                int scanIndex = GetFirstScanWithMassOverOrEqual(allowedInterval.Minimum);
                if (scanIndex < ArrayOfSortedMS2Scans.Length)
                {
                    var scanMass = MyScanPrecursorMasses[scanIndex];
                    while (scanMass <= allowedInterval.Maximum)
                    {
                        var scan = ArrayOfSortedMS2Scans[scanIndex];
                        yield return new ScanWithIndexAndNotchInfo(scan, allowedIntervalWithNotch.Notch, scanIndex);
                        scanIndex++;
                        if (scanIndex == ArrayOfSortedMS2Scans.Length)
                        {
                            break;
                        }

                        scanMass = MyScanPrecursorMasses[scanIndex];
                    }
                }
            }
        }

        private int GetFirstScanWithMassOverOrEqual(double minimum)
        {
            int index = Array.BinarySearch(MyScanPrecursorMasses, minimum);
            if (index < 0)
            {
                index = ~index;
            }

            // index of the first element that is larger than value
            return index;
        }
    }
}