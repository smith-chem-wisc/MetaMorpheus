using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.SingleCellMBR
{
    public class SingleCellMBREngine : MetaMorpheusEngine
    {
        private readonly PeptideSpectralMatch[] PeptideSpectralMatches;
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private readonly FlashLfqResults FlashLfqResults;
        private readonly SpectralLibrary SpectralLibrary;
        private readonly MassDiffAcceptor SearchMode;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;

        private readonly double[] MyScanPrecursorMasses;

        public SingleCellMBREngine(PeptideSpectralMatch[] allPeptides, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, FlashLfqResults flashLfqResults,
            List<Modification> variableModifications, List<Modification> fixedModifications, MassDiffAcceptor searchMode,
            CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, SpectralLibrary spectralLibrary,
            List<string> nestedIds)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            PeptideSpectralMatches = allPeptides;
            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            SearchMode = searchMode;
            FlashLfqResults = flashLfqResults;
            MyScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            SpectralLibrary = spectralLibrary;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Getting ms2 scans...");

            double peptidesSearched = 0;
            int oldPercentProgress = 0;

            // one lock for each MS2 scan; a scan can only be accessed by one thread at a time
            var myLocks = new object[PeptideSpectralMatches.Length];
            for (int i = 0; i < myLocks.Length; i++)
            {
                myLocks[i] = new object();
            }

            Status("Performing classic search...");

            if (PeptideSpectralMatches.Any())
            {
                int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
                int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
                Parallel.ForEach(threads, (i) =>
                {
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

                    for (; i < PeptideSpectralMatches.Length; i += maxThreadsPerFile)
                    {
                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops) { return; }

                        if (SpectralLibrary != null)
                        {
                            int[] newAAlocations = new int[PeptideSpectralMatches[i].BaseSequence.Length];
                        }

                        // clear fragments from the last peptide
                        foreach (var fragmentSet in targetFragmentsForEachDissociationType)
                        {
                            fragmentSet.Value.Clear();
                        }

                        // score each scan that has an acceptable precursor mass
                        foreach (ScanWithIndexAndNotchInfo scan in GetAcceptableScans(PeptideSpectralMatches[i].PeptideMonisotopicMass.Value, SearchMode))
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
                                foreach (var peptide in PeptideSpectralMatches[i].BestMatchingPeptides)
                                {
                                    peptide.Peptide.Fragment(dissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);
                                    // match theoretical target ions to spectrum
                                    List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan, peptideTheorProducts, CommonParameters);

                                    // calculate the peptide's score
                                    double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons);

                                    AddPeptideCandidateToPsm(scan, myLocks, thisScore, peptide.Peptide, matchedIons);
                                }
                            }

                        }

                        // report search progress (proteins searched so far out of total proteins in database)
                        peptidesSearched++;
                        var percentProgress = (int)((peptidesSearched / PeptideSpectralMatches.Length) * 100);

                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportProgress(new ProgressEventArgs(percentProgress, "Performing single cell MBR... ", NestedIds));
                        }
                    }
                });
            }

            foreach (PeptideSpectralMatch psm in PeptideSpectralMatches.Where(p => p != null))
            {
                psm.ResolveAllAmbiguities();
            }

            return new MetaMorpheusEngineResults(this);
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

        private void AddPeptideCandidateToPsm(ScanWithIndexAndNotchInfo scan, object[] myLocks, double thisScore, PeptideWithSetModifications peptide, List<MatchedFragmentIon> matchedIons)
        {
            bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;

            // this is thread-safe because even if the score improves from another thread writing to this PSM,
            // the lock combined with AddOrReplace method will ensure thread safety
            if (meetsScoreCutoff)
            {
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
        }
    }
}