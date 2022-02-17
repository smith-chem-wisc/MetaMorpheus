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
        private readonly List<Protein> Proteins;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        
        
        private readonly double[] MyScanPrecursorMasses;

        public SingleCellMBREngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, FlashLfqResults flashLfqResults,
            List<Modification> variableModifications, List<Modification> fixedModifications, List<Protein> proteinList, MassDiffAcceptor searchMode, 
            CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, SpectralLibrary spectralLibrary, 
            List<string> nestedIds, bool writeSpectralLibrary)
            : base(commonParameters, fileSpecificParameters, nestedIds)
            {
                PeptideSpectralMatches = globalPsms;
                ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
                FlashLfqResults = flashLfqResults;
                MyScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
                VariableModifications = variableModifications;
                FixedModifications = fixedModifications;
                SpectralLibrary = spectralLibrary;
                Proteins = proteinList;
            }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Getting ms2 scans...");

            double proteinsSearched = 0;
            int oldPercentProgress = 0;

            // one lock for each MS2 scan; a scan can only be accessed by one thread at a time
            var myLocks = new object[PeptideSpectralMatches.Length];
            for (int i = 0; i < myLocks.Length; i++)
            {
                myLocks[i] = new object();
            }

            Status("Performing classic search...");

            if (Proteins.Any())
            {
                int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
                int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
                Parallel.ForEach(threads, (i) =>
                {
                    var targetFragmentsForEachDissociationType = new Dictionary<DissociationType, List<Product>>();
                    var decoyFragmentsForEachDissociationType = new Dictionary<DissociationType, List<Product>>();

                    // check if we're supposed to autodetect dissociation type from the scan header or not
                    if (CommonParameters.DissociationType == DissociationType.Autodetect)
                    {
                        foreach (var item in GlobalVariables.AllSupportedDissociationTypes.Where(p => p.Value != DissociationType.Autodetect))
                        {
                            targetFragmentsForEachDissociationType.Add(item.Value, new List<Product>());
                            decoyFragmentsForEachDissociationType.Add(item.Value, new List<Product>());
                        }
                    }
                    else
                    {
                        targetFragmentsForEachDissociationType.Add(CommonParameters.DissociationType, new List<Product>());
                        decoyFragmentsForEachDissociationType.Add(CommonParameters.DissociationType, new List<Product>());
                    }

                    for (; i < Proteins.Count; i += maxThreadsPerFile)
                    {
                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops) { return; }

                        // digest each protein into peptides and search for each peptide in all spectra within precursor mass tolerance
                        foreach (PeptideWithSetModifications peptide in Proteins[i].Digest(CommonParameters.DigestionParams, FixedModifications, VariableModifications, null, null))
                        {
                            PeptideWithSetModifications reversedOnTheFlyDecoy = null;

                            if (SpectralLibrary != null)
                            {
                                int[] newAAlocations = new int[peptide.BaseSequence.Length];
                                reversedOnTheFlyDecoy = peptide.GetReverseDecoyFromTarget(newAAlocations);
                            }

                            // clear fragments from the last peptide
                            foreach (var fragmentSet in targetFragmentsForEachDissociationType)
                            {
                                fragmentSet.Value.Clear();
                                decoyFragmentsForEachDissociationType[fragmentSet.Key].Clear();
                            }

                            // score each scan that has an acceptable precursor mass
                            foreach (ScanWithIndexAndNotchInfo scan in GetAcceptableScans(peptide.MonoisotopicMass, SearchMode))
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
                                    peptide.Fragment(dissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);
                                }

                                // match theoretical target ions to spectrum
                                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan, peptideTheorProducts, CommonParameters,
                                        matchAllCharges: WriteSpectralLibrary);

                                // calculate the peptide's score
                                double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons, fragmentsCanHaveDifferentCharges: WriteSpectralLibrary);

                                AddPeptideCandidateToPsm(scan, myLocks, thisScore, peptide, matchedIons);

                                if (SpectralLibrary != null)
                                {
                                    DecoyScoreForSpectralLibrarySearch(scan, reversedOnTheFlyDecoy, decoyFragmentsForEachDissociationType, dissociationType, myLocks);
                                }
                            }
                        }

                        // report search progress (proteins searched so far out of total proteins in database)
                        proteinsSearched++;
                        var percentProgress = (int)((proteinsSearched / Proteins.Count) * 100);

                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportProgress(new ProgressEventArgs(percentProgress, "Performing classic search... ", NestedIds));
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
    }
}