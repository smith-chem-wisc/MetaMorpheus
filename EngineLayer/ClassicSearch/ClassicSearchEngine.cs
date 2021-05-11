using MassSpectrometry;
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
    public class ClassicSearchEngine : MetaMorpheusEngine
    {
        private readonly SpectralLibrary SpectralLibrary;
        private readonly MassDiffAcceptor SearchMode;
        private readonly List<Protein> Proteins;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        private readonly List<SilacLabel> SilacLabels;
        private readonly (SilacLabel StartLabel, SilacLabel EndLabel)? TurnoverLabels;
        private readonly PeptideSpectralMatch[] PeptideSpectralMatches;
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private readonly double[] MyScanPrecursorMasses;
        public Dictionary<PeptideWithSetModifications, PeptideWithSetModifications> decoyTargetPairs;

        public ClassicSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans,
            List<Modification> variableModifications, List<Modification> fixedModifications, List<SilacLabel> silacLabels, SilacLabel startLabel, SilacLabel endLabel,
            List<Protein> proteinList, MassDiffAcceptor searchMode, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            SpectralLibrary spectralLibrary, List<string> nestedIds)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            PeptideSpectralMatches = globalPsms;
            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            MyScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            SilacLabels = silacLabels;
            if (startLabel != null || endLabel != null) //else it's null
            {
                TurnoverLabels = (startLabel, endLabel);
            }
            Proteins = proteinList;
            SearchMode = searchMode;
            SpectralLibrary = spectralLibrary;
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
                    var fragmentsForDissociationTypes = new Dictionary<DissociationType, List<Product>>();

                    // check if we're supposed to autodetect dissociation type from the scan header or not
                    if (CommonParameters.DissociationType == DissociationType.Autodetect)
                    {
                        foreach (var item in GlobalVariables.AllSupportedDissociationTypes.Where(p => p.Value != DissociationType.Autodetect))
                        {
                            fragmentsForDissociationTypes.Add(item.Value, new List<Product>());
                        }
                    }
                    else
                    {
                        fragmentsForDissociationTypes.Add(CommonParameters.DissociationType, new List<Product>());
                    }

                    for (; i < Proteins.Count; i += maxThreadsPerFile)
                    {
                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops) { return; }

                        // digest each protein into peptides and search for each peptide in all spectra within precursor mass tolerance
                        foreach (PeptideWithSetModifications peptide in Proteins[i].Digest(CommonParameters.DigestionParams, FixedModifications, VariableModifications, SilacLabels, TurnoverLabels))
                        {

                            int[] newAAlocations = new int[peptide.BaseSequence.Length];
                            PeptideWithSetModifications decoy = DecoyOnTheFly.GetReverseDecoyFromTarget(peptide, newAAlocations);

                            // we need a function to get the original target sequence of a decoy peptide
                            lock (decoyTargetPairs)
                            {
                                if (!decoyTargetPairs.ContainsKey(decoy))
                                {
                                    decoyTargetPairs.Add(decoy, peptide);
                                }
                            }

                            foreach (var fragmentSet in fragmentsForDissociationTypes)
                            {
                                fragmentSet.Value.Clear();
                            }

                            foreach (ScanWithIndexAndNotchInfo scan in GetAcceptableScans(peptide.MonoisotopicMass, SearchMode))
                            {
                                //if (SpectralLibrary != null && !SpectralLibrary.ContainsSpectrum(peptide.FullSequence, scan.TheScan.PrecursorCharge))
                                //{
                                //    continue;
                                //}

                                var dissociationType = CommonParameters.DissociationType == DissociationType.Autodetect ?
                                    scan.TheScan.TheScan.DissociationType.Value : CommonParameters.DissociationType;

                                if (!fragmentsForDissociationTypes.TryGetValue(dissociationType, out var peptideTheorProducts))
                                {
                                    //TODO: print some kind of warning here. the scan header dissociation type was unknown
                                    continue;
                                }

                                var decoyPeptideTheorProducts = new List<Product>();
                                

                                // check if we've already generated theoretical fragments for this peptide+dissociation type
                                if (peptideTheorProducts.Count == 0)
                                {
                                    peptide.Fragment(dissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);
                                }

                                if (decoyPeptideTheorProducts.Count == 0)
                                {
                                    decoy.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, decoyPeptideTheorProducts);
                                }

                                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan, peptideTheorProducts, CommonParameters);
                                List<MatchedFragmentIon> decoyMatchedIons = MatchFragmentIons(scan.TheScan, decoyPeptideTheorProducts, CommonParameters);
                                double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons);
                                double decoyScore = CalculatePeptideScore(scan.TheScan.TheScan, decoyMatchedIons);
                                
                                //bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;
                                bool meetsScoreCutoff = Math.Max(thisScore, decoyScore) >= CommonParameters.ScoreCutoff;

                                // this is thread-safe because even if the score improves from another thread writing to this PSM,
                                // the lock combined with AddOrReplace method will ensure thread safety
                                if (meetsScoreCutoff)
                                {
                                    // valid hit (met the cutoff score); lock the scan to prevent other threads from accessing it
                                    lock (myLocks[scan.ScanIndex])
                                    {
                                        bool scoreImprovement = PeptideSpectralMatches[scan.ScanIndex] == null || (Math.Max(thisScore, decoyScore) - PeptideSpectralMatches[scan.ScanIndex].RunnerUpScore) > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

                                        if (scoreImprovement)
                                        {
                                            if (PeptideSpectralMatches[scan.ScanIndex] == null)
                                            {
                                                if (thisScore >= decoyScore)
                                                {
                                                    PeptideSpectralMatches[scan.ScanIndex] = new PeptideSpectralMatch(peptide, scan.Notch, thisScore, scan.ScanIndex, scan.TheScan, CommonParameters, matchedIons, 0);
                                                }
                                                else
                                                {
                                                    PeptideSpectralMatches[scan.ScanIndex] = new PeptideSpectralMatch(decoy, scan.Notch, decoyScore, scan.ScanIndex, scan.TheScan, CommonParameters, decoyMatchedIons, 0);

                                                }
                                            }
                                            else
                                            {
                                                if (thisScore >= decoyScore)
                                                {
                                                    PeptideSpectralMatches[scan.ScanIndex].AddOrReplace(peptide, thisScore, scan.Notch, CommonParameters.ReportAllAmbiguity, matchedIons, 0);
                                                }
                                                else
                                                {
                                                    PeptideSpectralMatches[scan.ScanIndex].AddOrReplace(decoy, decoyScore, scan.Notch, CommonParameters.ReportAllAmbiguity, decoyMatchedIons, 0);
                                                }
                                            }
                                        }
                                    }
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