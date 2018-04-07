using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.ClassicSearch
{
    public class ClassicSearchEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly MassDiffAcceptor searchMode;

        private readonly List<Protein> proteins;

        private readonly List<ModificationWithMass> fixedModifications;

        private readonly List<ModificationWithMass> variableModifications;

        private readonly PeptideSpectralMatch[] peptideSpectralMatches;

        private readonly Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans;

        private readonly double[] myScanPrecursorMasses;

        private readonly List<ProductType> lp;

        private readonly bool addCompIons;

        private readonly ICommonParameters commonParameters;

        private readonly List<DissociationType> dissociationTypes;
        private readonly Tolerance productMassTolerance;

        #endregion Private Fields

        #region Public Constructors

        public ClassicSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, List<ProductType> lp, MassDiffAcceptor searchMode, bool addCompIons, ICommonParameters CommonParameters, Tolerance productMassTolerance, List<string> nestedIds) : base(nestedIds)
        {
            this.peptideSpectralMatches = globalPsms;
            this.arrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            this.myScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.proteins = proteinList;
            this.searchMode = searchMode;
            this.lp = lp;
            this.addCompIons = addCompIons;
            this.dissociationTypes = DetermineDissociationType(lp);
            this.commonParameters = CommonParameters;
            this.productMassTolerance = productMassTolerance;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Getting ms2 scans...");
            
            double proteinsSearched = 0;
            int oldPercentProgress = 0;
            TerminusType terminusType = ProductTypeMethod.IdentifyTerminusType(lp);

            // one lock for each MS2 scan; a scan can only be accessed by one thread at a time
            var myLocks = new object[peptideSpectralMatches.Length];
            for (int i = 0; i < myLocks.Length; i++)
            {
                myLocks[i] = new object();
            }

            Status("Performing classic search...");
            
            if (proteins.Any())
            {
                Parallel.ForEach(Partitioner.Create(0, proteins.Count), partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        // digest each protein into peptides and search for each peptide in all spectra within precursor mass tolerance
                        foreach (var peptide in proteins[i].Digest(commonParameters.DigestionParams, fixedModifications, variableModifications))
                        {
                            var compactPeptide = peptide.CompactPeptide(terminusType);

                            var productMasses = compactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                            Array.Sort(productMasses);

                            foreach (ScanWithIndexAndNotchInfo scan in GetAcceptableScans(compactPeptide.MonoisotopicMassIncludingFixedMods, searchMode))
                            {
                                double scanPrecursorMass = scan.theScan.PrecursorMass;

                                var thisScore = CalculatePeptideScore(scan.theScan.TheScan, productMassTolerance, productMasses, scanPrecursorMass, dissociationTypes, addCompIons, 0);
                                bool meetsScoreCutoff = thisScore > commonParameters.ScoreCutoff;
                                bool scoreImprovement = peptideSpectralMatches[scan.scanIndex] == null || (thisScore - peptideSpectralMatches[scan.scanIndex].RunnerUpScore) > -PeptideSpectralMatch.tolForScoreDifferentiation;

                                // this is thread-safe because even if the score improves from another thread writing to this PSM,
                                // the lock combined with AddOrReplace method will ensure thread safety
                                if ((meetsScoreCutoff && scoreImprovement) || commonParameters.CalculateEValue)
                                {
                                    // valid hit (met the cutoff score); lock the scan to prevent other threads from accessing it
                                    lock (myLocks[scan.scanIndex])
                                    {
                                        if (peptideSpectralMatches[scan.scanIndex] == null)
                                        {
                                            peptideSpectralMatches[scan.scanIndex] = new PeptideSpectralMatch(compactPeptide, scan.notch, thisScore, scan.scanIndex, scan.theScan);
                                        }
                                        else
                                        {
                                            peptideSpectralMatches[scan.scanIndex].AddOrReplace(compactPeptide, thisScore, scan.notch, commonParameters.ReportAllAmbiguity);
                                        }

                                        if (commonParameters.CalculateEValue)
                                        {
                                            peptideSpectralMatches[scan.scanIndex].AddThisScoreToScoreDistribution(thisScore);
                                        }
                                    }
                                }
                            }
                        }

                        // report search progress (proteins searched so far out of total proteins in database)
                        proteinsSearched++;
                        var percentProgress = (int)((proteinsSearched / proteins.Count) * 100);

                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportProgress(new ProgressEventArgs(percentProgress, "Performing classic search... ", nestedIds));
                        }
                    }
                });
            }

            // remove peptides below the score cutoff that were stored to calculate expectation values
            if (commonParameters.CalculateEValue)
            {
                for (int i = 0; i < peptideSpectralMatches.Length; i++)
                {
                    if (peptideSpectralMatches[i] != null && peptideSpectralMatches[i].Score < commonParameters.ScoreCutoff)
                    {
                        peptideSpectralMatches[i] = null;
                    }
                }
            }

            return new MetaMorpheusEngineResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(double peptideMonoisotopicMass, MassDiffAcceptor searchMode)
        {
            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervals(peptideMonoisotopicMass).ToList())
            {
                DoubleRange allowedInterval = allowedIntervalWithNotch.allowedInterval;
                int scanIndex = GetFirstScanWithMassOverOrEqual(allowedInterval.Minimum);
                if (scanIndex < arrayOfSortedMS2Scans.Length)
                {
                    var scanMass = myScanPrecursorMasses[scanIndex];
                    while (scanMass <= allowedInterval.Maximum)
                    {
                        var theScan = arrayOfSortedMS2Scans[scanIndex];
                        yield return new ScanWithIndexAndNotchInfo(theScan, allowedIntervalWithNotch.Notch, scanIndex);
                        scanIndex++;
                        if (scanIndex == arrayOfSortedMS2Scans.Length)
                            break;
                        scanMass = myScanPrecursorMasses[scanIndex];
                    }
                }
            }
        }

        private int GetFirstScanWithMassOverOrEqual(double minimum)
        {
            int index = Array.BinarySearch(myScanPrecursorMasses, minimum);
            if (index < 0)
                index = ~index;

            // index of the first element that is larger than value
            return index;
        }

        #endregion Private Methods
    }
}