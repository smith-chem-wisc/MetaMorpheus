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

        private readonly List<Protein> proteinList;

        private readonly List<ModificationWithMass> fixedModifications;

        private readonly List<ModificationWithMass> variableModifications;

        private readonly Psm[] globalPsms;

        private readonly Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans;

        private readonly double[] myScanPrecursorMasses;

        private readonly List<ProductType> lp;

        private readonly bool addCompIons;

        private readonly ICommonParameters commonParameters;

        private readonly List<DissociationType> dissociationTypes;
        private readonly Tolerance productMassTolerance;

        #endregion Private Fields

        #region Public Constructors

        public ClassicSearchEngine(Psm[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, List<ProductType> lp, MassDiffAcceptor searchMode, bool addCompIons, ICommonParameters CommonParameters, Tolerance productMassTolerance, List<string> nestedIds) : base(nestedIds)
        {
            this.globalPsms = globalPsms;
            this.arrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            this.myScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.proteinList = proteinList;
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
            Status("In classic search engine!");

            int totalProteins = proteinList.Count;

            var observedPeptides = new HashSet<CompactPeptide>();

            Status("Getting ms2 scans...");

            var lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;
            TerminusType terminusType = ProductTypeMethod.IdentifyTerminusType(lp);
            Status("Starting classic search...");
            if (totalProteins > 0)
            {
                Parallel.ForEach(Partitioner.Create(0, totalProteins), partitionRange =>
                {
                    var psms = new Psm[arrayOfSortedMS2Scans.Length];
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var protein = proteinList[i];
                        var digestedList = protein.Digest(commonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                        foreach (var yyy in digestedList)
                        {
                            var correspondingCompactPeptide = yyy.CompactPeptide(terminusType);
                            if (!commonParameters.ConserveMemory)
                            {
                                var peptideWasObserved = observedPeptides.Contains(correspondingCompactPeptide);
                                if (peptideWasObserved)
                                    continue;
                                lock (observedPeptides)
                                {
                                    peptideWasObserved = observedPeptides.Contains(correspondingCompactPeptide);
                                    if (peptideWasObserved)
                                        continue;
                                    observedPeptides.Add(correspondingCompactPeptide);
                                }
                            }

                            var productMasses = correspondingCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                            Array.Sort(productMasses);

                            var searchMode = this.searchMode;
                            foreach (ScanWithIndexAndNotchInfo scanWithIndexAndNotchInfo in GetAcceptableScans(correspondingCompactPeptide.MonoisotopicMassIncludingFixedMods, searchMode).ToList())
                            {
                                if (scanWithIndexAndNotchInfo.theScan.OneBasedScanNumber != 10419)
                                    continue;
                                double thePrecursorMass = scanWithIndexAndNotchInfo.theScan.PrecursorMass;
                                var score = CalculatePeptideScore(scanWithIndexAndNotchInfo.theScan.TheScan, productMassTolerance, productMasses, thePrecursorMass, dissociationTypes, addCompIons, 0);

                                if (score > commonParameters.ScoreCutoff || commonParameters.CalculateEValue)
                                {
                                    if (psms[scanWithIndexAndNotchInfo.scanIndex] == null)
                                        psms[scanWithIndexAndNotchInfo.scanIndex] = new Psm(correspondingCompactPeptide, scanWithIndexAndNotchInfo.notch, score, scanWithIndexAndNotchInfo.scanIndex, scanWithIndexAndNotchInfo.theScan);
                                    else
                                        psms[scanWithIndexAndNotchInfo.scanIndex].AddOrReplace(correspondingCompactPeptide, score, scanWithIndexAndNotchInfo.notch, commonParameters.ReportAllAmbiguity);
                                    if (commonParameters.CalculateEValue)
                                        psms[scanWithIndexAndNotchInfo.scanIndex].UpdateAllScores(score);
                                }
                            }
                        }
                    }
                    lock (lockObject)
                    {
                        for (int i = 0; i < globalPsms.Length; i++)
                            if (psms[i] != null)
                            {
                                if (globalPsms[i] == null)
                                    globalPsms[i] = psms[i];
                                else
                                {
                                    if (commonParameters.CalculateEValue)
                                        globalPsms[i].SumAllScores(psms[i]);
                                    globalPsms[i].AddOrReplace(psms[i], commonParameters.ReportAllAmbiguity);
                                }
                            }
                        proteinsSeen += partitionRange.Item2 - partitionRange.Item1;
                        var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                        if (new_progress > old_progress)
                        {
                            ReportProgress(new ProgressEventArgs(new_progress, "Performing classic search...", nestedIds));
                            old_progress = new_progress;
                        }
                    }
                });
            }
            if (commonParameters.CalculateEValue)
                Parallel.ForEach(Partitioner.Create(0, globalPsms.Length), partitionRange =>
                 {
                     for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                         if (globalPsms[i] != null && globalPsms[i].Score < commonParameters.ScoreCutoff)
                             globalPsms[i] = null;
                 });
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