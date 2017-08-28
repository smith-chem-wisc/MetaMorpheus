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

        private readonly List<MassDiffAcceptor> searchModes;

        private readonly List<Protein> proteinList;

        private readonly List<ModificationWithMass> fixedModifications;

        private readonly List<ModificationWithMass> variableModifications;

        private readonly Psm[][] globalPsms;

        private readonly Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans;

        private readonly double[] myScanPrecursorMasses;

        private readonly List<ProductType> lp;

        private readonly bool addCompIons;

        private readonly CommonParameters CommonParameters;

        #endregion Private Fields

        #region Public Constructors

        public ClassicSearchEngine(Psm[][] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, List<ProductType> lp, List<MassDiffAcceptor> searchModes, bool addCompIons, CommonParameters CommonParameters, List<string> nestedIds) : base(nestedIds)
        {
            this.globalPsms = globalPsms;
            this.arrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            this.myScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.proteinList = proteinList;
            this.searchModes = searchModes;
            this.lp = lp;
            this.addCompIons = addCompIons;
            this.CommonParameters = CommonParameters;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("In classic search engine!", nestedIds);

            int totalProteins = proteinList.Count;

            var observedPeptides = new HashSet<CompactPeptide>();

            Status("Getting ms2 scans...", nestedIds);

            var lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;
            TerminusType terminusType = ProductTypeToTerminusType.IdentifyTerminusType(lp);
            Status("Starting classic search loop...", nestedIds);
            Parallel.ForEach(Partitioner.Create(0, totalProteins), partitionRange =>
            {
                var psms = new Psm[searchModes.Count][];
                for (int searchModeIndex = 0; searchModeIndex < searchModes.Count; searchModeIndex++)
                    psms[searchModeIndex] = new Psm[arrayOfSortedMS2Scans.Length];
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(CommonParameters.DigestionParams, fixedModifications).ToList();
                    foreach (var peptide in digestedList)
                    {
                        var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(CommonParameters.DigestionParams, variableModifications).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            var correspondingCompactPeptide = yyy.CompactPeptide(terminusType);
                            if (!CommonParameters.ConserveMemory)
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
                            double[] matchedIonMassesListPositiveIsMatch = new double[productMasses.Length];

                            for (int searchModeIndex = 0; searchModeIndex < searchModes.Count; searchModeIndex++)
                            {
                                var searchMode = searchModes[searchModeIndex];
                                foreach (ScanWithIndexAndNotchInfo scanWithIndexAndNotchInfo in GetAcceptableScans(correspondingCompactPeptide.MonoisotopicMassIncludingFixedMods, searchMode).ToList())
                                {
                                    double thePrecursorMass = scanWithIndexAndNotchInfo.theScan.PrecursorMass;
                                    var score = MatchIons(scanWithIndexAndNotchInfo.theScan.TheScan, CommonParameters.ProductMassTolerance, productMasses, matchedIonMassesListPositiveIsMatch, this.addCompIons, thePrecursorMass, lp);

                                    if (score > CommonParameters.ScoreCutoff)
                                    {
                                        if (psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex] == null)
                                            psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex] = new Psm(correspondingCompactPeptide, scanWithIndexAndNotchInfo.notch, score, scanWithIndexAndNotchInfo.scanIndex, scanWithIndexAndNotchInfo.theScan, CommonParameters.ExcelCompatible);
                                        else 
                                            psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex].AddOrReplace(correspondingCompactPeptide, score, scanWithIndexAndNotchInfo.notch, CommonParameters.ReportAllAmbiguity);
                                    }
                                }
                            }
                        }
                    }
                }
                lock (lockObject)
                {
                    for (int searchModeIndex = 0; searchModeIndex < searchModes.Count; searchModeIndex++)
                        for (int i = 0; i < globalPsms[searchModeIndex].Length; i++)
                            if (psms[searchModeIndex][i] != null)
                            {
                                if (globalPsms[searchModeIndex][i] == null)
                                    globalPsms[searchModeIndex][i] = psms[searchModeIndex][i];
                                else
                                {
                                    globalPsms[searchModeIndex][i].AddOrReplace(psms[searchModeIndex][i], CommonParameters.ReportAllAmbiguity);
                                }
                            }
                    proteinsSeen += partitionRange.Item2 - partitionRange.Item1;
                    var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In classic search loop", nestedIds));
                        old_progress = new_progress;
                    }
                }
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
                        yield return new ScanWithIndexAndNotchInfo(theScan, allowedIntervalWithNotch.notch, scanIndex);
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