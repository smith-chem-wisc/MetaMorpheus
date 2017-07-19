﻿using MzLibUtil;
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

        private const int max_mods_for_peptide = 3;

        private readonly int maximumMissedCleavages;
        private readonly int? minPeptideLength;
        private readonly int? maxPeptideLength;
        private readonly int maximumVariableModificationIsoforms;
        private readonly List<MassDiffAcceptor> searchModes;

        private readonly List<Protein> proteinList;

        private readonly Protease protease;

        private readonly List<ModificationWithMass> fixedModifications;

        private readonly List<ModificationWithMass> variableModifications;

        private readonly Tolerance productMassTolerance;

        private readonly Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans;

        private readonly double[] myScanPrecursorMasses;

        private readonly List<ProductType> lp;

        private readonly bool conserveMemory;

        private readonly bool addCompIons;

        #endregion Private Fields

        #region Public Constructors

        public ClassicSearchEngine(Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, Tolerance productMassTolerance, Protease protease, List<MassDiffAcceptor> searchModes, int maximumMissedCleavages, int? minPeptideLength, int? maxPeptideLength, int maximumVariableModificationIsoforms, List<ProductType> lp, List<string> nestedIds, bool conserveMemory, bool addCompIons) : base(nestedIds)
        {
            this.arrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            this.myScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.proteinList = proteinList;
            this.productMassTolerance = productMassTolerance;
            this.maximumMissedCleavages = maximumMissedCleavages;
            this.minPeptideLength = minPeptideLength;
            this.maxPeptideLength = maxPeptideLength;
            this.maximumVariableModificationIsoforms = maximumVariableModificationIsoforms;
            this.searchModes = searchModes;
            this.protease = protease;
            this.lp = lp;
            this.addCompIons = addCompIons;
            this.conserveMemory = conserveMemory;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("In classic search engine!", nestedIds);

            int totalProteins = proteinList.Count;

            //var observed_base_sequences = new HashSet<string>();
            var observed_sequences = new HashSet<string>();

            Status("Getting ms2 scans...", nestedIds);

            var outerPsms = new PsmParent[searchModes.Count][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                outerPsms[aede] = new PsmParent[arrayOfSortedMS2Scans.Length];

            var lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;

            Status("Starting classic search loop...", nestedIds);
            Parallel.ForEach(Partitioner.Create(0, totalProteins), partitionRange =>
            {
                var psms = new PsmParent[searchModes.Count][];
                for (int searchModeIndex = 0; searchModeIndex < searchModes.Count; searchModeIndex++)
                    psms[searchModeIndex] = new PsmParent[arrayOfSortedMS2Scans.Length];
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, InitiatorMethionineBehavior.Variable, fixedModifications).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length <= 1)
                            continue;

                        var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, maximumVariableModificationIsoforms, max_mods_for_peptide).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            if (!conserveMemory)
                            {
                                var hc = yyy.Sequence;
                                var observed = observed_sequences.Contains(hc);
                                if (observed)
                                    continue;
                                lock (observed_sequences)
                                {
                                    observed = observed_sequences.Contains(hc);
                                    if (observed)
                                        continue;
                                    observed_sequences.Add(hc);
                                }
                            }

                            var productMasses = yyy.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                            Array.Sort(productMasses);
                            double[] matchedIonMassesListPositiveIsMatch = new double[productMasses.Length];

                            for (int searchModeIndex = 0; searchModeIndex < searchModes.Count; searchModeIndex++)
                            {
                                var searchMode = searchModes[searchModeIndex];
                                foreach (ScanWithIndexAndNotchInfo scanWithIndexAndNotchInfo in GetAcceptableScans(yyy.MonoisotopicMass, searchMode).ToList())
                                {
                                    double thePrecursorMass=scanWithIndexAndNotchInfo.theScan.PrecursorMass;
                                    var score = PsmParent.MatchIons(scanWithIndexAndNotchInfo.theScan.TheScan, productMassTolerance, productMasses, matchedIonMassesListPositiveIsMatch, this.addCompIons, thePrecursorMass, lp);
                                    if (score > 1)
                                    {
                                        var psm = new PsmClassic(yyy, scanWithIndexAndNotchInfo.notch, score, scanWithIndexAndNotchInfo.scanIndex, scanWithIndexAndNotchInfo.theScan);
                                        var currentBestPsmList = psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex];
                                        if (currentBestPsmList == null)
                                            psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex] = psm;
                                        else
                                        {
                                            var singleIsPreferable = PsmClassic.FirstIsPreferable(psm, currentBestPsmList as PsmClassic, variableModifications);
                                            if (singleIsPreferable.HasValue && singleIsPreferable.Value)
                                                psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex] = psm;
                                            else if (!singleIsPreferable.HasValue)
                                                psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex].NumAmbiguous++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                lock (lockObject)
                {
                    for (int searchModeIndex = 0; searchModeIndex < searchModes.Count; searchModeIndex++)
                        for (int i = 0; i < outerPsms[searchModeIndex].Length; i++)
                            if (psms[searchModeIndex][i] != null)
                            {
                                if (outerPsms[searchModeIndex][i] == null)
                                    outerPsms[searchModeIndex][i] = psms[searchModeIndex][i];
                                else
                                {
                                    var firstIsPreferable = PsmClassic.FirstIsPreferable(psms[searchModeIndex][i] as PsmClassic, outerPsms[searchModeIndex][i] as PsmClassic, variableModifications);
                                    if (firstIsPreferable.HasValue && firstIsPreferable.Value)
                                        outerPsms[searchModeIndex][i] = psms[searchModeIndex][i];
                                    else if (!firstIsPreferable.HasValue)
                                        outerPsms[searchModeIndex][i].NumAmbiguous++;
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
            return new SearchResults(outerPsms, this);
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