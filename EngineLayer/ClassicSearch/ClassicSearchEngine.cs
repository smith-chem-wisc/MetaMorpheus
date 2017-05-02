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

        private const int max_mods_for_peptide = 3;

        private readonly int maximumMissedCleavages;
        private readonly int? minPeptideLength;
        private readonly int? maxPeptideLength;
        private readonly int maximumVariableModificationIsoforms;
        private readonly List<SearchMode> searchModes;

        private readonly List<Protein> proteinList;

        private readonly Protease protease;

        private readonly List<ModificationWithMass> fixedModifications;

        private readonly List<ModificationWithMass> variableModifications;

        private readonly Tolerance productMassTolerance;

        private readonly Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans;

        private readonly double[] myScanPrecursorMasses;

        private readonly List<ProductType> lp;
        private readonly List<string> nestedIds;

        private readonly bool conserveMemory;

        #endregion Private Fields

        #region Public Constructors

        public ClassicSearchEngine(Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, Tolerance productMassTolerance, Protease protease, List<SearchMode> searchModes, int maximumMissedCleavages, int? minPeptideLength, int? maxPeptideLength, int maximumVariableModificationIsoforms, List<ProductType> lp, List<string> nestedIds, bool conserveMemory)
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
            this.nestedIds = nestedIds;
            this.conserveMemory = conserveMemory;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("In classic search engine!", nestedIds);

            var searchResults = new ClassicSearchResults(this);

            int totalProteins = proteinList.Count;

            //var observed_base_sequences = new HashSet<string>();
            var observed_sequences = new HashSet<string>();

            Status("Getting ms2 scans...", nestedIds);

            var outerPsms = new List<PsmParent>[searchModes.Count][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                outerPsms[aede] = new List<PsmParent>[arrayOfSortedMS2Scans.Length];

            var lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;

            Status("Starting classic search loop...", nestedIds);
            Parallel.ForEach(Partitioner.Create(0, totalProteins), partitionRange =>
            {
                var psms = new List<PsmParent>[searchModes.Count][];
                for (int searchModeIndex = 0; searchModeIndex < searchModes.Count; searchModeIndex++)
                    psms[searchModeIndex] = new List<PsmParent>[arrayOfSortedMS2Scans.Length];
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, InitiatorMethionineBehavior.Variable, fixedModifications).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length <= 1)
                            continue;
                        //if (peptide.BaseLeucineSequence.Equals("GPK"))
                        //Console.WriteLine("Found GPK in protein " + protein.Accession + ". NumLocMods = " + peptide.NumKnownPossibleLocMods);
                        //if (peptide.NumKnownPossibleLocMods == 0 && !conserveMemory)
                        //{
                        //    var hc = peptide.BaseLeucineSequence;
                        //    if (peptide.BaseLeucineSequence.Equals("GPK"))
                        //        Console.WriteLine("In protein " + protein.Accession + " hc = " + hc);
                        //    var observed = observed_base_sequences.Contains(hc);
                        //    if (peptide.BaseLeucineSequence.Equals("GPK"))
                        //        Console.WriteLine("observed1 status of GPK in protein " + protein.Accession + " is " + observed);
                        //    if (observed)
                        //        continue;
                        //    lock (observed_base_sequences)
                        //    {
                        //        if (peptide.BaseLeucineSequence.Equals("GPK"))
                        //            Console.WriteLine("Locking in protein " + protein.Accession);
                        //        observed = observed_base_sequences.Contains(hc);
                        //        if (peptide.BaseLeucineSequence.Equals("GPK"))
                        //            Console.WriteLine("observed2 status of GPK in protein " + protein.Accession + " is " + observed);
                        //        if (observed)
                        //            continue;
                        //        if (peptide.BaseLeucineSequence.Equals("GPK"))
                        //            Console.WriteLine("Adding GPK in protein " + protein.Accession);
                        //        observed_base_sequences.Add(hc);
                        //    }
                        //}
                        //if (peptide.BaseLeucineSequence.Equals("GPK"))
                        //    Console.WriteLine("Ok, looking at it in protein " + protein.Accession);

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
                                    var scan = scanWithIndexAndNotchInfo.theScan;

                                    var score = PsmWithMultiplePossiblePeptides.MatchIons(scan.TheScan, productMassTolerance, productMasses, matchedIonMassesListPositiveIsMatch);
                                    var psm = new PsmClassic(yyy, scan, score, scanWithIndexAndNotchInfo.notch);
                                    if (psm.score > 1)
                                    {
                                        if (psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex] == null)
                                        {
                                            psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex] = new List<PsmParent> { psm };
                                            matchedIonMassesListPositiveIsMatch = new double[productMasses.Length];
                                        }
                                        else
                                        {
                                            PsmParent current_best_psm = psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex].First();
                                            if (PsmClassic.FirstIsPreferable(psm, current_best_psm as PsmClassic, variableModifications))
                                            {
                                                psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex] = new List<PsmParent> { psm };
                                                matchedIonMassesListPositiveIsMatch = new double[productMasses.Length];
                                            }
                                            else if (current_best_psm.score == psm.score)
                                            {
                                                psms[searchModeIndex][scanWithIndexAndNotchInfo.scanIndex].Add(psm);
                                                matchedIonMassesListPositiveIsMatch = new double[productMasses.Length];
                                            }
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
                                    if (PsmClassic.FirstIsPreferable(psms[searchModeIndex][i].First() as PsmClassic, outerPsms[searchModeIndex][i].First() as PsmClassic, variableModifications))
                                        outerPsms[searchModeIndex][i] = psms[searchModeIndex][i];
                                    else if (psms[searchModeIndex][i].First().score == outerPsms[searchModeIndex][i].First().score)
                                        outerPsms[searchModeIndex][i].AddRange(psms[searchModeIndex][i]);
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
            searchResults.OuterPsms = outerPsms;
            return searchResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(double peptideMonoisotopicMass, SearchMode searchMode)
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