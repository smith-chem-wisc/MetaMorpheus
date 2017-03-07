using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.ClassicSearch
{
    public class ClassicSearchEngine : MyEngine
    {

        #region Private Fields

        private const int max_mods_for_peptide = 3;

        private readonly int maximumMissedCleavages;
        private readonly int maximumVariableModificationIsoforms;
        private readonly List<SearchMode> searchModes;

        private readonly List<Protein> proteinList;

        private readonly Protease protease;

        private readonly List<ModificationWithMass> fixedModifications;

        private readonly List<ModificationWithMass> variableModifications;

        private readonly Tolerance productMassTolerance;

        private readonly LocalMS2Scan[] arrayOfSortedMS2Scans;

        private readonly double[] myScanPrecursorMasses;

        private readonly int myMsDataFileNumSpectra;
        private readonly string fileName;

        private readonly List<ProductType> lp;
        private readonly List<string> nestedIds;

        #endregion Private Fields

        #region Public Constructors

        public ClassicSearchEngine(LocalMS2Scan[] arrayOfSortedMS2Scans, int myMsDataFileNumSpectra, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, List<Protein> proteinList, Tolerance productMassTolerance, Protease protease, List<SearchMode> searchModes, int maximumMissedCleavages, int maximumVariableModificationIsoforms, string fileName, List<ProductType> lp, List<string> nestedIds)
        {
            this.arrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            this.myScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.MonoisotopicPrecursorMass).ToArray();
            this.myMsDataFileNumSpectra = myMsDataFileNumSpectra;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.proteinList = proteinList;
            this.productMassTolerance = productMassTolerance;
            this.maximumMissedCleavages = maximumMissedCleavages;
            this.maximumVariableModificationIsoforms = maximumVariableModificationIsoforms;
            this.searchModes = searchModes;
            this.protease = protease;
            this.fileName = fileName;
            this.lp = lp;
            this.nestedIds = nestedIds;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            Status("In classic search engine!", nestedIds);

            var searchResults = new ClassicSearchResults(this);

            int totalProteins = proteinList.Count;

            var level3_observed = new HashSet<string>();
            var level4_observed = new HashSet<string>();

            Status("Getting ms2 scans...", nestedIds);

            var outerPsms = new PsmClassic[searchModes.Count][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                outerPsms[aede] = new PsmClassic[myMsDataFileNumSpectra];

            var lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;

            Status("Starting classic search loop...", nestedIds);
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                var psms = new PsmClassic[searchModes.Count][];
                for (int aede = 0; aede < searchModes.Count; aede++)
                    psms[aede] = new PsmClassic[myMsDataFileNumSpectra];
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, InitiatorMethionineBehavior.Variable, fixedModifications).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length <= 1)
                            continue;

                        if (peptide.numLocMods == 0)
                        {
                            var hc = peptide.BaseLeucineSequence;
                            var observed = level3_observed.Contains(hc);
                            if (observed)
                                continue;
                            lock (level3_observed)
                            {
                                observed = level3_observed.Contains(hc);
                                if (observed)
                                    continue;
                                level3_observed.Add(hc);
                            }
                        }

                        var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, maximumVariableModificationIsoforms, max_mods_for_peptide).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            if (peptide.numLocMods > 0)
                            {
                                var hc = yyy.Sequence;
                                var observed = level4_observed.Contains(hc);
                                if (observed)
                                    continue;
                                lock (level4_observed)
                                {
                                    observed = level4_observed.Contains(hc);
                                    if (observed)
                                        continue;
                                    level4_observed.Add(hc);
                                }
                            }

                            var sortedProductMasses = yyy.FastSortedProductMasses(lp);
                            double[] matchedIonsArray = new double[sortedProductMasses.Length];

                            for (int aede = 0; aede < searchModes.Count; aede++)
                            {
                                var searchMode = searchModes[aede];
                                foreach (Tuple<LocalMS2Scan, int> theTuple in GetAcceptableScans(yyy.MonoisotopicMass, searchMode).ToList())
                                {
                                    var scan = theTuple.Item1;

                                    //if(scan.TheScan.DissociationType == MassSpectrometry.DissociationType.ETD)
                                    //{
                                    //
                                    //}

                                    var score = PsmWithMultiplePossiblePeptides.MatchIons(scan.TheScan, productMassTolerance, sortedProductMasses, matchedIonsArray);
                                    var psm = new PsmClassic(yyy, fileName, scan.RetentionTime, scan.MonoisotopicPrecursorIntensity, scan.MonoisotopicPrecursorMass, scan.OneBasedScanNumber, scan.PrecursorCharge, scan.NumPeaks, scan.TotalIonCurrent, scan.MonoisotopicPrecursorMZ, score, theTuple.Item2);
                                    if (psm.score > 1)
                                    {
                                        PsmClassic current_best_psm = psms[aede][scan.OneBasedScanNumber - 1];
                                        if (current_best_psm == null || PsmClassic.FirstIsPreferable(psm, current_best_psm))
                                        {
                                            psms[aede][scan.OneBasedScanNumber - 1] = psm;
                                            matchedIonsArray = new double[sortedProductMasses.Length];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                lock (lockObject)
                {
                    for (int aede = 0; aede < searchModes.Count; aede++)
                        for (int i = 0; i < outerPsms[aede].Length; i++)
                            if (psms[aede][i] != null)
                                if (outerPsms[aede][i] == null || PsmClassic.FirstIsPreferable(psms[aede][i], outerPsms[aede][i]))
                                    outerPsms[aede][i] = psms[aede][i];
                    proteinsSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In classic search loop"));
                        old_progress = new_progress;
                    }
                }
            });
            searchResults.OuterPsms = outerPsms;
            return searchResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private IEnumerable<Tuple<LocalMS2Scan, int>> GetAcceptableScans(double peptideMonoisotopicMass, SearchMode searchMode)
        {
            foreach (Tuple<DoubleRange, int> theTuple in searchMode.GetAllowedPrecursorMassIntervals(peptideMonoisotopicMass).ToList())
            {
                DoubleRange ye = theTuple.Item1;
                int scanIndex = GetFirstScanWithMassOverOrEqual(ye.Minimum);
                if (scanIndex < arrayOfSortedMS2Scans.Length)
                {
                    var scan = arrayOfSortedMS2Scans[scanIndex];
                    while (scan.MonoisotopicPrecursorMass <= ye.Maximum)
                    {
                        yield return new Tuple<LocalMS2Scan, int>(scan, theTuple.Item2);
                        scanIndex++;
                        if (scanIndex == arrayOfSortedMS2Scans.Length)
                            break;
                        scan = arrayOfSortedMS2Scans[scanIndex];
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