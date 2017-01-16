using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace InternalLogicEngineLayer
{
    public class ClassicSearchEngine : MyEngine
    {

        #region Private Fields

        private const int maximumMissedCleavages = 2;
        private const int maximumVariableModificationIsoforms = 4096;
        private const int max_mods_for_peptide = 3;
        private readonly List<SearchMode> searchModes;

        private readonly List<Protein> proteinList;

        private readonly Protease protease;

        private readonly List<MorpheusModification> fixedModifications;

        private readonly List<MorpheusModification> variableModifications;

        private readonly Tolerance productMassTolerance;

        private readonly LocalMs2Scan[] myMsDataFile;

        private readonly int spectraFileIndex;
        private readonly int myMsDataFileNumSpectra;

        #endregion Private Fields

        #region Public Constructors

        public ClassicSearchEngine(LocalMs2Scan[] myMsDataFile, int myMsDataFileNumSpectra, int spectraFileIndex, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<Protein> proteinList, Tolerance productMassTolerance, Protease protease, List<SearchMode> searchModes) : base(2)
        {
            this.myMsDataFile = myMsDataFile;
            this.myMsDataFileNumSpectra = myMsDataFileNumSpectra;
            this.spectraFileIndex = spectraFileIndex;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.proteinList = proteinList;
            this.productMassTolerance = productMassTolerance;
            this.protease = protease;
            this.searchModes = searchModes;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            status("In classic search engine!");

            var searchResults = new ClassicSearchResults(this);

            int totalProteins = proteinList.Count;

            var level3_observed = new HashSet<string>();
            var level4_observed = new HashSet<string>();

            var lp = new List<ProductType> { ProductType.b, ProductType.y };

            status("Getting ms2 scans...");

            var outerPsms = new ClassicSpectrumMatch[searchModes.Count][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                outerPsms[aede] = new ClassicSpectrumMatch[myMsDataFileNumSpectra];

            var lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;

            status("Starting classic search loop...");
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                var psms = new ClassicSpectrumMatch[searchModes.Count][];
                for (int aede = 0; aede < searchModes.Count; aede++)
                    psms[aede] = new ClassicSpectrumMatch[myMsDataFileNumSpectra];
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, InitiatorMethionineBehavior.Variable).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length == 1 || peptide.Length > byte.MaxValue - 2)
                            continue;

                        if (peptide.OneBasedPossibleLocalizedModifications.Count == 0)
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

                        peptide.SetFixedModifications(fixedModifications);

                        var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(variableModifications, maximumVariableModificationIsoforms, max_mods_for_peptide).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            if (peptide.OneBasedPossibleLocalizedModifications.Count > 0)
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
                                foreach (LocalMs2Scan scan in GetAcceptableScans(myMsDataFile, yyy.MonoisotopicMass, searchMode).ToList())
                                {
                                    var score = PSMwithTargetDecoyKnown.MatchIons(scan.theScan, productMassTolerance, sortedProductMasses, matchedIonsArray);
                                    var psm = new ClassicSpectrumMatch(score, yyy, scan.precursorMass, scan.monoisotopicPrecursorMZ, scan.OneBasedScanNumber, scan.RetentionTime, scan.monoisotopicPrecursorCharge, scan.NumPeaks, scan.TotalIonCurrent, scan.monoisotopicPrecursorIntensity, spectraFileIndex);
                                    if (psm.Score > 1)
                                    {
                                        ClassicSpectrumMatch current_best_psm = psms[aede][scan.OneBasedScanNumber - 1];
                                        if (current_best_psm == null || ClassicSpectrumMatch.FirstIsPreferable(psm, current_best_psm))
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
                                if (outerPsms[aede][i] == null || ClassicSpectrumMatch.FirstIsPreferable(psms[aede][i], outerPsms[aede][i]))
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
            searchResults.outerPsms = outerPsms;
            return searchResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private IEnumerable<LocalMs2Scan> GetAcceptableScans(LocalMs2Scan[] listOfSortedms2Scans, double peptideMonoisotopicMass, SearchMode searchMode)
        {
            foreach (DoubleRange ye in searchMode.GetAllowedPrecursorMassIntervals(peptideMonoisotopicMass).ToList())
            {
                int scanIndex = GetFirstScanWithMassOverOrEqual(listOfSortedms2Scans, ye.Minimum);
                if (scanIndex < listOfSortedms2Scans.Length)
                {
                    var scan = listOfSortedms2Scans[scanIndex];
                    while (scan.precursorMass <= ye.Maximum)
                    {
                        yield return scan;
                        scanIndex++;
                        if (scanIndex == listOfSortedms2Scans.Length)
                            break;
                        scan = listOfSortedms2Scans[scanIndex];
                    }
                }
            }
        }

        private int GetFirstScanWithMassOverOrEqual(LocalMs2Scan[] listOfSortedms2Scans, double minimum)
        {
            int index = Array.BinarySearch(listOfSortedms2Scans, minimum);
            if (index < 0)
                index = ~index;

            // index of the first element that is larger than value
            return index;
        }

        #endregion Private Methods

    }
}