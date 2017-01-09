using MetaMorpheus;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace IndexSearchAndAnalyze
{
    public class ClassicSearchEngine : MyEngine
    {
        public ClassicSearchEngine(ClassicSearchParams searchParams)
        {
            this.myParams = searchParams;
        }
        
        protected override MyResults RunSpecific()
        {
            var searchParams = (ClassicSearchParams)myParams;
            searchParams.allTasksParams.status("In classic search engine!");

            int totalProteins = searchParams.proteinList.Count;

            HashSet<string> level3_observed = new HashSet<string>();
            HashSet<string> level4_observed = new HashSet<string>();

            var lp = new List<ProductType>() { ProductType.b, ProductType.y };

            var listOfSortedms2Scans = searchParams.myMsDataFile.Where(b => b.MsnOrder == 2).Select(b => new LocalMs2Scan(b)).OrderBy(b => b.precursorMass).ToArray();

            var outerPsms = new ClassicSpectrumMatch[searchParams.myMsDataFile.NumSpectra];
            object lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;

            searchParams.allTasksParams.ReportProgress(new ProgressEventArgs(0, "Starting classic search loop"));
            Parallel.ForEach(Partitioner.Create(0, totalProteins), fff =>
            {
                var psms = new ClassicSpectrumMatch[searchParams.myMsDataFile.NumSpectra];
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var protein = searchParams.proteinList[i];
                    var digestedList = protein.Digest(searchParams.protease, 2, InitiatorMethionineBehavior.Variable).ToList();
                    foreach (var peptide in digestedList)
                    {
                        if (peptide.Length == 1 || peptide.Length > 252)
                            continue;

                        if (peptide.OneBasedPossibleLocalizedModifications.Count == 0)
                        {
                            lock (level3_observed)
                            {
                                var hc = peptide.BaseLeucineSequence;
                                var observed = level3_observed.Contains(hc);
                                if (observed)
                                    continue;
                                level3_observed.Add(hc);
                            }
                        }

                        peptide.SetFixedModifications(searchParams.fixedModifications);

                        var ListOfModifiedPeptides = peptide.GetPeptideWithSetModifications(searchParams.variableModifications, 4098, 3, searchParams.localizeableModifications).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            if (peptide.OneBasedPossibleLocalizedModifications.Count > 0)
                            {
                                lock (level4_observed)
                                {
                                    var hc = yyy.Sequence;
                                    var observed = level4_observed.Contains(hc);
                                    if (observed)
                                        continue;
                                    level4_observed.Add(hc);
                                }

                                var ps = new CompactPeptide(yyy, searchParams.variableModifications, searchParams.localizeableModifications);
                                var sortedProductMasses = yyy.FastSortedProductMasses(lp);
                                double[] matchedIonsArray = new double[sortedProductMasses.Length];
                                ps.MonoisotopicMass = (float)yyy.MonoisotopicMass;

                                foreach (LocalMs2Scan scan in GetAcceptableScans(listOfSortedms2Scans, yyy.MonoisotopicMass, searchParams.searchMode).ToList())
                                {
                                    var score = PSMwithPeptide.MatchIons(scan.b, searchParams.productMassTolerance, sortedProductMasses, matchedIonsArray);
                                    ClassicSpectrumMatch psm = new ClassicSpectrumMatch(score, ps, matchedIonsArray, scan.precursorMass, scan.monoisotopicPrecursorMZ, scan.OneBasedScanNumber, scan.RetentionTime, scan.monoisotopicPrecursorCharge, scan.NumPeaks, scan.TotalIonCurrent, scan.monoisotopicPrecursorIntensity, searchParams.spectraFileIndex);
                                    if (psm.score > 1)
                                    {
                                        ClassicSpectrumMatch current_best_psm = psms[scan.OneBasedScanNumber - 1];
                                        if (current_best_psm == null || ClassicSpectrumMatch.FirstIsPreferable(psm, current_best_psm))
                                        {
                                            psms[scan.OneBasedScanNumber - 1] = psm;
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
                    for (int i = 0; i < outerPsms.Length; i++)
                        if (psms[i] != null)
                            if (outerPsms[i] == null || ClassicSpectrumMatch.FirstIsPreferable(psms[i], outerPsms[i]))
                                outerPsms[i] = psms[i];
                    proteinsSeen += fff.Item2 - fff.Item1;
                    int new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                    if (new_progress > old_progress)
                    {
                        searchParams.allTasksParams.ReportProgress(new ProgressEventArgs(new_progress, "In classic search loop"));
                        old_progress = new_progress;
                    }
                }
            });

            List<NewPsm> newPsms = GetNewPsmsList(outerPsms);

            return new ClassicSearchResults(searchParams, newPsms);
        }

        private IEnumerable<LocalMs2Scan> GetAcceptableScans(LocalMs2Scan[] listOfSortedms2Scans, double peptideMonoisotopicMass, SearchMode searchMode)
        {
            IntervalSearchMode intervalSearchMode = searchMode as IntervalSearchMode;
            if (intervalSearchMode != null)
            {
                foreach (DoubleRange ye in intervalSearchMode.GetAllowedPrecursorMassIntervals(peptideMonoisotopicMass).ToList())
                {
                    int scanIndex = GetFirstScanWithMassOverOrEqual(listOfSortedms2Scans, ye.Minimum);
                    var scan = listOfSortedms2Scans[scanIndex];
                    while (scan.precursorMass <= ye.Maximum)
                    {
                        yield return scan;
                        scanIndex++;
                        scan = listOfSortedms2Scans[scanIndex];
                    }
                }
            }
            FunctionSearchMode functionSearchMode = searchMode as FunctionSearchMode;
            if (functionSearchMode != null)
            {
                foreach (var scan in listOfSortedms2Scans)
                {
                    if (functionSearchMode.Accepts(scan.precursorMass, peptideMonoisotopicMass))
                        yield return scan;
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

        private List<NewPsm> GetNewPsmsList(ClassicSpectrumMatch[] outerPsms)
        {
            List<NewPsm> newPsms = new List<NewPsm>(outerPsms.Length);
            foreach (var classicPSM in outerPsms)
            {
                if (classicPSM != null)
                {
                    newPsms.Add(new NewPsm(classicPSM.scanPrecursorMZ, classicPSM.scanNumber, classicPSM.scanRT, classicPSM.scanPrecursorCharge, classicPSM.scanExperimentalPeaks, classicPSM.TotalIonCurrent, classicPSM.scanPrecursorIntensity, classicPSM.spectraFileIndex, classicPSM.ps, classicPSM.score));
                }
            }
            return newPsms;
        }
    }
}