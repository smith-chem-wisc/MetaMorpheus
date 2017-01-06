using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
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
            searchParams.outputAction("In classic search method!");

            int totalProteins = searchParams.proteinList.Count;

            HashSet<string> level3_observed = new HashSet<string>();
            HashSet<string> level4_observed = new HashSet<string>();

            var lp = new List<ProductType>() { ProductType.b, ProductType.y };

            var listOfScans = searchParams.myMsDataFile.ToList();

            var outerPsms = new ClassicSpectrumMatch[searchParams.myMsDataFile.NumSpectra];
            object lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;

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
                                var huh = yyy.FastSortedProductMasses(lp);
                                double[] hehe = new double[huh.Length];
                                ps.MonoisotopicMass = (float)yyy.MonoisotopicMass;

                                foreach (IMsDataScan<IMzSpectrum<MzPeak>> scan in listOfScans)
                                {
                                    if (scan.MsnOrder == 2)
                                    {
                                        double selectedMZ;
                                        scan.TryGetSelectedIonGuessMonoisotopicMZ(out selectedMZ);
                                        int selectedCharge;
                                        scan.TryGetSelectedIonGuessChargeStateGuess(out selectedCharge);
                                        var precursorMass = selectedMZ.ToMass(selectedCharge);
                                        double precursorIntensity;
                                        scan.TryGetSelectedIonGuessMonoisotopicIntensity(out precursorIntensity);
                                        if (searchParams.searchMode.Accepts(precursorMass - yyy.MonoisotopicMass))
                                        {
                                            var score = PSMwithPeptide.MatchIons(scan, searchParams.productMassTolerance, huh, hehe);
                                            ClassicSpectrumMatch psm = new ClassicSpectrumMatch(score, ps, hehe, precursorMass, selectedMZ, scan.OneBasedScanNumber, scan.RetentionTime,selectedCharge, scan.MassSpectrum.xArray.Length, scan.TotalIonCurrent, precursorIntensity, searchParams.spectraFileIndex);
                                            if (psm.score > 1)
                                            {
                                                ClassicSpectrumMatch current_best_psm = psms[scan.OneBasedScanNumber - 1];
                                                if (current_best_psm == null || ClassicSpectrumMatch.FirstIsPreferable(psm, current_best_psm))
                                                {
                                                    psms[scan.OneBasedScanNumber - 1] = psm;
                                                    hehe = new double[huh.Length];
                                                }
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
                    for (int i = 0; i < outerPsms.Length; i++)
                        if (psms[i] != null)
                            if (outerPsms[i] == null || ClassicSpectrumMatch.FirstIsPreferable(psms[i], outerPsms[i]))
                                outerPsms[i] = psms[i];
                    proteinsSeen += fff.Item2 - fff.Item1;
                    int new_progress = (int)((double)proteinsSeen / (totalProteins) * 100);
                    if (new_progress > old_progress)
                    {
                        searchParams.progressAction(new_progress);
                        old_progress = new_progress;
                    }
                }
            });

            List<NewPsm> newPsms = GetNewPsmsList(outerPsms);

            return new ClassicSearchResults(searchParams, newPsms);
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