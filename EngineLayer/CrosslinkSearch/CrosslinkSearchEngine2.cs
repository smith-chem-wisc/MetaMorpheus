using Chemistry;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkSearchEngine2 : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly List<int>[] fragmentIndex;

        private readonly float[] keys;

        private readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;

        private readonly List<CompactPeptide> peptideIndex;

        private readonly List<ProductType> lp;

        //Crosslink parameters
        private readonly CrosslinkerTypeClass crosslinker;

        private readonly int CrosslinkSearchTopNum;
        private readonly bool CrosslinkSearchWithCrosslinkerMod;
        private readonly Tolerance XLprecusorMsTl;
        private readonly Tolerance XLBetaPrecusorMsTl;

        private readonly Dictionary<ModificationWithMass, ushort> modsDictionary;
        private readonly List<Protein> proteinList;
        private readonly List<ModificationWithMass> variableModifications;
        private readonly List<ModificationWithMass> fixedModifications;

        private readonly ICommonParameters CommonParameters;

        private readonly List<PsmCross> psmCross;
        private MassDiffAcceptor XLBetaSearchMode;
        private MassDiffAcceptor XLPrecusorSearchMode;

        #endregion Private Fields

        #region Public Constructors

        public CrosslinkSearchEngine2(List<PsmCross> psmCross, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, CrosslinkerTypeClass crosslinker, int CrosslinkSearchTopNum, bool CrosslinkSearchWithCrosslinkerMod, Tolerance XLprecusorMsTl, Tolerance XLBetaPrecusorMsTl, Dictionary<ModificationWithMass, ushort> modsDictionary, List<ProductType> lp, List<Protein> proteinList, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, ICommonParameters CommonParameters, List<string> nestedIds) : base(nestedIds)
        {
            this.psmCross = psmCross;
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.keys = keys;
            this.fragmentIndex = fragmentIndex;

            this.crosslinker = crosslinker;
            this.CrosslinkSearchTopNum = CrosslinkSearchTopNum;
            this.CrosslinkSearchWithCrosslinkerMod = CrosslinkSearchWithCrosslinkerMod;
            this.XLprecusorMsTl = XLprecusorMsTl;
            this.XLBetaPrecusorMsTl = XLBetaPrecusorMsTl;
            //AnalysisEngine parameters
            this.modsDictionary = modsDictionary;
            this.lp = lp;
            this.proteinList = proteinList;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.CommonParameters = CommonParameters;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            MassDiffAcceptor XLsearchMode = new OpenSearchMode();
            XLBetaSearchMode = new SinglePpmAroundZeroSearchMode(XLBetaPrecusorMsTl.Value);
            if (XLBetaPrecusorMsTl.ToString().Contains("Absolute"))
            {
                XLBetaSearchMode = new SingleAbsoluteAroundZeroSearchMode(XLBetaPrecusorMsTl.Value);
            }
            XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(XLprecusorMsTl.Value);
            if (XLBetaPrecusorMsTl.ToString().Contains("Absolute"))
            {
                XLPrecusorSearchMode = new SingleAbsoluteAroundZeroSearchMode(XLprecusorMsTl.Value);
            }

            Status("In crosslink search engine...");

            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;

            //var searchModesCount = searchModes.Count;
            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;

            //Crosslink data storage: ListOfAllScan < ListOfCurrentScanTopPsms > [searchModesCount]
            List<PsmCross>[] newPsmsTop = new List<PsmCross>[listOfSortedms2Scans.Length];

            //Find Top matched peptides and then find Crosslinked peptides from them.
            //Parallel.ForEach(Partitioner.Create(0, 1), fff =>
            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2ScansLength), fff =>
            {
                List<BestPeptideScoreNotch> bestPeptideScoreNotch = new List<BestPeptideScoreNotch>();
                double worstScores = new double();

                double[] fullPeptideScores = new double[peptideIndexCount];
                //Find the Top matched peptides
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var thisScan = listOfSortedms2Scans[i];
                    var thisScanprecursorMass = thisScan.PrecursorMass;
                    Array.Clear(fullPeptideScores, 0, peptideIndexCount);
                    CalculatePeptideScores(thisScan.TheScan, fullPeptideScores);
                    bestPeptideScoreNotch.Clear();
                    worstScores = 0;

                    for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                    {
                        var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
                        CompactPeptide candidatePeptide = peptideIndex[possibleWinningPeptideIndex];

                        if (consideredScore >= 1)
                        {
                            // Check if makes sense to add due to peptidescore!

                            //currentWorstScore to mark the current worst score and peptide for comparation and removal.
                            double currentWorstScore = worstScores;
                            //From all scored peptides to choose the Top Num ones

                            if (bestPeptideScoreNotch != null && bestPeptideScoreNotch.Count == CrosslinkSearchTopNum)
                            {
                                // Full! Need to compare with old worst match
                                if (Math.Abs(currentWorstScore - consideredScore) < 1e-9)
                                {
                                    // Score is same as the worst, need to see if accepts and if prefer the new one
                                    int notch = XLsearchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                    if (notch >= 0)
                                    {
                                        bestPeptideScoreNotch.RemoveAt(CrosslinkSearchTopNum - 1);
                                        bestPeptideScoreNotch.Add(new BestPeptideScoreNotch(candidatePeptide, consideredScore, notch));
                                        bestPeptideScoreNotch = bestPeptideScoreNotch.OrderByDescending(p => p.BestScore).ToList();
                                        worstScores = bestPeptideScoreNotch.Last().BestScore;
                                    }
                                }
                                else if (currentWorstScore < consideredScore)
                                {
                                    // Score is better than the worst, only make sure it is acceptable
                                    int notch = XLsearchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                    if (notch >= 0)
                                    {
                                        bestPeptideScoreNotch.RemoveAt(CrosslinkSearchTopNum - 1);
                                        bestPeptideScoreNotch.Add(new BestPeptideScoreNotch(candidatePeptide, consideredScore, notch));
                                        bestPeptideScoreNotch = bestPeptideScoreNotch.OrderByDescending(p => p.BestScore).ToList();
                                        worstScores = bestPeptideScoreNotch.Last().BestScore;
                                    }
                                }
                            }
                            // Did not exist! Only make sure that it is acceptable.
                            else
                            {
                                int notch = XLsearchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                if (notch >= 0)
                                {
                                    if (bestPeptideScoreNotch == null)
                                    {
                                        bestPeptideScoreNotch = new List<BestPeptideScoreNotch>();
                                    }
                                    bestPeptideScoreNotch.Add(new BestPeptideScoreNotch(candidatePeptide, consideredScore, notch));
                                    bestPeptideScoreNotch = bestPeptideScoreNotch.OrderByDescending(p => p.BestScore).ToList();
                                    worstScores = bestPeptideScoreNotch.Last().BestScore;
                                }
                            }
                        }
                    }

                    //Create parameters to store current Crosslinked psm data?

                    if (bestPeptideScoreNotch != null)
                    {
                        List<PsmCross> currentScanPsmParent = new List<PsmCross>();
                        foreach (var bestPeptide in bestPeptideScoreNotch)
                        {
                            if (bestPeptide != null)
                            {
                                currentScanPsmParent.Add(new PsmCross(bestPeptide.BestPeptide, bestPeptide.BestNotch, bestPeptide.BestScore, i, thisScan));
                            }
                        }
                        for (int itop = 0; itop < currentScanPsmParent.Count; itop++)
                        {
                            currentScanPsmParent[itop].XlRank = new int[] { itop, 0 };
                        }
                        newPsmsTop[i] = currentScanPsmParent;
                    }
                }
                lock (outputObject)
                {
                    scansSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)scansSeen / (listOfSortedms2ScansLength) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In Crosslink search loop", nestedIds));
                        old_progress = new_progress;
                    }
                }
            });

            //Start search the Beta PsmCross
            List<PsmCross> allAlphaPsms = new List<PsmCross>();
            List<PsmCross> AllCrossPsms = new List<PsmCross>();
            List<Ms2ScanWithSpecificMass> selectedMS2Scans = new List<Ms2ScanWithSpecificMass>();
            for (int i = 0; i < newPsmsTop.Length; i++)
            {
                if (newPsmsTop[i] != null)
                {
                    for (int j = 0; j < newPsmsTop[i].Count; j++)
                    {
                        if (XLPrecusorSearchMode.Accepts(newPsmsTop[i][j].ScanPrecursorMass, newPsmsTop[i][j].compactPeptide.MonoisotopicMassIncludingFixedMods) >= 0)
                        {
                            newPsmsTop[i][j].XLTotalScore = newPsmsTop[i][j].Score;
                            newPsmsTop[i][j].CrossType = PsmCrossType.Singe;
                            AllCrossPsms.Add(newPsmsTop[i][j]);
                        }
                        if (XLPrecusorSearchMode.Accepts(newPsmsTop[i][j].ScanPrecursorMass, newPsmsTop[i][j].compactPeptide.MonoisotopicMassIncludingFixedMods + 156.0786) >= 0)
                        {
                            newPsmsTop[i][j].XLTotalScore = newPsmsTop[i][j].Score;
                            newPsmsTop[i][j].CrossType = PsmCrossType.DeadEnd;
                            AllCrossPsms.Add(newPsmsTop[i][j]);
                        }
                        if (XLPrecusorSearchMode.Accepts(newPsmsTop[i][j].ScanPrecursorMass, newPsmsTop[i][j].compactPeptide.MonoisotopicMassIncludingFixedMods + 138.06808) >= 0)
                        {
                            newPsmsTop[i][j].XLTotalScore = newPsmsTop[i][j].Score;
                            newPsmsTop[i][j].CrossType = PsmCrossType.Loop;
                            AllCrossPsms.Add(newPsmsTop[i][j]);
                        }
                        if (newPsmsTop[i][j].ScanPrecursorMass - newPsmsTop[i][j].compactPeptide.MonoisotopicMassIncludingFixedMods - crosslinker.TotalMass > 500)
                        {
                            allAlphaPsms.Add(newPsmsTop[i][j]);
                        }
                    }
                }
            }
            for (int i = 0; i < allAlphaPsms.Count; i++)
            {
                selectedMS2Scans.Add(listOfSortedms2Scans[allAlphaPsms[i].ScanIndex]);
            }
            var xlPsmCross = ClassSearchTheBetaPeptide(selectedMS2Scans.ToArray(), allAlphaPsms.ToArray(), CommonParameters.ConserveMemory);
            AllCrossPsms.AddRange(xlPsmCross);
            var AllCrossPsmsGroupOrder = (from c in AllCrossPsms
                                          group c by c.ScanNumber into grp
                                          select grp.OrderByDescending(c => c.XLTotalScore).ToList()).ToList();
            foreach (var item in AllCrossPsmsGroupOrder)
            {
                if (item.Count > 1)
                {
                    item[0].DScore = item[0].XLTotalScore - item[1].XLTotalScore;
                }
                else { item[0].DScore = item[0].XLTotalScore; }
            }
            var AllCrossPsmsRe = AllCrossPsmsGroupOrder.Select(p => p.First()).ToList();

            psmCross.AddRange(AllCrossPsmsRe);
            return new CrosslinkSearchResults(psmCross, this);
        }

        #endregion Protected Methods

        #region Private Methods

        private void CalculatePeptideScores(IMsDataScan<IMzSpectrum<IMzPeak>> spectrum, double[] peptideScores)
        {
            for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
            {
                var theAdd = 1 + spectrum.MassSpectrum.YArray[i] / spectrum.TotalIonCurrent;
                var experimentalPeakInDaltons = spectrum.MassSpectrum.XArray[i] - Constants.protonMass;
                float closestPeak;
                var ipos = Array.BinarySearch(keys, (float)experimentalPeakInDaltons);
                if (ipos < 0)
                    ipos = ~ipos;

                if (ipos > 0)
                {
                    var downIpos = ipos - 1;
                    // Try down
                    while (downIpos >= 0)
                    {
                        closestPeak = keys[downIpos];
                        if (CommonParameters.ProductMassTolerance.Within(experimentalPeakInDaltons, closestPeak))
                        {
                            foreach (var heh in fragmentIndex[downIpos])
                                peptideScores[heh] += theAdd;
                        }
                        else
                            break;
                        downIpos--;
                    }
                }
                if (ipos < keys.Length)
                {
                    var upIpos = ipos;
                    // Try here and up
                    while (upIpos < keys.Length)
                    {
                        closestPeak = keys[upIpos];
                        if (CommonParameters.ProductMassTolerance.Within(experimentalPeakInDaltons, closestPeak))
                        {
                            foreach (var heh in fragmentIndex[upIpos])
                                peptideScores[heh] += theAdd;
                        }
                        else
                            break;
                        upIpos++;
                    }
                }
            }
        }

        private List<PsmCross> ClassSearchTheBetaPeptide(Ms2ScanWithSpecificMass[] selectedScan, PsmCross[] selectedPsmParent, bool conserveMemory)
        {
            double[] selectedScanPrecusor = selectedScan.Select(p => p.PrecursorMass).ToArray();
            double[] AlphaPeptidePrecusor = selectedPsmParent.Select(p => p.compactPeptide.MonoisotopicMassIncludingFixedMods).ToArray();
            double[] BetaPeptidePrecusor = selectedScanPrecusor.Zip(AlphaPeptidePrecusor, (one, two) => one - two - crosslinker.TotalMass).ToArray();
            //Sort the BetaPeptidePrecusor which is necessary for the GetAcceptableScans
            Array.Sort(BetaPeptidePrecusor.ToArray(), selectedPsmParent);
            Array.Sort(BetaPeptidePrecusor, selectedScan);

            Status("In xlclassic search engine!");

            int totalProteins = proteinList.Count;

            var observed_sequences = new HashSet<CompactPeptide>();

            Status("Getting ms2 scans...");

            var outerPsms = new PsmCross[selectedScan.Length];

            var lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;
            TerminusType terminusType = ProductTypeMethod.IdentifyTerminusType(lp);
            Status("Starting xlclassic search loop...");
            //Parallel.ForEach(Partitioner.Create(0, 1), partitionRange =>
            //Parallel.ForEach(Partitioner.Create(0, totalProteins), partitionRange =>
            Parallel.ForEach(Partitioner.Create(0, totalProteins),
                new ParallelOptions { MaxDegreeOfParallelism = 1 }, partitionRange =>
            {
                var psms = new PsmCross[selectedScan.Length];
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                    foreach (var yyy in digestedList)
                    {
                        var correspondingCompactPeptide = yyy.CompactPeptide(terminusType);
                        if (!conserveMemory)
                        {
                            var observed = observed_sequences.Contains(correspondingCompactPeptide);
                            if (observed)
                                continue;
                            lock (observed_sequences)
                            {
                                observed = observed_sequences.Contains(correspondingCompactPeptide);
                                if (observed)
                                    continue;
                                observed_sequences.Add(correspondingCompactPeptide);
                            }
                        }

                        var productMasses = correspondingCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                        Array.Sort(productMasses);

                        foreach (ScanWithIndexAndNotchInfo scanWithIndexAndNotchInfo in GetAcceptableScans(BetaPeptidePrecusor, yyy.MonoisotopicMass, XLBetaSearchMode, selectedScan).ToList())
                        {
                            var score = CalculatePeptideScore(scanWithIndexAndNotchInfo.theScan.TheScan, CommonParameters.ProductMassTolerance, productMasses, yyy.MonoisotopicMass, new List<DissociationType>(), false);

                            if (score > 1 && PsmCross.xlPosCal(correspondingCompactPeptide, crosslinker).Count != 0)
                            {
                                var psm = new PsmCross(correspondingCompactPeptide, scanWithIndexAndNotchInfo.notch, score, scanWithIndexAndNotchInfo.scanIndex, scanWithIndexAndNotchInfo.theScan);
                                PsmCross.XLCalculateTotalProductMassesMightHave(scanWithIndexAndNotchInfo.theScan, psm, crosslinker, lp, CommonParameters.ProductMassTolerance);
                                double currentBestPsmLocalScore = 0;
                                if (psms[scanWithIndexAndNotchInfo.scanIndex] == null)
                                {
                                    psms[scanWithIndexAndNotchInfo.scanIndex] = psm;
                                    currentBestPsmLocalScore = psm.XLBestScore;
                                }
                                else
                                {
                                    var psmLocalScore = psm.XLBestScore;
                                    if (currentBestPsmLocalScore < psmLocalScore)
                                    {
                                        psms[scanWithIndexAndNotchInfo.scanIndex] = psm;
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
                        {
                            double outerPsmsLocalScore = 0;
                            if (outerPsms[i] == null)
                            {
                                outerPsms[i] = psms[i];
                                outerPsmsLocalScore = outerPsms[i].XLBestScore;
                            }
                            else
                            {
                                var psmLocalScore = psms[i].XLBestScore;
                                if (outerPsmsLocalScore < psmLocalScore)
                                {
                                    outerPsms[i] = psms[i];
                                }
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

            List<PsmCross> newPsmsTop = new List<PsmCross>();
            for (int i = 0; i < outerPsms.Length; i++)
            {
                if (outerPsms[i] != null && PsmCross.xlPosCal(selectedPsmParent[i].compactPeptide, crosslinker).Count != 0)
                {
                    PsmCross.XLCalculateTotalProductMassesMightHave(selectedScan[i], selectedPsmParent[i], crosslinker, lp, CommonParameters.ProductMassTolerance);
                    selectedPsmParent[i].XLTotalScore = selectedPsmParent[i].XLBestScore + outerPsms[i].XLBestScore;
                    selectedPsmParent[i].BetaPsmCross = outerPsms[i];
                    selectedPsmParent[i].CrossType = PsmCrossType.Cross;
                    newPsmsTop.Add(selectedPsmParent[i]);
                }
            }
            return newPsmsTop;
        }

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(double[] myScanPrecursorMasses, double peptideMonoisotopicMass, MassDiffAcceptor searchMode, Ms2ScanWithSpecificMass[] selectedScan)
        {
            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervals(peptideMonoisotopicMass).ToList())
            {
                DoubleRange allowedInterval = allowedIntervalWithNotch.allowedInterval;
                int scanIndex = GetFirstScanWithMassOverOrEqual(myScanPrecursorMasses, allowedInterval.Minimum);
                if (scanIndex < selectedScan.Length)
                {
                    var scanMass = myScanPrecursorMasses[scanIndex];
                    while (scanMass <= allowedInterval.Maximum)
                    {
                        var theScan = selectedScan[scanIndex];
                        yield return new ScanWithIndexAndNotchInfo(theScan, allowedIntervalWithNotch.notch, scanIndex);
                        scanIndex++;
                        if (scanIndex == selectedScan.Length)
                            break;
                        scanMass = myScanPrecursorMasses[scanIndex];
                    }
                }
            }
        }

        private int GetFirstScanWithMassOverOrEqual(double[] myScanPrecursorMasses, double minimum)
        {
            //could not use binarySearch, cause the precursorMass is not sorted.
            int index = Array.BinarySearch(myScanPrecursorMasses, minimum);
            if (index < 0)
                index = ~index;

            // index of the first element that is larger than value
            return index;
        }

        #endregion Private Methods
    }
}