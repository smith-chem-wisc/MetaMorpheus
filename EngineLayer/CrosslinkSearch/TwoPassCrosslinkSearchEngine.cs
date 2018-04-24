using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.CrosslinkSearch
{
    public class TwoPassCrosslinkSearchEngine : MetaMorpheusEngine
    {
        #region Protected Fields

        protected const int fragmentBinsPerDalton = 1000;
        protected readonly List<int>[] fragmentIndex;
        protected readonly PeptideSpectralMatch[] globalPsms;
        protected readonly List<PsmCross> globalPsmsCross;
        protected readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;
        protected readonly List<CompactPeptide> peptideIndex;
        protected readonly List<ProductType> lp;
        protected readonly int currentPartition;
        protected readonly ICommonParameters CommonParameters;
        protected readonly bool addCompIons;
        protected readonly MassDiffAcceptor massDiffAcceptor;
        protected readonly List<DissociationType> dissociationTypes;

        #endregion Protected Fields

        #region Private Fields

        //Crosslink parameters
        private readonly CrosslinkerTypeClass crosslinker;

        private readonly bool CrosslinkSearchTop;
        private readonly int CrosslinkSearchTopNum;

        //private readonly bool CrosslinkSearchWithCrosslinkerMod;
        private readonly Tolerance XLPrecusorMsTl;

        //private readonly Tolerance XLBetaPrecusorMsTl;
        private readonly bool quench_H2O;

        private readonly bool quench_NH2;
        private readonly bool quench_Tris;
        private readonly bool charge_2_3;
        private readonly bool charge_2_3_PrimeFragment;
        private MassDiffAcceptor XLPrecusorSearchMode;

        #endregion Private Fields

        #region Public Constructors

        public TwoPassCrosslinkSearchEngine(List<PsmCross> globalPsmsCross, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition, ICommonParameters CommonParameters, bool addCompIons, Tolerance XLPrecusorMsTl, CrosslinkerTypeClass crosslinker, bool CrosslinkSearchTop, int CrosslinkSearchTopNum, bool quench_H2O, bool quench_NH2, bool quench_Tris, bool charge_2_3, bool charge_2_3_PrimeFragment, List<string> nestedIds) : base(nestedIds)
        {
            this.globalPsmsCross = globalPsmsCross;
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.fragmentIndex = fragmentIndex;
            this.lp = lp;
            this.currentPartition = currentPartition + 1;
            this.CommonParameters = CommonParameters;
            this.addCompIons = addCompIons;
            //Here use LowTheoreticalDiffAcceptor in practice doesn't work in 12/2/2017
            this.massDiffAcceptor = new OpenSearchMode();
            this.dissociationTypes = DetermineDissociationType(lp);
            this.XLPrecusorMsTl = XLPrecusorMsTl;
            XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(XLPrecusorMsTl.Value);
            //if (XLBetaPrecusorMsTl.ToString().Contains("Absolute"))
            //{
            //    XLPrecusorSearchMode = new SingleAbsoluteAroundZeroSearchMode(XLPrecusorMsTl.Value);
            //}
            this.crosslinker = crosslinker;
            this.CrosslinkSearchTop = CrosslinkSearchTop;
            this.CrosslinkSearchTopNum = CrosslinkSearchTopNum;
            this.quench_H2O = quench_H2O;
            this.quench_NH2 = quench_NH2;
            this.quench_Tris = quench_Tris;
            this.charge_2_3 = charge_2_3;
            this.charge_2_3_PrimeFragment = charge_2_3_PrimeFragment;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing crosslink search... " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;

            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile }, range =>
            {
                byte[] scoringTable = new byte[peptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = listOfSortedms2Scans[i];

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(scan);
                    List<BestPeptideScoreNotch> bestPeptideScoreNotchList = new List<BestPeptideScoreNotch>();
                    // get allowed theoretical masses from the known experimental mass
                    // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)
                    // this is just PRELIMINARY precursor-mass filtering
                    // additional checks are made later to ensure that the theoretical precursor mass is acceptable
                    var notches = massDiffAcceptor.GetAllowedPrecursorMassIntervals(scan.PrecursorMass);

                    double lowestMassPeptideToLookFor = Double.NegativeInfinity;
                    double highestMassPeptideToLookFor = Double.PositiveInfinity;

                    double largestMassDiff = notches.Max(p => p.allowedInterval.Maximum);
                    double smallestMassDiff = notches.Min(p => p.allowedInterval.Minimum);

                    if (!Double.IsInfinity(largestMassDiff))
                    {
                        double largestOppositeMassDiff = -1 * (notches.Max(p => p.allowedInterval.Maximum) - scan.PrecursorMass);
                        lowestMassPeptideToLookFor = scan.PrecursorMass + largestOppositeMassDiff;
                    }
                    if (!Double.IsNegativeInfinity(smallestMassDiff))
                    {
                        double smallestOppositeMassDiff = -1 * (notches.Min(p => p.allowedInterval.Minimum) - scan.PrecursorMass);
                        highestMassPeptideToLookFor = scan.PrecursorMass + smallestOppositeMassDiff;
                    }

                    // first-pass scoring
                    IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor);

                    // done with indexed scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        List<int> idsRankedByScore = new List<int>();
                        foreach (var id in idsOfPeptidesPossiblyObserved)
                        {
                            idsRankedByScore.Add(id);
                        }

                        idsRankedByScore = idsRankedByScore.OrderByDescending(p => scoringTable[p]).ToList();
                        idsRankedByScore = idsRankedByScore.Take(CrosslinkSearchTopNum).ToList();
                        if (CrosslinkSearchTop)
                        {
                            idsOfPeptidesPossiblyObserved = idsRankedByScore;
                        }
                        foreach (var id in idsOfPeptidesPossiblyObserved)
                        {
                            var peptide = peptideIndex[id];

                            double thePrecursorMass = scan.PrecursorMass;
                            int notch = massDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMassIncludingFixedMods);
                            BestPeptideScoreNotch bestPeptideScoreNotch = new BestPeptideScoreNotch(peptide, scoringTable[id], notch);
                            bestPeptideScoreNotchList.Add(bestPeptideScoreNotch);
                        }

                        var possiblePsmCross = FindCrosslinkedPeptide(scan, bestPeptideScoreNotchList, i);
                        if (possiblePsmCross != null)
                        {
                            lock(globalPsmsCross)
                                globalPsmsCross.Add(possiblePsmCross);
                        }
                    }
                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / listOfSortedms2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing crosslink search... " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));
                    }
                }
            });

            return new MetaMorpheusEngineResults(this);
        }

        protected List<int> GetBinsToSearch(Ms2ScanWithSpecificMass scan)
        {
            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();
            foreach (var peakMz in scan.TheScan.MassSpectrum.XArray)
            {
                // assume charge state 1 to calculate mass tolerance
                double experimentalFragmentMass = Chemistry.ClassExtensions.ToMass(peakMz, 1);

                // get theoretical fragment bins within mass tolerance
                int obsFragmentFloorMass = (int)Math.Floor((CommonParameters.ProductMassTolerance.GetMinimumValue(experimentalFragmentMass)) * fragmentBinsPerDalton);
                if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)
                    obsFragmentFloorMass = obsPreviousFragmentCeilingMz;
                int obsFragmentCeilingMass = (int)Math.Ceiling((CommonParameters.ProductMassTolerance.GetMaximumValue(experimentalFragmentMass)) * fragmentBinsPerDalton);
                obsPreviousFragmentCeilingMz = obsFragmentCeilingMass + 1;
                for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                    if (fragmentIndex[fragmentBin] != null)
                        binsToSearch.Add(fragmentBin);
            }
            return binsToSearch;
        }

        #endregion Protected Methods

        #region Private Methods

        private int BinarySearchBinForPrecursorIndex(List<int> peptideIdsInThisBin, double peptideMassToLookFor)
        {
            int m = 0;
            int l = 0;
            int r = peptideIdsInThisBin.Count - 1;

            // binary search in the fragment bin for precursor mass
            while (l <= r)
            {
                m = l + ((r - l) / 2);

                if (r - l < 2)
                    break;
                if (peptideIndex[peptideIdsInThisBin[m]].MonoisotopicMassIncludingFixedMods < peptideMassToLookFor)
                    l = m + 1;
                else
                    r = m - 1;
            }
            if (m > 0)
                m--;
            return m;
        }

        private void IndexedScoring(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = fragmentIndex[binsToSearch[i]];

                //get index for minimum monoisotopic allowed
                int lowestPeptideMassIndex = Double.IsInfinity(lowestMassPeptideToLookFor) ? 0 : BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, lowestMassPeptideToLookFor);

                // get index for highest mass allowed
                int highestPeptideMassIndex = peptideIdsInThisBin.Count - 1;

                if (!Double.IsInfinity(highestMassPeptideToLookFor))
                {
                    highestPeptideMassIndex = BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, highestMassPeptideToLookFor);

                    for (int j = highestPeptideMassIndex; j < peptideIdsInThisBin.Count; j++)
                    {
                        int nextId = peptideIdsInThisBin[j];
                        var nextPep = peptideIndex[nextId];
                        if (nextPep.MonoisotopicMassIncludingFixedMods < highestMassPeptideToLookFor)
                            highestPeptideMassIndex = j;
                        else
                            break;
                    }
                }

                // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                {
                    int id = peptideIdsInThisBin[j];
                    scoringTable[id]++;

                    // add possible search results to the hashset of id's
                    if (scoringTable[id] == byteScoreCutoff && massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMassIncludingFixedMods) >= 0)
                        idsOfPeptidesPossiblyObserved.Add(id);
                }
            }
        }

        //Targetting function: to find two peptides that in the Top matched peptides
        private PsmCross FindCrosslinkedPeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int i)
        {
            List<PsmCross> bestPsmCrossList = new List<PsmCross>();
            PsmCross bestPsmCross = null;
            string crosslinkerModSitesAll = new string((crosslinker.CrosslinkerModSites + crosslinker.CrosslinkerModSites2).ToCharArray().Distinct().ToArray());
            for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
            {
                //Single Peptide
                if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods) >= 0)
                {
                    var productMasses = theScanBestPeptide[ind].BestPeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                    Array.Sort(productMasses);
                    double score = CalculatePeptideScore(theScan.TheScan, CommonParameters.ProductMassTolerance, productMasses, theScan.PrecursorMass, dissociationTypes, addCompIons, 0);
                    var psmCrossSingle = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, score, i, theScan);
                    psmCrossSingle.XLTotalScore = psmCrossSingle.Score;
                    psmCrossSingle.CrossType = PsmCrossType.Singe;

                    bestPsmCrossList.Add(psmCrossSingle);
                }
                //Deadend Peptide
                else if (quench_Tris && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.DeadendMassTris) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan);              
                    var xlPos = PsmCross.XlPosCal(psmCrossEnd.compactPeptide, crosslinkerModSitesAll);
                    if (xlPos.Count() >= 1)
                    {
                        PsmCross.XlLocalization(theScan, psmCrossEnd, crosslinker.DeadendMassTris, crosslinker, lp, CommonParameters.ProductMassTolerance, false, false, xlPos);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.XLBestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndTris;
                        bestPsmCrossList.Add(psmCrossEnd);
                    }              
                }
                else if (quench_H2O  && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.DeadendMassH2O) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan);
                    var xlPos = PsmCross.XlPosCal(psmCrossEnd.compactPeptide, crosslinkerModSitesAll);
                    if (xlPos.Count() >= 1)
                    {
                        PsmCross.XlLocalization(theScan, psmCrossEnd, crosslinker.DeadendMassH2O, crosslinker, lp, CommonParameters.ProductMassTolerance, false, false, xlPos);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.XLBestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndH2O;
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                else if (quench_NH2  && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.DeadendMassNH2) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan);
                    var xlPos = PsmCross.XlPosCal(psmCrossEnd.compactPeptide, crosslinkerModSitesAll);
                    if (xlPos.Count() >= 1)
                    {
                        PsmCross.XlLocalization(theScan, psmCrossEnd, crosslinker.DeadendMassNH2, crosslinker, lp, CommonParameters.ProductMassTolerance, false, false, xlPos);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.XLBestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndNH2;
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                //loop peptide
                else if (crosslinker.LoopMass!=0 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.LoopMass) >= 0)
                {
                    var psmCrossLoop = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan);
                    var xlPos = PsmCross.XlPosCal(psmCrossLoop.compactPeptide, crosslinker.CrosslinkerModSites);
                    if (xlPos.Count() >= 2)
                    {
                        PsmCross.XlLocalizationForLoopCrosslink(theScan, psmCrossLoop, crosslinker.LoopMass, crosslinker, lp, CommonParameters.ProductMassTolerance, xlPos);
                        psmCrossLoop.XLTotalScore = psmCrossLoop.XLBestScore;
                        psmCrossLoop.CrossType = PsmCrossType.Loop;
                        bestPsmCrossList.Add(psmCrossLoop);
                    }                   
                }
                //Cross-linked peptide
                else if (theScan.PrecursorMass - theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods >= 200)
                {
                    var x = theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods;
                    for (int inx = ind; inx < theScanBestPeptide.Count; inx++)
                    {
                        var y = theScanBestPeptide[inx].BestPeptide.MonoisotopicMassIncludingFixedMods;
                        if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, x + y + crosslinker.TotalMass) >= 0)
                        {
                            var psmCrossAlpha = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan);
                            var psmCrossBeta = new PsmCross(theScanBestPeptide[inx].BestPeptide, theScanBestPeptide[inx].BestNotch, theScanBestPeptide[inx].BestScore, i, theScan);

                            PsmCross psmCross = null;
                            if (crosslinker.CrosslinkerModSites == crosslinker.CrosslinkerModSites2)
                            {
                                psmCross = XlSitesSame(theScan, psmCrossAlpha, psmCrossBeta, crosslinker, ind, inx);
                            }
                            else
                            {
                                psmCross = XlSitesTwoDiff(theScan, psmCrossAlpha, psmCrossBeta, crosslinker, ind, inx);
                            }
                            
                            if (psmCross!=null)
                            {
                                bestPsmCrossList.Add(psmCross);
                            }                            
                        }
                    }
                }
            }

            if (bestPsmCrossList.Count != 0)
            {
                bestPsmCross = bestPsmCrossList.OrderByDescending(p => p.XLTotalScore).First();
                if (bestPsmCrossList.Count > 1)
                {
                    bestPsmCross.DScore = Math.Abs(bestPsmCrossList.First().XLTotalScore - bestPsmCrossList[1].XLTotalScore);
                }
            }

            return bestPsmCross;
        }

        //If XlSites is of same types
        private PsmCross XlSitesSame(Ms2ScanWithSpecificMass theScan, PsmCross psmCrossAlpha, PsmCross psmCrossBeta,
            CrosslinkerTypeClass crosslinker, int ind, int inx)
        {
            PsmCross psmCross = null;
            var xlPosAlpha = PsmCross.XlPosCal(psmCrossAlpha.compactPeptide, crosslinker.CrosslinkerModSites);
            var xlPosBeta = PsmCross.XlPosCal(psmCrossBeta.compactPeptide, crosslinker.CrosslinkerModSites);
            if (xlPosAlpha.Count() >= 1 && xlPosBeta.Count() >= 1)
            {
                PsmCross.XlLocalization(theScan, psmCrossAlpha, psmCrossBeta.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosAlpha);
                PsmCross.XlLocalization(theScan, psmCrossBeta, psmCrossAlpha.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosBeta);

                if (psmCrossAlpha.XLBestScore < psmCrossBeta.XLBestScore)
                {
                    var swap = psmCrossAlpha;
                    psmCrossAlpha = psmCrossBeta;
                    psmCrossBeta = swap;
                }
                psmCrossAlpha.XlRank = new int[] { ind, inx };
                psmCrossAlpha.XLTotalScore = psmCrossAlpha.XLBestScore + psmCrossBeta.XLBestScore;
                psmCrossAlpha.XLQvalueTotalScore = Math.Sqrt(psmCrossAlpha.XLBestScore) * psmCrossBeta.XLBestScore;
                psmCrossAlpha.BetaPsmCross = psmCrossBeta;
                if (crosslinker.Cleavable)
                {
                    psmCrossAlpha.ParentIonMaxIntensityRanks.AddRange(psmCrossBeta.ParentIonMaxIntensityRanks);
                    psmCrossAlpha.ParentIonExistNum += psmCrossBeta.ParentIonExistNum;
                }
                psmCrossAlpha.CrossType = PsmCrossType.Cross;

                psmCross = psmCrossAlpha;
            }
            return psmCross;
        }

        //Deal with XlSites are of two different types
        private PsmCross XlSitesTwoDiff(Ms2ScanWithSpecificMass theScan, PsmCross psmCrossAlpha, PsmCross psmCrossBeta, 
            CrosslinkerTypeClass crosslinker, int ind, int inx)
        {
            PsmCross psmCross = null;

            var xlPosAlpha1 = PsmCross.XlPosCal(psmCrossAlpha.compactPeptide, crosslinker.CrosslinkerModSites);
            var xlPosBeta1 = PsmCross.XlPosCal(psmCrossBeta.compactPeptide, crosslinker.CrosslinkerModSites);
            var xlPosAlpha2 = PsmCross.XlPosCal(psmCrossAlpha.compactPeptide, crosslinker.CrosslinkerModSites2);
            var xlPosBeta2 = PsmCross.XlPosCal(psmCrossBeta.compactPeptide, crosslinker.CrosslinkerModSites2);

            if ((xlPosAlpha1.Count() >= 1 && xlPosBeta2.Count() >= 1) || (xlPosAlpha2.Count() >= 1 && xlPosBeta1.Count() >= 1))
            {

                if (xlPosAlpha2.Count() < 1 || xlPosBeta1.Count() < 1)
                {
                    PsmCross.XlLocalization(theScan, psmCrossAlpha, psmCrossBeta.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosAlpha1);
                    PsmCross.XlLocalization(theScan, psmCrossBeta, psmCrossAlpha.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosBeta2);
                }

                if (xlPosAlpha1.Count() < 1 || xlPosBeta2.Count() < 1)
                {
                    PsmCross.XlLocalization(theScan, psmCrossAlpha, psmCrossBeta.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosAlpha2);
                    PsmCross.XlLocalization(theScan, psmCrossBeta, psmCrossAlpha.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosBeta1);
                }

                if ((xlPosAlpha1.Count() >= 1 && xlPosBeta2.Count() >= 1) && (xlPosAlpha2.Count() >= 1 && xlPosBeta1.Count() >= 1))
                {
                    var psmCrossAlpha1 = psmCrossAlpha;
                    var psmCrossAlpha2 = psmCrossAlpha;
                    var psmCrossBeta1 = psmCrossBeta;
                    var psmCrossBeta2 = psmCrossBeta;

                    PsmCross.XlLocalization(theScan, psmCrossAlpha1, psmCrossBeta.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosAlpha1);
                    PsmCross.XlLocalization(theScan, psmCrossBeta2, psmCrossAlpha.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosBeta2);
                    PsmCross.XlLocalization(theScan, psmCrossAlpha2, psmCrossBeta.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosAlpha2);
                    PsmCross.XlLocalization(theScan, psmCrossBeta1, psmCrossAlpha.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, CommonParameters.ProductMassTolerance, charge_2_3, charge_2_3_PrimeFragment, xlPosBeta1);
                    if (psmCrossAlpha1.XLBestScore + psmCrossBeta2.XLBestScore > psmCrossAlpha2.XLBestScore + psmCrossBeta1.XLBestScore)
                    {
                        psmCrossAlpha = psmCrossAlpha1;
                        psmCrossBeta = psmCrossBeta2;
                    }
                    else
                    {
                        psmCrossAlpha = psmCrossAlpha2;
                        psmCrossBeta = psmCrossBeta1;
                    }
                }

                if (psmCrossAlpha.XLBestScore < psmCrossBeta.XLBestScore)
                {
                    var swap = psmCrossAlpha;
                    psmCrossAlpha = psmCrossBeta;
                    psmCrossBeta = swap;
                }
                psmCrossAlpha.XlRank = new int[] { ind, inx };
                psmCrossAlpha.XLTotalScore = psmCrossAlpha.XLBestScore + psmCrossBeta.XLBestScore;
                psmCrossAlpha.XLQvalueTotalScore = Math.Sqrt(psmCrossAlpha.XLBestScore) * psmCrossBeta.XLBestScore;
                psmCrossAlpha.BetaPsmCross = psmCrossBeta;
                if (crosslinker.Cleavable)
                {
                    psmCrossAlpha.ParentIonMaxIntensityRanks.AddRange(psmCrossBeta.ParentIonMaxIntensityRanks);
                    psmCrossAlpha.ParentIonExistNum += psmCrossBeta.ParentIonExistNum;
                }
                psmCrossAlpha.CrossType = PsmCrossType.Cross;

                psmCross = psmCrossAlpha;
            }
            return psmCross;
        }

        #endregion Private Methods
    }
}