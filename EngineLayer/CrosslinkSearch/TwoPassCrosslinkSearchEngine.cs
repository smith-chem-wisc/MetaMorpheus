using MassSpectrometry;
using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.CrosslinkSearch
{
    public class TwoPassCrosslinkSearchEngine : MetaMorpheusEngine
    {
        protected const int FragmentBinsPerDalton = 1000;
        protected readonly List<int>[] FragmentIndex;
        protected readonly PeptideSpectralMatch[] GlobalPsms;
        protected readonly List<PsmCross> GlobalPsmsCross;
        protected readonly Ms2ScanWithSpecificMass[] ListOfSortedms2Scans;
        protected readonly List<CompactPeptide> PeptideIndex;
        protected readonly List<ProductType> ProductTypes;
        protected readonly int CurrentPartition;
        protected readonly bool AddComplementaryIons;
        protected readonly MassDiffAcceptor MassDiffAcceptor;
        protected readonly List<DissociationType> DissociationTypes;

        //Crosslink parameters
        private bool _searchGlycan;
        private readonly CrosslinkerTypeClass Crosslinker;
        private readonly bool CrosslinkSearchTop;
        private readonly int CrosslinkSearchTopNum;
        private readonly Tolerance XLPrecusorMsTl;
        private readonly bool QuenchH2O;
        private readonly bool QuenchNH2;
        private readonly bool QuenchTris;
        private readonly bool Charge_2_3;
        private MassDiffAcceptor XLPrecusorSearchMode;

        public TwoPassCrosslinkSearchEngine(List<PsmCross> globalPsmsCross, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition, CommonParameters commonParameters, bool addCompIons, bool searchGlycan, bool searchGlycanBgYgIndex, Tolerance XLPrecusorMsTl, CrosslinkerTypeClass crosslinker, bool CrosslinkSearchTop, int CrosslinkSearchTopNum, bool quench_H2O, bool quench_NH2, bool quench_Tris, bool charge_2_3, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            this.GlobalPsmsCross = globalPsmsCross;
            this.ListOfSortedms2Scans = listOfSortedms2Scans;
            this.PeptideIndex = peptideIndex;
            this.FragmentIndex = fragmentIndex;
            this.ProductTypes = lp;
            this.CurrentPartition = currentPartition + 1;
            this.AddComplementaryIons = addCompIons;
            this._searchGlycan = searchGlycan;
            this.MassDiffAcceptor = new OpenSearchMode();
            this.DissociationTypes = DetermineDissociationType(lp);
            this.XLPrecusorMsTl = XLPrecusorMsTl;
            XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(XLPrecusorMsTl.Value);
            this.Crosslinker = crosslinker;
            this.CrosslinkSearchTop = CrosslinkSearchTop;
            this.CrosslinkSearchTopNum = CrosslinkSearchTopNum;
            this.QuenchH2O = quench_H2O;
            this.QuenchNH2 = quench_NH2;
            this.QuenchTris = quench_Tris;
            this.Charge_2_3 = charge_2_3;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing crosslink search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = (byte)commonParameters.ScoreCutoff;

            Parallel.ForEach(Partitioner.Create(0, ListOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (range, loopState) =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();
                byte[] scoringTableGly = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObservedGly = new List<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        loopState.Stop();
                        return;
                    }

                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    Array.Clear(scoringTableGly, 0, scoringTableGly.Length);
                    idsOfPeptidesPossiblyObservedGly.Clear();

                    var scan = ListOfSortedms2Scans[i];

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(scan);
                    List<BestPeptideScoreNotch> bestPeptideScoreNotchList = new List<BestPeptideScoreNotch>();
                    // get allowed theoretical masses from the known experimental mass
                    // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)
                    // this is just PRELIMINARY precursor-mass filtering
                    // additional checks are made later to ensure that the theoretical precursor mass is acceptable
                    var notches = MassDiffAcceptor.GetAllowedPrecursorMassIntervals(scan.PrecursorMass);

                    double lowestMassPeptideToLookFor = Double.NegativeInfinity;
                    double highestMassPeptideToLookFor = Double.PositiveInfinity;

                    double largestMassDiff = notches.Max(p => p.AllowedInterval.Maximum);
                    double smallestMassDiff = notches.Min(p => p.AllowedInterval.Minimum);

                    if (!Double.IsInfinity(largestMassDiff))
                    {
                        double largestOppositeMassDiff = -1 * (notches.Max(p => p.AllowedInterval.Maximum) - scan.PrecursorMass);
                        lowestMassPeptideToLookFor = scan.PrecursorMass + largestOppositeMassDiff;
                    }
                    if (!Double.IsNegativeInfinity(smallestMassDiff))
                    {
                        double smallestOppositeMassDiff = -1 * (notches.Min(p => p.AllowedInterval.Minimum) - scan.PrecursorMass);
                        highestMassPeptideToLookFor = scan.PrecursorMass + smallestOppositeMassDiff;
                    }

                    // first-pass scoring
                    IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor, FragmentIndex);
                    
                    // done with indexed scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        if (idsOfPeptidesPossiblyObservedGly.Any())
                        {
                            foreach (var iG in idsOfPeptidesPossiblyObservedGly)
                            {
                                if (!idsOfPeptidesPossiblyObserved.Contains(iG))
                                { idsOfPeptidesPossiblyObserved.Add(iG); }
                                scoringTable[iG] += scoringTableGly[iG];
                            }
                        }
                        List<int> idsRankedByScore = idsOfPeptidesPossiblyObserved.Where(p => scoringTable[p] > byteScoreCutoff).OrderByDescending(p => scoringTable[p]).ToList();                      
                        if (CrosslinkSearchTop)
                        {
                            idsRankedByScore = idsRankedByScore.Take(CrosslinkSearchTopNum).ToList();
                        }
                        foreach (var id in idsRankedByScore)
                        {
                            var peptide = PeptideIndex[id];

                            double thePrecursorMass = scan.PrecursorMass;
                            int notch = MassDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMassIncludingFixedMods);
                            BestPeptideScoreNotch bestPeptideScoreNotch = new BestPeptideScoreNotch(peptide, scoringTable[id], notch);
                            bestPeptideScoreNotchList.Add(bestPeptideScoreNotch);
                        }

                        {
                            var possiblePsmCross = FindCrosslinkedPeptide(scan, bestPeptideScoreNotchList, i);
                            if (possiblePsmCross != null)
                            {
                                lock (GlobalPsmsCross)
                                    GlobalPsmsCross.Add(possiblePsmCross);
                            }
                        }                      
                    }
                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedms2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing crosslink search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));
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
                int obsFragmentFloorMass = (int)Math.Floor((commonParameters.ProductMassTolerance.GetMinimumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);
                if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)
                    obsFragmentFloorMass = obsPreviousFragmentCeilingMz;
                int obsFragmentCeilingMass = (int)Math.Ceiling((commonParameters.ProductMassTolerance.GetMaximumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);
                obsPreviousFragmentCeilingMz = obsFragmentCeilingMass + 1;
                for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                    if (FragmentIndex[fragmentBin] != null)
                        binsToSearch.Add(fragmentBin);
            }
            return binsToSearch;
        }

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
                if (PeptideIndex[peptideIdsInThisBin[m]].MonoisotopicMassIncludingFixedMods < peptideMassToLookFor)
                    l = m + 1;
                else
                    r = m - 1;
            }
            if (m > 0)
                m--;
            return m;
        }

        private void IndexedScoring(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor, List<int>[] theFragmentIndex)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = theFragmentIndex[binsToSearch[i]];

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
                        var nextPep = PeptideIndex[nextId];
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
                    if (scoringTable[id] == byteScoreCutoff && MassDiffAcceptor.Accepts(scanPrecursorMass, PeptideIndex[id].MonoisotopicMassIncludingFixedMods) >= 0)
                        idsOfPeptidesPossiblyObserved.Add(id);
                }

            }
        }

        //Targetting function: to find two peptides that in the Top matched peptides
        private PsmCross FindCrosslinkedPeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int i)
        {
            List<PsmCross> bestPsmCrossList = new List<PsmCross>();
            PsmCross bestPsmCross = null;
            string crosslinkerModSitesAll = new string((Crosslinker.CrosslinkerModSites + Crosslinker.CrosslinkerModSites2).ToCharArray().Distinct().ToArray());
            for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
            {
                //Single Peptide
                if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods) >= 0)
                {
                    var psmCrossSingle = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams);
                    var pmmh = psmCrossSingle.GetTheoreticalFragmentIons(ProductTypes);
                    psmCrossSingle.MatchedIons = MatchFragmentIons(theScan.TheScan.MassSpectrum, pmmh, commonParameters);
                    psmCrossSingle.BestScore = CalculatePeptideScore(theScan.TheScan, psmCrossSingle.MatchedIons, 0);               
                    psmCrossSingle.XLTotalScore = psmCrossSingle.BestScore;
                    psmCrossSingle.CrossType = PsmCrossType.Singe;
                    psmCrossSingle.XlRank = new List<int> { ind };

                    bestPsmCrossList.Add(psmCrossSingle);
                }
                //Deadend Peptide
                else if (QuenchTris && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + Crosslinker.DeadendMassTris) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams);
                    var xlPos = psmCrossEnd.XlPosCal(crosslinkerModSitesAll);
                    if (xlPos.Count() >= 1)
                    {
                        var pmmhList = psmCrossEnd.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPos, Crosslinker.DeadendMassTris);
                        psmCrossEnd.GetBestMatch(theScan, pmmhList, commonParameters);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.BestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndTris;
                        psmCrossEnd.XlRank = new List<int> { ind };
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                else if (QuenchH2O && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + Crosslinker.DeadendMassH2O) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams);
                    var xlPos = psmCrossEnd.XlPosCal(crosslinkerModSitesAll);
                    if (xlPos.Count() >= 1)
                    {
                        var pmmhList = psmCrossEnd.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPos, Crosslinker.DeadendMassH2O);
                        psmCrossEnd.GetBestMatch(theScan, pmmhList, commonParameters);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.BestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndH2O;
                        psmCrossEnd.XlRank = new List<int> { ind };
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                else if (QuenchNH2 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + Crosslinker.DeadendMassNH2) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams);
                    var xlPos = psmCrossEnd.XlPosCal(crosslinkerModSitesAll);
                    if (xlPos.Count() >= 1)
                    {
                        var pmmhList = psmCrossEnd.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPos, Crosslinker.DeadendMassNH2);
                        psmCrossEnd.GetBestMatch(theScan, pmmhList, commonParameters);
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndNH2;
                        psmCrossEnd.XlRank = new List<int> { ind };
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                //loop peptide
                else if (Crosslinker.LoopMass != 0 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + Crosslinker.LoopMass) >= 0)
                {
                    var psmCrossLoop = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams);
                    var xlPos = psmCrossLoop.XlPosCal(Crosslinker.CrosslinkerModSites);
                    if (xlPos.Count() >= 2)
                    {
                        var pmmhList = psmCrossLoop.XlLoopGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPos, Crosslinker.LoopMass);
                        psmCrossLoop.GetBestMatch(theScan, pmmhList, commonParameters);
                        psmCrossLoop.XLTotalScore = psmCrossLoop.BestScore;
                        psmCrossLoop.CrossType = PsmCrossType.Loop;
                        psmCrossLoop.XlRank = new List<int> { ind };
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
                        if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, x + y + Crosslinker.TotalMass) >= 0)
                        {
                            var psmCrossAlpha = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams);
                            var psmCrossBeta = new PsmCross(theScanBestPeptide[inx].BestPeptide, theScanBestPeptide[inx].BestNotch, theScanBestPeptide[inx].BestScore, i, theScan, commonParameters.DigestionParams);

                            PsmCross psmCross = null;
                            if (Crosslinker.CrosslinkerModSites == Crosslinker.CrosslinkerModSites2)
                            {
                                psmCross = XlSitesSame(theScan, psmCrossAlpha, psmCrossBeta, Crosslinker, ind, inx, x, y);
                            }
                            else
                            {
                                psmCross = XlSitesTwoDiff(theScan, psmCrossAlpha, psmCrossBeta, Crosslinker, ind, inx, x, y);
                            }

                            if (psmCross != null)
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
            CrosslinkerTypeClass crosslinker, int ind, int inx, double x, double y)
        {
            PsmCross psmCross = null;
            var xlPosAlpha = psmCrossAlpha.XlPosCal( crosslinker.CrosslinkerModSites);
            var xlPosBeta = psmCrossBeta.XlPosCal(crosslinker.CrosslinkerModSites);
            if (xlPosAlpha.Count() >= 1 && xlPosBeta.Count() >= 1)
            {
                var pmmhListAlpha = psmCrossAlpha.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosAlpha, y);
                psmCrossAlpha.GetBestMatch(theScan, pmmhListAlpha, commonParameters);
                var pmmhListBeta = psmCrossBeta.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosBeta, x);
                psmCrossBeta.GetBestMatch(theScan, pmmhListBeta, commonParameters);

                if (psmCrossAlpha.BestScore < psmCrossBeta.BestScore)
                {
                    var swap = psmCrossAlpha;
                    psmCrossAlpha = psmCrossBeta;
                    psmCrossBeta = swap;
                }
                psmCrossAlpha.XlRank = new List<int> { ind, inx };
                psmCrossAlpha.XLTotalScore = psmCrossAlpha.BestScore + psmCrossBeta.BestScore;
                psmCrossAlpha.XLQvalueTotalScore = Math.Sqrt(psmCrossAlpha.BestScore) * psmCrossBeta.BestScore;
                psmCrossAlpha.BetaPsmCross = psmCrossBeta;
                if (crosslinker.Cleavable)
                {        
                    psmCrossAlpha.ParentIonMaxIntensityRanks = psmCrossAlpha.MatchedIons.Where(p => p.TheoreticalFragmentIon.ProductType == ProductType.None).Select(p => p.IntensityRank).ToList();
                    psmCrossAlpha.ParentIonExistNum = psmCrossAlpha.ParentIonMaxIntensityRanks.Count;
                }
                psmCrossAlpha.CrossType = PsmCrossType.Cross;

                psmCross = psmCrossAlpha;
            }
            return psmCross;
        }

        //Deal with XlSites are of two different types
        private PsmCross XlSitesTwoDiff(Ms2ScanWithSpecificMass theScan, PsmCross psmCrossAlpha, PsmCross psmCrossBeta,
            CrosslinkerTypeClass crosslinker, int ind, int inx, double x, double y)
        {
            PsmCross psmCross = null;

            var xlPosAlpha1 = psmCrossAlpha.XlPosCal(crosslinker.CrosslinkerModSites);
            var xlPosBeta1 = psmCrossBeta.XlPosCal(crosslinker.CrosslinkerModSites);
            var xlPosAlpha2 = psmCrossAlpha.XlPosCal(crosslinker.CrosslinkerModSites2);
            var xlPosBeta2 = psmCrossBeta.XlPosCal(crosslinker.CrosslinkerModSites2);

            if ((xlPosAlpha1.Count() >= 1 && xlPosBeta2.Count() >= 1) || (xlPosAlpha2.Count() >= 1 && xlPosBeta1.Count() >= 1))
            {
                if (xlPosAlpha2.Count() < 1 || xlPosBeta1.Count() < 1)
                {
                    var pmmhListAlpha = psmCrossAlpha.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosAlpha1, y);
                    psmCrossAlpha.GetBestMatch(theScan, pmmhListAlpha, commonParameters);
                    var pmmhListBeta = psmCrossBeta.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosBeta2, x);
                    psmCrossBeta.GetBestMatch(theScan, pmmhListBeta, commonParameters);
                }

                if (xlPosAlpha1.Count() < 1 || xlPosBeta2.Count() < 1)
                {
                    var pmmhListAlpha = psmCrossAlpha.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosAlpha2, y);
                    psmCrossAlpha.GetBestMatch(theScan, pmmhListAlpha, commonParameters);
                    var pmmhListBeta = psmCrossBeta.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosBeta1, x);
                    psmCrossBeta.GetBestMatch(theScan, pmmhListBeta, commonParameters);

                }

                if ((xlPosAlpha1.Count() >= 1 && xlPosBeta2.Count() >= 1) && (xlPosAlpha2.Count() >= 1 && xlPosBeta1.Count() >= 1))
                {
                    var psmCrossAlpha1 = psmCrossAlpha;
                    var psmCrossAlpha2 = psmCrossAlpha;
                    var psmCrossBeta1 = psmCrossBeta;
                    var psmCrossBeta2 = psmCrossBeta;

                    var pmmhListAlpha1 = psmCrossAlpha.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosAlpha1, y);
                    psmCrossAlpha1.GetBestMatch(theScan, pmmhListAlpha1, commonParameters);
                    var pmmhListBeta1 = psmCrossBeta.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosBeta1, x);
                    psmCrossBeta1.GetBestMatch(theScan, pmmhListBeta1, commonParameters);
                    var pmmhListAlpha2 = psmCrossAlpha.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosAlpha2, y);
                    psmCrossAlpha2.GetBestMatch(theScan, pmmhListAlpha2, commonParameters);
                    var pmmhListBeta2 = psmCrossBeta.XlGetTheoreticalFramentIons(ProductTypes, Charge_2_3, Crosslinker, xlPosBeta2, x);
                    psmCrossBeta2.GetBestMatch(theScan, pmmhListBeta2, commonParameters);

                    if (psmCrossAlpha1.BestScore + psmCrossBeta2.BestScore > psmCrossAlpha2.BestScore + psmCrossBeta1.BestScore)
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

                if (psmCrossAlpha.BestScore < psmCrossBeta.BestScore)
                {
                    var swap = psmCrossAlpha;
                    psmCrossAlpha = psmCrossBeta;
                    psmCrossBeta = swap;
                }
                psmCrossAlpha.XlRank = new List<int> { ind, inx };
                psmCrossAlpha.XLTotalScore = psmCrossAlpha.BestScore + psmCrossBeta.BestScore;
                psmCrossAlpha.XLQvalueTotalScore = Math.Sqrt(psmCrossAlpha.BestScore) * psmCrossBeta.BestScore;
                psmCrossAlpha.BetaPsmCross = psmCrossBeta;
                if (crosslinker.Cleavable)
                {
                    psmCrossAlpha.ParentIonMaxIntensityRanks = psmCrossAlpha.MatchedIons.Where(p => p.TheoreticalFragmentIon.ProductType == ProductType.None).Select(p => p.IntensityRank).ToList();
                    psmCrossAlpha.ParentIonExistNum = psmCrossAlpha.ParentIonMaxIntensityRanks.Count;
                }
                psmCrossAlpha.CrossType = PsmCrossType.Cross;

                psmCross = psmCrossAlpha;
            }
            return psmCross;
        }
    }
}