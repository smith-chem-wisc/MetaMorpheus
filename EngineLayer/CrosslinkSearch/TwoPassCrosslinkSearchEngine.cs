using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.Fragmentation;
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
        protected readonly List<CrosslinkSpectralMatch> GlobalPsmsCross;
        protected readonly Ms2ScanWithSpecificMass[] ListOfSortedms2Scans;
        protected readonly List<PeptideWithSetModifications> PeptideIndex;
        protected readonly List<ProductType> ProductTypes;
        protected readonly int CurrentPartition;
        protected readonly MassDiffAcceptor MassDiffAcceptor;
        protected readonly DissociationType DissociationType;

        // crosslinker molecule
        private readonly Crosslinker Crosslinker;

        private readonly bool CrosslinkSearchTop;
        private readonly int CrosslinkSearchTopNum;

        //private readonly bool CrosslinkSearchWithCrosslinkerMod;
        private readonly Tolerance XLPrecusorMsTl;

        //private readonly Tolerance XLBetaPrecusorMsTl;
        private readonly bool QuenchH2O;

        private readonly bool QuenchNH2;
        private readonly bool QuenchTris;
        private readonly bool Charge_2_3;
        private readonly bool Charge_2_3_PrimeFragment;
        private MassDiffAcceptor XLPrecusorSearchMode;

        public TwoPassCrosslinkSearchEngine(List<CrosslinkSpectralMatch> globalPsmsCross, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition, CommonParameters commonParameters, bool addCompIons, Tolerance XLPrecusorMsTl, Crosslinker crosslinker, bool CrosslinkSearchTop, int CrosslinkSearchTopNum, bool quench_H2O, bool quench_NH2, bool quench_Tris, bool charge_2_3, bool charge_2_3_PrimeFragment, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            this.GlobalPsmsCross = globalPsmsCross;
            this.ListOfSortedms2Scans = listOfSortedms2Scans;
            this.PeptideIndex = peptideIndex;
            this.FragmentIndex = fragmentIndex;
            this.ProductTypes = lp;
            this.CurrentPartition = currentPartition + 1;
            //Here use LowTheoreticalDiffAcceptor in practice doesn't work in 12/2/2017
            this.MassDiffAcceptor = new OpenSearchMode();
            this.DissociationType = commonParameters.DissociationType;
            this.XLPrecusorMsTl = XLPrecusorMsTl;
            XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(XLPrecusorMsTl.Value);
            //if (XLBetaPrecusorMsTl.ToString().Contains("Absolute"))
            //{
            //    XLPrecusorSearchMode = new SingleAbsoluteAroundZeroSearchMode(XLPrecusorMsTl.Value);
            //}
            this.Crosslinker = crosslinker;
            this.CrosslinkSearchTop = CrosslinkSearchTop;
            this.CrosslinkSearchTopNum = CrosslinkSearchTopNum;
            this.QuenchH2O = quench_H2O;
            this.QuenchNH2 = quench_NH2;
            this.QuenchTris = quench_Tris;
            this.Charge_2_3 = charge_2_3;
            this.Charge_2_3_PrimeFragment = charge_2_3_PrimeFragment;
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
                    var scan = ListOfSortedms2Scans[i];

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = ModernSearchEngine.GetBinsToSearch(scan, commonParameters, FragmentIndex);
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
                    ModernSearchEngine.IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor, FragmentIndex, PeptideIndex, XLPrecusorSearchMode, 0);

                    // done with indexed scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        List<int> idsRankedByScore = new List<int>();
                        foreach (var id in idsOfPeptidesPossiblyObserved)
                        {
                            idsRankedByScore.Add(id);
                        }

                        idsRankedByScore = idsRankedByScore.OrderByDescending(p => scoringTable[p]).ToList();

                        // take top N hits for this scan
                        idsRankedByScore = idsRankedByScore.Take(CrosslinkSearchTopNum).ToList();

                        if (CrosslinkSearchTop)
                        {
                            idsOfPeptidesPossiblyObserved = idsRankedByScore;
                        }
                        foreach (var id in idsOfPeptidesPossiblyObserved)
                        {
                            var peptide = PeptideIndex[id];

                            double thePrecursorMass = scan.PrecursorMass;
                            int notch = MassDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMass);
                            BestPeptideScoreNotch bestPeptideScoreNotch = new BestPeptideScoreNotch(peptide, scoringTable[id], notch);
                            bestPeptideScoreNotchList.Add(bestPeptideScoreNotch);
                        }

                        // combine individual peptide hits with crosslinker mass to find best crosslink PSM hit
                        CrosslinkSpectralMatch possiblePsmCross = FindCrosslinkedPeptide(scan, bestPeptideScoreNotchList, i);
                        if (possiblePsmCross != null)
                        {
                            lock (GlobalPsmsCross)
                            {
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

        //Targetting function: to find two peptides that in the Top matched peptides
        private CrosslinkSpectralMatch FindCrosslinkedPeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int i)
        {
            PeptideWithSetModifications peptide = theScanBestPeptide.First().BestPeptide;
            var products = peptide.Fragment(DissociationType, FragmentationTerminus.Both).ToList();
            var matchedFragmentIons = MatchFragmentIons(theScan.TheScan.MassSpectrum, products, commonParameters);
            
            List<CrosslinkSpectralMatch> bestPsmCrossList = new List<CrosslinkSpectralMatch>();
            CrosslinkSpectralMatch bestPsmCross = null;
            string crosslinkerModSitesAll = new string((Crosslinker.CrosslinkerModSites + Crosslinker.CrosslinkerModSites2).ToCharArray().Distinct().ToArray());
            for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
            {
                var xlPos = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinkerModSitesAll, theScanBestPeptide[ind].BestPeptide);

                //Single Peptide
                if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass) >= 0)
                {
                    var psmCrossSingle = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams, matchedFragmentIons);
                    psmCrossSingle.BestScore = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons, 0);
                    psmCrossSingle.XLTotalScore = psmCrossSingle.BestScore;
                    psmCrossSingle.CrossType = PsmCrossType.Single;
                    psmCrossSingle.XlRank = new List<int> { ind };

                    bestPsmCrossList.Add(psmCrossSingle);
                }
                //Deadend Peptide
                else if (QuenchTris && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass + Crosslinker.DeadendMassTris) >= 0)
                {
                    var psmCrossEnd = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams, matchedFragmentIons);
                    
                    if (xlPos.Count() >= 1)
                    {
                        var pmmhList = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPos, Crosslinker.DeadendMassTris, theScanBestPeptide[ind].BestPeptide);
                        psmCrossEnd.GetBestMatch(theScan, pmmhList, commonParameters);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.BestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndTris;
                        psmCrossEnd.XlRank = new List<int> { ind };
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                else if (QuenchH2O && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass + Crosslinker.DeadendMassH2O) >= 0)
                {
                    var psmCrossEnd = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams, matchedFragmentIons);
                    if (xlPos.Count() >= 1)
                    {
                        var pmmhList = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPos, Crosslinker.DeadendMassH2O, theScanBestPeptide[ind].BestPeptide);
                        psmCrossEnd.GetBestMatch(theScan, pmmhList, commonParameters);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.BestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndH2O;
                        psmCrossEnd.XlRank = new List<int> { ind };
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                else if (QuenchNH2 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass + Crosslinker.DeadendMassNH2) >= 0)
                {
                    var psmCrossEnd = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams, matchedFragmentIons);
                    if (xlPos.Count() >= 1)
                    {
                        var pmmhList = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPos, Crosslinker.DeadendMassNH2, theScanBestPeptide[ind].BestPeptide);
                        psmCrossEnd.GetBestMatch(theScan, pmmhList, commonParameters);
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndNH2;
                        psmCrossEnd.XlRank = new List<int> { ind };
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                //loop peptide
                else if (Crosslinker.LoopMass != 0 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass + Crosslinker.LoopMass) >= 0)
                {
                    var psmCrossLoop = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams, matchedFragmentIons);
                    if (xlPos.Count() >= 2)
                    {
                        var pmmhList = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPos, Crosslinker.LoopMass, theScanBestPeptide[ind].BestPeptide);
                        psmCrossLoop.GetBestMatch(theScan, pmmhList, commonParameters);
                        psmCrossLoop.XLTotalScore = psmCrossLoop.BestScore;
                        psmCrossLoop.CrossType = PsmCrossType.Loop;
                        psmCrossLoop.XlRank = new List<int> { ind };
                        bestPsmCrossList.Add(psmCrossLoop);
                    }
                }
                //Cross-linked peptide
                else if (theScan.PrecursorMass - theScanBestPeptide[ind].BestPeptide.MonoisotopicMass >= 200)
                {
                    double alphaPeptideMass = theScanBestPeptide[ind].BestPeptide.MonoisotopicMass;
                    
                    for (int inx = ind; inx < theScanBestPeptide.Count; inx++)
                    {
                        var y = theScanBestPeptide[inx].BestPeptide.MonoisotopicMass;
                        if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, alphaPeptideMass + y + Crosslinker.TotalMass) >= 0)
                        {
                            var psmCrossAlpha = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, commonParameters.DigestionParams, matchedFragmentIons);

                            // score beta peptide
                            PeptideWithSetModifications betaPeptide = theScanBestPeptide[inx].BestPeptide;
                            var betaProducts = betaPeptide.Fragment(DissociationType, FragmentationTerminus.Both).ToList();
                            var betaMatchedIons = MatchFragmentIons(theScan.TheScan.MassSpectrum, betaProducts, commonParameters);

                            var psmCrossBeta = new CrosslinkSpectralMatch(theScanBestPeptide[inx].BestPeptide, theScanBestPeptide[inx].BestNotch, theScanBestPeptide[inx].BestScore, i, theScan, commonParameters.DigestionParams, betaMatchedIons);

                            CrosslinkSpectralMatch psmCross = null;
                            if (Crosslinker.CrosslinkerModSites == Crosslinker.CrosslinkerModSites2)
                            {
                                psmCross = XlSitesSame(theScan, psmCrossAlpha, psmCrossBeta, Crosslinker, ind, inx, alphaPeptideMass, y);
                            }
                            else
                            {
                                psmCross = XlSitesTwoDiff(theScan, psmCrossAlpha, psmCrossBeta, Crosslinker, ind, inx, alphaPeptideMass, y);
                            }

                            if (psmCross != null)
                            {
                                psmCross.BetaPeptide.ResolveAllAmbiguities();
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

            bestPsmCross.ResolveAllAmbiguities();
            return bestPsmCross;
        }

        /// <summary>
        /// Used if crosslink sites are of the same residue type (e.g., K to K)
        /// </summary>
        private CrosslinkSpectralMatch XlSitesSame(Ms2ScanWithSpecificMass theScan, CrosslinkSpectralMatch psmCrossAlpha, CrosslinkSpectralMatch psmCrossBeta,
            Crosslinker crosslinker, int ind, int inx, double x, double y)
        {
            CrosslinkSpectralMatch psmCross = null;
            var xlPosAlpha = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
            var xlPosBeta = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites, psmCrossBeta.BestMatchingPeptideWithSetMods.First().Pwsm);

            if (xlPosAlpha.Count() >= 1 && xlPosBeta.Count() >= 1)
            {
                var pmmhListAlpha = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosAlpha, y, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
                psmCrossAlpha.GetBestMatch(theScan, pmmhListAlpha, commonParameters);

                var pmmhListBeta = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosAlpha, y, psmCrossBeta.BestMatchingPeptideWithSetMods.First().Pwsm);
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
                psmCrossAlpha.BetaPeptide = psmCrossBeta;

                if (crosslinker.Cleavable)
                {
                    //TODO: re-enable intensity ranks
                    //psmCrossAlpha.ParentIonMaxIntensityRanks = psmCrossAlpha.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).Select(p => p.IntensityRank).ToList();
                    psmCrossAlpha.ParentIonMaxIntensityRanks = new List<int>();

                    psmCrossAlpha.ParentIonExistNum = psmCrossAlpha.ParentIonMaxIntensityRanks.Count;
                }
                psmCrossAlpha.CrossType = PsmCrossType.Cross;

                psmCross = psmCrossAlpha;
            }
            return psmCross;
        }

        /// <summary>
        /// Used if crosslink sites are two different residue type (e.g., K to R)
        /// </summary>
        private CrosslinkSpectralMatch XlSitesTwoDiff(Ms2ScanWithSpecificMass theScan, CrosslinkSpectralMatch psmCrossAlpha, CrosslinkSpectralMatch psmCrossBeta,
            Crosslinker crosslinker, int ind, int inx, double x, double y)
        {
            CrosslinkSpectralMatch psmCross = null;
            
            var xlPosAlpha1 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
            var xlPosBeta1 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites, psmCrossBeta.BestMatchingPeptideWithSetMods.First().Pwsm);
            var xlPosAlpha2 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites2, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
            var xlPosBeta2 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites2, psmCrossBeta.BestMatchingPeptideWithSetMods.First().Pwsm);

            if ((xlPosAlpha1.Count() >= 1 && xlPosBeta2.Count() >= 1) || (xlPosAlpha2.Count() >= 1 && xlPosBeta1.Count() >= 1))
            {
                if (xlPosAlpha2.Count() < 1 || xlPosBeta1.Count() < 1)
                {
                    var pmmhListAlpha = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosAlpha2, y, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
                    psmCrossAlpha.GetBestMatch(theScan, pmmhListAlpha, commonParameters);
                    var pmmhListBeta = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosBeta2, x, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
                    psmCrossBeta.GetBestMatch(theScan, pmmhListBeta, commonParameters);
                }

                if (xlPosAlpha1.Count() < 1 || xlPosBeta2.Count() < 1)
                {
                    var pmmhListAlpha = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosAlpha2, y, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
                    psmCrossAlpha.GetBestMatch(theScan, pmmhListAlpha, commonParameters);
                    var pmmhListBeta = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosBeta1, x, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
                    psmCrossBeta.GetBestMatch(theScan, pmmhListBeta, commonParameters);
                }

                if ((xlPosAlpha1.Count() >= 1 && xlPosBeta2.Count() >= 1) && (xlPosAlpha2.Count() >= 1 && xlPosBeta1.Count() >= 1))
                {
                    var psmCrossAlpha1 = psmCrossAlpha;
                    var psmCrossAlpha2 = psmCrossAlpha;
                    var psmCrossBeta1 = psmCrossBeta;
                    var psmCrossBeta2 = psmCrossBeta;
                    
                    var pmmhListAlpha1 = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosAlpha1, y, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
                    psmCrossAlpha1.GetBestMatch(theScan, pmmhListAlpha1, commonParameters);
                    var pmmhListBeta1 = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosBeta1, x, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
                    psmCrossBeta1.GetBestMatch(theScan, pmmhListBeta1, commonParameters);
                    var pmmhListAlpha2 = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosAlpha2, y, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
                    psmCrossAlpha2.GetBestMatch(theScan, pmmhListAlpha2, commonParameters);
                    var pmmhListBeta2 = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType, Charge_2_3, Crosslinker, xlPosBeta2, x, psmCrossAlpha.BestMatchingPeptideWithSetMods.First().Pwsm);
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
                psmCrossAlpha.BetaPeptide = psmCrossBeta;
                if (crosslinker.Cleavable)
                {
                    // TODO: re-enable intensity rank
                    //psmCrossAlpha.ParentIonMaxIntensityRanks = psmCrossAlpha.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).Select(p => p.IntensityRank).ToList();
                    psmCrossAlpha.ParentIonMaxIntensityRanks = new List<int>();
                    psmCrossAlpha.ParentIonExistNum = psmCrossAlpha.ParentIonMaxIntensityRanks.Count;
                }
                psmCrossAlpha.CrossType = PsmCrossType.Cross;

                psmCross = psmCrossAlpha;
            }
            return psmCross;
        }
    }
}