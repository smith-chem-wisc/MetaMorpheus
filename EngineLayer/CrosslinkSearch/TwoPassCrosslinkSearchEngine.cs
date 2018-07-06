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
        private readonly CrosslinkerTypeClass Crosslinker;

        private readonly bool CrosslinkSearchTop;
        private readonly int CrosslinkSearchTopNum;
        //private readonly bool CrosslinkSearchWithCrosslinkerMod;

        private readonly Tolerance XLPrecusorMsTl;

        //private readonly Tolerance XLBetaPrecusorMsTl;
        private MassDiffAcceptor XLPrecusorSearchMode;

        private readonly bool QuenchH2O;
        private readonly bool QuenchNH2;
        private readonly bool QuenchTris;
        private readonly bool Charge_2_3;
        private readonly bool Charge_2_3_PrimeFragment;

        public TwoPassCrosslinkSearchEngine(List<PsmCross> globalPsmsCross, Ms2ScanWithSpecificMass[] listOfSortedms2Scans,
                List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex,
                List<ProductType> productTypes, int currentPartition, CommonParameters commonParameters, bool addCompIons,
                Tolerance xlPrecusorMsTl, CrosslinkerTypeClass crosslinker, bool crosslinkSearchTop, int crosslinkSearchTopNum,
                bool quench_H2O, bool quench_NH2, bool quench_Tris, bool charge_2_3, bool charge_2_3_PrimeFragment,
                List<string> nestedIds)
            : base(commonParameters, nestedIds)
        {
            GlobalPsmsCross = globalPsmsCross;
            ListOfSortedms2Scans = listOfSortedms2Scans;
            PeptideIndex = peptideIndex;
            FragmentIndex = fragmentIndex;
            ProductTypes = productTypes;
            CurrentPartition = currentPartition + 1;
            AddComplementaryIons = addCompIons;
            //Here use LowTheoreticalDiffAcceptor in practice doesn't work in 12/2/2017
            MassDiffAcceptor = new OpenSearchMode();
            DissociationTypes = DetermineDissociationType(productTypes);
            XLPrecusorMsTl = xlPrecusorMsTl;
            XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(xlPrecusorMsTl.Value);
            //if (XLBetaPrecusorMsTl.ToString().Contains("Absolute"))
            //{
            //    XLPrecusorSearchMode = new SingleAbsoluteAroundZeroSearchMode(XLPrecusorMsTl.Value);
            //}
            Crosslinker = crosslinker;
            CrosslinkSearchTop = crosslinkSearchTop;
            CrosslinkSearchTopNum = crosslinkSearchTopNum;
            QuenchH2O = quench_H2O;
            QuenchNH2 = quench_NH2;
            QuenchTris = quench_Tris;
            Charge_2_3 = charge_2_3;
            Charge_2_3_PrimeFragment = charge_2_3_PrimeFragment;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing crosslink search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)base.CommonParameters.ScoreCutoff;

            Parallel.ForEach(Partitioner.Create(0, ListOfSortedms2Scans.Length),
                new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                (range, loopState) =>
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
                    List<int> allBinsToSearch = GetBinsToSearch(scan);
                    List<BestPeptideScoreNotch> bestPeptideScoreNotchList = new List<BestPeptideScoreNotch>();
                    // get allowed theoretical masses from the known experimental mass
                    // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)
                    // this is just PRELIMINARY precursor-mass filtering
                    // additional checks are made later to ensure that the theoretical precursor mass is acceptable
                    var notches = MassDiffAcceptor.GetAllowedPrecursorMassIntervals(scan.PrecursorMass);

                    double lowestMassPeptideToLookFor = double.NegativeInfinity;
                    double highestMassPeptideToLookFor = double.PositiveInfinity;

                    double largestMassDiff = notches.Max(p => p.AllowedInterval.Maximum);
                    double smallestMassDiff = notches.Min(p => p.AllowedInterval.Minimum);

                    if (!double.IsInfinity(largestMassDiff))
                    {
                        double largestOppositeMassDiff = -1 * (notches.Max(p => p.AllowedInterval.Maximum) - scan.PrecursorMass);
                        lowestMassPeptideToLookFor = scan.PrecursorMass + largestOppositeMassDiff;
                    }
                    if (!double.IsNegativeInfinity(smallestMassDiff))
                    {
                        double smallestOppositeMassDiff = -1 * (notches.Min(p => p.AllowedInterval.Minimum) - scan.PrecursorMass);
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
                            var peptide = PeptideIndex[id];

                            double thePrecursorMass = scan.PrecursorMass;
                            int notch = MassDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMassIncludingFixedMods);
                            BestPeptideScoreNotch bestPeptideScoreNotch = new BestPeptideScoreNotch(peptide, scoringTable[id], notch);
                            bestPeptideScoreNotchList.Add(bestPeptideScoreNotch);
                        }

                        var possiblePsmCross = FindCrosslinkedPeptide(scan, bestPeptideScoreNotchList, i);
                        if (possiblePsmCross != null)
                        {
                            lock (GlobalPsmsCross)
                                GlobalPsmsCross.Add(possiblePsmCross);
                        }
                    }
                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedms2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing crosslink search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
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
                int obsFragmentFloorMass = (int)Math.Floor((CommonParameters.ProductMassTolerance.GetMinimumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);
                if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)
                {
                    obsFragmentFloorMass = obsPreviousFragmentCeilingMz;
                }
                int obsFragmentCeilingMass = (int)Math.Ceiling((CommonParameters.ProductMassTolerance.GetMaximumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);
                obsPreviousFragmentCeilingMz = obsFragmentCeilingMass + 1;
                for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                {
                    if (FragmentIndex[fragmentBin] != null)
                    {
                        binsToSearch.Add(fragmentBin);
                    }
                }
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
                {
                    break;
                }
                if (PeptideIndex[peptideIdsInThisBin[m]].MonoisotopicMassIncludingFixedMods < peptideMassToLookFor)
                {
                    l = m + 1;
                }
                else
                {
                    r = m - 1;
                }
            }
            if (m > 0)
            {
                m--;
            }
            return m;
        }

        private void IndexedScoring(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = FragmentIndex[binsToSearch[i]];

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
                        {
                            highestPeptideMassIndex = j;
                        }
                        else
                        {
                            break;
                        }
                    }
                }

                // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                {
                    int id = peptideIdsInThisBin[j];
                    scoringTable[id]++;

                    // add possible search results to the hashset of id's
                    if (scoringTable[id] == byteScoreCutoff && MassDiffAcceptor.Accepts(scanPrecursorMass, PeptideIndex[id].MonoisotopicMassIncludingFixedMods) >= 0)
                    {
                        idsOfPeptidesPossiblyObserved.Add(id);
                    }
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
                    var productMasses = theScanBestPeptide[ind].BestPeptide.ProductMassesMightHaveDuplicatesAndNaNs(ProductTypes);
                    Array.Sort(productMasses);
                    double score = CalculatePeptideScoreOld(theScan.TheScan, CommonParameters.ProductMassTolerance, productMasses, theScan.PrecursorMass, DissociationTypes, AddComplementaryIons, 0);
                    var psmCrossSingle = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, score, i, theScan, CommonParameters.DigestionParams);
                    psmCrossSingle.XLTotalScore = psmCrossSingle.Score;
                    psmCrossSingle.CrossType = PsmCrossType.Singe;

                    bestPsmCrossList.Add(psmCrossSingle);
                }
                //Deadend Peptide
                else if (QuenchTris && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + Crosslinker.DeadendMassTris) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, CommonParameters.DigestionParams);
                    var xlPos = PsmCross.XlPosCal(psmCrossEnd.CompactPeptide, crosslinkerModSitesAll);
                    if (xlPos.Count >= 1)
                    {
                        PsmCross.XlLocalization(theScan, psmCrossEnd, Crosslinker.DeadendMassTris, Crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, false, false, xlPos);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.XLBestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndTris;
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                else if (QuenchH2O && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + Crosslinker.DeadendMassH2O) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, CommonParameters.DigestionParams);
                    var xlPos = PsmCross.XlPosCal(psmCrossEnd.CompactPeptide, crosslinkerModSitesAll);
                    if (xlPos.Count >= 1)
                    {
                        PsmCross.XlLocalization(theScan, psmCrossEnd, Crosslinker.DeadendMassH2O, Crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, false, false, xlPos);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.XLBestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndH2O;
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                else if (QuenchNH2 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + Crosslinker.DeadendMassNH2) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, CommonParameters.DigestionParams);
                    var xlPos = PsmCross.XlPosCal(psmCrossEnd.CompactPeptide, crosslinkerModSitesAll);
                    if (xlPos.Count >= 1)
                    {
                        PsmCross.XlLocalization(theScan, psmCrossEnd, Crosslinker.DeadendMassNH2, Crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, false, false, xlPos);
                        psmCrossEnd.XLTotalScore = psmCrossEnd.XLBestScore;
                        psmCrossEnd.CrossType = PsmCrossType.DeadEndNH2;
                        bestPsmCrossList.Add(psmCrossEnd);
                    }
                }
                //loop peptide
                else if (Crosslinker.LoopMass != 0 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + Crosslinker.LoopMass) >= 0)
                {
                    var psmCrossLoop = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, CommonParameters.DigestionParams);
                    var xlPos = PsmCross.XlPosCal(psmCrossLoop.CompactPeptide, Crosslinker.CrosslinkerModSites);
                    if (xlPos.Count >= 2)
                    {
                        PsmCross.XlLocalizationForLoopCrosslink(theScan, psmCrossLoop, Crosslinker.LoopMass, Crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, xlPos);
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
                        if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, x + y + Crosslinker.TotalMass) >= 0)
                        {
                            var psmCrossAlpha = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan, base.CommonParameters.DigestionParams);
                            var psmCrossBeta = new PsmCross(theScanBestPeptide[inx].BestPeptide, theScanBestPeptide[inx].BestNotch, theScanBestPeptide[inx].BestScore, i, theScan, base.CommonParameters.DigestionParams);

                            PsmCross psmCross = null;
                            if (Crosslinker.CrosslinkerModSites == Crosslinker.CrosslinkerModSites2)
                            {
                                psmCross = XlSitesSame(theScan, psmCrossAlpha, psmCrossBeta, Crosslinker, ind, inx);
                            }
                            else
                            {
                                psmCross = XlSitesTwoDiff(theScan, psmCrossAlpha, psmCrossBeta, Crosslinker, ind, inx);
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
            CrosslinkerTypeClass crosslinker, int ind, int inx)
        {
            PsmCross psmCross = null;
            var xlPosAlpha = PsmCross.XlPosCal(psmCrossAlpha.CompactPeptide, crosslinker.CrosslinkerModSites);
            var xlPosBeta = PsmCross.XlPosCal(psmCrossBeta.CompactPeptide, crosslinker.CrosslinkerModSites);
            if (xlPosAlpha.Count >= 1 && xlPosBeta.Count >= 1)
            {
                PsmCross.XlLocalization(theScan, psmCrossAlpha, psmCrossBeta.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosAlpha);
                PsmCross.XlLocalization(theScan, psmCrossBeta, psmCrossAlpha.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosBeta);

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

            var xlPosAlpha1 = PsmCross.XlPosCal(psmCrossAlpha.CompactPeptide, crosslinker.CrosslinkerModSites);
            var xlPosBeta1 = PsmCross.XlPosCal(psmCrossBeta.CompactPeptide, crosslinker.CrosslinkerModSites);
            var xlPosAlpha2 = PsmCross.XlPosCal(psmCrossAlpha.CompactPeptide, crosslinker.CrosslinkerModSites2);
            var xlPosBeta2 = PsmCross.XlPosCal(psmCrossBeta.CompactPeptide, crosslinker.CrosslinkerModSites2);

            if ((xlPosAlpha1.Count >= 1 && xlPosBeta2.Count >= 1) || (xlPosAlpha2.Count >= 1 && xlPosBeta1.Count >= 1))
            {
                if (xlPosAlpha2.Count < 1 || xlPosBeta1.Count < 1)
                {
                    PsmCross.XlLocalization(theScan, psmCrossAlpha, psmCrossBeta.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosAlpha1);
                    PsmCross.XlLocalization(theScan, psmCrossBeta, psmCrossAlpha.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosBeta2);
                }

                if (xlPosAlpha1.Count < 1 || xlPosBeta2.Count < 1)
                {
                    PsmCross.XlLocalization(theScan, psmCrossAlpha, psmCrossBeta.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosAlpha2);
                    PsmCross.XlLocalization(theScan, psmCrossBeta, psmCrossAlpha.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosBeta1);
                }

                if ((xlPosAlpha1.Count >= 1 && xlPosBeta2.Count >= 1) && (xlPosAlpha2.Count >= 1 && xlPosBeta1.Count >= 1))
                {
                    var psmCrossAlpha1 = psmCrossAlpha;
                    var psmCrossAlpha2 = psmCrossAlpha;
                    var psmCrossBeta1 = psmCrossBeta;
                    var psmCrossBeta2 = psmCrossBeta;

                    PsmCross.XlLocalization(theScan, psmCrossAlpha1, psmCrossBeta.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosAlpha1);
                    PsmCross.XlLocalization(theScan, psmCrossBeta2, psmCrossAlpha.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosBeta2);
                    PsmCross.XlLocalization(theScan, psmCrossAlpha2, psmCrossBeta.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosAlpha2);
                    PsmCross.XlLocalization(theScan, psmCrossBeta1, psmCrossAlpha.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, ProductTypes, CommonParameters.ProductMassTolerance, Charge_2_3, Charge_2_3_PrimeFragment, xlPosBeta1);
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
    }
}