using Chemistry;
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
    public class CrosslinkSearchEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly List<int>[] fragmentIndex;

        private readonly Tolerance fragmentTolerance;

        private readonly float[] keys;

        private readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;

        private readonly List<CompactPeptide> peptideIndex;

        //Crosslink parameters
        private readonly CrosslinkerTypeClass crosslinker;

        private readonly int CrosslinkSearchTopNum;
        private readonly bool CrosslinkSearchWithCrosslinkerMod;
        private readonly Tolerance XLprecusorMsTl;
        private readonly Tolerance XLBetaPrecusorMsTl;
        private readonly List<ProductType> lp;
        private readonly Dictionary<ModificationWithMass, ushort> modsDictionary;
        private readonly List<PsmCross> psmCross;

        private MassDiffAcceptor XLPrecusorSearchMode;

        #endregion Private Fields

        #region Public Constructors

        public CrosslinkSearchEngine(List<PsmCross> psmCross, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, Tolerance fragmentTolerance, CrosslinkerTypeClass crosslinker, int CrosslinkSearchTopNum, bool CrosslinkSearchWithCrosslinkerMod, Tolerance XLprecusorMsTl, Tolerance XLBetaPrecusorMsTl, List<ProductType> lp, Dictionary<ModificationWithMass, ushort> modsDictionary, List<string> nestedIds) : base(nestedIds)
        {
            this.psmCross = psmCross;
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.keys = keys;
            this.fragmentIndex = fragmentIndex;
            this.fragmentTolerance = fragmentTolerance;
            this.crosslinker = crosslinker;
            this.CrosslinkSearchTopNum = CrosslinkSearchTopNum;
            this.CrosslinkSearchWithCrosslinkerMod = CrosslinkSearchWithCrosslinkerMod;
            this.XLprecusorMsTl = XLprecusorMsTl;
            this.XLBetaPrecusorMsTl = XLBetaPrecusorMsTl;
            //AnalysisEngine parameters
            this.lp = lp;
            this.modsDictionary = modsDictionary;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            MassDiffAcceptor XLsearchMode = new OpenSearchMode();
            XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(XLprecusorMsTl.Value);
            if (XLBetaPrecusorMsTl.ToString().Contains("Absolute"))
            {
                XLPrecusorSearchMode = new SingleAbsoluteAroundZeroSearchMode(XLprecusorMsTl.Value);
            }

            Status("In crosslink search engine...");

            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;

            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;

            //Crosslink data storage
            List<PsmCross> newPsmsCross = new List<PsmCross>();

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
                        //Function that find the two crosslinked peptide
                        var possiblePsmCross = FindCrosslinkedPeptide(thisScan, bestPeptideScoreNotch, i);

                        if (possiblePsmCross != null)
                        {
                            newPsmsCross.Add(possiblePsmCross);
                        }
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
            psmCross.AddRange(newPsmsCross);
            return new CrosslinkSearchResults(newPsmsCross, this);
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
                        if (fragmentTolerance.Within(experimentalPeakInDaltons, closestPeak))
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
                        if (fragmentTolerance.Within(experimentalPeakInDaltons, closestPeak))
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

        //Targetting function: to find two peptides that in the Top matched peptides
        private PsmCross FindCrosslinkedPeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int i)
        {
            List<PsmCross> bestPsmCrossList = new List<PsmCross>();
            PsmCross bestPsmCross = null;
            for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
            {
                //Single Peptide
                if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods) >= 0)
                {
                    var psmCrossSingle = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan);
                    psmCrossSingle.XLTotalScore = psmCrossSingle.Score;
                    psmCrossSingle.CrossType = PsmCrossType.Singe;
                    //if (bestPsmCrossList.Count >= 2)
                    //{
                    //    bestPsmCrossList.OrderByDescending(p => p.XLTotalScore);
                    //    bestPsmCrossList.RemoveAt(1);
                    //    bestPsmCrossList.Add(psmCrossSingle);
                    //}
                    //else
                    //{
                    //    bestPsmCrossList.Add(psmCrossSingle);
                    //}
                    bestPsmCrossList.Add(psmCrossSingle);
                }
                //Deadend Peptide
                else if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + 156.0786) >= 0)
                {
                    var psmCrossEnd = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan);
                    //The Score need to recaculate.
                    psmCrossEnd.XLTotalScore = psmCrossEnd.Score;
                    psmCrossEnd.CrossType = PsmCrossType.DeadEnd;
                    //if (bestPsmCrossList.Count >= 2)
                    //{
                    //    bestPsmCrossList.OrderByDescending(p => p.XLTotalScore);
                    //    bestPsmCrossList.RemoveAt(1);
                    //    bestPsmCrossList.Add(psmCrossEnd);
                    //}
                    //else
                    //{
                    //    bestPsmCrossList.Add(psmCrossEnd);
                    //}
                    bestPsmCrossList.Add(psmCrossEnd);
                }
                //loop peptide
                else if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods + 138.06808) >= 0)
                {
                    var psmCrossLoop = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan);
                    psmCrossLoop.XLTotalScore = psmCrossLoop.Score;
                    psmCrossLoop.CrossType = PsmCrossType.Loop;
                    //if (bestPsmCrossList.Count >= 2)
                    //{
                    //    bestPsmCrossList.OrderByDescending(p => p.XLTotalScore);
                    //    bestPsmCrossList.RemoveAt(1);
                    //    bestPsmCrossList.Add(psmCrossLoop);
                    //}
                    //else
                    //{
                    //    bestPsmCrossList.Add(psmCrossLoop);
                    //}
                    bestPsmCrossList.Add(psmCrossLoop);
                }
                //Cross-linked peptide
                else if (theScan.PrecursorMass - theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods >= 500)
                {
                    var x = theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods;
                    for (int inx = ind; inx < theScanBestPeptide.Count; inx++)
                    {
                        var y = theScanBestPeptide[inx].BestPeptide.MonoisotopicMassIncludingFixedMods;
                        if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, x + y + crosslinker.TotalMass) >= 0 && PsmCross.xlPosCal(theScanBestPeptide[ind].BestPeptide, crosslinker).Count != 0 && PsmCross.xlPosCal(theScanBestPeptide[inx].BestPeptide, crosslinker).Count != 0)
                        {
                            var psmCrossAlpha = new PsmCross(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, theScanBestPeptide[ind].BestScore, i, theScan);
                            var psmCrossBeta = new PsmCross(theScanBestPeptide[inx].BestPeptide, theScanBestPeptide[inx].BestNotch, theScanBestPeptide[inx].BestScore, i, theScan);
                            psmCrossAlpha.XlRank = new int[] { ind, inx };
                            PsmCross.XLCalculateTotalProductMassesMightHave(theScan, psmCrossAlpha, crosslinker, lp, fragmentTolerance);
                            PsmCross.XLCalculateTotalProductMassesMightHave(theScan, psmCrossBeta, crosslinker, lp, fragmentTolerance);
                            psmCrossAlpha.XLTotalScore = psmCrossAlpha.XLBestScore + psmCrossBeta.XLBestScore;
                            psmCrossAlpha.CrossType = PsmCrossType.Cross;
                            psmCrossAlpha.BetaPsmCross = psmCrossBeta;
                            //if (bestPsmCrossList.Count >= 2)
                            //{
                            //    bestPsmCrossList.OrderByDescending(p => p.XLTotalScore);
                            //    bestPsmCrossList.RemoveAt(1);
                            //    bestPsmCrossList.Add(psmCrossAlpha);
                            //}
                            //else
                            //{
                            //    bestPsmCrossList.Add(psmCrossAlpha);
                            //}
                            bestPsmCrossList.Add(psmCrossAlpha);
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
                else { bestPsmCross.DScore = bestPsmCross.XLTotalScore; }
            }

            return bestPsmCross;
        }

        #endregion Private Methods
    }
}