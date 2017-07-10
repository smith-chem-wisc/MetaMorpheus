using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading.Tasks;
using System.Linq;
using Proteomics;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkSearchEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private const double tolInDaForPreferringHavingMods = 0.03;

        private readonly List<int>[] fragmentIndex;

        private readonly Tolerance fragmentTolerance;

        private readonly float[] keys;

        private readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;

        private readonly List<CompactPeptide> peptideIndex;

        private readonly MassDiffAcceptor searchMode;

        //Crosslink parameters
        private CrosslinkerTypeClass crosslinker;
        private readonly int CrosslinkSearchTopNum;
        private readonly bool CrosslinkSearchWithCrosslinkerMod;
        private readonly Tolerance XLprecusorMsTl;

        private readonly List<ProductType> lp;
        private readonly Dictionary<ModificationWithMass, ushort> modsDictionary;

        #endregion Private Fields

        #region Public Constructors

        public CrosslinkSearchEngine(Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, Tolerance fragmentTolerance, MassDiffAcceptor searchMode, CrosslinkerTypeClass crosslinker, int CrosslinkSearchTopNum, bool CrosslinkSearchWithCrosslinkerMod, Tolerance XLprecusorMsTl, List<ProductType> lp, Dictionary<ModificationWithMass, ushort> modsDictionary, List<string> nestedIds) : base(nestedIds)
        {
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.keys = keys;
            this.fragmentIndex = fragmentIndex;
            this.fragmentTolerance = fragmentTolerance;
            this.searchMode = searchMode;
            this.crosslinker = crosslinker;
            this.CrosslinkSearchTopNum = CrosslinkSearchTopNum;
            this.CrosslinkSearchWithCrosslinkerMod = CrosslinkSearchWithCrosslinkerMod;
            this.XLprecusorMsTl = XLprecusorMsTl;
            //AnalysisEngine parameters
            this.lp = lp;
            this.modsDictionary = modsDictionary;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            MassDiffAcceptor XLsearchMode = new OpenSearchMode();

            Status("In crosslink search engine...", nestedIds);

            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;

            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;

            //Crosslink data storage
            List<Tuple<PsmCross, PsmCross>> newPsmsTopTuple = new List<Tuple<PsmCross, PsmCross>>();

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
                            #region
                            if (bestPeptideScoreNotch != null && bestPeptideScoreNotch.Count == CrosslinkSearchTopNum)
                            {
                                // Full! Need to compare with old worst match
                                if (Math.Abs(currentWorstScore - consideredScore) < 1e-9)
                                {
                                    // Score is same as the worst, need to see if accepts and if prefer the new one
                                    int notch = XLsearchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                    if (notch >= 0 && FirstIsPreferableWithoutScore(candidatePeptide, bestPeptideScoreNotch.Last().BestPeptide, thisScanprecursorMass))
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
                            #endregion

                        }
                    }

                    //Create parameters to store current Crosslinked psm data?

                    if (bestPeptideScoreNotch != null)
                    {
                        //Function that find the two crosslinked peptide
                        var crosslinkPeptidePairList = FindCrosslinkedPeptide(thisScan, bestPeptideScoreNotch);

                        if (crosslinkPeptidePairList != null)
                        {
                            foreach (var crosslinkPeptidePair in crosslinkPeptidePairList)
                            {
                                var psmCross1 = new PsmCross(crosslinkPeptidePair.Item1.BestPeptide, crosslinkPeptidePair.Item1.BestNotch, crosslinkPeptidePair.Item1.BestScore, i, thisScan);
                                var psmCross2 = new PsmCross(crosslinkPeptidePair.Item2.BestPeptide, crosslinkPeptidePair.Item2.BestNotch, crosslinkPeptidePair.Item2.BestScore, i, thisScan);
                                if (psmCross1.CompactPeptide.BaseSequence.Contains((byte)crosslinker.CrosslinkerModSite.First()) && psmCross2.CompactPeptide.BaseSequence.Contains((byte)crosslinker.CrosslinkerModSite.First()))
                                {
                                    XLCalculateTotalProductMassesMightHave(thisScan, psmCross1);
                                    XLCalculateTotalProductMassesMightHave(thisScan, psmCross2);
                                    var currentTuplePair = new Tuple<PsmCross, PsmCross>(psmCross1, psmCross2);
                                    newPsmsTopTuple.Add(currentTuplePair);
                                }
                            }
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
            return new CrosslinkSearchResults(newPsmsTopTuple, this);
        }

        #endregion Protected Methods

        #region Private Methods

        // Want this to return false more!! So less computation is done. So second is preferable more often.
        private static bool FirstIsPreferableWithoutScore(CompactPeptide first, CompactPeptide second, double pm)
        {
            if (Math.Abs(first.MonoisotopicMassIncludingFixedMods - pm) < tolInDaForPreferringHavingMods && Math.Abs(second.MonoisotopicMassIncludingFixedMods - pm) > tolInDaForPreferringHavingMods)
                return true;
            if (Math.Abs(first.MonoisotopicMassIncludingFixedMods - pm) > tolInDaForPreferringHavingMods && Math.Abs(second.MonoisotopicMassIncludingFixedMods - pm) < tolInDaForPreferringHavingMods)
                return false;

            if (first.varMod1Type == 0 && second.varMod1Type > 0)
                return true;
            if (first.varMod1Type > 0 && second.varMod1Type == 0)
                return false;
            if (first.varMod2Type == 0 && second.varMod2Type > 0)
                return true;
            if (first.varMod2Type > 0 && second.varMod2Type == 0)
                return false;
            if (first.varMod3Type == 0 && second.varMod3Type > 0)
                return true;
            if (first.varMod3Type > 0 && second.varMod3Type == 0)
                return false;

            return false;
        }

        private void CalculatePeptideScores(IMsDataScan<IMzSpectrum<IMzPeak>> spectrum, double[] peptideScores)
        {
            foreach (var experimentalPeak in spectrum.MassSpectrum)
            {
                var theAdd = 1 + experimentalPeak.Intensity / spectrum.TotalIonCurrent;
                var experimentalPeakInDaltons = experimentalPeak.Mz - Constants.protonMass;
                float closestPeak = float.NaN;
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
        private List<Tuple<BestPeptideScoreNotch, BestPeptideScoreNotch>> FindCrosslinkedPeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide)
        {
            var bestPeptideScoreNotchList = new List<Tuple<BestPeptideScoreNotch, BestPeptideScoreNotch>>();

            for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
            {

                var x = theScanBestPeptide[ind].BestPeptide.MonoisotopicMassIncludingFixedMods;
                for (int inx = ind; inx < theScanBestPeptide.Count; inx++)
                {
                    var y = theScanBestPeptide[inx].BestPeptide.MonoisotopicMassIncludingFixedMods;
                    if (crosslinker.Cleavable == true && CrosslinkSearchWithCrosslinkerMod == true)
                    {
                        if (XLprecusorMsTl.Within(theScan.PrecursorMass, x + y + crosslinker.TotalMass - crosslinker.CleaveMassShort * 2)
                            || XLprecusorMsTl.Within(theScan.PrecursorMass, x + y + crosslinker.TotalMass + crosslinker.CleaveMassShort * 2)
                            || XLprecusorMsTl.Within(theScan.PrecursorMass, x + y))
                        {
                            Tuple<BestPeptideScoreNotch, BestPeptideScoreNotch> BestPeptideScoreNotchPair = new Tuple<BestPeptideScoreNotch, BestPeptideScoreNotch>(theScanBestPeptide[ind], theScanBestPeptide[inx]);
                            bestPeptideScoreNotchList.Add(BestPeptideScoreNotchPair);
                            return bestPeptideScoreNotchList;
                        }
                    }
                    else
                    {
                        if (XLprecusorMsTl.Within(theScan.PrecursorMass, x + y + crosslinker.TotalMass))
                        {
                            Tuple<BestPeptideScoreNotch, BestPeptideScoreNotch> BestPeptideScoreNotchPair = new Tuple<BestPeptideScoreNotch, BestPeptideScoreNotch>(theScanBestPeptide[ind], theScanBestPeptide[inx]);
                            bestPeptideScoreNotchList.Add(BestPeptideScoreNotchPair);
                            return bestPeptideScoreNotchList;
                        }
                    }
                }
            }
            return bestPeptideScoreNotchList;
        }

        private void XLCalculateTotalProductMassesMightHave(Ms2ScanWithSpecificMass theScan, PsmCross psmCross)
        {
            var modMass = theScan.PrecursorMass - psmCross.CompactPeptide.MonoisotopicMassIncludingFixedMods - crosslinker.TotalMass;
            int length = psmCross.CompactPeptide.BaseSequence.Length;
            var pmmh = psmCross.ProductMassesMightHaveDuplicatesAndNaNs(lp, modsDictionary);
            PsmCross.ProductMassesMightHave pmmhNew = new PsmCross.ProductMassesMightHave();
            int pos = -1;
            List<PsmCross.ProductMassesMightHave> pmmhList = new List<PsmCross.ProductMassesMightHave>();
            for (int ipos = 0; ipos < psmCross.CompactPeptide.BaseSequence.Length; ipos++)
            {
                if (psmCross.CompactPeptide.BaseSequence[ipos] == (byte)crosslinker.CrosslinkerModSite.First())
                {
                    pos = ipos;
                    PsmCross.ProductMassesMightHave pmmhCurr = new PsmCross.ProductMassesMightHave();
                    List<double> x = new List<double>();
                    List<string> y = new List<string>();
                    for (int i = 0; i < pmmh.ProductMz.Length; i++)
                    {
                        var cr = pmmh.ProductName[i][0];
                        var nm = Int32.Parse(System.Text.RegularExpressions.Regex.Match(pmmh.ProductName[i], @"\d+").Value);
                        if (crosslinker.Cleavable)
                        {
                            x.Add(theScan.PrecursorMass - modMass - crosslinker.CleaveMassLong);
                            y.Add("PepS");
                            x.Add(theScan.PrecursorMass - modMass - crosslinker.CleaveMassShort);
                            y.Add("PepL");
                        }
                        if ((cr == 'b' || cr == 'c') && nm < pos + 1)
                        {
                            x.Add(pmmh.ProductMz[i]);
                            y.Add(pmmh.ProductName[i]);
                        }
                        if ((cr == 'y' || cr == 'z') && nm < length - pos)
                        {
                            x.Add(pmmh.ProductMz[i]);
                            y.Add(pmmh.ProductName[i]);
                        }
                        if (cr == 'b' && nm >= pos + 1)
                        {
                            x.Add(pmmh.ProductMz[i] + modMass + crosslinker.TotalMass);
                            y.Add("t1b" + nm.ToString());

                            x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 2);
                            y.Add("t2b" + nm.ToString());

                            if (crosslinker.Cleavable)
                            {
                                x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                                y.Add("sb" + nm.ToString());

                                x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                                y.Add("lb" + nm.ToString());
                            }
                        }

                        if (cr == 'c' && nm >= pos)
                        {
                            x.Add(pmmh.ProductMz[i] + modMass + crosslinker.TotalMass);
                            y.Add("t1c" + nm.ToString());

                            x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 2);
                            y.Add("t2c" + nm.ToString());

                            if (crosslinker.Cleavable)
                            {
                                x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                                y.Add("sc" + nm.ToString());

                                x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                                y.Add("lc" + nm.ToString());
                            }
                        }

                        if (cr == 'y' && (nm >= length - pos))
                        {
                            x.Add(pmmh.ProductMz[i] + modMass + crosslinker.TotalMass);
                            y.Add("t1y" + nm.ToString());

                            x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 2);
                            y.Add("t2y" + nm.ToString());

                            if (crosslinker.Cleavable)
                            {
                                x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                                y.Add("sy" + nm.ToString());

                                x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                                y.Add("ly" + nm.ToString());
                            }
                        }

                        if (cr == 'z' && (nm >= length - pos))
                        {
                            x.Add(pmmh.ProductMz[i] + modMass + crosslinker.TotalMass);
                            y.Add("t1z" + nm.ToString());

                            x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 2);
                            y.Add("t2z" + nm.ToString());

                            if (crosslinker.Cleavable)
                            {
                                x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                                y.Add("sz" + nm.ToString());

                                x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                                y.Add("lz" + nm.ToString());
                            }
                        }
                    }
                    pmmhCurr.ProductMz = x.ToArray();
                    pmmhCurr.ProductName = y.ToArray();
                    Array.Sort(pmmhCurr.ProductMz, pmmhCurr.ProductName);
                    pmmhList.Add(pmmhCurr);
                }
            }

            //If the peptide did contain the crosslink amino acid
            //if (pos != -1)
            //{
            List<double> scoreList = new List<double>();
            List<MatchedIonInfo> miil = new List<MatchedIonInfo>();
            foreach (var pmm in pmmhList)
            {
                var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(pmm.ProductMz.Length);
                double pmmScore = PsmParent.MatchIons(theScan.TheScan, fragmentTolerance, pmm.ProductMz, matchedIonMassesListPositiveIsMatch.MatchedIonMz);
                miil.Add(matchedIonMassesListPositiveIsMatch);
                scoreList.Add(pmmScore);
            }

            pmmhNew = pmmhList[scoreList.IndexOf(scoreList.Max())];
            //psmCross.pmmh = pmmhNew;
            psmCross.XLBestScore = scoreList.Max();
            psmCross.matchedIonInfo = miil[scoreList.IndexOf(scoreList.Max())];
            //}

        }

        #endregion Private Methods
    }
}

