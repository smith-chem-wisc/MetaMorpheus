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

        private const int max_mods_for_peptide = 3;
        private readonly List<int>[] fragmentIndex;

        private readonly Tolerance fragmentTolerance;

        private readonly float[] keys;

        private readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;

        private readonly List<CompactPeptide> peptideIndex;

        private readonly List<ProductType> lp;
        //Crosslink parameters

        private readonly int CrosslinkSearchTopNum;
        private readonly bool CrosslinkSearchWithCrosslinkerMod;
        private readonly Tolerance XLprecusorMsTl;
        private readonly Tolerance XLBetaPrecusorMsTl;
        private readonly Dictionary<ModificationWithMass, ushort> modsDictionary;
        private readonly List<Protein> proteinList;
        private readonly Protease protease;
        private readonly int maximumMissedCleavages;
        private readonly int? minPeptideLength;
        private readonly int? maxPeptideLength;
        private readonly List<ModificationWithMass> variableModifications;
        private readonly List<ModificationWithMass> fixedModifications;
        private readonly int maxModIsoforms;
        private CrosslinkerTypeClass crosslinker;
        private MassDiffAcceptor XLBetaSearchMode;

        #endregion Private Fields

        //AnalysisEngine parameters
        //private readonly Dictionary<ModificationWithMass, ushort> modsDictionary;

        #region Public Constructors

        public CrosslinkSearchEngine2(Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, Tolerance fragmentTolerance, CrosslinkerTypeClass crosslinker, int CrosslinkSearchTopNum, bool CrosslinkSearchWithCrosslinkerMod, Tolerance XLprecusorMsTl, Tolerance XLBetaPrecusorMsTl, Dictionary<ModificationWithMass, ushort> modsDictionary, List<ProductType> lp, List<Protein> proteinList, Protease protease, int maximumMissedCleavages, int? minPeptideLength, int? maxPeptideLength, List<ModificationWithMass> variableModifications, List<ModificationWithMass> fixedModifications, int maxModIsoforms, List<string> nestedIds) : base(nestedIds)
        {
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
            this.modsDictionary = modsDictionary;
            this.lp = lp;
            this.proteinList = proteinList;
            this.protease = protease;
            this.maximumMissedCleavages = maximumMissedCleavages;
            this.minPeptideLength = minPeptideLength;
            this.maxPeptideLength = maxPeptideLength;
            this.variableModifications = variableModifications;
            this.fixedModifications = fixedModifications;
            this.maxModIsoforms = maxModIsoforms;
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

            Status("In crosslink search engine...", nestedIds);

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
                            currentScanPsmParent[itop].TopPosition = new int[] { itop, 0 };
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
            List<Ms2ScanWithSpecificMass> selectedMS2Scans = new List<Ms2ScanWithSpecificMass>();
            for (int i = 0; i < newPsmsTop.Length; i++)
            {
                if (newPsmsTop[i] != null)
                {
                    for (int j = 0; j < newPsmsTop[i].Count; j++)
                    {
                        if (newPsmsTop[i][j].ScanPrecursorMass - newPsmsTop[i][j].CompactPeptide.MonoisotopicMassIncludingFixedMods - crosslinker.TotalMass > 100 && PsmCross.XlPosCal(newPsmsTop[i][j].CompactPeptide, crosslinker) != null)
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
            var xlPsmCross = ClassSearchTheBetaPeptide(selectedMS2Scans.ToArray(), allAlphaPsms.ToArray(), true);
            var xlPsmCrossRe = (from c in xlPsmCross
                                group c by c.Item1.ScanNumber into grp
                                select grp.OrderByDescending(c => c.Item1.XLBestScore + c.Item2.XLBestScore).FirstOrDefault()).ToList();
            return new CrosslinkSearchResults(xlPsmCrossRe, this);
        }

        #endregion Protected Methods

        #region Private Methods

        private void CalculatePeptideScores(IMsDataScan<IMzSpectrum<IMzPeak>> spectrum, double[] peptideScores)
        {
            for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
            {
                var theAdd = 1 + spectrum.MassSpectrum[i].Intensity / spectrum.TotalIonCurrent;
                var experimentalPeakInDaltons = spectrum.MassSpectrum[i].Mz - Constants.protonMass;
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

        private void XLCalculateTotalProductMassesMightHave(Ms2ScanWithSpecificMass theScan, PsmCross psmCross)
        {
            var modMass = theScan.PrecursorMass - psmCross.CompactPeptide.MonoisotopicMassIncludingFixedMods - crosslinker.TotalMass;
            int length = psmCross.CompactPeptide.NTerminalMasses.Length;
            var pmmh = psmCross.ProductMassesMightHaveDuplicatesAndNaNs(lp);
            ProductMassesMightHave pmmhTop = new ProductMassesMightHave();
            //int pos = -1;
            List<ProductMassesMightHave> pmmhList = new List<ProductMassesMightHave>();

            var linkPos = PsmCross.XlPosCal(psmCross.CompactPeptide, crosslinker);
            linkPos.Add(0);
            foreach (var ipos in linkPos)
            {
                //pos = ipos;
                ProductMassesMightHave pmmhCurr = new ProductMassesMightHave()
                {
                    Xlpos = ipos
                };
                List<double> x = new List<double>();
                List<string> y = new List<string>();
                if (crosslinker.Cleavable)
                {
                    x.Add(theScan.PrecursorMass - modMass - crosslinker.CleaveMassLong);
                    y.Add("PepS");
                    x.Add(theScan.PrecursorMass - modMass - crosslinker.CleaveMassShort);
                    y.Add("PepL");
                }
                for (int i = 0; i < pmmh.ProductMz.Length; i++)
                {
                    var cr = pmmh.ProductName[i][0];
                    var nm = Int32.Parse(System.Text.RegularExpressions.Regex.Match(pmmh.ProductName[i], @"\d+").Value);
                    if ((cr == 'b' || cr == 'c') && nm < ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i]);
                        y.Add(pmmh.ProductName[i]);
                    }
                    if ((cr == 'y' || cr == 'z') && nm < length - ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i]);
                        y.Add(pmmh.ProductName[i]);
                    }
                    if (cr == 'b' && nm >= ipos + 1)
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

                    if (cr == 'c' && nm >= ipos)
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

                    if (cr == 'y' && (nm >= length - ipos + 1))
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

                    if (cr == 'z' && (nm >= length - ipos + 1))
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

            List<double> scoreList = new List<double>();
            List<MatchedIonInfo> miil = new List<MatchedIonInfo>();
            foreach (var pmm in pmmhList)
            {
                var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(pmm.ProductMz.Length);
                double pmmScore = PsmCross.XLMatchIons(theScan.TheScan, fragmentTolerance, pmm.ProductMz, pmm.ProductName, matchedIonMassesListPositiveIsMatch);
                miil.Add(matchedIonMassesListPositiveIsMatch);
                scoreList.Add(pmmScore);
            }

            psmCross.XLBestScore = scoreList.Max();
            psmCross.MatchedIonInfo = miil[scoreList.IndexOf(scoreList.Max())];
            psmCross.Xlpos = pmmhList[scoreList.IndexOf(scoreList.Max())].Xlpos + 1;
        }

        private List<Tuple<PsmCross, PsmCross>> ClassSearchTheBetaPeptide(Ms2ScanWithSpecificMass[] selectedScan, PsmCross[] selectedPsmParent, bool conserveMemory = false)
        {
            double[] selectedScanPrecusor = selectedScan.Select(p => p.PrecursorMass).ToArray();
            double[] AlphaPeptidePrecusor = selectedPsmParent.Select(p => (double)p.CompactPeptide.MonoisotopicMassIncludingFixedMods).ToArray();
            double[] BetaPeptidePrecusor = selectedScanPrecusor.Zip(AlphaPeptidePrecusor, (one, two) => one - two - crosslinker.TotalMass).ToArray();
            //Sort the BetaPeptidePrecusor which is necessary for the GetAcceptableScans
            Array.Sort(BetaPeptidePrecusor.ToArray(), selectedPsmParent);
            Array.Sort(BetaPeptidePrecusor, selectedScan);

            Status("In xlclassic search engine!", nestedIds);

            int totalProteins = proteinList.Count;

            var observed_sequences = new HashSet<CompactPeptide>();

            Status("Getting ms2 scans...", nestedIds);

            var outerPsms = new PsmCross[selectedScan.Length];

            var lockObject = new object();
            int proteinsSeen = 0;
            int old_progress = 0;

            Status("Starting xlclassic search loop...", nestedIds);
            //Parallel.ForEach(Partitioner.Create(0, 1), partitionRange =>
            Parallel.ForEach(Partitioner.Create(0, totalProteins), partitionRange =>
            {
                var psms = new PsmCross[selectedScan.Length];
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    var protein = proteinList[i];
                    var digestedList = protein.Digest(protease, maximumMissedCleavages, minPeptideLength, maxPeptideLength, InitiatorMethionineBehavior.Variable, fixedModifications).ToList();
                    foreach (var peptide in digestedList)
                    {
                        var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, maxModIsoforms, max_mods_for_peptide).ToList();
                        foreach (var yyy in ListOfModifiedPeptides)
                        {
                            var correspondingCompactPeptide = yyy.CompactPeptide(ProductTypeToTerminusType.IdentifyTerminusType(lp));
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
                            //var productNames = correspondingCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp).ProductName;
                            Array.Sort(productMasses);
                            var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(productMasses.Length);

                            foreach (ScanWithIndexAndNotchInfo scanWithIndexAndNotchInfo in GetAcceptableScans(BetaPeptidePrecusor, yyy.MonoisotopicMass, XLBetaSearchMode, selectedScan).ToList())
                            {
                                var score = MatchIons(scanWithIndexAndNotchInfo.theScan.TheScan, fragmentTolerance, productMasses, matchedIonMassesListPositiveIsMatch.MatchedIonMz, true, yyy.MonoisotopicMass, lp);
                                if (score > 1)
                                {
                                    var psm = new PsmCross(correspondingCompactPeptide, scanWithIndexAndNotchInfo.notch, score, scanWithIndexAndNotchInfo.scanIndex, scanWithIndexAndNotchInfo.theScan);
                                    XLCalculateTotalProductMassesMightHave(scanWithIndexAndNotchInfo.theScan, psm);
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

            List<Tuple<PsmCross, PsmCross>> newPsmsTopTuple = new List<Tuple<PsmCross, PsmCross>>();
            for (int i = 0; i < outerPsms.Length; i++)
            {
                if (outerPsms[i] != null)
                {
                    XLCalculateTotalProductMassesMightHave(selectedScan[i], selectedPsmParent[i]);
                    selectedPsmParent[i].XLTotalScore = selectedPsmParent[i].XLBestScore + outerPsms[i].XLBestScore;
                    newPsmsTopTuple.Add(new Tuple<PsmCross, PsmCross>(selectedPsmParent[i], outerPsms[i]));
                }
            }
            return newPsmsTopTuple;
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