using EngineLayer.ModernSearch;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MassSpectrometry;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkSearchEngine : ModernSearchEngine
    {
        protected readonly List<CrosslinkSpectralMatch>[] GlobalCsms;

        // crosslinker molecule
        private readonly Crosslinker Crosslinker;

        private readonly bool CrosslinkSearchTopN;
        private readonly int TopN;
        private readonly bool QuenchH2O;
        private readonly bool QuenchNH2;
        private readonly bool QuenchTris;
        private MassDiffAcceptor XLPrecusorSearchMode;
        private Modification TrisDeadEnd;
        private Modification H2ODeadEnd;
        private Modification NH2DeadEnd;
        private Modification Loop;
        private readonly char[] AllCrosslinkerSites;
        private readonly List<int>[] SecondFragmentIndex;

        public CrosslinkSearchEngine(List<CrosslinkSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters, Crosslinker crosslinker, bool CrosslinkSearchTop, int CrosslinkSearchTopNum,
            bool quench_H2O, bool quench_NH2, bool quench_Tris, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.Crosslinker = crosslinker;
            this.CrosslinkSearchTopN = CrosslinkSearchTop;
            this.TopN = CrosslinkSearchTopNum;
            this.QuenchH2O = quench_H2O;
            this.QuenchNH2 = quench_NH2;
            this.QuenchTris = quench_Tris;
            SecondFragmentIndex = secondFragmentIndex;
            if (CommonParameters.ChildScanDissociationType!=DissociationType.Unknown && DissociationTypeGenerateSameTypeOfIons(CommonParameters.DissociationType, CommonParameters.ChildScanDissociationType))
            {
                SecondFragmentIndex = FragmentIndex;
            }
            GenerateCrosslinkModifications(crosslinker);
            AllCrosslinkerSites = Crosslinker.CrosslinkerModSites.ToCharArray().Concat(Crosslinker.CrosslinkerModSites2.ToCharArray()).Distinct().ToArray();

            if (commonParameters.PrecursorMassTolerance is PpmTolerance)
            {
                XLPrecusorSearchMode = new SinglePpmAroundZeroSearchMode(commonParameters.PrecursorMassTolerance.Value);
            }
            else
            {
                XLPrecusorSearchMode = new SingleAbsoluteAroundZeroSearchMode(commonParameters.PrecursorMassTolerance.Value);
            }
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing crosslink search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (scanIndex) =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();
                byte[] secondScoringTable = new byte[PeptideIndex.Count];
                List<int> childIdsOfPeptidesPossiblyObserved = new List<int>();

                for (; scanIndex < ListOfSortedMs2Scans.Length; scanIndex += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = ListOfSortedMs2Scans[scanIndex];

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(scan, FragmentIndex, CommonParameters.DissociationType);
                    List<BestPeptideScoreNotch> bestPeptideScoreNotchList = new List<BestPeptideScoreNotch>();

                    // first-pass scoring
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, Double.PositiveInfinity, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);

                    //child scan first-pass scoring
                    if (scan.ChildScans != null && CommonParameters.ChildScanDissociationType != DissociationType.LowCID)
                    {
                        Array.Clear(secondScoringTable, 0, secondScoringTable.Length);
                        childIdsOfPeptidesPossiblyObserved.Clear();

                        List<int> childBinsToSearch = new List<int>();
     
                        foreach (var aChildScan in scan.ChildScans)
                        {
                            var x = GetBinsToSearch(aChildScan, SecondFragmentIndex, CommonParameters.ChildScanDissociationType);
                            childBinsToSearch.AddRange(x);
                        }

                        IndexedScoring(SecondFragmentIndex, childBinsToSearch, secondScoringTable, byteScoreCutoff, childIdsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, Double.PositiveInfinity, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.ChildScanDissociationType);

                        foreach (var childId in childIdsOfPeptidesPossiblyObserved)
                        {
                            if (!idsOfPeptidesPossiblyObserved.Contains(childId))
                            {
                                idsOfPeptidesPossiblyObserved.Add(childId);
                            }
                            scoringTable[childId] = (byte)(scoringTable[childId] + secondScoringTable[childId]);
                        }
                    }

                    // done with indexed scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        if (CrosslinkSearchTopN)
                        {
                            // take top N hits for this scan
                            idsOfPeptidesPossiblyObserved = idsOfPeptidesPossiblyObserved.OrderByDescending(p => scoringTable[p]).Take(TopN).ToList();
                        }

                        foreach (var id in idsOfPeptidesPossiblyObserved)
                        {
                            PeptideWithSetModifications peptide = PeptideIndex[id];

                            int notch = MassDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMass);
                            bestPeptideScoreNotchList.Add(new BestPeptideScoreNotch(peptide, scoringTable[id], notch));
                        }

                        // combine individual peptide hits with crosslinker mass to find best crosslink PSM hit
                        var csms = FindCrosslinkedPeptide(scan, bestPeptideScoreNotchList, scanIndex);

                        if (csms == null || csms.Count == 0 || csms.Where(p => p != null).Count() == 0)
                        {
                            progress++;
                            continue;
                        }

                        if (GlobalCsms[scanIndex] == null)
                        {
                            GlobalCsms[scanIndex] = new List<CrosslinkSpectralMatch>();
                        }

                        GlobalCsms[scanIndex].AddRange(csms.Where(p => p != null).OrderByDescending(p => p.XLTotalScore));
                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing crosslink search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
                    }
                }
            });

            return new MetaMorpheusEngineResults(this);
        }

        /// <summary>
        /// 
        /// </summary>
        private void GenerateCrosslinkModifications(Crosslinker crosslinker)
        {
            ModificationMotif.TryGetMotif("X", out var motif);
            TrisDeadEnd = new Modification(_originalId: "Tris Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker.DeadendMassTris);
            H2ODeadEnd = new Modification(_originalId: "H2O Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker.DeadendMassH2O);
            NH2DeadEnd = new Modification(_originalId: "NH2 Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker.DeadendMassNH2);
            Loop = new Modification(_originalId: "Loop", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker.LoopMass);
        }

        /// <summary>
        /// 
        /// </summary>
        private List<CrosslinkSpectralMatch> FindCrosslinkedPeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int scanIndex)
        {
            List<CrosslinkSpectralMatch> possibleMatches = new List<CrosslinkSpectralMatch>();

            for (int alphaIndex = 0; alphaIndex < theScanBestPeptide.Count; alphaIndex++)
            {
                PeptideWithSetModifications bestPeptide = theScanBestPeptide[alphaIndex].BestPeptide;

                //Single Peptide
                if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass) >= 0)
                {
                    List<Product> products = bestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                    var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
                    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                    var psmCrossSingle = new CrosslinkSpectralMatch(bestPeptide, theScanBestPeptide[alphaIndex].BestNotch, score, scanIndex, theScan, CommonParameters.DigestionParams, matchedFragmentIons)
                    {
                        CrossType = PsmCrossType.Single,
                        XlRank = new List<int> { alphaIndex }
                    };

                    possibleMatches.Add(psmCrossSingle);
                }
                // Deadend Peptide
                else if (QuenchTris && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass + Crosslinker.DeadendMassTris) >= 0)
                {
                    List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);

                    if (possibleCrosslinkLocations.Any())
                    {
                        // tris deadend
                        possibleMatches.Add(LocalizeDeadEndSite(bestPeptide, theScan, CommonParameters, possibleCrosslinkLocations, TrisDeadEnd, theScanBestPeptide[alphaIndex].BestNotch, scanIndex, alphaIndex));
                    }
                }
                else if (QuenchH2O && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass + Crosslinker.DeadendMassH2O) >= 0)
                {
                    List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);

                    if (possibleCrosslinkLocations.Any())
                    {
                        // H2O deadend
                        possibleMatches.Add(LocalizeDeadEndSite(bestPeptide, theScan, CommonParameters, possibleCrosslinkLocations, H2ODeadEnd, theScanBestPeptide[alphaIndex].BestNotch, scanIndex, alphaIndex));
                    }
                }
                else if (QuenchNH2 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass + Crosslinker.DeadendMassNH2) >= 0)
                {
                    List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);

                    if (possibleCrosslinkLocations.Any())
                    {
                        // NH2 deadend
                        possibleMatches.Add(LocalizeDeadEndSite(bestPeptide, theScan, CommonParameters, possibleCrosslinkLocations, NH2DeadEnd, theScanBestPeptide[alphaIndex].BestNotch, scanIndex, alphaIndex));
                    }
                }
                // loop peptide
                else if (Crosslinker.LoopMass != 0 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass + Crosslinker.LoopMass) >= 0)
                {
                    List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);

                    if (possibleCrosslinkLocations.Count >= 2)
                    {
                        possibleMatches.Add(LocalizeLoopSites(bestPeptide, theScan, CommonParameters, possibleCrosslinkLocations, Loop, theScanBestPeptide[alphaIndex].BestNotch, scanIndex, alphaIndex));
                    }
                }
                // Cross-linked peptide
                else if (theScan.PrecursorMass - bestPeptide.MonoisotopicMass >= (CommonParameters.DigestionParams.MinPeptideLength * 50))
                {
                    List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);
                    if (!possibleCrosslinkLocations.Any())
                    {
                        continue;
                    }

                    PeptideWithSetModifications alphaPeptide = bestPeptide;

                    for (int betaIndex = 0; betaIndex < theScanBestPeptide.Count; betaIndex++)
                    {
                        PeptideWithSetModifications betaPeptide = theScanBestPeptide[betaIndex].BestPeptide;

                        if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, alphaPeptide.MonoisotopicMass + betaPeptide.MonoisotopicMass + Crosslinker.TotalMass) >= 0)
                        {
                            List<int> possibleBetaCrosslinkSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, betaPeptide);

                            if (!possibleBetaCrosslinkSites.Any())
                            {
                                continue;
                            }

                            CrosslinkSpectralMatch csm = LocalizeCrosslinkSites(theScan, theScanBestPeptide[alphaIndex], theScanBestPeptide[betaIndex], Crosslinker, alphaIndex, betaIndex);

                            possibleMatches.Add(csm);
                        }
                    }
                }
            }

            return possibleMatches;
        }

        /// <summary>
        /// Localizes the crosslink position on the alpha and beta peptides
        /// </summary>
        private CrosslinkSpectralMatch LocalizeCrosslinkSites(Ms2ScanWithSpecificMass theScan, BestPeptideScoreNotch alphaPeptide, BestPeptideScoreNotch betaPeptide, Crosslinker crosslinker, int ind, int inx)
        {
            CrosslinkSpectralMatch localizedCrosslinkedSpectralMatch = null;

            List<Tuple<List<int>, List<int>>> pairs = new List<Tuple<List<int>, List<int>>>();

            if (crosslinker.CrosslinkerModSites.Equals(crosslinker.CrosslinkerModSites2))
            {
                List<int> possibleAlphaXlSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), alphaPeptide.BestPeptide);
                List<int> possibleBetaXlSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), betaPeptide.BestPeptide);

                pairs.Add(new Tuple<List<int>, List<int>>(possibleAlphaXlSites, possibleBetaXlSites));
            }
            else
            {
                List<int> possibleAlphaXlSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), alphaPeptide.BestPeptide);
                List<int> possibleBetaXlSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites2.ToCharArray(), betaPeptide.BestPeptide);

                pairs.Add(new Tuple<List<int>, List<int>>(possibleAlphaXlSites, possibleBetaXlSites));

                possibleAlphaXlSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites2.ToCharArray(), alphaPeptide.BestPeptide);
                possibleBetaXlSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), betaPeptide.BestPeptide);

                pairs.Add(new Tuple<List<int>, List<int>>(possibleAlphaXlSites, possibleBetaXlSites));
            }

            foreach (var pair in pairs)
            {
                List<int> possibleAlphaXlSites = pair.Item1;
                List<int> possibleBetaXlSites = pair.Item2;

                if (possibleAlphaXlSites.Any() && possibleBetaXlSites.Any())
                {
                    int bestAlphaSite = 0;
                    int bestBetaSite = 0;
                    List<MatchedFragmentIon> bestMatchedAlphaIons = new List<MatchedFragmentIon>();
                    List<MatchedFragmentIon> bestMatchedBetaIons = new List<MatchedFragmentIon>();
                    Dictionary<int, List<MatchedFragmentIon>> bestMatchedChildAlphaIons = new Dictionary<int, List<MatchedFragmentIon>>();
                    Dictionary<int, List<MatchedFragmentIon>> bestMatchedChildBetaIons = new Dictionary<int, List<MatchedFragmentIon>>();
                    double bestAlphaLocalizedScore = 0;
                    double bestBetaLocalizedScore = 0;

                    var fragmentsForEachAlphaLocalizedPossibility = CrosslinkedPeptide.XlGetTheoreticalFragments(CommonParameters.DissociationType,
                        Crosslinker, possibleAlphaXlSites, betaPeptide.BestPeptide.MonoisotopicMass, alphaPeptide.BestPeptide).ToList();

                    foreach (int possibleSite in possibleAlphaXlSites)
                    {
                        foreach (var setOfFragments in fragmentsForEachAlphaLocalizedPossibility.Where(v => v.Item1 == possibleSite))
                        {
                            var matchedChildAlphaIons = new Dictionary<int, List<MatchedFragmentIon>>();
                            var matchedIons = MatchFragmentIons(theScan, setOfFragments.Item2, CommonParameters);
                            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

                            // search child scans (MS2+MS3)
                            foreach (Ms2ScanWithSpecificMass childScan in theScan.ChildScans)
                            {
                                var matchedChildIons = ScoreChildScan(theScan, childScan, possibleSite, alphaPeptide, betaPeptide);
                                
                                if (matchedChildIons == null)
                                {
                                    continue;
                                }

                                matchedChildAlphaIons.Add(childScan.OneBasedScanNumber, matchedChildIons);
                                //double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);

                                //score += childScore;
                            }

                            if (score > bestAlphaLocalizedScore)
                            {
                                bestAlphaLocalizedScore = score;
                                bestAlphaSite = possibleSite;
                                bestMatchedAlphaIons = matchedIons;
                                bestMatchedChildAlphaIons = matchedChildAlphaIons;
                            }
                        }
                    }

                    var fragmentsForEachBetaLocalizedPossibility = CrosslinkedPeptide.XlGetTheoreticalFragments(CommonParameters.DissociationType,
                        Crosslinker, possibleBetaXlSites, alphaPeptide.BestPeptide.MonoisotopicMass, betaPeptide.BestPeptide).ToList();

                    var alphaMz = new HashSet<double>(bestMatchedAlphaIons.Select(p => p.Mz));

                    foreach (int possibleSite in possibleBetaXlSites)
                    {
                        foreach (var setOfFragments in fragmentsForEachBetaLocalizedPossibility.Where(v => v.Item1 == possibleSite))
                        {
                            var matchedIons = MatchFragmentIons(theScan, setOfFragments.Item2, CommonParameters);
                            var matchedChildBetaIons = new Dictionary<int, List<MatchedFragmentIon>>();

                            // remove any matched beta ions that also matched to the alpha peptide
                            matchedIons.RemoveAll(p => alphaMz.Contains(p.Mz));

                            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

                            // search child scans (MS2+MS3)
                            foreach (Ms2ScanWithSpecificMass childScan in theScan.ChildScans)
                            {
                                var matchedChildIons = ScoreChildScan(theScan, childScan, possibleSite, betaPeptide, alphaPeptide);

                                if (matchedChildIons == null)
                                {
                                    continue;
                                }

                                matchedChildBetaIons.Add(childScan.OneBasedScanNumber, matchedChildIons);
                                //double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);

                                //score += childScore;
                            }

                            if (score > bestBetaLocalizedScore)
                            {
                                bestBetaLocalizedScore = score;
                                bestBetaSite = possibleSite;
                                bestMatchedBetaIons = matchedIons;
                                bestMatchedChildBetaIons = matchedChildBetaIons;
                            }
                        }
                    }

                    if (bestAlphaLocalizedScore < CommonParameters.ScoreCutoff ||
                        bestBetaLocalizedScore < CommonParameters.ScoreCutoff)
                    {
                        return null;
                    }

                    var localizedAlpha = new CrosslinkSpectralMatch(alphaPeptide.BestPeptide, alphaPeptide.BestNotch, bestAlphaLocalizedScore, 0, theScan, alphaPeptide.BestPeptide.DigestionParams, bestMatchedAlphaIons);
                    var localizedBeta = new CrosslinkSpectralMatch(betaPeptide.BestPeptide, betaPeptide.BestNotch, bestBetaLocalizedScore, 0, theScan, betaPeptide.BestPeptide.DigestionParams, bestMatchedBetaIons);

                    localizedAlpha.XlRank = new List<int> { ind, inx };
                    localizedAlpha.XLTotalScore = localizedAlpha.Score + localizedBeta.Score;
                    localizedAlpha.BetaPeptide = localizedBeta;

                    localizedAlpha.ChildMatchedFragmentIons = bestMatchedChildAlphaIons;
                    localizedBeta.ChildMatchedFragmentIons = bestMatchedChildBetaIons;

                    if (crosslinker.Cleavable)
                    {
                        //TODO: re-enable intensity ranks
                        //psmCrossAlpha.ParentIonMaxIntensityRanks = psmCrossAlpha.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).Select(p => p.IntensityRank).ToList();
                        //localizedAlpha.ParentIonMaxIntensityRanks = new List<int>();

                        //localizedAlpha.ParentIonExistNum = psmCrossAlpha.ParentIonMaxIntensityRanks.Count;
                    }

                    localizedAlpha.CrossType = PsmCrossType.Cross;
                    localizedCrosslinkedSpectralMatch = localizedAlpha;
                    localizedCrosslinkedSpectralMatch.LinkPositions = new List<int> { bestAlphaSite };
                    localizedCrosslinkedSpectralMatch.BetaPeptide.LinkPositions = new List<int> { bestBetaSite };
                }
            }

            return localizedCrosslinkedSpectralMatch;
        }

        private List<MatchedFragmentIon> ScoreChildScan(Ms2ScanWithSpecificMass parentScan, Ms2ScanWithSpecificMass childScan, int possibleSite, BestPeptideScoreNotch mainPeptide, BestPeptideScoreNotch otherPeptide)
        {
            bool shortMassAlphaMs3 = XLPrecusorSearchMode.Accepts(childScan.PrecursorMass, mainPeptide.BestPeptide.MonoisotopicMass + Crosslinker.CleaveMassShort) >= 0;
            bool longMassAlphaMs3 = XLPrecusorSearchMode.Accepts(childScan.PrecursorMass, mainPeptide.BestPeptide.MonoisotopicMass + Crosslinker.CleaveMassLong) >= 0;

            List<Product> childProducts;

            if (Crosslinker.Cleavable && (shortMassAlphaMs3 || longMassAlphaMs3))
            {
                double massToLocalize = shortMassAlphaMs3 ? Crosslinker.CleaveMassShort : Crosslinker.CleaveMassLong;
                if (mainPeptide.BestPeptide.AllModsOneIsNterminus.TryGetValue(possibleSite + 1, out var existingMod))
                {
                    massToLocalize += existingMod.MonoisotopicMass.Value;
                }

                Dictionary<int, Modification> mod = new Dictionary<int, Modification> { { possibleSite + 1, new Modification(_monoisotopicMass: massToLocalize) } };

                foreach (var otherExistingMod in mainPeptide.BestPeptide.AllModsOneIsNterminus.Where(p => p.Key != possibleSite + 1))
                {
                    mod.Add(otherExistingMod.Key, otherExistingMod.Value);
                }

                var peptideWithMod = new PeptideWithSetModifications(mainPeptide.BestPeptide.Protein, mainPeptide.BestPeptide.DigestionParams,
                    mainPeptide.BestPeptide.OneBasedStartResidueInProtein, mainPeptide.BestPeptide.OneBasedEndResidueInProtein,
                    mainPeptide.BestPeptide.CleavageSpecificityForFdrCategory, mainPeptide.BestPeptide.PeptideDescription,
                    mainPeptide.BestPeptide.MissedCleavages, mod, mainPeptide.BestPeptide.NumFixedMods);

                childProducts = peptideWithMod.Fragment(CommonParameters.ChildScanDissociationType, FragmentationTerminus.Both).ToList();
            }
            else if (Math.Abs(childScan.PrecursorMass - parentScan.PrecursorMass) < 0.01 && CommonParameters.DissociationType != CommonParameters.ChildScanDissociationType)
            {
                // same species got fragmented twice, the second time with a different dissociation type
                childProducts = CrosslinkedPeptide.XlGetTheoreticalFragments(CommonParameters.ChildScanDissociationType,
                    Crosslinker, new List<int> { possibleSite }, otherPeptide.BestPeptide.MonoisotopicMass, mainPeptide.BestPeptide).First().Item2;
            }
            else
            {
                return null;
            }

            var matchedChildIons = MatchFragmentIons(childScan, childProducts, CommonParameters);

            return matchedChildIons;
        }

        /// <summary>
        /// Localizes the deadend mod to a residue
        /// </summary>
        private CrosslinkSpectralMatch LocalizeDeadEndSite(PeptideWithSetModifications originalPeptide, Ms2ScanWithSpecificMass theScan, CommonParameters commonParameters,
            List<int> possiblePositions, Modification deadEndMod, int notch, int scanIndex, int peptideIndex)
        {
            double bestScore = 0;
            List<MatchedFragmentIon> bestMatchingFragments = new List<MatchedFragmentIon>();
            PeptideWithSetModifications bestLocalizedPeptide = null;
            int bestPosition = 0;

            foreach (int location in possiblePositions)
            {
                Dictionary<int, Modification> mods = originalPeptide.AllModsOneIsNterminus.ToDictionary(p => p.Key, p => p.Value);
                if (mods.ContainsKey(location + 1))
                {
                    var alreadyAnnotatedMod = mods[location + 1];
                    double combinedMass = mods[location + 1].MonoisotopicMass.Value + deadEndMod.MonoisotopicMass.Value;
                    Modification combinedMod = new Modification(_originalId: alreadyAnnotatedMod.OriginalId + "+" + deadEndMod.OriginalId, _modificationType: "Crosslink", _target: alreadyAnnotatedMod.Target, _locationRestriction: "Anywhere.", _monoisotopicMass: combinedMass);
                    mods[location + 1] = combinedMod;
                }
                else
                {
                    mods.Add(location + 1, deadEndMod);
                }

                var localizedPeptide = new PeptideWithSetModifications(originalPeptide.Protein, originalPeptide.DigestionParams, originalPeptide.OneBasedStartResidueInProtein,
                    originalPeptide.OneBasedEndResidueInProtein, originalPeptide.CleavageSpecificityForFdrCategory, originalPeptide.PeptideDescription, originalPeptide.MissedCleavages, mods, originalPeptide.NumFixedMods);

                var products = localizedPeptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                var matchedFragmentIons = MatchFragmentIons(theScan, products, commonParameters);

                double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                if (score > bestScore)
                {
                    bestMatchingFragments = matchedFragmentIons;
                    bestScore = score;
                    bestLocalizedPeptide = localizedPeptide;
                    bestPosition = location;
                }
            }

            if (bestScore < commonParameters.ScoreCutoff)
            {
                return null;
            }

            var csm = new CrosslinkSpectralMatch(bestLocalizedPeptide, notch, bestScore, scanIndex, theScan, originalPeptide.DigestionParams, bestMatchingFragments);

            if (deadEndMod == TrisDeadEnd)
            {
                csm.CrossType = PsmCrossType.DeadEndTris;
            }
            else if (deadEndMod == H2ODeadEnd)
            {
                csm.CrossType = PsmCrossType.DeadEndH2O;
            }
            else if (deadEndMod == NH2DeadEnd)
            {
                csm.CrossType = PsmCrossType.DeadEndNH2;
            }

            csm.LinkPositions = new List<int> { bestPosition };
            csm.XlRank = new List<int> { peptideIndex };

            return csm;
        }

        /// <summary>
        /// Localizes the loop to a begin and end residue
        /// </summary>
        private CrosslinkSpectralMatch LocalizeLoopSites(PeptideWithSetModifications originalPeptide, Ms2ScanWithSpecificMass theScan, CommonParameters commonParameters,
            List<int> possiblePositions, Modification loopMod, int notch, int scanIndex, int peptideIndex)
        {
            var possibleFragmentSets = CrosslinkedPeptide.XlLoopGetTheoreticalFragments(commonParameters.DissociationType, Loop, possiblePositions, originalPeptide);
            double bestScore = 0;
            Tuple<int, int> bestModPositionSites = null;
            List<MatchedFragmentIon> bestMatchingFragments = new List<MatchedFragmentIon>();

            foreach (var setOfPositions in possibleFragmentSets)
            {
                var matchedFragmentIons = MatchFragmentIons(theScan, setOfPositions.Value, commonParameters);

                double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                if (score > bestScore)
                {
                    bestMatchingFragments = matchedFragmentIons;
                    bestScore = score;
                    bestModPositionSites = setOfPositions.Key;
                }
            }

            if (bestScore < commonParameters.ScoreCutoff)
            {
                return null;
            }

            var csm = new CrosslinkSpectralMatch(originalPeptide, notch, bestScore, scanIndex, theScan, originalPeptide.DigestionParams, bestMatchingFragments)
            {
                CrossType = PsmCrossType.Loop,
                XlRank = new List<int> { peptideIndex },
                LinkPositions = new List<int> { bestModPositionSites.Item1, bestModPositionSites.Item2 }
            };

            return csm;
        }

        //TO DO: A better method can be implemented in mzLib.
        public static bool DissociationTypeGenerateSameTypeOfIons(DissociationType d, DissociationType childD)
        {
            if (d == childD)
            {
                return true;
            }
            if (d == DissociationType.CID && childD == DissociationType.HCD)
            {
                return true;
            }
            if (d == DissociationType.HCD && childD == DissociationType.CID)
            {
                return true;
            }
            if (d == DissociationType.ETD && childD == DissociationType.ECD)
            {
                return true;
            }
            if (d == DissociationType.ECD && childD == DissociationType.ETD)
            {
                return true;
            }
            return false;
        }
    }
}