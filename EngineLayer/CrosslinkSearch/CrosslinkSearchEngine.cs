using EngineLayer.ModernSearch;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkSearchEngine : ModernSearchEngine
    {
        protected readonly CrosslinkSpectralMatch[] GlobalCsms;

        private OpenSearchType OpenSearchType;
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
        private char[] AllCrosslinkerSites;

        public CrosslinkSearchEngine(CrosslinkSpectralMatch[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, int currentPartition, CommonParameters commonParameters, OpenSearchType openSearchType, Crosslinker crosslinker, bool CrosslinkSearchTop, int CrosslinkSearchTopNum,
            bool quench_H2O, bool quench_NH2, bool quench_Tris, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.OpenSearchType = openSearchType;
            this.Crosslinker = crosslinker;
            this.CrosslinkSearchTopN = CrosslinkSearchTop;
            this.TopN = CrosslinkSearchTopNum;
            this.QuenchH2O = quench_H2O;
            this.QuenchNH2 = quench_NH2;
            this.QuenchTris = quench_Tris;
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

            if (OpenSearchType == OpenSearchType.NGlyco)
            {
                //groupedGlycans = GlobalVariables.NGlycans.GroupBy(p => p.Mass).ToDictionary(p => p.Key, p => p.ToList());
                Glycans = GlobalVariables.NGlycans.ToArray();
            }

            if (OpenSearchType == OpenSearchType.OGlyco)
            {
                var OGlycanBoxes = Glycan.BuildGlycanBoxes(GlobalVariables.OGlycans.ToList(), 3);
                GroupedOGlycanBoxes = OGlycanBoxes.GroupBy(p => p.Mass).ToDictionary(p=>p.Key, p=>p.ToList());
            }
        }

        //public Dictionary<double, List<Glycan>> groupedGlycans { get; }
        private Glycan[] Glycans {get;}
        public Dictionary<double, List<GlycanBox>> GroupedOGlycanBoxes { get; }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing crosslink search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = (byte)commonParameters.ScoreCutoff;

            Parallel.ForEach(Partitioner.Create(0, ListOfSortedMs2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, (range, loopState) =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();

                for (int scanIndex = range.Item1; scanIndex < range.Item2; scanIndex++)
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
                    var scan = ListOfSortedMs2Scans[scanIndex];

                    

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(scan);
                    List<BestPeptideScoreNotch> bestPeptideScoreNotchList = new List<BestPeptideScoreNotch>();

                    // first-pass scoring
                    IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, Double.PositiveInfinity, PeptideIndex, MassDiffAcceptor, 0);

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

                        CrosslinkSpectralMatch csm;
                        if (OpenSearchType == OpenSearchType.NGlyco)
                        {
                            //Glycopeptide Search
                            csm = FindNGlycopeptide(scan, bestPeptideScoreNotchList, scanIndex);
                        }
                        else if (OpenSearchType == OpenSearchType.Crosslink)
                        {                            
                            // combine individual peptide hits with crosslinker mass to find best crosslink PSM hit
                            csm = FindCrosslinkedPeptide(scan, bestPeptideScoreNotchList, scanIndex);
                        }
                        else
                        {
                            csm = FindOGlycopeptide(scan, bestPeptideScoreNotchList, scanIndex);
                        }


                        if (csm == null)
                        {
                            progress++;
                            continue;
                        }

                        // this scan might already have a hit from a different database partition; check to see if the score improves
                        if (GlobalCsms[scanIndex] == null || GlobalCsms[scanIndex].XLTotalScore < csm.XLTotalScore)
                        {
                            GlobalCsms[scanIndex] = csm;
                        }
                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing crosslink search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));
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
        private CrosslinkSpectralMatch FindCrosslinkedPeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int scanIndex)
        {
            List<CrosslinkSpectralMatch> possibleMatches = new List<CrosslinkSpectralMatch>();

            for (int alphaIndex = 0; alphaIndex < theScanBestPeptide.Count; alphaIndex++)
            {
                PeptideWithSetModifications bestPeptide = theScanBestPeptide[alphaIndex].BestPeptide;

                //Single Peptide
                if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass) >= 0)
                {
                    List<Product> products = bestPeptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                    var matchedFragmentIons = MatchFragmentIons(theScan, products, commonParameters);
                    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons, 0);

                    var psmCrossSingle = new CrosslinkSpectralMatch(bestPeptide, theScanBestPeptide[alphaIndex].BestNotch, score, scanIndex, theScan, commonParameters.DigestionParams, matchedFragmentIons);
                    psmCrossSingle.CrossType = PsmCrossType.Single;
                    psmCrossSingle.XlRank = new List<int> { alphaIndex };

                    possibleMatches.Add(psmCrossSingle);
                }
                // Deadend Peptide
                else if (QuenchTris && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass + Crosslinker.DeadendMassTris) >= 0)
                {
                    List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);

                    if (possibleCrosslinkLocations.Any())
                    {
                        // tris deadend
                        possibleMatches.Add(LocalizeDeadEndSite(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, TrisDeadEnd, theScanBestPeptide[alphaIndex].BestNotch, scanIndex, alphaIndex));
                    }
                }
                else if (QuenchH2O && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass + Crosslinker.DeadendMassH2O) >= 0)
                {
                    List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);

                    if (possibleCrosslinkLocations.Any())
                    {
                        // H2O deadend
                        possibleMatches.Add(LocalizeDeadEndSite(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, H2ODeadEnd, theScanBestPeptide[alphaIndex].BestNotch, scanIndex, alphaIndex));
                    }
                }
                else if (QuenchNH2 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass + Crosslinker.DeadendMassNH2) >= 0)
                {
                    List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);

                    if (possibleCrosslinkLocations.Any())
                    {
                        // NH2 deadend
                        possibleMatches.Add(LocalizeDeadEndSite(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, NH2DeadEnd, theScanBestPeptide[alphaIndex].BestNotch, scanIndex, alphaIndex));
                    }
                }
                // loop peptide
                else if (Crosslinker.LoopMass != 0 && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, bestPeptide.MonoisotopicMass + Crosslinker.LoopMass) >= 0)
                {
                    List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, bestPeptide);

                    if (possibleCrosslinkLocations.Count >= 2)
                    {
                        possibleMatches.Add(LocalizeLoopSites(bestPeptide, theScan, commonParameters, possibleCrosslinkLocations, Loop, theScanBestPeptide[alphaIndex].BestNotch, scanIndex, alphaIndex));
                    }
                }
                // Cross-linked peptide
                else if (theScan.PrecursorMass - bestPeptide.MonoisotopicMass >= (commonParameters.DigestionParams.MinPeptideLength * 50))
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

            // get the best match for this spectrum
            // bestPsmCross will be null if there are no valid hits
            possibleMatches.RemoveAll(v => v == null);
            possibleMatches = possibleMatches.OrderByDescending(p => p.XLTotalScore).ToList();
            var bestPsmCross = possibleMatches.FirstOrDefault();

            // resolve ambiguities
            if (bestPsmCross != null)
            {
                bestPsmCross.ResolveAllAmbiguities();

                if (bestPsmCross.BetaPeptide != null)
                {
                    bestPsmCross.BetaPeptide.ResolveAllAmbiguities();
                }
            }

            // calculate delta score
            if (possibleMatches.Count > 1)
            {
                bestPsmCross.DeltaScore = possibleMatches[0].XLTotalScore - possibleMatches[1].XLTotalScore;
            }

            return bestPsmCross;
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
                    double bestAlphaLocalizedScore = 0;
                    double bestBetaLocalizedScore = 0;

                    var fragmentsForEachAlphaLocalizedPossibility = CrosslinkedPeptide.XlGetTheoreticalFragments(commonParameters.DissociationType,
                        Crosslinker, possibleAlphaXlSites, betaPeptide.BestPeptide.MonoisotopicMass, alphaPeptide.BestPeptide).ToList();

                    foreach (int possibleSite in possibleAlphaXlSites)
                    {
                        foreach (var setOfFragments in fragmentsForEachAlphaLocalizedPossibility.Where(v => v.Item1 == possibleSite))
                        {
                            var matchedIons = MatchFragmentIons(theScan, setOfFragments.Item2, commonParameters);
                            double score = CalculatePeptideScore(theScan.TheScan, matchedIons, 0);

                            if (score > bestAlphaLocalizedScore)
                            {
                                bestAlphaLocalizedScore = score;
                                bestAlphaSite = possibleSite;
                                bestMatchedAlphaIons = matchedIons;
                            }
                        }
                    }

                    var fragmentsForEachBetaLocalizedPossibility = CrosslinkedPeptide.XlGetTheoreticalFragments(commonParameters.DissociationType,
                        Crosslinker, possibleBetaXlSites, alphaPeptide.BestPeptide.MonoisotopicMass, betaPeptide.BestPeptide).ToList();

                    var alphaMz = new HashSet<double>(bestMatchedAlphaIons.Select(p => p.Mz));

                    foreach (int possibleSite in possibleBetaXlSites)
                    {
                        foreach (var setOfFragments in fragmentsForEachBetaLocalizedPossibility.Where(v => v.Item1 == possibleSite))
                        {
                            var matchedIons = MatchFragmentIons(theScan, setOfFragments.Item2, commonParameters);

                            // remove any matched beta ions that also matched to the alpha peptide
                            matchedIons.RemoveAll(p => alphaMz.Contains(p.Mz));

                            double score = CalculatePeptideScore(theScan.TheScan, matchedIons, 0);

                            if (score > bestBetaLocalizedScore)
                            {
                                bestBetaLocalizedScore = score;
                                bestBetaSite = possibleSite;
                                bestMatchedBetaIons = matchedIons;
                            }
                        }
                    }

                    if (bestAlphaLocalizedScore < commonParameters.ScoreCutoff ||
                        bestBetaLocalizedScore < commonParameters.ScoreCutoff)
                    {
                        return null;
                    }

                    var localizedAlpha = new CrosslinkSpectralMatch(alphaPeptide.BestPeptide, alphaPeptide.BestNotch, bestAlphaLocalizedScore, 0, theScan, alphaPeptide.BestPeptide.DigestionParams, bestMatchedAlphaIons);
                    var localizedBeta = new CrosslinkSpectralMatch(betaPeptide.BestPeptide, betaPeptide.BestNotch, bestBetaLocalizedScore, 0, theScan, betaPeptide.BestPeptide.DigestionParams, bestMatchedBetaIons);

                    localizedAlpha.XlRank = new List<int> { ind, inx };
                    localizedAlpha.XLTotalScore = localizedAlpha.Score + localizedBeta.Score;
                    localizedAlpha.BetaPeptide = localizedBeta;

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

                double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons, 0);

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

                double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons, 0);

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

            var csm = new CrosslinkSpectralMatch(originalPeptide, notch, bestScore, scanIndex, theScan, originalPeptide.DigestionParams, bestMatchingFragments);
            csm.CrossType = PsmCrossType.Loop;
            csm.XlRank = new List<int> { peptideIndex };
            csm.LinkPositions = new List<int> { bestModPositionSites.Item1, bestModPositionSites.Item2 };

            return csm;
        }

        private CrosslinkSpectralMatch FindNGlycopeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int scanIndex)
        {
            
            List<CrosslinkSpectralMatch> possibleMatches = new List<CrosslinkSpectralMatch>();
            CrosslinkSpectralMatch bestPsmCross = null;

            if (GlycoPeptides.ScanOxoniumIonFilter(theScan, commonParameters.DissociationType))
            {
                var TopGlycans = GlycoPeptides.MatchBestGlycan(theScan, Glycans, commonParameters).Where(p => p != null && p.Item2>=2).OrderByDescending(p => p.Item2).Take(100).OrderBy(p => p.Item3).ToArray();

                if (TopGlycans.Count() >0 )
                {
                    for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
                    {
                        List<int> modPos = CrosslinkSpectralMatch.GetPossibleModSites(theScanBestPeptide[ind].BestPeptide, new string[] { "Nxt", "Nxs" });
                        if (modPos.Count < 1)
                        {
                            continue;
                        }

                        var possibleGlycanMassLow = theScan.PrecursorMass * (1 - 10E-5) - theScanBestPeptide[ind].BestPeptide.MonoisotopicMass;
                        if (possibleGlycanMassLow < 200 || possibleGlycanMassLow >= TopGlycans.Last().Item3)
                        {
                            continue;
                        }

                        int iDLow = GlycoPeptides.BinarySearchGetIndex(TopGlycans.Select(p => p.Item3).ToArray(), possibleGlycanMassLow);

                        while (iDLow < TopGlycans.Count() && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass + TopGlycans[iDLow].Item3) >= 0)
                        {
                            var fragmentsForEachGlycanLocalizedPossibility = GlycoPeptides.NGlyGetTheoreticalFragments(theScan, commonParameters.DissociationType, modPos, theScanBestPeptide[ind].BestPeptide, Glycans[TopGlycans[iDLow].Item1]).ToList();

                            double bestLocalizedScore = 0;
                            int bestSite = 0;
                            List<MatchedFragmentIon> bestMatchedIons = new List<MatchedFragmentIon>();
                            foreach (int possibleSite in modPos)
                            {
                                var setOfFragments = fragmentsForEachGlycanLocalizedPossibility.Where(v => v.Item1 == possibleSite).FirstOrDefault();

                                var matchedIons = MatchFragmentIons(theScan, setOfFragments.Item2, commonParameters);

                                double score = CalculatePeptideScore(theScan.TheScan, matchedIons, 0);

                                if (score > bestLocalizedScore)
                                {
                                    bestLocalizedScore = score;
                                    bestSite = possibleSite;
                                    bestMatchedIons = matchedIons;
                                }

                            }
                            var psmCross = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, bestLocalizedScore, scanIndex, theScan, commonParameters.DigestionParams, bestMatchedIons);
                            psmCross.Glycan = new List<Glycan> { Glycans[TopGlycans[iDLow].Item1] };
                            //psmCross.XLTotalScore = xLTotalScore;
                            psmCross.XlRank = new List<int> { ind };
                            psmCross.LinkPositions = new List<int> { bestSite }; //TO DO: ambiguity modification site
                            possibleMatches.Add(psmCross);

                            iDLow++;
                        }
                    }
                }
                
            }
            else
            {
                for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
                {
                    if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass) >= 0)
                    {
                        List<Product> products = theScanBestPeptide[ind].BestPeptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                        var matchedFragmentIons = MatchFragmentIons(theScan, products, commonParameters);
                        double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons, 0);

                        var psmCrossSingle = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, score, scanIndex, theScan, commonParameters.DigestionParams, matchedFragmentIons);
                        psmCrossSingle.CrossType = PsmCrossType.Single;
                        psmCrossSingle.XlRank = new List<int> { ind };

                        possibleMatches.Add(psmCrossSingle);
                    }
                }
            }

            if (possibleMatches.Count != 0)
            {
                possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
                bestPsmCross = possibleMatches.First();
                bestPsmCross.ResolveAllAmbiguities();
                if (possibleMatches.Count > 1)
                {
                    //This DeltaScore will be 0 if there are more than one glycan matched.
                    bestPsmCross.DeltaScore = Math.Abs(possibleMatches.First().Score - possibleMatches[1].Score);
                    for (int iPsm = 1; iPsm < possibleMatches.Count; iPsm++)
                    {
                        //TO DO: What if there are more than one peptide
                        if (possibleMatches[iPsm].Score == bestPsmCross.Score && possibleMatches[iPsm].Glycan!=null)
                        {
                             bestPsmCross.Glycan.Add(possibleMatches[iPsm].Glycan.First());                     
                        }
                    }
                }
            }

            return bestPsmCross;
        }

        private CrosslinkSpectralMatch FindOGlycopeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int scanIndex)
        {
            List<CrosslinkSpectralMatch> possibleMatches = new List<CrosslinkSpectralMatch>();
            CrosslinkSpectralMatch bestPsmCross = null;
            for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
            {
                if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass) >= 0)
                {
                    List<Product> products = theScanBestPeptide[ind].BestPeptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                    var matchedFragmentIons = MatchFragmentIons(theScan, products, commonParameters);
                    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons, 0);

                    var psmCrossSingle = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, score, scanIndex, theScan, commonParameters.DigestionParams, matchedFragmentIons);
                    psmCrossSingle.CrossType = PsmCrossType.Single;
                    psmCrossSingle.XlRank = new List<int> { ind };

                    possibleMatches.Add(psmCrossSingle);
                }
                //TO DO: add if the scan contains diagnostic ions
                else if ((theScan.PrecursorMass - theScanBestPeptide[ind].BestPeptide.MonoisotopicMass >= 200) && GlycoPeptides.ScanOxoniumIonFilter(theScan, commonParameters.DissociationType))
                {
                    //Using glycanBoxes
                    var possibleGlycanMass = theScan.PrecursorMass - theScanBestPeptide[ind].BestPeptide.MonoisotopicMass;
                    var iD = Array.BinarySearch(GroupedOGlycanBoxes.Keys.ToArray(), possibleGlycanMass);
                    if (iD < 0) { iD = ~iD; }
                    if (iD < GroupedOGlycanBoxes.Count && iD > 0)
                    {
                        if ((GroupedOGlycanBoxes.Keys.ElementAt(iD) - possibleGlycanMass) > (possibleGlycanMass - GroupedOGlycanBoxes.Keys.ElementAt(iD - 1))) { iD = iD - 1; }
                    }
                    if (iD >= GroupedOGlycanBoxes.Count)
                    {
                        iD = iD - 1;
                    }
                    if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass + GroupedOGlycanBoxes.Keys.ElementAt(iD)) >= 0)
                    {
                        List<int> modPos = CrosslinkSpectralMatch.GetPossibleModSites(theScanBestPeptide[ind].BestPeptide, new string[] { "S", "T" });
                        if (modPos.Count > 0)
                        {
                            foreach (var glycanBox in GroupedOGlycanBoxes.Values.ElementAt(iD))
                            {
                                if (glycanBox.NumberOfGlycans > modPos.Count)
                                {
                                    continue;
                                }

                                var fragmentsForEachGlycanLocalizedPossibility = GlycoPeptides.OGlyGetTheoreticalFragments(commonParameters.DissociationType, modPos, theScanBestPeptide[ind].BestPeptide, glycanBox).ToList();
                                double bestLocalizedScore = 0;
                                //double xLTotalScore = 0;

                                List<MatchedFragmentIon> bestMatchedIons = new List<MatchedFragmentIon>();

                                foreach (var setOfFragments in fragmentsForEachGlycanLocalizedPossibility)
                                {
                                    var matchedIons = MatchFragmentIons(theScan, setOfFragments.Item2.Item2, commonParameters);
                                    double score = CalculatePeptideScore(theScan.TheScan, matchedIons, 0);

                                    if (score > bestLocalizedScore)
                                    {
                                        bestLocalizedScore = score;
                                        bestMatchedIons = matchedIons;
                                        //xLTotalScore = CalculatePeptideScore(theScan.TheScan, matchedIons, 0);
                                    }
                                }


                                var psmCross = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, bestLocalizedScore, scanIndex, theScan, commonParameters.DigestionParams, bestMatchedIons);
                                psmCross.glycanBoxes = new List<GlycanBox> { glycanBox };
                                //psmCross.XLTotalScore = xLTotalScore;
                                psmCross.XlRank = new List<int> { ind };
                                possibleMatches.Add(psmCross);

                            }
                        }

                    }
                }
            }

            if (possibleMatches.Count != 0)
            {
                possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
                bestPsmCross = possibleMatches.First();
                bestPsmCross.ResolveAllAmbiguities();
                if (possibleMatches.Count > 1)
                {
                    //This DeltaScore will be 0 if there are more than one glycan matched.
                    bestPsmCross.DeltaScore = Math.Abs(possibleMatches.First().Score - possibleMatches[1].Score);
                    for (int iPsm = 1; iPsm < possibleMatches.Count; iPsm++)
                    {
                        //TO DO: What if there are more than one peptide
                        if (possibleMatches[iPsm].Score == bestPsmCross.Score && possibleMatches[iPsm].Glycan != null)
                        {
                            bestPsmCross.Glycan.Add(possibleMatches[iPsm].Glycan.First());
                        }
                    }
                }
            }

            return bestPsmCross;
        }

    }
}