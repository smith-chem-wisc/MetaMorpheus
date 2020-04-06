using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkSearchEngine : ModernSearchEngine
    {
        public static readonly double ToleranceForMassDifferentiation = 1e-9;
        protected readonly List<CrosslinkSpectralMatch>[] GlobalCsms;

        // crosslinker molecule
        private readonly Crosslinker Crosslinker;

        private readonly int TopN;
        private readonly bool CleaveAtCrosslinkSite;
        private readonly bool QuenchH2O;
        private readonly bool QuenchNH2;
        private readonly bool QuenchTris;
        private readonly MassDiffAcceptor XLPrecusorSearchMode;
        private readonly MassDiffAcceptor XLProductSearchMode;
        private Modification TrisDeadEnd;
        private Modification H2ODeadEnd;
        private Modification NH2DeadEnd;
        private Modification Loop;
        private readonly char[] AllCrosslinkerSites;
        private readonly List<int>[] SecondFragmentIndex;
        private readonly double[] massTable;

        public CrosslinkSearchEngine(List<CrosslinkSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, Crosslinker crosslinker, int CrosslinkSearchTopNum,
            bool CleaveAtCrosslinkSite, bool quench_H2O, bool quench_NH2, bool quench_Tris, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, fileSpecificParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.Crosslinker = crosslinker;
            this.TopN = CrosslinkSearchTopNum;
            this.CleaveAtCrosslinkSite = CleaveAtCrosslinkSite;
            this.QuenchH2O = quench_H2O;
            this.QuenchNH2 = quench_NH2;
            this.QuenchTris = quench_Tris;

            this.massTable = peptideIndex.Select(p => p.MonoisotopicMass).ToArray();

            SecondFragmentIndex = secondFragmentIndex;
            if (CommonParameters.MS2ChildScanDissociationType!=DissociationType.Unknown && DissociationTypeGenerateSameTypeOfIons(CommonParameters.DissociationType, CommonParameters.MS2ChildScanDissociationType))
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

            if (commonParameters.ProductMassTolerance is PpmTolerance)
            {
                XLProductSearchMode = new SinglePpmAroundZeroSearchMode(commonParameters.ProductMassTolerance.Value);
            }
            else
            {
                XLProductSearchMode = new SingleAbsoluteAroundZeroSearchMode(commonParameters.ProductMassTolerance.Value);
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

                List<CrosslinkSpectralMatch> csms = new List<CrosslinkSpectralMatch>();
                HashSet<Tuple<int, int>> seenPair = new HashSet<Tuple<int, int>>();
                List<int> idsOfPeptidesTopN = new List<int>();
                byte scoreAtTopN = 0;
                int peptideCount = 0;

                for (; scanIndex < ListOfSortedMs2Scans.Length; scanIndex += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    csms.Clear();
                    seenPair.Clear();
                    idsOfPeptidesTopN.Clear();

                    var scan = ListOfSortedMs2Scans[scanIndex];

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(scan, FragmentIndex, CommonParameters.DissociationType);

                    //Limit the high bound limitation, here assume it is possible to has max 3 Da shift. This allows for correcting precursor in the future.
                    var high_bound_limitation = scan.PrecursorMass + 4;

                    // first-pass scoring
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);

                    //child scan first-pass scoring
                    if (scan.ChildScans != null && CommonParameters.MS2ChildScanDissociationType != DissociationType.Unknown && CommonParameters.MS2ChildScanDissociationType != DissociationType.LowCID)
                    {
                        Array.Clear(secondScoringTable, 0, secondScoringTable.Length);
                        childIdsOfPeptidesPossiblyObserved.Clear();

                        List<int> childBinsToSearch = new List<int>();
     
                        foreach (var aChildScan in scan.ChildScans)
                        {
                            var x = GetBinsToSearch(aChildScan, SecondFragmentIndex, CommonParameters.MS2ChildScanDissociationType);
                            childBinsToSearch.AddRange(x);
                        }

                        IndexedScoring(SecondFragmentIndex, childBinsToSearch, secondScoringTable, byteScoreCutoff, childIdsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.MS2ChildScanDissociationType);

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
                        scoreAtTopN = 0;
                        peptideCount = 0;

                        foreach (int id in idsOfPeptidesPossiblyObserved.OrderByDescending(p => scoringTable[p]))
                        {
                            peptideCount++;
                            // Whenever the count exceeds the TopN that we want to keep, we removed everything with a score lower than the score of the TopN-th peptide in the ids list
                            if (peptideCount == TopN)
                            {
                                scoreAtTopN = scoringTable[id];
                            }

                            if (scoringTable[id] < scoreAtTopN)
                            {
                                break;
                            }
                            idsOfPeptidesTopN.Add(id);
                        }
                        
                        //peptide candidates in idsOfPeptidesTopN are treated as alpha peptides. Then the mass of the beta peptides are calculated and searched from massTable.
                        FindCrosslinkedPeptide(scan, idsOfPeptidesTopN, scoringTable, massTable, byteScoreCutoff, scanIndex, csms, seenPair);

                        if (csms == null || csms.Count == 0 || csms.Where(p => p != null).Count() == 0)
                        {
                            progress++;
                            continue;
                        }

                        if (GlobalCsms[scanIndex] == null)
                        {
                            GlobalCsms[scanIndex] = new List<CrosslinkSpectralMatch>();
                        }
                        //scans are sorted here so that the first one in the list is the top score. The others are ignored for now. We sort it here because this step is parallelized.
                        GlobalCsms[scanIndex].AddRange(csms.Where(p => p != null));
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
        /// Deadend and loop crosslink are changed to a type of modification and searched as a modified peptide.
        /// </summary>
        /// <param name="crosslinker"></param>
        private void GenerateCrosslinkModifications(Crosslinker crosslinker)
        {
            ModificationMotif.TryGetMotif("X", out var motif);
            TrisDeadEnd = new Modification(_originalId: "Tris Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker.DeadendMassTris);
            H2ODeadEnd = new Modification(_originalId: "H2O Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker.DeadendMassH2O);
            NH2DeadEnd = new Modification(_originalId: "NH2 Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker.DeadendMassNH2);
            Loop = new Modification(_originalId: "Loop", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: Crosslinker.LoopMass);
        }


        private void FindCrosslinkedPeptide(Ms2ScanWithSpecificMass scan, List<int> idsOfPeptidesPossiblyObserved, byte[] scoringTable, double[] massTable, byte byteScoreCutoff, int scanIndex, List<CrosslinkSpectralMatch> possibleMatches, HashSet<Tuple<int, int>> seenPair)
        {
            //The code here is to generate intensity ranks of signature ions for cleavable crosslink.
            int rank = 1;
            double[] experimentFragmentMasses = null;
            double[] experimentFragmentIntensities = null;
            int[] intensityRanks = null;
            if (Crosslinker.Cleavable)
            {
                experimentFragmentMasses = scan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray();
                experimentFragmentIntensities = scan.ExperimentalFragments.Select(p => p.Peaks.Sum(q => q.intensity)).ToArray();
                intensityRanks = CrosslinkSpectralMatch.GenerateIntensityRanks(experimentFragmentIntensities);
            }

            List<Product> products = new List<Product>();
            foreach (var id in idsOfPeptidesPossiblyObserved)
            {
                List<int> possibleCrosslinkLocations = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, PeptideIndex[id], CommonParameters.DigestionParams.InitiatorMethionineBehavior, CleaveAtCrosslinkSite);
                
                if (XLPrecusorSearchMode.Accepts(scan.PrecursorMass, massTable[id]) >= 0)
                {
                    PeptideIndex[id].Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, products);
                    var matchedFragmentIons = MatchFragmentIons(scan, products, CommonParameters);
                    double score = CalculatePeptideScore(scan.TheScan, matchedFragmentIons);

                    if (score > byteScoreCutoff)
                    {
                        var psmCrossSingle = new CrosslinkSpectralMatch(PeptideIndex[id], 0, score, scanIndex, scan, CommonParameters, matchedFragmentIons)
                        {
                            CrossType = PsmCrossType.Single,
                        };

                        possibleMatches.Add(psmCrossSingle);
                    }
                }
                else if (QuenchTris && XLPrecusorSearchMode.Accepts(scan.PrecursorMass, massTable[id] + Crosslinker.DeadendMassTris) >= 0)
                {
                    if (possibleCrosslinkLocations != null)
                    {
                        // tris deadend
                        possibleMatches.Add(LocalizeDeadEndSite(PeptideIndex[id], scan, CommonParameters, possibleCrosslinkLocations, TrisDeadEnd, 0, scanIndex));
                    }
                }
                else if (QuenchH2O && XLPrecusorSearchMode.Accepts(scan.PrecursorMass, massTable[id] + Crosslinker.DeadendMassH2O) >= 0)
                {
                    if (possibleCrosslinkLocations != null)
                    {
                        // H2O deadend
                        possibleMatches.Add(LocalizeDeadEndSite(PeptideIndex[id], scan, CommonParameters, possibleCrosslinkLocations, H2ODeadEnd, 0, scanIndex));
                    }
                }
                else if (QuenchNH2 && XLPrecusorSearchMode.Accepts(scan.PrecursorMass, massTable[id] + Crosslinker.DeadendMassNH2) >= 0)
                {
                    if (possibleCrosslinkLocations != null)
                    {
                        // NH2 deadend
                        possibleMatches.Add(LocalizeDeadEndSite(PeptideIndex[id], scan, CommonParameters, possibleCrosslinkLocations, NH2DeadEnd, 0, scanIndex));
                    }
                }
                else if (Crosslinker.LoopMass != 0 && XLPrecusorSearchMode.Accepts(scan.PrecursorMass, massTable[id] + Crosslinker.LoopMass) >= 0)
                {
                    //TO THINK: Is there any cases that Loop Mass equals dead-end mass.
                    if (possibleCrosslinkLocations != null && possibleCrosslinkLocations.Count >= 2)
                    {
                        possibleMatches.Add(LocalizeLoopSites(PeptideIndex[id], scan, CommonParameters, possibleCrosslinkLocations, 0, scanIndex));
                    }
                }
                else if (scan.PrecursorMass - massTable[id] >= (CommonParameters.DigestionParams.MinPeptideLength * 50))
                {
                   if (possibleCrosslinkLocations == null)
                    {
                        continue;
                    }

                    double betaMass = scan.PrecursorMass - massTable[id] - Crosslinker.TotalMass;

                    double betaMassLow = XLPrecusorSearchMode.GetAllowedPrecursorMassIntervalsFromObservedMass(betaMass).First().AllowedInterval.Minimum;

                    double betaMassHigh = XLPrecusorSearchMode.GetAllowedPrecursorMassIntervalsFromObservedMass(betaMass).First().AllowedInterval.Maximum;

                    int betaMassLowIndex = BinarySearchGetIndex(massTable, betaMassLow);

                    while (betaMassLowIndex < massTable.Length && massTable[betaMassLowIndex] <= betaMassHigh)
                    {
                        if (scoringTable[betaMassLowIndex] >= byteScoreCutoff)
                        {
                            var key = new Tuple<int, int>(betaMassLowIndex, id);
                            if (!seenPair.Contains(key))
                            {
                                seenPair.Add(key);

                                List<int> possibleBetaCrosslinkSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(AllCrosslinkerSites, PeptideIndex[betaMassLowIndex], CommonParameters.DigestionParams.InitiatorMethionineBehavior, CleaveAtCrosslinkSite);

                                if (possibleBetaCrosslinkSites ==null)
                                {
                                    continue;
                                }

                                CrosslinkSpectralMatch csm = LocalizeCrosslinkSites(scan, id, betaMassLowIndex, Crosslinker, experimentFragmentMasses, intensityRanks);

                                if (csm!=null)
                                {
                                    csm.XlRank = rank;
                                }

                                possibleMatches.Add(csm);
                            }
                        }

                        betaMassLowIndex++;
                    }
                }

                rank++;
            }
        }
        

        /// <summary>
        /// Localizes the crosslink position on the alpha and beta peptides
        /// </summary>
        private CrosslinkSpectralMatch LocalizeCrosslinkSites(Ms2ScanWithSpecificMass theScan, int alphaIndex, int betaIndex, Crosslinker crosslinker, double[] experimentFragmentMasses, int[] intensityRanks)
        {
            CrosslinkSpectralMatch localizedCrosslinkedSpectralMatch = null;

            //The crosslink can crosslink same or different amino acid. Pairs are potential crosslink sites for alpha or beta. 
            List<Tuple<List<int>, List<int>>> pairs = new List<Tuple<List<int>, List<int>>>();

            if (crosslinker.CrosslinkerModSites.Equals(crosslinker.CrosslinkerModSites2))
            {
                List<int> possibleAlphaXlSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), PeptideIndex[alphaIndex], CommonParameters.DigestionParams.InitiatorMethionineBehavior, CleaveAtCrosslinkSite);
                List<int> possibleBetaXlSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), PeptideIndex[betaIndex], CommonParameters.DigestionParams.InitiatorMethionineBehavior, CleaveAtCrosslinkSite);

                pairs.Add(new Tuple<List<int>, List<int>>(possibleAlphaXlSites, possibleBetaXlSites));
            }
            else
            {
                List<int> possibleAlphaXlSites =CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), PeptideIndex[alphaIndex], CommonParameters.DigestionParams.InitiatorMethionineBehavior, CleaveAtCrosslinkSite);
                List<int> possibleBetaXlSites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites2.ToCharArray(), PeptideIndex[betaIndex], CommonParameters.DigestionParams.InitiatorMethionineBehavior, CleaveAtCrosslinkSite);

                pairs.Add(new Tuple<List<int>, List<int>>(possibleAlphaXlSites, possibleBetaXlSites));

                List<int> possibleAlphaXlSites2 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites2.ToCharArray(), PeptideIndex[alphaIndex], CommonParameters.DigestionParams.InitiatorMethionineBehavior, CleaveAtCrosslinkSite);
                List<int> possibleBetaXlSites2 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), PeptideIndex[betaIndex], CommonParameters.DigestionParams.InitiatorMethionineBehavior, CleaveAtCrosslinkSite);

                pairs.Add(new Tuple<List<int>, List<int>>(possibleAlphaXlSites2, possibleBetaXlSites2));
            }

            foreach (var pair in pairs)
            {
                if (pair.Item1!= null && pair.Item2!= null)
                {
                    int bestAlphaSite = 0;
                    int bestBetaSite = 0;
                    List<MatchedFragmentIon> bestMatchedAlphaIons = new List<MatchedFragmentIon>();
                    List<MatchedFragmentIon> bestMatchedBetaIons = new List<MatchedFragmentIon>();
                    Dictionary<int, List<MatchedFragmentIon>> bestMatchedChildAlphaIons = null;
                    Dictionary<int, List<MatchedFragmentIon>> bestMatchedChildBetaIons = null;
                    double bestAlphaLocalizedScore = 0;
                    double bestBetaLocalizedScore = 0;
                    double bestMS3AlphaScore = 0;
                    double bestMS3BetaScore = 0;

                    var fragmentsForEachAlphaLocalizedPossibility = CrosslinkedPeptide.XlGetTheoreticalFragments(CommonParameters.DissociationType,
                        Crosslinker, pair.Item1, PeptideIndex[betaIndex].MonoisotopicMass, PeptideIndex[alphaIndex]).ToList();

                    foreach (int possibleSite in pair.Item1)
                    {
                        foreach (var setOfFragments in fragmentsForEachAlphaLocalizedPossibility.Where(v => v.Item1 == possibleSite))
                        {
                            Dictionary<int, List<MatchedFragmentIon>> matchedChildAlphaIons = null;
                            var matchedIons = MatchFragmentIons(theScan, setOfFragments.Item2, CommonParameters);
                            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);
                            double ms3score = 0;

                            // search child scans (MS2+MS3)
                            foreach (Ms2ScanWithSpecificMass childScan in theScan.ChildScans)
                            {
                                var matchedChildIons = ScoreChildScan(theScan, childScan, possibleSite, PeptideIndex[alphaIndex], PeptideIndex[betaIndex]);
                                
                                if (matchedChildIons == null)
                                {
                                    continue;
                                }

                                if (matchedChildAlphaIons == null)
                                {
                                    matchedChildAlphaIons = new Dictionary<int, List<MatchedFragmentIon>>();
                                }

                                matchedChildAlphaIons.Add(childScan.OneBasedScanNumber, matchedChildIons);

                                double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);

                                if (childScan.TheScan.MsnOrder == 2 && CommonParameters.MS2ChildScanDissociationType!= DissociationType.LowCID)
                                {
                                    //Note that for  MS2(HCD)-(MS2)ETD type of data, add the childScore to score will bias the  alpha score more than beta score. 
                                    //We add 1/10 of childScore here, but better scoring is needed.
                                    score += childScore/10;
                                }
                                else
                                {
                                    ms3score += childScore;
                                }
                            }

                            if (score > bestAlphaLocalizedScore)
                            {
                                bestAlphaLocalizedScore = score;
                                bestMS3AlphaScore = ms3score;
                                bestAlphaSite = possibleSite;
                                bestMatchedAlphaIons = matchedIons;
                                bestMatchedChildAlphaIons = matchedChildAlphaIons;
                            }
                        }
                    }

                    var fragmentsForEachBetaLocalizedPossibility = CrosslinkedPeptide.XlGetTheoreticalFragments(CommonParameters.DissociationType,
                        Crosslinker, pair.Item2, PeptideIndex[alphaIndex].MonoisotopicMass, PeptideIndex[betaIndex]).ToList();

                    foreach (int possibleSite in pair.Item2)
                    {
                        foreach (var setOfFragments in fragmentsForEachBetaLocalizedPossibility.Where(v => v.Item1 == possibleSite))
                        {
                            var matchedIons = MatchFragmentIons(theScan, setOfFragments.Item2, CommonParameters);
                            Dictionary<int, List<MatchedFragmentIon>> matchedChildBetaIons = null;

                            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);
                            double ms3score = 0;

                            // search child scans (MS2+MS3)
                            foreach (Ms2ScanWithSpecificMass childScan in theScan.ChildScans)
                            {
                                var matchedChildIons = ScoreChildScan(theScan, childScan, possibleSite, PeptideIndex[betaIndex], PeptideIndex[alphaIndex]);

                                if (matchedChildIons == null)
                                {
                                    continue;
                                }

                                if (matchedChildBetaIons == null)
                                {
                                    matchedChildBetaIons = new Dictionary<int, List<MatchedFragmentIon>>();
                                }

                                matchedChildBetaIons.Add(childScan.OneBasedScanNumber, matchedChildIons);

                                double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);

                                if (childScan.TheScan.MsnOrder == 2 && CommonParameters.MS2ChildScanDissociationType != DissociationType.LowCID)
                                {
                                    //Note that for  MS2(HCD)-(MS2)ETD type of data, add the childScore to score will bias the  alpha score more than beta score. 
                                    //We add 1/10 of childScore here, but better scoring is needed.
                                    score += childScore/10;
                                }
                                else
                                {
                                    ms3score += childScore;
                                }
                            }

                            if (score > bestBetaLocalizedScore)
                            {
                                bestBetaLocalizedScore = score;
                                bestMS3BetaScore = ms3score;
                                bestBetaSite = possibleSite;
                                bestMatchedBetaIons = matchedIons;
                                bestMatchedChildBetaIons = matchedChildBetaIons;
                            }
                        }
                    }

                    //Remove any matched beta ions that also matched to the alpha peptide. The higher score one is alpha peptide.
                    if (PeptideIndex[alphaIndex].FullSequence != PeptideIndex[betaIndex].FullSequence)
                    {
                        if (bestAlphaLocalizedScore < bestBetaLocalizedScore)
                        {
                            var betaMz = new HashSet<double>(bestMatchedBetaIons.Select(p => p.Mz));
                            bestMatchedAlphaIons.RemoveAll(p => betaMz.Contains(p.Mz));
                            if ((int)bestAlphaLocalizedScore > bestMatchedAlphaIons.Count())
                            {
                                bestAlphaLocalizedScore = CalculatePeptideScore(theScan.TheScan, bestMatchedAlphaIons);
                            }
                        }
                        else
                        {
                            var alphaMz = new HashSet<double>(bestMatchedAlphaIons.Select(p => p.Mz));
                            bestMatchedBetaIons.RemoveAll(p => alphaMz.Contains(p.Mz));
                            if ((int)bestBetaLocalizedScore > bestMatchedBetaIons.Count())
                            {
                                bestBetaLocalizedScore = CalculatePeptideScore(theScan.TheScan, bestMatchedBetaIons);
                            }
                        }
                    }

                    if (bestAlphaLocalizedScore < CommonParameters.ScoreCutoff ||
                        bestBetaLocalizedScore < CommonParameters.ScoreCutoff)
                    {
                        return null;
                    }

                    var localizedAlpha = new CrosslinkSpectralMatch(PeptideIndex[alphaIndex], 0, bestAlphaLocalizedScore, 0, theScan, CommonParameters, bestMatchedAlphaIons);
                    var localizedBeta = new CrosslinkSpectralMatch(PeptideIndex[betaIndex], 0, bestBetaLocalizedScore, 0, theScan, CommonParameters, bestMatchedBetaIons);

                    localizedAlpha.ChildMatchedFragmentIons = bestMatchedChildAlphaIons;
                    localizedBeta.ChildMatchedFragmentIons = bestMatchedChildBetaIons;

                    localizedAlpha.LinkPositions = new List<int> { bestAlphaSite };
                    localizedBeta.LinkPositions = new List<int> { bestBetaSite };

                    if (crosslinker.Cleavable)
                    {
                        //cleavable crosslink parent ion information: intensity ranks

                        var alphaM = bestMatchedAlphaIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M);

                        localizedAlpha.ParentIonMaxIntensityRanks = new List<int>();

                        foreach (var am in alphaM)
                        {
                            var ind = BinarySearchGetIndex(experimentFragmentMasses, am.NeutralTheoreticalProduct.NeutralMass);
                            if (ind == experimentFragmentMasses.Length)
                            {
                                ind--;
                            }
                            localizedAlpha.ParentIonMaxIntensityRanks.Add(intensityRanks[ind]);
                        }
          
                        var betaM = bestMatchedBetaIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M);

                        localizedBeta.ParentIonMaxIntensityRanks = new List<int>();

                        foreach (var bm in betaM)
                        {
                            var ind = BinarySearchGetIndex(experimentFragmentMasses, bm.NeutralTheoreticalProduct.NeutralMass);
                            if (ind == experimentFragmentMasses.Length)
                            {
                                ind--;
                            }
                            localizedBeta.ParentIonMaxIntensityRanks.Add(intensityRanks[ind]);
                        }
                    }

                    localizedAlpha.MS3ChildScore = bestMS3AlphaScore;
                    localizedBeta.MS3ChildScore = bestMS3BetaScore;

                    //Decide which is alpha and which is beta.
                    if (bestAlphaLocalizedScore < bestBetaLocalizedScore)
                    {
                        var x = localizedAlpha;
                        localizedAlpha = localizedBeta;
                        localizedBeta = x;
                    }

                    localizedAlpha.BetaPeptide = localizedBeta;

                    //I think this is the only place where XLTotalScore is set
                    localizedAlpha.XLTotalScore = localizedAlpha.Score + localizedBeta.Score;        

                    localizedAlpha.CrossType = PsmCrossType.Cross;

                    localizedCrosslinkedSpectralMatch = localizedAlpha;
                }
            }

            return localizedCrosslinkedSpectralMatch;
        }

        private List<MatchedFragmentIon> ScoreChildScan(Ms2ScanWithSpecificMass parentScan, Ms2ScanWithSpecificMass childScan, int possibleSite, PeptideWithSetModifications mainPeptide, PeptideWithSetModifications otherPeptide)
        {
            bool shortMassAlphaMs3 = XLProductSearchMode.Accepts(childScan.PrecursorMass, mainPeptide.MonoisotopicMass + Crosslinker.CleaveMassShort) >= 0;
            bool longMassAlphaMs3 = XLProductSearchMode.Accepts(childScan.PrecursorMass, mainPeptide.MonoisotopicMass + Crosslinker.CleaveMassLong) >= 0;

            List<Product> childProducts = new List<Product>();

            //There are two situations now. 1) The childScan is MS3 scan and the crosslinker is cleavable. So the precursor mass of the MS3 scan must be same as signature ions.
            //2) The childScan is MS2 or MS3, but the precursor of the ChildScan is same. It is weird that the MS3 has same precursor mass as its parent scan, but it happens in some data.
            if (Crosslinker.Cleavable && childScan.TheScan.MsnOrder == 3 && (shortMassAlphaMs3 || longMassAlphaMs3))
            {
                double massToLocalize = shortMassAlphaMs3 ? Crosslinker.CleaveMassShort : Crosslinker.CleaveMassLong;
                if (mainPeptide.AllModsOneIsNterminus.TryGetValue(possibleSite + 1, out var existingMod))
                {
                    massToLocalize += existingMod.MonoisotopicMass.Value;
                }

                Dictionary<int, Modification> mod = new Dictionary<int, Modification> { { possibleSite + 1, new Modification(_monoisotopicMass: massToLocalize) } };

                foreach (var otherExistingMod in mainPeptide.AllModsOneIsNterminus.Where(p => p.Key != possibleSite + 1))
                {
                    mod.Add(otherExistingMod.Key, otherExistingMod.Value);
                }

                var peptideWithMod = new PeptideWithSetModifications(mainPeptide.Protein, mainPeptide.DigestionParams,
                    mainPeptide.OneBasedStartResidueInProtein, mainPeptide.OneBasedEndResidueInProtein,
                    mainPeptide.CleavageSpecificityForFdrCategory, mainPeptide.PeptideDescription,
                    mainPeptide.MissedCleavages, mod, mainPeptide.NumFixedMods);
          
                    peptideWithMod.Fragment(CommonParameters.MS3ChildScanDissociationType, FragmentationTerminus.Both, childProducts);              
            }
            else if (Math.Abs(childScan.PrecursorMass - parentScan.PrecursorMass) < 0.001)
            {
                if (childScan.TheScan.MsnOrder == 2)
                {
                    childProducts = CrosslinkedPeptide.XlGetTheoreticalFragments(CommonParameters.MS2ChildScanDissociationType,
                        Crosslinker, new List<int> { possibleSite }, otherPeptide.MonoisotopicMass, mainPeptide).First().Item2;
                }
                else
                {
                    //It tried to corver the situation that the MS3 scan have the same precursor mass as the parent MS2.  
                    //This is so rare that I only see it in one data. May need unit test in the future.
                    childProducts = CrosslinkedPeptide.XlGetTheoreticalFragments(CommonParameters.MS3ChildScanDissociationType,
                        Crosslinker, new List<int> { possibleSite }, otherPeptide.MonoisotopicMass, mainPeptide).First().Item2;
                }
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
            List<int> possiblePositions, Modification deadEndMod, int notch, int scanIndex)
        {
            double bestScore = 0;
            List<MatchedFragmentIon> bestMatchingFragments = new List<MatchedFragmentIon>();
            PeptideWithSetModifications bestLocalizedPeptide = null;
            int bestPosition = 0;
            List<Product> products = new List<Product>();

            foreach (int location in possiblePositions)
            {
                Dictionary<int, Modification> mods = originalPeptide.AllModsOneIsNterminus.ToDictionary(p => p.Key, p => p.Value);

                mods.Add(location + 1, deadEndMod);
 
                var localizedPeptide = new PeptideWithSetModifications(originalPeptide.Protein, originalPeptide.DigestionParams, originalPeptide.OneBasedStartResidueInProtein,
                    originalPeptide.OneBasedEndResidueInProtein, originalPeptide.CleavageSpecificityForFdrCategory, originalPeptide.PeptideDescription, originalPeptide.MissedCleavages, mods, originalPeptide.NumFixedMods);

                localizedPeptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both, products);
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

            var csm = new CrosslinkSpectralMatch(bestLocalizedPeptide, notch, bestScore, scanIndex, theScan, commonParameters, bestMatchingFragments);

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

            return csm;
        }

        /// <summary>
        /// Localizes the loop to a begin and end residue
        /// </summary>
        private CrosslinkSpectralMatch LocalizeLoopSites(PeptideWithSetModifications originalPeptide, Ms2ScanWithSpecificMass theScan, CommonParameters commonParameters,
            List<int> possiblePositions, int notch, int scanIndex)
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

            var csm = new CrosslinkSpectralMatch(originalPeptide, notch, bestScore, scanIndex, theScan, commonParameters, bestMatchingFragments)
            {
                CrossType = PsmCrossType.Loop,
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

        public static int BinarySearchGetIndex(double[] massArray, double targetMass)
        {

            // BinarySearch Returns:
            //     The index of the specified value in the specified array, if value is found; otherwise,
            //     a negative number. If value is not found and value is less than one or more elements
            //     in array, the negative number returned is the bitwise complement of the index
            //     of the first element that is larger than value. If value is not found and value
            //     is greater than all elements in array, the negative number returned is the bitwise
            //     complement of (the index of the last element plus 1). If this method is called
            //     with a non-sorted array, the return value can be incorrect and a negative number
            //     could be returned, even if value is present in array.
            var iD = Array.BinarySearch(massArray, targetMass);
            if (iD < 0) { iD = ~iD; }
            else
            {
                while (iD - 1 >= 0 && massArray[iD - 1] - targetMass >= -ToleranceForMassDifferentiation)
                {
                    iD--;
                }
            }
            return iD;
        }
    }
}