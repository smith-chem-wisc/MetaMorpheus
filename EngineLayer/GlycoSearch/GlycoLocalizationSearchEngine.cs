using EngineLayer.ModernSearch;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using EngineLayer;
using MassSpectrometry;

namespace EngineLayer.GlycoSearch
{
    public class GlycoLocalizationSearchEngine : ModernSearchEngine
    {
        protected readonly List<List<GlycoSpectralMatch>>[] GlobalCsms;

        private bool IsOGlycoSearch;
        // crosslinker molecule
        private readonly bool GlycoSearchTopN;
        private readonly int TopN;
        private readonly int _maxOGlycanNum;

        private readonly bool SearchGlycan182;
        private readonly Tolerance PrecusorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;

        private readonly List<int>[] SecondFragmentIndex;
        private readonly List<double>[] Precursors;

        public GlycoLocalizationSearchEngine(List<List<GlycoSpectralMatch>>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters,
             bool isOGlycoSearch, bool glycoSearchTop, int glycoSearchTopNum, bool searchGlycan182, int maxOGlycanNum, List<double>[] _precursors, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.IsOGlycoSearch = isOGlycoSearch;
            this.GlycoSearchTopN = glycoSearchTop;
            this.TopN = glycoSearchTopNum;
            this._maxOGlycanNum = maxOGlycanNum;

            SecondFragmentIndex = secondFragmentIndex;
            Precursors = _precursors;

            PrecusorSearchMode = commonParameters.PrecursorMassTolerance;
            ProductSearchMode = new SingleAbsoluteAroundZeroSearchMode(20); //For Oxinium ion only

            SearchGlycan182 = searchGlycan182;
            if (!isOGlycoSearch)
            {
                var NGlycans = Glycan.LoadGlycan(GlobalVariables.NGlycanLocation);

                if (SearchGlycan182)
                {
                    var NGlycans182 = Glycan.LoadKindGlycan(GlobalVariables.NGlycanLocation_182, NGlycans);
                    Glycans = NGlycans182.OrderBy(p => p.Mass).ToArray();
                }
                else
                {
                    Glycans = NGlycans.OrderBy(p => p.Mass).ToArray();
                }

                DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
            }
            else
            {
                GlycanBox.GlobalOGlycans = Glycan.LoadGlycan(GlobalVariables.OGlycanLocation).ToArray();
                GlycanBox.GlobalOGlycanModifications = GlycanBox.BuildGlobalOGlycanModifications(GlycanBox.GlobalOGlycans);
                OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(_maxOGlycanNum).OrderBy(p => p.Mass).ToArray();

                //ModBox.SelectedModifications = new Modification[5];
                //ModBox.SelectedModifications[0] = GlobalVariables.AllModsKnownDictionary["Oxidation on M"];
                //ModBox.SelectedModifications[1] = GlobalVariables.AllModsKnownDictionary["Hydroxylation on P"];
                //ModBox.SelectedModifications[2] = GlobalVariables.AllModsKnownDictionary["Hydroxylation on K"];
                //ModBox.SelectedModifications[3] = GlobalVariables.AllModsKnownDictionary["Glucosylgalactosyl on K"];
                //ModBox.SelectedModifications[4] = GlobalVariables.AllModsKnownDictionary["Galactosyl on K"];

                //ModBoxes = ModBox.BuildModBoxes(_maxOGlycanNum).Where(p => !p.MotifNeeded.ContainsKey("K") || (p.MotifNeeded.ContainsKey("K") && p.MotifNeeded["K"].Count <= 3)).OrderBy(p => p.Mass).ToArray();
            }
        }

        private Glycan[] Glycans { get; }
        private Glycan[] DecoyGlycans { get; }
        private GlycanBox[] OGlycanBoxes { get; }
        private ModBox[] ModBoxes { get; }

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
                    idsOfPeptidesTopN.Clear();

                    var scan = ListOfSortedMs2Scans[scanIndex];
                    var precursors = Precursors[scanIndex];

                    // get fragment bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(scan, FragmentIndex, CommonParameters.DissociationType);
                    List<int> childBinsToSearch = null;

                    //TO DO: limit the high bound limitation

                    // first-pass scoring
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, Double.PositiveInfinity, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);

                    //child scan first-pass scoring
                    if (scan.ChildScans != null && CommonParameters.ChildScanDissociationType != DissociationType.LowCID)
                    {
                        Array.Clear(secondScoringTable, 0, secondScoringTable.Length);
                        childIdsOfPeptidesPossiblyObserved.Clear();

                        childBinsToSearch = new List<int>();

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
                        scoreAtTopN = 0;
                        peptideCount = 0;
                        foreach (int id in idsOfPeptidesPossiblyObserved.OrderByDescending(p => scoringTable[p]))
                        {
                            peptideCount++;
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

                        List<List<GlycoSpectralMatch>> gsms;
                        if (IsOGlycoSearch == false)
                        {
                            //gsms = FindNGlycopeptide(scan, idsOfPeptidesTopN, scanIndex);
                            gsms = new List<List<GlycoSpectralMatch>>();
                        }
                        else
                        {
                            //gsms = FindOGlycopeptide(scan, idsOfPeptidesTopN, scanIndex);
                            gsms = FindOGlycopeptideHash(scan, precursors, idsOfPeptidesTopN, scanIndex, allBinsToSearch, childBinsToSearch, (int)byteScoreCutoff);
                            //gsms = FindModPep(scan, idsOfPeptidesTopN, scanIndex, (int)byteScoreCutoff);
                        }


                        if (gsms.Count == 0)
                        {
                            progress++;
                            continue;
                        }

                        if (GlobalCsms[scanIndex] == null)
                        {
                            GlobalCsms[scanIndex] = new List<List<GlycoSpectralMatch>>();
                        }

                        GlobalCsms[scanIndex].AddRange(gsms);
                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing glyco search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
                    }
                }
            });

            return new MetaMorpheusEngineResults(this);
        }

        private List<GlycoSpectralMatch> FindNGlycopeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            if (theScan.OxiniumIonNum < 2)
            {
                return possibleMatches;
            }

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                //Considering coisolation, it doesn't mean it must from a glycopeptide even the scan contains oxinium ions.
                //if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass) >= 0)
                //{
                //    List<Product> products = theScanBestPeptide[ind].BestPeptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                //    var matchedFragmentIons = MatchFragmentIons(theScan, products, commonParameters);
                //    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                //    var psmCrossSingle = new CrosslinkSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, score, scanIndex, theScan, commonParameters.DigestionParams, matchedFragmentIons);
                //    psmCrossSingle.CrossType = PsmCrossType.Single;
                //    psmCrossSingle.XlRank = new List<int> { ind };

                //    possibleMatches.Add(psmCrossSingle);
                //}

                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];

                if (theScan.OxiniumIonNum < 2)
                {
                    continue;
                }

                List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxt", "Nxs" });
                if (modPos.Count < 1)
                {
                    continue;
                }

                var possibleGlycanMassLow = theScan.PrecursorMass * (1 - 1E-5) - theScanBestPeptide.MonoisotopicMass;
                if (possibleGlycanMassLow < 200 || possibleGlycanMassLow > Glycans.Last().Mass)
                {
                    continue;
                }


                int iDLow = GlycoPeptides.BinarySearchGetIndex(Glycans.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);

                while (iDLow < Glycans.Count() && PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + (double)Glycans[iDLow].Mass / 1E5))
                {
                    double bestLocalizedScore = 0;
                    int bestSite = 0;
                    List<MatchedFragmentIon> bestMatchedIons = new List<MatchedFragmentIon>();
                    PeptideWithSetModifications peptideWithSetModifications = theScanBestPeptide;
                    foreach (int possibleSite in modPos)
                    {
                        var testPeptide = GlycoPeptides.GenerateGlycopeptide(possibleSite, theScanBestPeptide, Glycans[iDLow]);

                        List<Product> theoreticalProducts = testPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both).Where(p => p.ProductType != ProductType.M).ToList();
                        theoreticalProducts.AddRange(GlycoPeptides.GetGlycanYIons(theScan.PrecursorMass, Glycans[iDLow]));

                        var matchedIons = MatchOriginFragmentIons(theScan, theoreticalProducts, CommonParameters);

                        if (!GlycoPeptides.ScanTrimannosylCoreFilter(matchedIons, Glycans[iDLow]))
                        {
                            continue;
                        }

                        double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

                        if (score > bestLocalizedScore)
                        {
                            peptideWithSetModifications = testPeptide;
                            bestLocalizedScore = score;
                            bestSite = possibleSite;
                            bestMatchedIons = matchedIons;
                        }

                    }

                    var psmCross = new GlycoSpectralMatch(peptideWithSetModifications, 0, bestLocalizedScore, scanIndex, theScan, CommonParameters.DigestionParams, bestMatchedIons);
                    psmCross.Glycan = new List<Glycan> { Glycans[iDLow] };
                    psmCross.GlycanScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.Annotation.Contains('M')).ToList());
                    psmCross.DiagnosticIonScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.Annotation.Contains('D')).ToList());
                    psmCross.PeptideScore = psmCross.TotalScore - psmCross.GlycanScore - psmCross.DiagnosticIonScore;
                    psmCross.Rank = ind;
                    psmCross.localizations = new List<int[]> { new int[] { bestSite - 1 } }; //TO DO: ambiguity modification site
                    possibleMatches.Add(psmCross);

                    iDLow++;
                }
            }

            //if (possibleMatches.Count != 0)
            //{
            //    possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
            //    bestPsmCross = possibleMatches.First();
            //    bestPsmCross.ResolveAllAmbiguities();
            //    bestPsmCross.DeltaScore = bestPsmCross.XLTotalScore;
            //    if (possibleMatches.Count > 1)
            //    {
            //        //This DeltaScore will be 0 if there are more than one glycan matched.
            //        bestPsmCross.DeltaScore = Math.Abs(possibleMatches.First().Score - possibleMatches[1].Score);

            //        ////TO DO: Try to find other plausible glycans
            //        //for (int iPsm = 1; iPsm < possibleMatches.Count; iPsm++)
            //        //{
            //        //    possibleMatches[iPsm].ResolveAllAmbiguities();
            //        //    if (possibleMatches[iPsm].Score == bestPsmCross.Score && possibleMatches[iPsm].Glycan != null)
            //        //    {
            //        //        bestPsmCross.Glycan.Add(possibleMatches[iPsm].Glycan.First());
            //        //    }
            //        //}
            //    }
            //}

            if (possibleMatches.Count != 0)
            {
                possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
            }
            return possibleMatches;
        }

        private List<GlycoSpectralMatch> FindOGlycopeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];

                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass))
                {
                    List<Product> products = theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                    var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
                    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                    var psmCrossSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, score, scanIndex, theScan, CommonParameters.DigestionParams, matchedFragmentIons);
                    psmCrossSingle.Rank = ind;
                    psmCrossSingle.ResolveAllAmbiguities();

                    possibleMatches.Add(psmCrossSingle);
                }
                //TO DO: add if the scan contains diagnostic ions, the rule to check diagnostic ions
                else if ((theScan.PrecursorMass - theScanBestPeptide.MonoisotopicMass >= 100) && GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType) >= 1)
                {
                    //Using glycanBoxes
                    var possibleGlycanMassLow = theScan.PrecursorMass * (1 - 1E-5) - theScanBestPeptide.MonoisotopicMass;
                    if (possibleGlycanMassLow < 100 || possibleGlycanMassLow > OGlycanBoxes.Last().Mass)
                    {
                        continue;
                    }

                    int iDLow = GlycoPeptides.BinarySearchGetIndex(OGlycanBoxes.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);

                    while (iDLow < OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + (double)OGlycanBoxes[iDLow].Mass / 1E5)))
                    {
                        List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" });
                        if (modPos.Count >= OGlycanBoxes[iDLow].NumberOfGlycans)
                        {
                            var permutateModPositions = GlycoPeptides.GetPermutations(modPos, OGlycanBoxes[iDLow].GlycanIds);

                            double bestLocalizedScore = 0;
                            PeptideWithSetModifications[] bestPeptideWithMod = new PeptideWithSetModifications[1];
                            var bestMatchedIons = new List<MatchedFragmentIon>();
                            var bestChildMatchedIons = new Dictionary<int, List<MatchedFragmentIon>>();
                            List<int[]> localization = new List<int[]>();

                            foreach (var theModPositions in permutateModPositions)
                            {
                                var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(theModPositions.ToArray(), theScanBestPeptide, OGlycanBoxes[iDLow]);

                                var fragmentsForEachGlycoPeptide = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.DissociationType, theScanBestPeptide, peptideWithMod);

                                var matchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);

                                //TO DO: May need a new score function
                                double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

                                var allMatchedChildIons = new Dictionary<int, List<MatchedFragmentIon>>();

                                foreach (var childScan in theScan.ChildScans)
                                {
                                    var childFragments = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.ChildScanDissociationType, theScanBestPeptide, peptideWithMod);

                                    var matchedChildIons = MatchFragmentIons(childScan, childFragments, CommonParameters);

                                    if (matchedChildIons == null)
                                    {
                                        continue;
                                    }

                                    allMatchedChildIons.Add(childScan.OneBasedScanNumber, matchedChildIons);
                                    double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);

                                    //TO DO:may think a different way to use childScore
                                    score += childScore;
                                }

                                if (bestLocalizedScore < score)
                                {
                                    bestPeptideWithMod[0] = peptideWithMod;
                                    bestLocalizedScore = score;
                                    bestMatchedIons = matchedIons;
                                    bestChildMatchedIons = allMatchedChildIons;
                                    localization.Clear();
                                    localization.Add(theModPositions.ToArray());
                                }
                                else if (bestLocalizedScore == score)
                                {
                                    localization.Add(theModPositions.ToArray());
                                }

                            }

                            //TO CHECK: In theory, the bestPeptideWithMod[0] shouldn't be null, because the score always >= index-score >0. 
                            if (bestPeptideWithMod[0] != null)
                            {
                                var psmGlyco = new GlycoSpectralMatch(bestPeptideWithMod[0], 0, bestLocalizedScore, scanIndex, theScan, CommonParameters.DigestionParams, bestMatchedIons);
                                psmGlyco.glycanBoxes = new List<GlycanBox> { OGlycanBoxes[iDLow] };
                                psmGlyco.Rank = ind;
                                psmGlyco.ChildMatchedFragmentIons = bestChildMatchedIons;
                                psmGlyco.localizations = localization;
                                possibleMatches.Add(psmGlyco);
                            }
                        }

                        iDLow++;
                    }
                }
            }

            if (possibleMatches.Count != 0)
            {
                possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
            }

            return possibleMatches;
        }

        private List<List<GlycoSpectralMatch>> FindOGlycopeptideHash(Ms2ScanWithSpecificMass theScan, List<double> precursors, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, List<int> allBinsToSearch, List<int> childBinsToSearch, int scoreCutOff)
        {
            List<List<GlycoSpectralMatch>> glycoSpectralMatches = new List<List<GlycoSpectralMatch>>();

            foreach (var precursor in precursors)
            {
                List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

                for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
                {
                    var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];

                    if (PrecusorSearchMode.Within(precursor, theScanBestPeptide.MonoisotopicMass))
                    {
                        List<Product> products = theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                        var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
                        double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                        var psmCrossSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, score, scanIndex, theScan, CommonParameters.DigestionParams, matchedFragmentIons);
                        psmCrossSingle.Rank = ind;
                        psmCrossSingle.ResolveAllAmbiguities();

                        possibleMatches.Add(psmCrossSingle);
                    }
                    //TO DO: add if the scan contains diagnostic ions, the rule to check diagnostic ions
                    else if ((precursor - theScanBestPeptide.MonoisotopicMass >= 100) && GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType) >= 1)
                    {
                        //Using glycanBoxes
                        var possibleGlycanMassLow = precursor * (1 - 1E-5) - theScanBestPeptide.MonoisotopicMass;
                        if (possibleGlycanMassLow < 100 || possibleGlycanMassLow > OGlycanBoxes.Last().Mass)
                        {
                            continue;
                        }

                        List<Tuple<int, int[]>> glycanBoxId_localization = new List<Tuple<int, int[]>>();

                        int iDLow = GlycoPeptides.BinarySearchGetIndex(OGlycanBoxes.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);

                        while (iDLow < OGlycanBoxes.Count() && (PrecusorSearchMode.Within(precursor, theScanBestPeptide.MonoisotopicMass + (double)OGlycanBoxes[iDLow].Mass / 1E5)))
                        {
                            List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" });
                            if (modPos.Count >= OGlycanBoxes[iDLow].NumberOfGlycans)
                            {
                                var permutateModPositions = GlycoPeptides.GetPermutations(modPos, OGlycanBoxes[iDLow].GlycanIds);

                                foreach (var theModPositions in permutateModPositions)
                                {
                                    glycanBoxId_localization.Add(new Tuple<int, int[]>(iDLow, theModPositions));
                                }
                            }

                            iDLow++;
                        }

                        //TO DO: Consider the situation that there is no child Scan. Consider if the father scan EThcD Scan.
                        if (glycanBoxId_localization.Count > 0 && theScan.ChildScans.Count > 0)
                        {
                            HashSet<int> allPeaksForLocalization = new HashSet<int>(childBinsToSearch);

                            List<Product> products = theScanBestPeptide.Fragment(CommonParameters.ChildScanDissociationType, FragmentationTerminus.Both).ToList();

                            List<Tuple<int, int[]>> localizationCandidates = new List<Tuple<int, int[]>>();
                            int BestLocalizaionScore = 0;

                            for (int i = 0; i < glycanBoxId_localization.Count; i++)
                            {
                                var fragmentHash = GlycoPeptides.GetFragmentHash(products, glycanBoxId_localization[i], OGlycanBoxes, FragmentBinsPerDalton);
                                int currentLocalizationScore = allPeaksForLocalization.Intersect(fragmentHash).Count();
                                if (currentLocalizationScore > BestLocalizaionScore)
                                {
                                    localizationCandidates.Clear();
                                    localizationCandidates.Add(glycanBoxId_localization[i]);
                                }
                                else if (BestLocalizaionScore > 0 && currentLocalizationScore == BestLocalizaionScore)
                                {
                                    localizationCandidates.Add(glycanBoxId_localization[i]);
                                }
                            }

                            if (localizationCandidates.Count > 0)
                            {
                                var psmGlyco = CreateGsm(theScan, scanIndex, ind, theScanBestPeptide, localizationCandidates, CommonParameters);

                                possibleMatches.Add(psmGlyco);
                            }
                        }
                    }

                    if (possibleMatches.Count != 0)
                    {
                        possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
                    }
                }

                if (possibleMatches.Count != 0)
                {
                    glycoSpectralMatches.Add(possibleMatches);
                }
            }

            return glycoSpectralMatches;
        }

        private GlycoSpectralMatch CreateGsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, List<Tuple<int, int[]>> glycanBox_localization, CommonParameters commonParameters)
        {
            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(glycanBox_localization.First().Item2, peptide, OGlycanBoxes[glycanBox_localization.First().Item1]);

            var fragmentsForEachGlycoPeptide = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.DissociationType, peptide, peptideWithMod);

            var matchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);

            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

            var allMatchedChildIons = new Dictionary<int, List<MatchedFragmentIon>>();

            foreach (var childScan in theScan.ChildScans)
            {
                var childFragments = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.ChildScanDissociationType, peptide, peptideWithMod);

                var matchedChildIons = MatchFragmentIons(childScan, childFragments, CommonParameters);

                if (matchedChildIons == null)
                {
                    continue;
                }

                allMatchedChildIons.Add(childScan.OneBasedScanNumber, matchedChildIons);
                double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);

                //TO DO:may think a different way to use childScore
                score += childScore;
            }

            var psmGlyco = new GlycoSpectralMatch(peptideWithMod, 0, score, scanIndex, theScan, CommonParameters.DigestionParams, matchedIons);

            psmGlyco.Rank = rank;
            //TO DO
            psmGlyco.ChildMatchedFragmentIons = allMatchedChildIons;
            psmGlyco.glycanBoxes = new List<GlycanBox> { OGlycanBoxes[glycanBox_localization.First().Item1] };
            //psmGlyco.localizations = glycanBox_localization.Values.ToList();

            return psmGlyco;
        }

        private List<GlycoSpectralMatch> FindModPep(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];

                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass))
                {
                    List<Product> products = theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                    var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
                    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                    var psmCrossSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, score, scanIndex, theScan, CommonParameters.DigestionParams, matchedFragmentIons);
                    psmCrossSingle.Rank = ind;
                    psmCrossSingle.ResolveAllAmbiguities();

                    possibleMatches.Add(psmCrossSingle);
                }
                else if (theScan.PrecursorMass - theScanBestPeptide.MonoisotopicMass >= 10)
                {
                    //Using glycanBoxes
                    var possibleGlycanMassLow = theScan.PrecursorMass * (1 - 1E-5) - theScanBestPeptide.MonoisotopicMass;
                    if (possibleGlycanMassLow < ModBoxes.First().Mass || possibleGlycanMassLow > ModBoxes.Last().Mass)
                    {
                        continue;
                    }

                    Dictionary<ModBox, List<Tuple<int, int>[]>> glycanBox_localization = new Dictionary<ModBox, List<Tuple<int, int>[]>>();
                    double bestLocalizedScore = 0;

                    int iDLow = GlycoPeptides.BinarySearchGetIndex(ModBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow);

                    while (iDLow < ModBoxes.Count() && PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + ModBoxes[iDLow].Mass))
                    {
                        var permutateModPositions = ModBox.GetPossibleModSites(theScanBestPeptide, ModBoxes[iDLow]);

                        foreach (var theModPositions in permutateModPositions)
                        {
                            var peptideWithMod = ModBox.GetTheoreticalPeptide(theModPositions.ToArray(), theScanBestPeptide, ModBoxes[iDLow]);

                            var fragmentsForEachGlycoPeptide = peptideWithMod.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both).ToList();

                            var matchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);

                            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

                            if (score >= scoreCutOff && bestLocalizedScore < score)
                            {
                                glycanBox_localization.Clear();

                                bestLocalizedScore = score;
                                glycanBox_localization.Add(ModBoxes[iDLow], new List<Tuple<int, int>[]> { theModPositions });

                            }
                            else if (score >= scoreCutOff && bestLocalizedScore == score)
                            {
                                if (glycanBox_localization.ContainsKey(ModBoxes[iDLow]))
                                {
                                    glycanBox_localization[ModBoxes[iDLow]].Add(theModPositions);
                                }
                                else
                                {
                                    glycanBox_localization.Add(ModBoxes[iDLow], new List<Tuple<int, int>[]> { theModPositions });
                                }
                            }

                        }

                        iDLow++;
                    }

                    //In theory, the peptide_localization shouldn't be null, but it is possible that the real score is smaller than indexed score.
                    if (glycanBox_localization.Count > 0)
                    {
                        var psmGlyco = CreateGsm(theScan, scanIndex, ind, theScanBestPeptide, glycanBox_localization, CommonParameters);

                        possibleMatches.Add(psmGlyco);
                    }
                }
            }

            if (possibleMatches.Count != 0)
            {
                possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
            }

            return possibleMatches;
        }

        private GlycoSpectralMatch CreateGsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, Dictionary<ModBox, List<Tuple<int, int>[]>> glycanBox_localization, CommonParameters commonParameters)
        {
            var peptideWithMod = ModBox.GetTheoreticalPeptide(glycanBox_localization.First().Value.First(), peptide, glycanBox_localization.First().Key);

            var fragmentsForEachGlycoPeptide = peptideWithMod.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both).ToList();

            var matchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);

            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

            var psmGlyco = new GlycoSpectralMatch(peptideWithMod, 0, score, scanIndex, theScan, CommonParameters.DigestionParams, matchedIons);

            psmGlyco.Rank = rank;
            //TO DO
            //psmGlyco.localizations = glycanBox_localization.First().Value;

            return psmGlyco;
        }

    }
}