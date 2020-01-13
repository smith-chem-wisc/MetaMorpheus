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
    public class GlycoSearchEngine : ModernSearchEngine
    {
        protected readonly List<GlycoSpectralMatch>[] GlobalCsms;

        private bool IsOGlycoSearch;
        private readonly int TopN;
        private readonly int _maxOGlycanNum;
        private readonly int _glycanDatabaseIndex;

        private readonly Tolerance PrecusorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;

        private readonly List<int>[] SecondFragmentIndex;

        public GlycoSearchEngine(List<GlycoSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters,
             int glycanDatabaseIndex, bool isOGlycoSearch, int glycoSearchTopNum, int maxOGlycanNum, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.IsOGlycoSearch = isOGlycoSearch;
            this.TopN = glycoSearchTopNum;
            this._maxOGlycanNum = maxOGlycanNum;
            this._glycanDatabaseIndex = glycanDatabaseIndex;

            SecondFragmentIndex = secondFragmentIndex;
            PrecusorSearchMode = commonParameters.PrecursorMassTolerance;
            ProductSearchMode = new SinglePpmAroundZeroSearchMode(20); //For Oxinium ion only


            if (!isOGlycoSearch)
            {
                Glycans = GlycanDatabase.LoadGlycan(GlobalVariables.GlycanLocations.ElementAt(_glycanDatabaseIndex), !isOGlycoSearch).OrderBy(p=>p.Mass).ToArray();          
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
            }
            else
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.GlycanLocations.ElementAt(_glycanDatabaseIndex), !isOGlycoSearch).ToArray();
                GlycanBox.GlobalOGlycanModifications = GlycanBox.BuildGlobalOGlycanModifications(GlycanBox.GlobalOGlycans);
                OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(_maxOGlycanNum).OrderBy(p => p.Mass).ToArray();

                //This is for Hydroxyproline 
                ModBox.SelectedModifications = new Modification[4];      
                ModBox.SelectedModifications[0] = GlobalVariables.AllModsKnownDictionary["Hydroxylation on P"];
                ModBox.SelectedModifications[1] = GlobalVariables.AllModsKnownDictionary["Hydroxylation on K"];
                ModBox.SelectedModifications[2] = GlobalVariables.AllModsKnownDictionary["Glucosylgalactosyl on K"];
                ModBox.SelectedModifications[3] = GlobalVariables.AllModsKnownDictionary["Galactosyl on K"];
                ModBoxes = ModBox.BuildModBoxes(_maxOGlycanNum).Where(p => !p.MotifNeeded.ContainsKey("K") || (p.MotifNeeded.ContainsKey("K") && p.MotifNeeded["K"].Count <= 3)).OrderBy(p => p.Mass).ToArray();
            }
        }

        private Glycan[] Glycans {get;}
        //private Glycan[] DecoyGlycans { get; }
        private GlycanBox[] OGlycanBoxes { get; }
        private ModBox[] ModBoxes { get;  }

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

                        List<GlycoSpectralMatch> gsms;
                        if (IsOGlycoSearch == false)
                        {
                            gsms = FindNGlycopeptide(scan, idsOfPeptidesTopN, scanIndex);
                        }
                        else
                        {
                            //gsms = FindOGlycopeptideHash(scan, idsOfPeptidesTopN, scanIndex, allBinsToSearch, childBinsToSearch, (int)byteScoreCutoff);
                            //gsms = FindOGlycopeptideHashLocal(scan, idsOfPeptidesTopN, scanIndex, allBinsToSearch, childBinsToSearch, (int)byteScoreCutoff);                            
                            gsms = FindModPepHash(scan, idsOfPeptidesTopN, scanIndex, allBinsToSearch, (int)byteScoreCutoff);
                        }


                        if (gsms.Count == 0)
                        {
                            progress++;
                            continue;
                        }

                        if (GlobalCsms[scanIndex] == null)
                        {
                            GlobalCsms[scanIndex] = new List<GlycoSpectralMatch>();
                        }

                        GlobalCsms[scanIndex].AddRange(gsms.Where(p => p != null).OrderByDescending(p => p.TotalScore));
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

        private List<GlycoSpectralMatch> FindOGlycopeptideHash(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, List<int> allBinsToSearch, List<int> childBinsToSearch, int scoreCutOff)
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
                else if ((theScan.PrecursorMass - theScanBestPeptide.MonoisotopicMass >= 100))
                {
                    //Using glycanBoxes
                    var possibleGlycanMassLow = theScan.PrecursorMass * (1 - 1E-5) - theScanBestPeptide.MonoisotopicMass;
                    if (possibleGlycanMassLow < 100 || possibleGlycanMassLow > OGlycanBoxes.Last().Mass)
                    {
                        continue;
                    }

                    List<Tuple<int, int[]>> glycanBoxId_localization = new List<Tuple<int, int[]>>();

                    int iDLow = GlycoPeptides.BinarySearchGetIndex(OGlycanBoxes.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);

                    while (iDLow < OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + (double)OGlycanBoxes[iDLow].Mass / 1E5)))
                    {
                        List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" });
                        if (modPos.Count >= OGlycanBoxes[iDLow].NumberOfMods)
                        {
                            var permutateModPositions = GlycoPeptides.GetPermutations(modPos, OGlycanBoxes[iDLow].ModIds);

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

            return possibleMatches;
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

        private List<GlycoSpectralMatch> FindOGlycopeptideHashLocal(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, List<int> allBinsToSearch, List<int> childBinsToSearch, int scoreCutOff)
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
                else if ((theScan.PrecursorMass - theScanBestPeptide.MonoisotopicMass >= 100))
                {
                    //Filter by glycanBoxes mass difference.
                    var possibleGlycanMassLow = theScan.PrecursorMass * (1 - 1E-5) - theScanBestPeptide.MonoisotopicMass;

                    if (possibleGlycanMassLow < OGlycanBoxes.First().Mass || possibleGlycanMassLow > OGlycanBoxes.Last().Mass)
                    {
                        continue;
                    }

                    //Filter by OxoniumIon
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);

                    //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    if (oxoniumIonIntensities[9] > 0)
                    {
                        continue;
                    }

                    HashSet<int> allPeaksForLocalization = new HashSet<int>(childBinsToSearch);

                    List<Product> products = theScanBestPeptide.Fragment(CommonParameters.ChildScanDissociationType, FragmentationTerminus.Both).ToList();

                    int[] modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" }).ToArray();

                    double bestLocalizedScore = 1;

                    List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();
                    List<int> ids = new List<int>();

                    int iDLow = GlycoPeptides.BinarySearchGetIndex(OGlycanBoxes.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);

                    while (iDLow < OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + (double)OGlycanBoxes[iDLow].Mass / 1E5)))
                    {
                        if (modPos.Length >= OGlycanBoxes[iDLow].NumberOfMods && GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, OGlycanBoxes[iDLow]))
                        {
                            var boxes = GlycanBox.BuildChildOGlycanBoxes(OGlycanBoxes[iDLow].NumberOfMods, OGlycanBoxes[iDLow].ModIds).ToArray();

                            LocalizationGraph localizationGraph = new LocalizationGraph(modPos.Length, boxes.Length);
                            localizationGraph.Localization(modPos, OGlycanBoxes[iDLow], boxes, allPeaksForLocalization, products, theScanBestPeptide.Length);

                            double currentLocalizationScore = localizationGraph.array.Last().Last().maxCost;
                            if (currentLocalizationScore > bestLocalizedScore)
                            {
                                ids.Clear();
                                localizationGraphs.Clear();
                                ids.Add(iDLow);
                                localizationGraphs.Add(localizationGraph);
                            }
                            else if (bestLocalizedScore > 0 && currentLocalizationScore == bestLocalizedScore)
                            {
                                ids.Add(iDLow);
                                localizationGraphs.Add(localizationGraph);
                            }
                        }

                        iDLow++;
                    }

                    //In theory, the peptide_localization shouldn't be null, but it is possible that the real score is smaller than indexed score.
                    if (localizationGraphs.Count > 0)
                    {
                        List<Tuple<int, Tuple<int, int>[]>> localizationCandidates = new List<Tuple<int, Tuple<int, int>[]>>();
                        for (int i = 0; i < localizationGraphs.Count; i++)
                        {
                            var boxes = GlycanBox.BuildChildOGlycanBoxes(OGlycanBoxes[ids[i]].NumberOfMods, OGlycanBoxes[ids[i]].ModIds).ToArray();
                            var allPaths = LocalizationGraph.GetAllPaths(localizationGraphs[i].array, boxes);
                            var local = LocalizationGraph.GetLocalizedPeptide(localizationGraphs[i].array, modPos, boxes, allPaths.First());
                            localizationCandidates.Add(new Tuple<int, Tuple<int, int>[]>(ids[i], local));
                        }

                        var psmGlyco = CreateGsm(theScan, scanIndex, ind, theScanBestPeptide, localizationCandidates, CommonParameters);

                        possibleMatches.Add(psmGlyco);
                    }
                }

                if (possibleMatches.Count != 0)
                {
                    possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
                }
            }

            return possibleMatches;
        }

        private GlycoSpectralMatch CreateGsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, List<Tuple<int, Tuple<int, int>[]>> glycanBox_localization, CommonParameters commonParameters)
        {
            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(glycanBox_localization.First().Item2, peptide);

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

        private List<GlycoSpectralMatch> FindModPepHash(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, List<int> allBinsToSearch, int scoreCutOff)
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

                    List<Product> products = theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                    HashSet<int> allPeaksForLocalization = new HashSet<int>(allBinsToSearch);

                    double bestLocalizedScore = 1;

                    List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();
                    List<int> ids = new List<int>();

                    int iDLow = GlycoPeptides.BinarySearchGetIndex(ModBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow);

                    while (iDLow < ModBoxes.Count() && PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + ModBoxes[iDLow].Mass))
                    {
                        int[] modPos = ModBox.GetAllPossibleModSites(theScanBestPeptide, ModBoxes[iDLow]);

                        if (modPos==null)
                        {
                            iDLow++;
                            continue;
                        }

                        var boxes = ModBox.BuildChildModBoxes(ModBoxes[iDLow].NumberOfMods, ModBoxes[iDLow].ModIds).ToArray();

                        LocalizationGraph localizationGraph = new LocalizationGraph(modPos.Length, boxes.Length);
                        localizationGraph.LocalizationMod(modPos, ModBoxes[iDLow], boxes, allPeaksForLocalization, products, theScanBestPeptide);

                        double currentLocalizationScore = localizationGraph.array.Last().Last().maxCost;
                        if (currentLocalizationScore > bestLocalizedScore)
                        {
                            ids.Clear();
                            localizationGraphs.Clear();
                            ids.Add(iDLow);
                            localizationGraphs.Add(localizationGraph);
                        }
                        else if (bestLocalizedScore > 0 && currentLocalizationScore == bestLocalizedScore)
                        {
                            ids.Add(iDLow);
                            localizationGraphs.Add(localizationGraph);        
                        }

                        iDLow++;
                    }

                    //In theory, the peptide_localization shouldn't be null, but it is possible that the real score is smaller than indexed score.
                    if (localizationGraphs.Count > 0)
                    {
                        List<Tuple<int, Tuple<int, int>[]>> localizationCandidates = new List<Tuple<int, Tuple<int, int>[]>>();
                        for (int i = 0; i < localizationGraphs.Count; i++)
                        {
                            int[] modPos = ModBox.GetAllPossibleModSites(theScanBestPeptide, ModBoxes[ids[i]]);
                            var boxes = ModBox.BuildChildModBoxes(ModBoxes[ids[i]].NumberOfMods, ModBoxes[ids[i]].ModIds).ToArray();
                            var allPaths = LocalizationGraph.GetAllPaths(localizationGraphs[i].array, boxes);
                            var local = LocalizationGraph.GetLocalizedPeptide(localizationGraphs[i].array, modPos, boxes, allPaths.First());
                            localizationCandidates.Add(new Tuple<int, Tuple<int, int>[]>(iDLow, local));
                        }

                        var psmGlyco = CreateGsmMod(theScan, scanIndex, ind, theScanBestPeptide, localizationCandidates, CommonParameters);

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

        private GlycoSpectralMatch CreateGsmMod(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, List<Tuple<int, Tuple<int, int>[]>> localizationCandidates, CommonParameters commonParameters)
        {
            var peptideWithMod = ModBox.GetTheoreticalPeptide(localizationCandidates.First().Item2, peptide, ModBoxes[localizationCandidates.First().Item1]);

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