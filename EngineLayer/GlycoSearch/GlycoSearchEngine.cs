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
        private readonly int OxoniumIon204Index = 9; //Check Glycan.AllOxoniumIons
        protected readonly List<GlycoSpectralMatch>[] GlobalCsms;

        private GlycoSearchType GlycoSearchType;
        private readonly int TopN;
        private readonly int _maxOGlycanNum;
        private readonly bool OxoniumIonFilter; //To filt Oxonium Ion before searching a spectrum as glycopeptides. If we filter spectrum, it must contain oxonium ions such as 204 (HexNAc). 
        private readonly string _oglycanDatabase;
        private readonly string _nglycanDatabase;

        private readonly Tolerance PrecusorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;

        private readonly List<int>[] SecondFragmentIndex;

        public GlycoSearchEngine(List<GlycoSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
             string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, int glycoSearchTopNum, int maxOGlycanNum, bool oxoniumIonFilter, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, fileSpecificParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.GlycoSearchType = glycoSearchType;
            this.TopN = glycoSearchTopNum;
            this._maxOGlycanNum = maxOGlycanNum;
            this.OxoniumIonFilter = oxoniumIonFilter;
            this._oglycanDatabase = oglycanDatabase;
            this._nglycanDatabase = nglycanDatabase;

            SecondFragmentIndex = secondFragmentIndex;
            PrecusorSearchMode = commonParameters.PrecursorMassTolerance;
            ProductSearchMode = new SinglePpmAroundZeroSearchMode(20); //For Oxonium ion only


            if (glycoSearchType == GlycoSearchType.OGlycanSearch)
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.GlobalOGlycanModifications = GlycanBox.BuildGlobalOGlycanModifications(GlycanBox.GlobalOGlycans);
                GlycanBox.OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(_maxOGlycanNum, false).OrderBy(p => p.Mass).ToArray();
            }
            else if (glycoSearchType == GlycoSearchType.NGlycanSearch)
            {
                NGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false).OrderBy(p => p.Mass).ToArray();
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
            }
            else if (glycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.GlobalOGlycanModifications = GlycanBox.BuildGlobalOGlycanModifications(GlycanBox.GlobalOGlycans);
                GlycanBox.OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(_maxOGlycanNum, false).OrderBy(p => p.Mass).ToArray();

                NGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false).OrderBy(p => p.Mass).ToArray();
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
            }

        }

        private Glycan[] NGlycans { get; }
        //private Glycan[] DecoyGlycans { get; }

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
                  
                    //Limit the high bound limitation, here assume it is possible to has max 3 Da shift. This allows for correcting precursor in the future.
                    var high_bound_limitation = scan.PrecursorMass + 1;

                    // first-pass scoring
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);

                    //child scan first-pass scoring
                    //List<int> childBinsToSearch = null;
                    //if (scan.ChildScans != null && scan.ChildScans.Count > 0 && CommonParameters.MS2ChildScanDissociationType != DissociationType.LowCID)
                    //{
                    //    Array.Clear(secondScoringTable, 0, secondScoringTable.Length);
                    //    childIdsOfPeptidesPossiblyObserved.Clear();

                    //    childBinsToSearch = new List<int>();

                    //    foreach (var aChildScan in scan.ChildScans)
                    //    {
                    //        var x = GetBinsToSearch(aChildScan, SecondFragmentIndex, CommonParameters.MS2ChildScanDissociationType);
                    //        childBinsToSearch.AddRange(x);
                    //    }

                    //    IndexedScoring(SecondFragmentIndex, childBinsToSearch, secondScoringTable, byteScoreCutoff, childIdsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.MS2ChildScanDissociationType);

                    //    foreach (var childId in childIdsOfPeptidesPossiblyObserved)
                    //    {
                    //        if (!idsOfPeptidesPossiblyObserved.Contains(childId))
                    //        {
                    //            idsOfPeptidesPossiblyObserved.Add(childId);
                    //        }
                    //        scoringTable[childId] = (byte)(scoringTable[childId] + secondScoringTable[childId]);
                    //    }
                    //}

                    // done with indexed scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        scoreAtTopN = 0;
                        peptideCount = 0;
                        foreach (int id in idsOfPeptidesPossiblyObserved.OrderByDescending(p => scoringTable[p]))
                        {
                            if (scoringTable[id] < (int)byteScoreCutoff)
                            {
                                continue;
                            }
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

                        if (GlycoSearchType == GlycoSearchType.OGlycanSearch)
                        {
                            gsms = FindOGlycopeptideHashLocal(scan, idsOfPeptidesTopN, scanIndex, (int)byteScoreCutoff);    
                        }
                        else if(GlycoSearchType == GlycoSearchType.NGlycanSearch)
                        {                     
                            gsms = FindNGlycopeptide(scan, idsOfPeptidesTopN, scanIndex, (int)byteScoreCutoff);
                        }
                        else
                        {
                            gsms = Find_N_O_Glycopeptide(scan, idsOfPeptidesTopN, scanIndex, (int)byteScoreCutoff);
                        }


                        if (gsms.Count == 0)
                        {
                            progress++;
                            continue;
                        }

                        if (GlobalCsms[scanIndex] == null)
                        {
                            GlobalCsms[scanIndex] = new List<GlycoSpectralMatch>();
                            GlobalCsms[scanIndex].AddRange(gsms.Where(p => p != null).OrderByDescending(p => p.Score).Take(10));
                        }
                        else
                        {
                            gsms.AddRange(GlobalCsms[scanIndex]);
                            GlobalCsms[scanIndex].Clear();
                            GlobalCsms[scanIndex].AddRange(gsms.Where(p => p != null).OrderByDescending(p => p.Score).Take(10));
                        }
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

        private List<GlycoSpectralMatch> FindNGlycopeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];

                //Considering coisolation, it doesn't mean it must from a glycopeptide even the scan contains oxonium ions.
                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass))
                {
                    List<Product> products = new List<Product>();
                    theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, products);
                    var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
                    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                    var psmSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, score, scanIndex, theScan, CommonParameters, matchedFragmentIons);
                    psmSingle.Rank = ind;
                    possibleMatches.Add(psmSingle);
                }
                else
                {
                    //Filter by glycanBoxes mass difference.
                    var possibleGlycanMassLow = PrecusorSearchMode.GetMinimumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    var possibleGlycanMassHigh = PrecusorSearchMode.GetMaximumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    if (possibleGlycanMassHigh < (double)NGlycans.First().Mass/1E5 || possibleGlycanMassLow > (double)NGlycans.Last().Mass/1E5)
                    {
                        continue;
                    }

                    //Filter by OxoniumIon
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);

                    //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    if (OxoniumIonFilter && oxoniumIonIntensities[OxoniumIon204Index] == 0)
                    {
                        continue;
                    }

                    List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxt", "Nxs" });
                    if (modPos.Count < 1)
                    {
                        continue;
                    }

                    int iDLow = GlycoPeptides.BinarySearchGetIndex(NGlycans.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);      

                    while (iDLow < NGlycans.Length && PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + (double)NGlycans[iDLow].Mass / 1E5))
                    {
                        double bestLocalizedScore = scoreCutOff;
                        int bestSite = 0;
                        List<MatchedFragmentIon> bestMatchedIons = new List<MatchedFragmentIon>();
                        PeptideWithSetModifications peptideWithSetModifications = theScanBestPeptide;
                        foreach (int possibleSite in modPos)
                        {
                            var testPeptide = GlycoPeptides.GenerateGlycopeptide(possibleSite, theScanBestPeptide, NGlycans[iDLow]);

                            List<Product> theoreticalProducts = new List<Product>();
                            testPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, theoreticalProducts);
                            theoreticalProducts = theoreticalProducts.Where(p => p.ProductType != ProductType.M).ToList();
                            theoreticalProducts.AddRange(GlycoPeptides.GetGlycanYIons(theScan.PrecursorMass, NGlycans[iDLow]));

                            //TO DO: the current MatchFragmentIons only match one charge states.
                            var matchedIons = MatchFragmentIons(theScan, theoreticalProducts, CommonParameters);

                            if (!GlycoPeptides.ScanTrimannosylCoreFilter(matchedIons, NGlycans[iDLow]))
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

                        var psmCross = new GlycoSpectralMatch(peptideWithSetModifications, 0, bestLocalizedScore, scanIndex, theScan, CommonParameters, bestMatchedIons);
                        psmCross.NGlycan = new List<Glycan> { NGlycans[iDLow] };
                        psmCross.GlycanScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList());
                        psmCross.DiagnosticIonScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());
                        psmCross.PeptideScore = psmCross.Score - psmCross.GlycanScore - psmCross.DiagnosticIonScore;
                        psmCross.Rank = ind;
                        psmCross.NGlycanLocalizations = new List<int> { bestSite - 1 }; //TO DO: ambiguity modification site
                        possibleMatches.Add(psmCross);

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
        private List<GlycoSpectralMatch> FindOGlycopeptideHashLocal(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];

                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass))
                {
                    List<Product> products = new List<Product>(); 
                    theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, products);
                    var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
                    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                    if (score > scoreCutOff)
                    {
                        var psmCrossSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, score, scanIndex, theScan, CommonParameters, matchedFragmentIons);
                        psmCrossSingle.Rank = ind;

                        possibleMatches.Add(psmCrossSingle);
                    }

                }
                else if (theScan.PrecursorMass - theScanBestPeptide.MonoisotopicMass >= 100) //Filter out unknow non-glycan modifications.
                {
                    //Filter by glycanBoxes mass difference.
                    var possibleGlycanMassLow = PrecusorSearchMode.GetMinimumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    var possibleGlycanMassHigh = PrecusorSearchMode.GetMaximumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    if (possibleGlycanMassHigh < GlycanBox.OGlycanBoxes.First().Mass || possibleGlycanMassLow > GlycanBox.OGlycanBoxes.Last().Mass)
                    {
                        continue;
                    }

                    //Filter by OxoniumIon
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);

                    //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    if (OxoniumIonFilter && oxoniumIonIntensities[9] == 0)
                    {
                        continue;
                    }


                    int iDLow = GlycoPeptides.BinarySearchGetIndex(GlycanBox.OGlycanBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow);

                    int[] modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" }).OrderBy(p=>p).ToArray();

                    //Localization for O-glycopeptides only works on ETD related dissociationtype
                    //No localization can be done with MS2-HCD spectrum
                    if ((theScan.ChildScans.Count == 0 || !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType)) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType))
                    {
                        while(iDLow < GlycanBox.OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBox.OGlycanBoxes[iDLow].Mass)))
                        {
                            if (modPos.Length >= GlycanBox.OGlycanBoxes[iDLow].NumberOfMods && GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, GlycanBox.OGlycanBoxes[iDLow]))
                            {
                                List<Product> hcdProducts = new List<Product>();
                                theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, hcdProducts);
                                var hcdMatchedFragmentIons = MatchFragmentIons(theScan, hcdProducts, CommonParameters);
                                double hcdScore = CalculatePeptideScore(theScan.TheScan, hcdMatchedFragmentIons);

                                var psmCrossSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, hcdScore, scanIndex, theScan, CommonParameters, hcdMatchedFragmentIons);
                                psmCrossSingle.Rank = ind;

                                possibleMatches.Add(psmCrossSingle);
                                break;
                            }

                            iDLow++;
                            
                        }

                        continue;
                    }

                    var localizationScan = theScan;

                    if (theScan.ChildScans.Count > 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType))
                    {
                        localizationScan = theScan.ChildScans.First();
                    }

                    List<Product> products = new List<Product>();
                    theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);                    

                    double bestLocalizedScore = 0;

                    List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();

                    while (iDLow < GlycanBox.OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBox.OGlycanBoxes[iDLow].Mass)))
                    {
                        if (modPos.Length >= GlycanBox.OGlycanBoxes[iDLow].NumberOfMods && GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, GlycanBox.OGlycanBoxes[iDLow]))
                        {                            
                            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, GlycanBox.OGlycanBoxes[iDLow], GlycanBox.OGlycanBoxes[iDLow].ChildGlycanBoxes, iDLow);
                            LocalizationGraph.LocalizeOGlycan(localizationGraph, localizationScan, CommonParameters.ProductMassTolerance, products);

                            double currentLocalizationScore = localizationGraph.TotalScore;
                            if (currentLocalizationScore > bestLocalizedScore)
                            {
                                bestLocalizedScore = currentLocalizationScore;
                                localizationGraphs.Clear();
                                localizationGraphs.Add(localizationGraph);
                            }
                            else if (bestLocalizedScore > 0 && (currentLocalizationScore <= bestLocalizedScore + 0.00000001 && currentLocalizationScore >= bestLocalizedScore - 0.00000001))
                            {
                                localizationGraphs.Add(localizationGraph);
                            }
                        }

                        iDLow++;
                    }

                    //In theory, the peptide_localization shouldn't be null, but it is possible that the real score is smaller than indexed score.
                    if (localizationGraphs.Count > 0)
                    {
                        var firstPath = LocalizationGraph.GetFirstPath(localizationGraphs[0].array, localizationGraphs[0].ChildModBoxes);
                        var localizationCandidate = LocalizationGraph.GetLocalizedPath(localizationGraphs[0], firstPath);

                        var psmGlyco = CreateGsm(theScan, scanIndex, ind, theScanBestPeptide, localizationCandidate, oxoniumIonIntensities, localizationGraphs);

                        if (psmGlyco.Score > scoreCutOff)
                        {
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

        //List<Tuple<int, Tuple<int, int>[]>> <glycanBoxId, <mod site, glycan id>>
        private GlycoSpectralMatch CreateGsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, Route localization, double[] oxoniumIonIntensities, List<LocalizationGraph> localizationGraphs)
        {
            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(localization, peptide);

            var fragmentsForEachGlycoPeptide = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.DissociationType, peptide, peptideWithMod);

            var matchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);

            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

            var DiagnosticIonScore = CalculatePeptideScore(theScan.TheScan, matchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());

            var GlycanScore = CalculatePeptideScore(theScan.TheScan, matchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList());

            var PeptideScore = score - DiagnosticIonScore;

            var p = theScan.TheScan.MassSpectrum.Size * CommonParameters.ProductMassTolerance.GetRange(1000).Width / theScan.TheScan.MassSpectrum.Range.Width;

            int n = fragmentsForEachGlycoPeptide.Where(p => p.ProductType == ProductType.c || p.ProductType == ProductType.zDot).Count();

            var allMatchedChildIons = new Dictionary<int, List<MatchedFragmentIon>>();

            foreach (var childScan in theScan.ChildScans)
            {
                var childFragments = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.MS2ChildScanDissociationType, peptide, peptideWithMod);

                var matchedChildIons = MatchFragmentIons(childScan, childFragments, CommonParameters);

                n += childFragments.Where(p => p.ProductType == ProductType.c || p.ProductType == ProductType.zDot).Count();

                if (matchedChildIons == null)
                {
                    continue;
                }

                allMatchedChildIons.Add(childScan.OneBasedScanNumber, matchedChildIons);
                double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);

                double childDiagnosticIonScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());
                double childGlycanScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList());

                DiagnosticIonScore += childDiagnosticIonScore;
                GlycanScore += childGlycanScore;

                PeptideScore += childScore - childDiagnosticIonScore;
                //TO THINK:may think a different way to use childScore
                score += childScore;

                p += childScan.TheScan.MassSpectrum.Size * CommonParameters.ProductMassTolerance.GetRange(1000).Width / childScan.TheScan.MassSpectrum.Range.Width;

            }

            var psmGlyco = new GlycoSpectralMatch(peptideWithMod, 0, PeptideScore, scanIndex, theScan, CommonParameters, matchedIons);
            
            //TO DO: This p is from childScan p, it works for HCD-pd-EThcD, which may not work for other type.
            psmGlyco.ScanInfo_p = p;

            psmGlyco.Thero_n = n;

            psmGlyco.Rank = rank;

            psmGlyco.DiagnosticIonScore = DiagnosticIonScore;

            psmGlyco.GlycanScore = GlycanScore;

            psmGlyco.ChildMatchedFragmentIons = allMatchedChildIons;

            psmGlyco.LocalizationGraphs = localizationGraphs;

            if (oxoniumIonIntensities[5] <= 0.00000001)
            {
                psmGlyco.R138vs144 = 100000000;
            }
            else
            {
                psmGlyco.R138vs144 = oxoniumIonIntensities[4] / oxoniumIonIntensities[5];
            }     

            return psmGlyco;
        }

        private List<GlycoSpectralMatch> Find_N_O_Glycopeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];

                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass))
                {
                    List<Product> products = new List<Product>();
                    theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, products);
                    var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
                    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                    if (score > scoreCutOff)
                    {
                        var psmCrossSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, score, scanIndex, theScan, CommonParameters, matchedFragmentIons);
                        psmCrossSingle.Rank = ind;

                        possibleMatches.Add(psmCrossSingle);
                    }

                }
                else if (theScan.PrecursorMass - theScanBestPeptide.MonoisotopicMass >= 100) //Filter out unknow non-glycan modifications.
                {
                    //Filter by glycanBoxes mass difference.
                    var possibleGlycanMassLow = PrecusorSearchMode.GetMinimumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    var possibleGlycanMassHigh = PrecusorSearchMode.GetMaximumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    if (possibleGlycanMassHigh < GlycanBox.OGlycanBoxes.First().Mass)
                    {
                        continue;
                    }

                    //Filter by OxoniumIon
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);

                    //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    if (OxoniumIonFilter && oxoniumIonIntensities[OxoniumIon204Index] == 0)
                    {
                        continue;
                    }

                    //Find N-Glycan 
                    {
                        List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxt", "Nxs" });
                        if (modPos.Count < 1)
                        {
                            continue;
                        }

                        int iDLow = GlycoPeptides.BinarySearchGetIndex(NGlycans.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);

                        while (iDLow < NGlycans.Length && PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + (double)NGlycans[iDLow].Mass / 1E5))
                        {
                            double bestLocalizedScore = scoreCutOff;
                            int bestSite = 0;
                            List<MatchedFragmentIon> bestMatchedIons = new List<MatchedFragmentIon>();
                            PeptideWithSetModifications peptideWithSetModifications = theScanBestPeptide;
                            foreach (int possibleSite in modPos)
                            {
                                var testPeptide = GlycoPeptides.GenerateGlycopeptide(possibleSite, theScanBestPeptide, NGlycans[iDLow]);

                                List<Product> theoreticalProducts = new List<Product>();
                                testPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, theoreticalProducts);
                                theoreticalProducts = theoreticalProducts.Where(p => p.ProductType != ProductType.M).ToList();
                                theoreticalProducts.AddRange(GlycoPeptides.GetGlycanYIons(theScan.PrecursorMass, NGlycans[iDLow]));

                                //TO DO: the current MatchFragmentIons only match one charge states.
                                var matchedIons = MatchFragmentIons(theScan, theoreticalProducts, CommonParameters);

                                if (!GlycoPeptides.ScanTrimannosylCoreFilter(matchedIons, NGlycans[iDLow]))
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

                            var psmCross = new GlycoSpectralMatch(peptideWithSetModifications, 0, bestLocalizedScore, scanIndex, theScan, CommonParameters, bestMatchedIons);
                            psmCross.NGlycan = new List<Glycan> { NGlycans[iDLow] };
                            psmCross.GlycanScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList());
                            psmCross.DiagnosticIonScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());
                            psmCross.PeptideScore = psmCross.Score - psmCross.GlycanScore - psmCross.DiagnosticIonScore;
                            psmCross.Rank = ind;
                            psmCross.NGlycanLocalizations = new List<int> { bestSite - 1 }; //TO DO: ambiguity modification site
                            possibleMatches.Add(psmCross);

                            iDLow++;
                        }
                    }


                    //Find O-Glycan
                    {
                        int iDLow = GlycoPeptides.BinarySearchGetIndex(GlycanBox.OGlycanBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow);

                        int[] modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" }).OrderBy(p => p).ToArray();

                        //Localization for O-glycopeptides only works on ETD related dissociationtype
                        //No localization can be done with MS2-HCD spectrum
                        if ((theScan.ChildScans.Count == 0 || !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType)) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType))
                        {
                            while (iDLow < GlycanBox.OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBox.OGlycanBoxes[iDLow].Mass)))
                            {
                                if (modPos.Length >= GlycanBox.OGlycanBoxes[iDLow].NumberOfMods && GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, GlycanBox.OGlycanBoxes[iDLow]))
                                {
                                    List<Product> hcdProducts = new List<Product>();
                                    theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, hcdProducts);
                                    var hcdMatchedFragmentIons = MatchFragmentIons(theScan, hcdProducts, CommonParameters);
                                    double hcdScore = CalculatePeptideScore(theScan.TheScan, hcdMatchedFragmentIons);

                                    var psmCrossSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, hcdScore, scanIndex, theScan, CommonParameters, hcdMatchedFragmentIons);
                                    psmCrossSingle.Rank = ind;

                                    possibleMatches.Add(psmCrossSingle);
                                    break;
                                }

                                iDLow++;

                            }

                            continue;
                        }

                        var localizationScan = theScan;

                        if (theScan.ChildScans.Count > 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType))
                        {
                            localizationScan = theScan.ChildScans.First();
                        }

                        List<Product> products = new List<Product>();
                        theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

                        double bestLocalizedScore = 0;

                        List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();

                        while (iDLow < GlycanBox.OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBox.OGlycanBoxes[iDLow].Mass)))
                        {
                            if (modPos.Length >= GlycanBox.OGlycanBoxes[iDLow].NumberOfMods && GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, GlycanBox.OGlycanBoxes[iDLow]))
                            {
                                LocalizationGraph localizationGraph = new LocalizationGraph(modPos, GlycanBox.OGlycanBoxes[iDLow], GlycanBox.OGlycanBoxes[iDLow].ChildGlycanBoxes, iDLow);
                                LocalizationGraph.LocalizeOGlycan(localizationGraph, localizationScan, CommonParameters.ProductMassTolerance, products);

                                double currentLocalizationScore = localizationGraph.TotalScore;
                                if (currentLocalizationScore > bestLocalizedScore)
                                {
                                    bestLocalizedScore = currentLocalizationScore;
                                    localizationGraphs.Clear();
                                    localizationGraphs.Add(localizationGraph);
                                }
                                else if (bestLocalizedScore > 0 && (currentLocalizationScore <= bestLocalizedScore + 0.00000001 && currentLocalizationScore >= bestLocalizedScore - 0.00000001))
                                {
                                    localizationGraphs.Add(localizationGraph);
                                }
                            }

                            iDLow++;
                        }

                        //In theory, the peptide_localization shouldn't be null, but it is possible that the real score is smaller than indexed score.
                        if (localizationGraphs.Count > 0)
                        {
                            var firstPath = LocalizationGraph.GetFirstPath(localizationGraphs[0].array, localizationGraphs[0].ChildModBoxes);
                            var localizationCandidate = LocalizationGraph.GetLocalizedPath(localizationGraphs[0], firstPath);

                            var psmGlyco = CreateGsm(theScan, scanIndex, ind, theScanBestPeptide, localizationCandidate, oxoniumIonIntensities, localizationGraphs);

                            if (psmGlyco.Score > scoreCutOff)
                            {
                                possibleMatches.Add(psmGlyco);
                            }
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

    }
}