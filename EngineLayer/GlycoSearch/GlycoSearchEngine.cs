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

namespace EngineLayer.GlycoSearch
{
    public class GlycoSearchEngine : ModernSearchEngine
    {
        protected readonly GlycoSpectralMatch[] GlobalCsms;

        private OpenSearchType OpenSearchType;
        // crosslinker molecule
        private readonly Crosslinker Crosslinker;
        private readonly bool CrosslinkSearchTopN;
        private readonly int TopN;

        private readonly bool SearchGlycan182;
        private readonly MassDiffAcceptor XLPrecusorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;

        private readonly char[] AllCrosslinkerSites;

        public GlycoSearchEngine(GlycoSpectralMatch[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, int currentPartition, CommonParameters commonParameters, MassDiffAcceptor massDiffAcceptor, OpenSearchType openSearchType, Crosslinker crosslinker, bool CrosslinkSearchTop, int CrosslinkSearchTopNum, bool searchGlycan182, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.OpenSearchType = openSearchType;
            this.Crosslinker = crosslinker;
            this.CrosslinkSearchTopN = CrosslinkSearchTop;
            this.TopN = CrosslinkSearchTopNum;

            AllCrosslinkerSites = Crosslinker.CrosslinkerModSites.ToCharArray().Concat(Crosslinker.CrosslinkerModSites2.ToCharArray()).Distinct().ToArray();

            XLPrecusorSearchMode = massDiffAcceptor;
            ProductSearchMode = new SingleAbsoluteAroundZeroSearchMode(20); //For Oxinium ion only

            SearchGlycan182 = searchGlycan182;
            if (OpenSearchType == OpenSearchType.NGlyco)
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
                
                //groupedGlycans = NGlycans.GroupBy(p => p.Mass).ToDictionary(p => p.Key, p => p.ToList());
                DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
            }

            if (OpenSearchType == OpenSearchType.OGlyco)
            {
                var OGlycans = Glycan.LoadGlycan(GlobalVariables.OGlycanLocation);
                OGlycanBoxes = Glycan.BuildGlycanBoxes(OGlycans.ToList(), 3).ToArray();
                //GroupedOGlycanBoxes = OGlycanBoxes.GroupBy(p => p.Mass).ToDictionary(p=>p.Key, p=>p.ToList());
            }
        }

        //public Dictionary<double, List<Glycan>> groupedGlycans { get; }
        private Glycan[] Glycans {get;}
        private Glycan[] DecoyGlycans { get; }
        //public Dictionary<int, List<GlycanBox>> GroupedOGlycanBoxes { get; }
        private GlycanBox[] OGlycanBoxes { get; }

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

                    //TO DO: limit the high bound limitation

                    // first-pass scoring
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, Double.PositiveInfinity, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);

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

                        GlycoSpectralMatch csm;
                        if (OpenSearchType == OpenSearchType.NGlyco)
                        {
                            //Glycopeptide Search
                            csm = FindNGlycopeptide(scan, bestPeptideScoreNotchList, scanIndex);
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
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing crosslink search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
                    }
                }
            });

            return new MetaMorpheusEngineResults(this);
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

        private GlycoSpectralMatch FindNGlycopeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int scanIndex)
        {
            GlycoSpectralMatch bestPsmCross = null;

            if (theScan.OxiniumIonNum < 2)
            {
                return bestPsmCross;
            }

            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
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

                if (theScan.OxiniumIonNum < 2)
                {
                    continue;
                }

                List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide[ind].BestPeptide, new string[] { "Nxt", "Nxs" });
                if (modPos.Count < 1)
                {
                    continue;
                }

                var possibleGlycanMassLow = theScan.PrecursorMass * (1 - 1E-5) - theScanBestPeptide[ind].BestPeptide.MonoisotopicMass;
                if (possibleGlycanMassLow < 200 || possibleGlycanMassLow > Glycans.Last().Mass)
                {
                    continue;
                }


                int iDLow = GlycoPeptides.BinarySearchGetIndex(Glycans.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);

                while (iDLow < Glycans.Count() && XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass + (double)Glycans[iDLow].Mass / 1E5) >= 0)
                {
                    double bestLocalizedScore = 0;
                    int bestSite = 0;
                    List<MatchedFragmentIon> bestMatchedIons = new List<MatchedFragmentIon>();
                    PeptideWithSetModifications peptideWithSetModifications = theScanBestPeptide[ind].BestPeptide;
                    foreach (int possibleSite in modPos)
                    {
                        var testPeptide = GlycoPeptides.GenerateGlycopeptide(possibleSite, theScanBestPeptide[ind].BestPeptide, Glycans[iDLow]);

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

                    var psmCross = new GlycoSpectralMatch(peptideWithSetModifications, theScanBestPeptide[ind].BestNotch, bestLocalizedScore, scanIndex, theScan, CommonParameters.DigestionParams, bestMatchedIons);
                    psmCross.Glycan = new List<Glycan> { Glycans[iDLow] };
                    psmCross.GlycanScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.Annotation.Contains('M')).ToList());
                    psmCross.DiagnosticIonScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.Annotation.Contains('D')).ToList());
                    psmCross.PeptideScore = psmCross.XLTotalScore - psmCross.GlycanScore - psmCross.DiagnosticIonScore;
                    psmCross.XlRank = new List<int> { ind };
                    psmCross.LinkPositions = new List<int> { bestSite - 1 }; //TO DO: ambiguity modification site
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
                foreach (var pm in possibleMatches)
                {
                    pm.ResolveAllAmbiguities();
                }
                bestPsmCross = possibleMatches.First();
            }
            return bestPsmCross;
        }

        private GlycoSpectralMatch FindOGlycopeptide(Ms2ScanWithSpecificMass theScan, List<BestPeptideScoreNotch> theScanBestPeptide, int scanIndex)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();
            GlycoSpectralMatch bestPsmCross = null;
            for (int ind = 0; ind < theScanBestPeptide.Count; ind++)
            {
                if (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass) >= 0)
                {
                    List<Product> products = theScanBestPeptide[ind].BestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                    var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
                    double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

                    var psmCrossSingle = new GlycoSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, score, scanIndex, theScan, CommonParameters.DigestionParams, matchedFragmentIons);
                    psmCrossSingle.CrossType = PsmCrossType.Single;
                    psmCrossSingle.XlRank = new List<int> { ind };
                    psmCrossSingle.ResolveAllAmbiguities();

                    possibleMatches.Add(psmCrossSingle);
                }
                //TO DO: add if the scan contains diagnostic ions
                else if ((theScan.PrecursorMass - theScanBestPeptide[ind].BestPeptide.MonoisotopicMass >= 200) && GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType)>=1)
                {
                    //Using glycanBoxes
                    var possibleGlycanMassLow = theScan.PrecursorMass * (1 - 1E-5) - theScanBestPeptide[ind].BestPeptide.MonoisotopicMass;
                    if (possibleGlycanMassLow < 200 || possibleGlycanMassLow > OGlycanBoxes.Last().Mass)
                    {
                        continue;
                    }

                    int iDLow = GlycoPeptides.BinarySearchGetIndex(OGlycanBoxes.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);

                    while(iDLow < OGlycanBoxes.Count() && (XLPrecusorSearchMode.Accepts(theScan.PrecursorMass, theScanBestPeptide[ind].BestPeptide.MonoisotopicMass + (double)OGlycanBoxes[iDLow].Mass/1E5) >= 0))
                    {
                        List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide[ind].BestPeptide, new string[] { "S", "T" });
                        if (modPos.Count >= OGlycanBoxes[iDLow].NumberOfGlycans)
                        {
                            

                            var fragmentsForEachGlycanLocalizedPossibility = GlycoPeptides.OGlyGetTheoreticalFragmentsUnlocalize(CommonParameters.DissociationType, modPos, theScanBestPeptide[ind].BestPeptide, OGlycanBoxes[iDLow]);
                            var bestMatchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycanLocalizedPossibility, CommonParameters);
                            double bestLocalizedScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons);

                            //var fragmentsForEachGlycanLocalizedPossibility = GlycoPeptides.OGlyGetTheoreticalFragments(commonParameters.DissociationType, modPos, theScanBestPeptide[ind].BestPeptide, OGlycanBoxes[iDLow]).ToList();
                            //double bestLocalizedScore = 0;
                            //List<MatchedFragmentIon> bestMatchedIons = new List<MatchedFragmentIon>();
                            //foreach (var setOfFragments in fragmentsForEachGlycanLocalizedPossibility)
                            //{
                            //    var matchedIons = MatchFragmentIons(theScan, setOfFragments.Item2.Item2, commonParameters);
                            //    double score = CalculatePeptideScore(theScan.TheScan, matchedIons, 0);

                            //    if (score > bestLocalizedScore)
                            //    {
                            //        bestLocalizedScore = score;
                            //        bestMatchedIons = matchedIons;
                            //    }
                            //}

                            var psmCross = new GlycoSpectralMatch(theScanBestPeptide[ind].BestPeptide, theScanBestPeptide[ind].BestNotch, bestLocalizedScore, scanIndex, theScan, CommonParameters.DigestionParams, bestMatchedIons);
                            psmCross.glycanBoxes = new List<GlycanBox> { OGlycanBoxes[iDLow] };
                            psmCross.XlRank = new List<int> { ind };
                            possibleMatches.Add(psmCross);

                        }
                        iDLow++;
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
                }
            }

            return bestPsmCross;
        }

    }
}