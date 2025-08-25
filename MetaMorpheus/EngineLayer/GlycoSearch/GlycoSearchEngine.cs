using EngineLayer.ModernSearch;
using MzLibUtil;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using EngineLayer;
using MassSpectrometry;
using EngineLayer.SpectrumMatch;

namespace EngineLayer.GlycoSearch
{
    public class GlycoSearchEngine : ModernSearchEngine
    {
        public static readonly double ToleranceForMassDifferentiation = 1e-9;
        private readonly int OxoniumIon204Index = 9;               // Check Glycan.AllOxoniumIons
        protected readonly List<GlycoSpectralMatch>[] GlobalGsms;  // Why don't we call it GlobalGsms?

        private GlycoSearchType GlycoSearchType;
        private readonly int TopN;              // DDA top Peak number.
        private readonly int _maxOGlycanNum;
        private readonly bool OxoniumIonFilter; // To filt Oxonium Ion before searching a spectrum as glycopeptides. If we filter spectrum, it must contain oxonium ions such as 204 (HexNAc). 
        private readonly string _oglycanDatabase;
        private readonly string _nglycanDatabase;

        private readonly Tolerance PrecusorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;

        private readonly List<int>[] SecondFragmentIndex;

        // The constructor for GlycoSearchEngine, we can load the parameter for the searhcing like mode, topN, maxOGlycanNum, oxoniumIonFilter, datsbase, etc.
        public GlycoSearchEngine(List<GlycoSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
             string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, int glycoSearchTopNum, int maxOGlycanNum, bool oxoniumIonFilter, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, fileSpecificParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalGsms = globalCsms;
            this.GlycoSearchType = glycoSearchType;
            this.TopN = glycoSearchTopNum;
            this._maxOGlycanNum = maxOGlycanNum;
            this.OxoniumIonFilter = oxoniumIonFilter;
            this._oglycanDatabase = oglycanDatabase;
            this._nglycanDatabase = nglycanDatabase;

            SecondFragmentIndex = secondFragmentIndex;
            PrecusorSearchMode = commonParameters.PrecursorMassTolerance;
            ProductSearchMode = new SinglePpmAroundZeroSearchMode(20); //For Oxonium ion only


            if (glycoSearchType == GlycoSearchType.OGlycanSearch) //if we do the O-glycan search, we need to load the O-glycan database and generate the glycoBox.
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == _oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(_maxOGlycanNum, false).OrderBy(p => p.Mass).ToArray(); //generate glycan box for O-glycan search
            }
            else if (glycoSearchType == GlycoSearchType.NGlycanSearch) //because the there is only one glycan in N-glycanpeptide, so we don't need to build the n-glycanBox here.
            {
                NGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false).OrderBy(p => p.Mass).ToArray();
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
            }
            else if (glycoSearchType == GlycoSearchType.N_O_GlycanSearch) //search both N-glycan and O-glycan is still not tested and build completely yet.
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == _oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(_maxOGlycanNum, false).OrderBy(p => p.Mass).ToArray();

                NGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false).OrderBy(p => p.Mass).ToArray();
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
            }

        }

        private Glycan[] NGlycans { get; }
        //private Glycan[] DecoyGlycans { get; }

        /// <summary>
        /// Run the glycoSearchEngine, the main function for the glycoSearchEngine.
        /// Four steps:
        /// (1) run a modern search engine to get the peptide candidates.
        /// (2) match the peptide candidates with the precursor mass.
        /// (3) use the mass shift to generate the route for the glycan localization.
        /// (4) evaluate the highest score for the glycan localization and generate the glycoSpectralMatch.
        /// </summary>
        /// <returns> SearchResult </returns>
        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing crosslink search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;  // MaxThreads = deafult is 7.
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray(); // We can do the parallel search on different threads
            Parallel.ForEach(threads, (scanIndex) =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();

                byte[] secondScoringTable = new byte[PeptideIndex.Count]; // We didn't use that right now.
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

                    // filtering the peptides candidate with the cufoff and limit the topN peptides.
                    if (idsOfPeptidesPossiblyObserved.Any()) 
                    {
                        scoreAtTopN = 0;
                        peptideCount = 0;
                        foreach (int id in idsOfPeptidesPossiblyObserved.OrderByDescending(p => scoringTable[p])) //from the higest score to the lowest score
                        {
                            if (scoringTable[id] < (int)byteScoreCutoff) //if the score is lower than the cutoff, we can skip this peptide.
                            {
                                continue;
                            }
                            peptideCount++;
                            if (peptideCount == TopN)
                            {
                                scoreAtTopN = scoringTable[id]; //ScoreAtTopN = The score of the last peptide in the TopN list.
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
                            gsms = FindOGlycopeptideHashLocal(scan, idsOfPeptidesTopN, scanIndex, (int)byteScoreCutoff); // Use the peptide candidate and the scan to generate the gsms.
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

                        if (GlobalGsms[scanIndex] == null)
                        {
                            GlobalGsms[scanIndex] = new List<GlycoSpectralMatch>(); //the first one finished task, create teh new gsms list.
                        }
                        else
                        {
                            gsms.AddRange(GlobalGsms[scanIndex]);
                            GlobalGsms[scanIndex].Clear();                       
                        }

                        Add2GlobalGsms(ref gsms, scanIndex);

                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing glyco search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
                    }   //percentProgress = 100, "Performing glyco search...1/1", NestedIds = 3.
                }
            });

            return new MetaMorpheusEngineResults(this); //Storage the result information into the result class.
        }

        private void Add2GlobalGsms(ref List<GlycoSpectralMatch> gsms, int scanIndex)
        {
            //keep top 10 candidates.
            double preScore = 0;
            int gsmsCount = 1;
            string preString = "";

            foreach (var gsm in gsms.Where(p => p != null).OrderByDescending(p => p.Score).ThenBy(c => c.FullSequence))
            {
                if (gsmsCount <= 10) 
                {
                    gsm.ResolveAllAmbiguities(); //Try to resolve any case that have the same sequence in the PSM.

                    if (gsmsCount == 1) //If the gsms number is 1, we don't need to check the score and sequence.
                    {
                        preScore = gsm.Score;
                        preString = gsm.FullSequence;

                        GlobalGsms[scanIndex].Add(gsm);
                        gsmsCount++;
                    }
                    else 
                    {
                        if (gsm.Score - preScore < ToleranceForMassDifferentiation && 
                        gsm.Score - preScore > -ToleranceForMassDifferentiation)
                        {
                            string currentString = gsm.FullSequence;

                            if (preString == currentString) //If peptides have the same sequence and their score is almost the same
                            {
                                foreach (SpectralMatchHypothesis bestMatchPeptide in gsm.BestMatchingBioPolymersWithSetMods) // We should add tje new ProteinMatch to the gsm. 
                                {                                                                                                               // Because the indentical sequence may from the different protein.
                                    GlobalGsms[scanIndex].Last().AddProteinMatch(bestMatchPeptide);
                                }
                            }
                            else
                            {
                                preString = currentString;
                                GlobalGsms[scanIndex].Add(gsm);
                                gsmsCount++;
                            }
                        }
                    }
                }
                else
                {
                    break;
                }
            }
        }

        //For FindOGlycan, generate the gsms for O-glycan search
        private GlycoSpectralMatch CreateGsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, Route localization, double[] oxoniumIonIntensities, List<LocalizationGraph> localizationGraphs)
        {
            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(localization, peptide);

            var fragmentsForEachGlycoPeptide = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.DissociationType, CommonParameters.CustomIons, peptide, peptideWithMod);

            var matchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);

            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

            var DiagnosticIonScore = CalculatePeptideScore(theScan.TheScan, matchedIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());

            var GlycanScore = CalculatePeptideScore(theScan.TheScan, matchedIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList());

            var PeptideScore = score - DiagnosticIonScore;

            var p = theScan.TheScan.MassSpectrum.Size * CommonParameters.ProductMassTolerance.GetRange(1000).Width / theScan.TheScan.MassSpectrum.Range.Width;

            int n = fragmentsForEachGlycoPeptide.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).Count();

            var allMatchedChildIons = new Dictionary<int, List<MatchedFragmentIon>>();

            foreach (var childScan in theScan.ChildScans)
            {
                var childFragments = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.MS2ChildScanDissociationType, CommonParameters.CustomIons, peptide, peptideWithMod);

                var matchedChildIons = MatchFragmentIons(childScan, childFragments, CommonParameters);

                n += childFragments.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).Count();

                if (matchedChildIons == null)
                {
                    continue;
                }

                allMatchedChildIons.Add(childScan.OneBasedScanNumber, matchedChildIons);
                double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);

                double childDiagnosticIonScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());
                double childGlycanScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList());

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
                psmGlyco.R138vs144 = oxoniumIonIntensities[4] / oxoniumIonIntensities[5]; // if the ratio is high, that means the glycan is more likely to be N-glycan. Oppsitely, ration is small means close to O-glycan.
            }

            return psmGlyco;
        }

        /// <summary>
        /// If the peptide mass is perfectly match with the precursor mass, we can directly generate the gsms for the peptide. Store the gsms into the possibleMatches.
        /// </summary>
        /// <param name="theScan"></param>
        /// <param name="scanIndex"></param>
        /// <param name="scoreCutOff"></param>
        /// <param name="theScanBestPeptide"> The peptide candidate </param>
        /// <param name="ind"></param>
        /// <param name="possibleMatches"> The space to store the gsms </param>
        private void FindSingle(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications theScanBestPeptide, int ind, ref List<GlycoSpectralMatch> possibleMatches)
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

        /// <summary>
        /// Match the mass of the peptide candidate with the precursor mass. Try to generate the Gsms for the Scan. Gsms will be stored in the possibleMatches.
        /// </summary>
        /// <param name="theScan"></param>
        /// <param name="scanIndex"></param>
        /// <param name="scoreCutOff"></param>
        /// <param name="theScanBestPeptide"> peptide candidate </param>
        /// <param name="ind"></param>
        /// <param name="possibleGlycanMassLow"> The precursor mass </param>
        /// <param name="oxoniumIonIntensities"></param>
        /// <param name="possibleMatches"> The space to store the gsms </param>
        private void FindOGlycan(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications theScanBestPeptide, int ind, double possibleGlycanMassLow, double[] oxoniumIonIntensities, ref List<GlycoSpectralMatch> possibleMatches)
        {
            // The glycanBoxes will be filtered by the oxonium ions. If the oxonium ions don't make sense, we will remove the glycanBox.


            int iDLow = GlycoPeptides.BinarySearchGetIndex(GlycanBox.OGlycanBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow); // try to find the index that closet match to the "possibleGlycanMassLow" within the glycanBox

            SortedDictionary<int, string> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" }); //list all of the possible glycoslation site/postition

            var localizationScan = theScan;
            List<Product> products = new List<Product>(); // product list for the theoretical fragment ions

            //For HCD-pd-ETD or CD-pd-EThcD type of data, we generate the different rpoducts.
            if (theScan.ChildScans.Count > 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType, CommonParameters.CustomIons))
            {
                localizationScan = theScan.ChildScans.First();
                theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
            }

            //For ETD type of data
            if (theScan.ChildScans.Count == 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType, CommonParameters.CustomIons))
            {
                theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
            }

            //Localization for O-glycopeptides only works on ETD related dissociationtype
            //No localization can be done with MS2-HCD spectrum
            //TO THINK: there is a special situation. The HCD only scan from  HCD-pd-EThcD data can be a glycopeptide, but there is no ETD, so there is no localization. What to do with this?
            bool is_HCD_only_data = !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType, CommonParameters.CustomIons) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType, CommonParameters.CustomIons);
            if (is_HCD_only_data) // In the HCD, there is no Y  ion, so we don't need to consider the modification here.
            {
                theScanBestPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            }

            double bestLocalizedScore = 0;

            List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>(); // if we also have ETD, then we will search the localization

            while (iDLow < GlycanBox.OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBox.OGlycanBoxes[iDLow].Mass))) // verify the glycan mass is invaild (within the range and match with mass shift)
            {
                if (OxoniumIonFilter && !GlycoPeptides.DiagonsticFilter(oxoniumIonIntensities, GlycanBox.OGlycanBoxes[iDLow])) // if the filter is turned on, we need to check does the oxoiums make sense.
                {
                    iDLow++; // if the oxonium ions don't make sense (there is no 204, or without their diagnostic ion), we can skip this glycan.
                    continue;
                }
                if (GraphCheck(modPos, GlycanBox.OGlycanBoxes[iDLow])) // the glycosite number should be larger than the possible glycan number.
                {
                    LocalizationGraph localizationGraph = new LocalizationGraph(modPos, GlycanBox.OGlycanBoxes[iDLow], GlycanBox.OGlycanBoxes[iDLow].ChildGlycanBoxes, iDLow);
                    LocalizationGraph.LocalizeOGlycan(localizationGraph, localizationScan, CommonParameters.ProductMassTolerance, products); //create the localization graph with the glycan mass and the possible glycosite.

                    double currentLocalizationScore = localizationGraph.TotalScore;
                    if (currentLocalizationScore > bestLocalizedScore) //Try to find the best glycanBox with the highest score.
                    {
                        bestLocalizedScore = currentLocalizationScore;
                        localizationGraphs.Clear();
                        localizationGraphs.Add(localizationGraph); // we only keep the best glycanBox and its localizationgraph.
                    }
                    else if ((is_HCD_only_data || bestLocalizedScore > 0) && (currentLocalizationScore <= bestLocalizedScore + 0.00000001 && currentLocalizationScore >= bestLocalizedScore - 0.00000001))
                    {
                        localizationGraphs.Add(localizationGraph);
                    }
                }

                iDLow++;
            }

            //In theory, the peptide_localization shouldn't be null, but it is possible that the real score is smaller than indexed score.
            if (localizationGraphs.Count > 0)
            {
                var firstPath = LocalizationGraph.GetFirstPath(localizationGraphs[0].array, localizationGraphs[0].ChildModBoxes); //Get the first path from the localization graph.
                var localizationCandidate = LocalizationGraph.GetLocalizedPath(localizationGraphs[0], firstPath); //Get the route of the localization from the first path inforation

                var psmGlyco = CreateGsm(theScan, scanIndex, ind, theScanBestPeptide, localizationCandidate, oxoniumIonIntensities, localizationGraphs); //Create the glycoSpectralMatch

                if (psmGlyco.Score > scoreCutOff)
                {
                    possibleMatches.Add(psmGlyco);
                }
            }
        }

        private void FindNGlycan(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications theScanBestPeptide, int ind, double possibleGlycanMassLow, double[] oxoniumIonIntensities, ref List<GlycoSpectralMatch> possibleMatches)
        {
            List<int> modPos_Nxs = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxs" }).Select(p => p.Key).ToList();
            List<int> modPos_Nxt = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxt" }).Select(p => p.Key).ToList();
            if (modPos_Nxs.Count < 1 && modPos_Nxt.Count < 1) // if there is no possible glycosylation site, we can skip this peptide.
            {
                return;
            }

            int iDLow = GlycoPeptides.BinarySearchGetIndex(NGlycans.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);
            while (iDLow < NGlycans.Length && PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + (double)NGlycans[iDLow].Mass / 1E5))
            {
                double bestLocalizedScore = scoreCutOff;
                int bestSite = 0;
                List<MatchedFragmentIon> bestMatchedIons = new List<MatchedFragmentIon>();
                PeptideWithSetModifications[] peptideWithSetModifications = new PeptideWithSetModifications[1];

                // Get the correct modification position based on the glycan target type
                List<int> modPos = NGlycans[iDLow].Target.ToString() == "Nxs" ? modPos_Nxs : modPos_Nxt;
                if (modPos.Count < 1)
                {
                    iDLow++;
                    continue; // if there is no possible glycosylation site, we can skip this glycan.
                }

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
                        peptideWithSetModifications[0] = testPeptide;
                        bestLocalizedScore = score;
                        bestSite = possibleSite;
                        bestMatchedIons = matchedIons;
                    }

                }

                if (peptideWithSetModifications[0] != null)
                {
                    var psmGlyco = new GlycoSpectralMatch(peptideWithSetModifications[0], 0, bestLocalizedScore, scanIndex, theScan, CommonParameters, bestMatchedIons);
                    psmGlyco.NGlycan = new List<Glycan> { NGlycans[iDLow] };
                    psmGlyco.GlycanScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList());
                    psmGlyco.DiagnosticIonScore = CalculatePeptideScore(theScan.TheScan, bestMatchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());
                    psmGlyco.PeptideScore = psmGlyco.Score - psmGlyco.GlycanScore - psmGlyco.DiagnosticIonScore;
                    psmGlyco.Rank = ind;
                    psmGlyco.NGlycanLocalizations = new List<int> { bestSite - 1 }; //TO DO: ambiguity modification site

                    if (oxoniumIonIntensities[5] <= 0.00000001)
                    {
                        psmGlyco.R138vs144 = 100000000;
                    }
                    else
                    {
                        psmGlyco.R138vs144 = oxoniumIonIntensities[4] / oxoniumIonIntensities[5];
                    }

                    possibleMatches.Add(psmGlyco);
                }

                iDLow++;
            }
        }

        // Conduct the search and generate the gsms for N-glycan search
        private List<GlycoSpectralMatch> FindNGlycopeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];

                //Considering coisolation, it doesn't mean it must from a glycopeptide even the scan contains oxonium ions.
                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass))
                {
                    FindSingle(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, ref possibleMatches);
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
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode);

                    //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    if (OxoniumIonFilter && oxoniumIonIntensities[OxoniumIon204Index] == 0)
                    {
                        continue;
                    }

                    //Find N-Glycan 
                    FindNGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, possibleGlycanMassLow, oxoniumIonIntensities, ref possibleMatches);

                }             
            }

            if (possibleMatches.Count != 0)
            {
                possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
            }
            return possibleMatches;
        }

        
        // Match the mass of the peptide candiate with the precursor mass, then try to generate the gsms object as output
        /// <summary>
        /// This is a general function for gsm generating. It was operated after the Modern Search.
        /// Two Step:
        /// (1) Match the mass of the peptide candiate with the precursor mass, then decide to go to which function to generate the gsms object.
        /// (2) Catch the gsms object and store it into the possibleMatches then return.
        /// </summary>
        /// <param name="theScan"> The MS2 Scan </param>
        /// <param name="idsOfPeptidesPossiblyObserved"> The peptide candidate from the modern Search </param>
        /// <param name="scanIndex"></param>
        /// <param name="scoreCutOff"></param>
        /// <returns> The Gsms collection.</returns>
        private List<GlycoSpectralMatch> FindOGlycopeptideHashLocal(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();


            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]]; // Get the peptide from the candidate list.

                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass)) // If the peptide mass is indentical to the precursor mass (or within the tolerance), we can directly search the glycopeptide.
                {
                    FindSingle(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, ref possibleMatches);
                }
                else if (theScan.PrecursorMass - theScanBestPeptide.MonoisotopicMass >= 100) //If not, we need to consider the glycan mass difference.
                {
                    //Filter by glycanBoxes mass difference.
                    var possibleGlycanMassLow = PrecusorSearchMode.GetMinimumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    var possibleGlycanMassHigh = PrecusorSearchMode.GetMaximumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    if (possibleGlycanMassHigh < GlycanBox.OGlycanBoxes.First().Mass || possibleGlycanMassLow > GlycanBox.OGlycanBoxes.Last().Mass)
                    {
                        continue; // if the glycan mass difference is out of the range of the glycan box, we can skip this peptide.
                    }

                    //Filter by OxoniumIon
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode);

                    //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    if (OxoniumIonFilter && oxoniumIonIntensities[9] == 0)
                    {
                        continue;
                    }

                    //Find O-Glycan
                    FindOGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, possibleGlycanMassLow, oxoniumIonIntensities, ref possibleMatches);
                }

                if (possibleMatches.Count != 0)
                {
                    possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
                }
            }

            return possibleMatches;
        }
        private List<GlycoSpectralMatch> Find_N_O_Glycopeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];

                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass))
                {
                    FindSingle(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, ref possibleMatches);

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
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode);

                    //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    if (OxoniumIonFilter && oxoniumIonIntensities[OxoniumIon204Index] == 0)
                    {
                        continue;
                    }

                    //Find N-Glycan 
                    FindNGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, possibleGlycanMassLow, oxoniumIonIntensities, ref possibleMatches);

                    //Find O-Glycan
                    FindOGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, possibleGlycanMassLow, oxoniumIonIntensities, ref possibleMatches);

                }

                if (possibleMatches.Count != 0)
                {
                    possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
                }
            }

            return possibleMatches;
        }

        /// <summary>
        /// Valid the Graph created by this modPos and glycanBox.
        /// Check if the motif in peptide is sufficient to cover the motif in glycanBox.
        /// </summary>
        /// <param name="modPos"></param>
        /// <param name="glycanBox"></param>
        /// <returns></returns>
        private bool GraphCheck(SortedDictionary<int, string> modPos, GlycanBox glycanBox)
        {
            // If the motifs number is less than the glycanBox, we can skip this graph.
            if (modPos.Count < glycanBox.NumberOfMods)
                return false;

            // Calculate the motif in glycanBox.
            var motifInBox = new Dictionary<string, int>();
            foreach (var modId in glycanBox.ModIds)
            {
                var motif = GlycanBox.GlobalOGlycans[modId].Target.ToString();

                if (!motifInBox.ContainsKey(motif))
                {
                    motifInBox[motif] = 0;
                }
                motifInBox[motif]++;
            }

            // Calculate the motif in peptide.
            var motifInPeptide = new Dictionary<string, int>();
            var modPos_motif = modPos.Values.ToArray();
            foreach (var motif in modPos_motif)
            {
                if (!motifInPeptide.ContainsKey(motif))
                {
                    motifInPeptide[motif] = 0;
                }
                motifInPeptide[motif]++;
            }

            // Check if the motif in peptide is sufficient to cover the motif in glycanBox.
            foreach (var motif in motifInBox)
            {
                if (!motifInPeptide.ContainsKey(motif.Key) || motifInPeptide[motif.Key] < motif.Value)
                {
                    return false; 
                }
            }

            return true;
        }


    }
}