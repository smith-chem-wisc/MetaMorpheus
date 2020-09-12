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
using Chemistry;

namespace EngineLayer.GlycoSearch
{
    public class GlycoSearchEngine : ModernSearchEngine
    {
        public static readonly double ToleranceForMassDifferentiation = 1e-9;
        private readonly int OxoniumIon204Index = 9; //Check Glycan.AllOxoniumIons
        protected readonly List<GlycoSpectralMatch>[] GlobalCsms;

        private GlycoSearchType GlycoSearchType;
        private readonly int TopN;
        private readonly int _maxOGlycanNum;
        private readonly bool OxoniumIonFilter; //To filt Oxonium Ion before searching a spectrum as glycopeptides. If we filter spectrum, it must contain oxonium ions such as 204 (HexNAc). 
        private readonly bool IndexChildScan; //Whether to index the child scan or not is a question.
        private readonly string _oglycanDatabase;
        private readonly string _nglycanDatabase;

        private readonly Tolerance PrecusorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;

        private readonly List<int>[] SecondFragmentIndex;

        public GlycoSearchEngine(List<GlycoSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
             string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, int glycoSearchTopNum, int maxOGlycanNum, bool oxoniumIonFilter, bool indexChildScan, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, fileSpecificParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.GlycoSearchType = glycoSearchType;
            this.TopN = glycoSearchTopNum;
            this._maxOGlycanNum = maxOGlycanNum;
            this.OxoniumIonFilter = oxoniumIonFilter;
            this.IndexChildScan = indexChildScan;
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
                //NGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);

            }
            else if (glycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.GlobalOGlycanModifications = GlycanBox.BuildGlobalOGlycanModifications(GlycanBox.GlobalOGlycans);
                GlycanBox.OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(_maxOGlycanNum, false).OrderBy(p => p.Mass).ToArray();

                NGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false).OrderBy(p => p.Mass).ToArray();
                //TO THINK: Glycan Decoy database.
                //NGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);

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
                    List<int> allBinsToSearch = GetGlycanBinsToSearch(scan, FragmentIndex, (GlycoSearchType == GlycoSearchType.NGlycanSearch || GlycoSearchType == GlycoSearchType.N_O_GlycanSearch));

                    //SecondFragmentIndex will only be used for HCD-trig-ETD data for N-Glycopeptide. 
                    //If SecondFragmentIndex == null, we only use the b/y FragmentIndex for child scans.
                    if (IndexChildScan && SecondFragmentIndex == null && scan.ChildScans != null && scan.ChildScans.Count > 0)
                    {
                        List<int> childBinsToSearch = new List<int>();      

                        foreach (var aChildScan in scan.ChildScans)
                        {
                            var x = GetGlycanBinsToSearch(aChildScan, FragmentIndex, (GlycoSearchType == GlycoSearchType.NGlycanSearch || GlycoSearchType == GlycoSearchType.N_O_GlycanSearch));
                            childBinsToSearch.AddRange(x);
                        }

                        allBinsToSearch.AddRange(childBinsToSearch);
                    }
                  
                    //Limit the high bound limitation, here assume it is possible to has max 3 Da shift. This allows for correcting precursor in the future.
                    var high_bound_limitation = scan.PrecursorMass + 1;

                    // first-pass scoring
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);

                    //child scan first - pass scoring
                    if (IndexChildScan && SecondFragmentIndex != null && scan.ChildScans != null && scan.ChildScans.Count > 0)
                    {
                        //SecondFragmentIndex will only be used for HCD-trig-ETD data for N-Glycopeptide. 
                        List<int> childBinsToSearch = null;

                        Array.Clear(secondScoringTable, 0, secondScoringTable.Length);
                        childIdsOfPeptidesPossiblyObserved.Clear();

                        childBinsToSearch = new List<int>();

                        foreach (var aChildScan in scan.ChildScans)
                        {
                            var x = GetGlycanBinsToSearch(aChildScan, SecondFragmentIndex, (GlycoSearchType == GlycoSearchType.NGlycanSearch || GlycoSearchType == GlycoSearchType.N_O_GlycanSearch));
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
                        }
                        else
                        {
                            gsms.AddRange(GlobalCsms[scanIndex]);
                            GlobalCsms[scanIndex].Clear();                       
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
                    }
                }
            });

            return new MetaMorpheusEngineResults(this);
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
                    gsm.ResolveAllAmbiguities();

                    if (gsmsCount == 1)
                    {
                        preScore = gsm.Score;
                        preString = gsm.FullSequence;

                        GlobalCsms[scanIndex].Add(gsm);
                        gsmsCount++;
                    }
                    else
                    {
                        if (gsm.Score - preScore < ToleranceForMassDifferentiation &&
                        gsm.Score - preScore > -ToleranceForMassDifferentiation)
                        {
                            string currentString = gsm.FullSequence;

                            if (preString == currentString)
                            {
                                foreach (var bestMatchPeptide in gsm.BestMatchingPeptides)
                                {
                                    GlobalCsms[scanIndex].Last().AddProteinMatch(bestMatchPeptide, gsm.PeptidesToMatchingFragments[bestMatchPeptide.Peptide]);

                                }
                            }
                            else
                            {
                                preString = currentString;
                                GlobalCsms[scanIndex].Add(gsm);
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

        /// <summary>
        /// Deprecated. (The function is in IndexingEngine, which is now used in Crosslink and Glyco Search)
        /// </summary>
        private List<int> GetGlycanBinsToSearch(Ms2ScanWithSpecificMass scan, List<int>[] FragmentIndex, bool AddGlycanShiftBin)
        {
            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();

            foreach (var envelope in scan.ExperimentalFragments)
            {
                // assume charge state 1 to calculate mass tolerance
                double experimentalFragmentMass = envelope.MonoisotopicMass;

                // get theoretical fragment bins within mass tolerance
                int obsFragmentFloorMass = (int)Math.Floor((CommonParameters.ProductMassTolerance.GetMinimumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);
                int obsFragmentCeilingMass = (int)Math.Ceiling((CommonParameters.ProductMassTolerance.GetMaximumValue(experimentalFragmentMass)) * FragmentBinsPerDalton);

                // prevents double-counting peaks close in m/z and lower-bound out of range exceptions
                if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)
                {
                    obsFragmentFloorMass = obsPreviousFragmentCeilingMz;
                }
                obsPreviousFragmentCeilingMz = obsFragmentCeilingMass + 1;

                // prevent upper-bound index out of bounds errors;
                // lower-bound is handled by the previous "if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)" statement
                if (obsFragmentCeilingMass >= FragmentIndex.Length)
                {
                    obsFragmentCeilingMass = FragmentIndex.Length - 1;

                    if (obsFragmentFloorMass >= FragmentIndex.Length)
                    {
                        obsFragmentFloorMass = FragmentIndex.Length - 1;
                    }
                }

                // search mass bins within a tolerance
                for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                {
                    if (FragmentIndex[fragmentBin] != null)
                    {
                        binsToSearch.Add(fragmentBin);
                    }
                }

                // add glycan add complementary ions
                if (AddGlycanShiftBin)
                {
                    //Here, we assume the fragment mass shifts contains a HexNAc 203.07937 or Cross-ring mass 83.03819. 
                    //Note this is for HCD or sceHCD spectra of N-Glycan.
                    //The fragment theory of N-Glycopeptide can be find in pGlyco2 paper.
                    double[] glycanMassShifts = new double[2] { 203.07937, 83.03819 };

                    foreach (var protonMassShift in glycanMassShifts)
                    {  
                        int compFragmentFloorMass = obsFragmentFloorMass - (int)Math.Round(protonMassShift * FragmentBinsPerDalton);
                        int compFragmentCeilingMass = obsFragmentCeilingMass - (int)Math.Round(protonMassShift * FragmentBinsPerDalton);

                        // prevent index out of bounds errors
                        if (compFragmentCeilingMass >= FragmentIndex.Length)
                        {
                            compFragmentCeilingMass = FragmentIndex.Length - 1;

                            if (compFragmentFloorMass >= FragmentIndex.Length)
                            {
                                compFragmentFloorMass = FragmentIndex.Length - 1;
                            }
                        }
                        if (compFragmentFloorMass < 0)
                        {
                            compFragmentFloorMass = 0;
                        }

                        for (int fragmentBin = compFragmentFloorMass; fragmentBin <= compFragmentCeilingMass; fragmentBin++)
                        {
                            if (FragmentIndex[fragmentBin] != null)
                            {
                                binsToSearch.Add(fragmentBin);
                            }
                        }
                    }
                }
            }

            return binsToSearch;
        }

        private void FindSingle(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications theScanBestPeptide, int ind, ref List<GlycoSpectralMatch> possibleMatches)
        {
            List<Product> products = new List<Product>();
            theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, products);
            var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
            double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

            var allMatchedChildIons = new Dictionary<int, List<MatchedFragmentIon>>();

            foreach (var childScan in theScan.ChildScans)
            {
                //People always use CID with low res. This special code works for Nic Scott's data. (High-HCD, High-EThcD, Low-CID, High-secHCD)
                if (childScan.TheScan.DissociationType == DissociationType.CID)
                {
                    continue;
                }

                List<Product> childFragments = new List<Product>();
                theScanBestPeptide.Fragment(CommonParameters.MS2ChildScanDissociationType, FragmentationTerminus.Both, childFragments);

                var matchedChildIons = MatchFragmentIons(childScan, childFragments, CommonParameters);

                if (matchedChildIons == null)
                {
                    continue;
                }

                allMatchedChildIons.Add(childScan.OneBasedScanNumber, matchedChildIons);
                double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);

                //TO THINK:may think a different way to use childScore
                score += childScore;
            }

            if (score > scoreCutOff)
            {
                var psmCrossSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, score, scanIndex, theScan, CommonParameters, matchedFragmentIons);
                psmCrossSingle.Rank = ind;
                if (allMatchedChildIons.Count() >0)
                {
                    psmCrossSingle.ChildMatchedFragmentIons = allMatchedChildIons;
                }
                possibleMatches.Add(psmCrossSingle);
            }
        }

        //For FindOGlycan
        private GlycoSpectralMatch CreateOGsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, Route localization, double[] oxoniumIonIntensities, List<LocalizationGraph> localizationGraphs)
        {
            var peptideWithMod = GlycoPeptides.OGlyGetTheoreticalPeptide(localization, peptide);

            var fragmentsForEachGlycoPeptide = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.DissociationType, peptide, peptideWithMod);

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
                //People always use CID with low res. This special code works for Nic Scott's data. (High-HCD, High-EThcD, Low-CID, High-secHCD)
                if (childScan.TheScan.DissociationType == DissociationType.CID)
                {
                    continue;
                }
                var childFragments = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.MS2ChildScanDissociationType, peptide, peptideWithMod);

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

            psmGlyco.OxoniumIonIntensity = oxoniumIonIntensities;

            var product_nn = GlycoPeptides.GetIndicatorYIon(peptide.MonoisotopicMass, "NN");
            psmGlyco.PepNN = GlycoPeptides.MatchIndicatorYIon(theScan, product_nn, CommonParameters);

            var product_nh = GlycoPeptides.GetIndicatorYIon(peptide.MonoisotopicMass, "NH");
            psmGlyco.PepNH = GlycoPeptides.MatchIndicatorYIon(theScan, product_nh, CommonParameters);

            //psmGlyco.LongestconcatenatedYion

            return psmGlyco;
        }
        
        private void FindOGlycan(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications theScanBestPeptide, int ind, double possibleGlycanMassLow, double[] oxoniumIonIntensities, ref List<GlycoSpectralMatch> possibleMatches)
        {
            int iDLow = GlycoPeptides.BinarySearchGetIndex(GlycanBox.OGlycanBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow);

            int[] modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" }).OrderBy(p => p).ToArray();

            //Localization for O-glycopeptides only works on ETD related dissociationtype
            //No localization can be done with MS2-HCD spectrum
            bool is_HCD_only_data = !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType);
            if (is_HCD_only_data)
            {
                List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();

                while (iDLow < GlycanBox.OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBox.OGlycanBoxes[iDLow].Mass)))
                {
                    if (modPos.Length >= GlycanBox.OGlycanBoxes[iDLow].NumberOfMods && GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, GlycanBox.OGlycanBoxes[iDLow]))
                    {
                        //Construct the localizationGraph, but didn't run the localization. The purpose is to use the GlycanBox Infomation.
                        LocalizationGraph localizationGraph = new LocalizationGraph(modPos, GlycanBox.OGlycanBoxes[iDLow], GlycanBox.OGlycanBoxes[iDLow].ChildGlycanBoxes, iDLow);
                        localizationGraphs.Add(localizationGraph);
                    }

                    iDLow++;
                }
                if (localizationGraphs.Count > 0)
                {   
                    var route = LocalizationGraph.GetAnyOnePath(localizationGraphs[0]);

                    var psmGlyco = CreateOGsm(theScan, scanIndex, ind, theScanBestPeptide, route, oxoniumIonIntensities, localizationGraphs);

                    if (psmGlyco.Score > scoreCutOff)
                    {
                        possibleMatches.Add(psmGlyco);
                    }
                }
            }
            else
            {
                var localizationScan = theScan;

                //We need to consider which scan should be used for localization.               
                //For ETD_only or ETD-pd-HCD type of data. Just keep the first scan as localization scan.
                //For ETD-pd-ETD type of data. There is no implimentation yet. Is there such type of data exist? Is there another type of data.               
                //For HCD-pd-ETD or CD-pd-EThcD type of data. We should consider the first childScan as localization scan.
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
                    var firstPath = LocalizationGraph.GetFirstPath(localizationGraphs[0].array, localizationGraphs[0].ChildModBoxes);
                    var localizationCandidate = LocalizationGraph.GetLocalizedPath(localizationGraphs[0], firstPath);

                    var psmGlyco = CreateOGsm(theScan, scanIndex, ind, theScanBestPeptide, localizationCandidate, oxoniumIonIntensities, localizationGraphs);

                    if (psmGlyco.Score > scoreCutOff)
                    {
                        possibleMatches.Add(psmGlyco);
                    }
                }
            }
           
        }

        //For Find NGlycan. Create NglycoPeptideSpectrumMatch and handle N-Glycan localization. 
        private GlycoSpectralMatch CreateNGsm(Ms2ScanWithSpecificMass theScan, double[] oxoniumIonIntensities, List<int> modPos, PeptideWithSetModifications[] bestPeptides, 
            List<MatchedFragmentIon>[] bestMatchedIons, Dictionary<int, List<MatchedFragmentIon>>[] bestChildMathcedIons, double[][] bestScores, int theo_n, int scanIndex, int iDLow, int ind)
        {
            //if (!GlycoPeptides.ScanTrimannosylCoreFilter(matchedIons, NGlycans[iDLow]))
            //{
            //    continue;
            //}

            int bestModIndex = 0;

            //There doesn't contain Level2 for N-Glycopeptide localization Level, since we only consider one glycan on one peptide.
            LocalizationLevel localizationLevel = LocalizationLevel.Level1b;
            Dictionary<int, double> siteProbability = new Dictionary<int, double>();

            //Assign localization level and localization site probability. 
            if (modPos.Count > 1)
            {
                bestModIndex = Array.IndexOf(bestScores[0], bestScores[0].Max());

                //Calculate the site specific localization probability.
                var p = theScan.TheScan.MassSpectrum.Size * CommonParameters.ProductMassTolerance.GetRange(1000).Width / theScan.TheScan.MassSpectrum.Range.Width;        

                foreach (var childScan in theScan.ChildScans)
                {
                    p += childScan.TheScan.MassSpectrum.Size * CommonParameters.ProductMassTolerance.GetRange(1000).Width / childScan.TheScan.MassSpectrum.Range.Width;             
                }

                double[] Ps = new double[modPos.Count];
                for (int i = 0; i < modPos.Count; i++)
                {                  
                    double k = bestScores[1][i]; 

                    //To understand the math, ref to "phosphoRS" papar.    
                    var cp = 1 / (1 - MathNet.Numerics.Distributions.Binomial.CDF(p, theo_n, k) + MathNet.Numerics.Distributions.Binomial.PMF(p, theo_n, (int)k));

                    Ps[i] = cp;
                }
                
                for (int i = 0; i < modPos.Count; i++)
                {
                    siteProbability.Add(modPos[i], Ps[i]/Ps.Sum());
                }

                //If the best and the second best localization score are the same
                var x = bestScores[1].OrderByDescending(p => p).ToList();
                if (x[0] > x[1] && siteProbability.ElementAt(bestModIndex).Value >= 0.75)
                {
                    localizationLevel = LocalizationLevel.Level1;
                }
            }
            else
            {
                siteProbability.Add(modPos[0], 1); //If there is only one N-Glycosite, the site probability for the site will be 1. 
            }

            var psmGlyco = new GlycoSpectralMatch(bestPeptides[bestModIndex], 0, bestScores[0][bestModIndex], scanIndex, theScan, CommonParameters, bestMatchedIons[bestModIndex]);
            psmGlyco.ModPos = modPos;
            psmGlyco.ChildMatchedFragmentIons = bestChildMathcedIons[bestModIndex];
            psmGlyco.NGlycan = new List<Glycan> { NGlycans[iDLow] };
            psmGlyco.GlycanScore = bestScores[2][bestModIndex];
            psmGlyco.DiagnosticIonScore = bestScores[3][bestModIndex];
            psmGlyco.PeptideScore = psmGlyco.Score - psmGlyco.GlycanScore;
            psmGlyco.Rank = ind;
            psmGlyco.LocalizationLevel = localizationLevel;
            psmGlyco.NGlycoSiteSpeciLocalProb = siteProbability;

            psmGlyco.OxoniumIonIntensity = oxoniumIonIntensities;

            var product_nn = GlycoPeptides.GetIndicatorYIon(bestPeptides[bestModIndex].MonoisotopicMass - psmGlyco.NGlycan.First().Mass/1E5, "NN");
            psmGlyco.PepNN = GlycoPeptides.MatchIndicatorYIon(theScan, product_nn, CommonParameters);

            var product_nh = GlycoPeptides.GetIndicatorYIon(bestPeptides[bestModIndex].MonoisotopicMass - psmGlyco.NGlycan.First().Mass / 1E5, "NH");
            psmGlyco.PepNH = GlycoPeptides.MatchIndicatorYIon(theScan, product_nh, CommonParameters);

            return psmGlyco;
        }

        private void FindNGlycan(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications theScanBestPeptide, int ind, double possibleGlycanMassLow, double[] oxoniumIonIntensities, ref List<GlycoSpectralMatch> possibleMatches)
        {
            List<int> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxt", "Nxs" });
            if (modPos.Count < 1)
            {
                return;
            }

            int iDLow = GlycoPeptides.BinarySearchGetIndex(NGlycans.Select(p => (double)p.Mass / 1E5).ToArray(), possibleGlycanMassLow);

            while (iDLow < NGlycans.Length && PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + (double)NGlycans[iDLow].Mass / 1E5))
            {
                if (!GlycoPeptides.NGlyOxoniumIonsAnalysis(oxoniumIonIntensities, NGlycans[iDLow]))
                {
                    iDLow++;
                    continue;
                }

                PeptideWithSetModifications[] bestPeptides = new PeptideWithSetModifications[modPos.Count];
                List<MatchedFragmentIon>[] bestMatchedIons = new List<MatchedFragmentIon>[modPos.Count];
                Dictionary<int, List<MatchedFragmentIon>>[] bestChildMathcedIons = new Dictionary<int, List<MatchedFragmentIon>>[modPos.Count];
                double[][] bestScores = new double[4][];  //allScores, localizationScores,YScores, OxoScores
                for (int i = 0; i < 4; i++)
                {
                    bestScores[i] = new double[modPos.Count];
                }
                
                int theo_n = 0; // For site specific probability calculation.

                for (int i = 0; i < modPos.Count; i++)
                {
                    var possibleSite = modPos[i];

                    var testPeptide = GlycoPeptides.GenerateNGlycopeptide(possibleSite, theScanBestPeptide, NGlycans[iDLow]);

                    List<Product> theoreticalProducts = GlycoPeptides.NGlyGetTheoreticalFragments(CommonParameters.DissociationType, theScanBestPeptide, testPeptide, NGlycans[iDLow]);

                    theo_n += theoreticalProducts.Where(p => p.ProductType != ProductType.M && p.ProductType != ProductType.D).Count();

                    //TO DO: the current MatchFragmentIons only match one charge states. For NGlycopeptide, the Y-ions always contain more than one charge state.
                    //var matchedIons = GlycoPeptides.GlyMatchOriginFragmentIons(theScan, theoreticalProducts, CommonParameters);
                    var matchedIons = MatchFragmentIons(theScan, theoreticalProducts, CommonParameters);


                    var matchedLocalizationIons = matchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M && p.NeutralTheoreticalProduct.ProductType != ProductType.D).ToList(); //TO DO: Select matched localzation ions
                    double localizationScore = CalculatePeptideScore(theScan.TheScan, matchedLocalizationIons);

                    var matchedYIons = matchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M ).ToList(); 
                    double YScore = CalculatePeptideScore(theScan.TheScan, matchedYIons);

                    double score = localizationScore + YScore;

                    var matchedOxoIons = matchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList();
                    double OxoScore = CalculatePeptideScore(theScan.TheScan, matchedOxoIons);

                    var allMatchedChildIons = new Dictionary<int, List<MatchedFragmentIon>>();
                    foreach (var childScan in theScan.ChildScans)
                    {
                        List<Product> theoreticalChildProducts = GlycoPeptides.NGlyGetTheoreticalFragments(CommonParameters.MS2ChildScanDissociationType, theScanBestPeptide, testPeptide, NGlycans[iDLow]);

                        theo_n += theoreticalChildProducts.Where(p => p.ProductType != ProductType.M && p.ProductType != ProductType.D).Count();

                        //var matchedChildIons = GlycoPeptides.GlyMatchOriginFragmentIons(childScan, theoreticalChildProducts, CommonParameters);
                        var matchedChildIons = MatchFragmentIons(childScan, theoreticalChildProducts, CommonParameters);
                     
                        allMatchedChildIons.Add(childScan.OneBasedScanNumber, matchedChildIons);

                        if (matchedChildIons == null)
                        {
                            continue;
                        }

                        var childMatchedLocalizationIons = matchedChildIons.Where(p=>p.NeutralTheoreticalProduct.ProductType!= ProductType.M && p.NeutralTheoreticalProduct.ProductType != ProductType.D).ToList(); //TO DO: Select child matched localzation ions
                        double childLocalizationScore = CalculatePeptideScore(childScan.TheScan, childMatchedLocalizationIons);


                        var childMatchYIons = matchedChildIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList();
                        double childYScore = CalculatePeptideScore(childScan.TheScan, childMatchYIons);

                        double childScore = childLocalizationScore + childYScore;

                        var childMatchedOxoIons = matchedChildIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList();
                        double childOxoScore = CalculatePeptideScore(childScan.TheScan, childMatchedOxoIons);

                        //TO THINK:may think a different way to use childScore
                        score += childScore;
                        localizationScore += childLocalizationScore;
                        YScore += childYScore;
                        OxoScore += childOxoScore;
                    }

                    bestPeptides[i] = testPeptide;
                    bestMatchedIons[i] = matchedIons;
                    bestChildMathcedIons[i] = allMatchedChildIons;
                    bestScores[0][i] = score;
                    bestScores[1][i] = localizationScore;
                    bestScores[2][i] = YScore;
                    bestScores[3][i] = OxoScore;

                }

                if (bestScores[1].Max() > scoreCutOff && bestScores[2].Max() > 0)
                {
                    var psmGlyco = CreateNGsm(theScan, oxoniumIonIntensities, modPos, bestPeptides, bestMatchedIons, bestChildMathcedIons, bestScores, theo_n, scanIndex, iDLow, ind);

                    possibleMatches.Add(psmGlyco);
                }

                iDLow++;
            }
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
                    ////The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);

                    if (OxoniumIonFilter)
                    {
                        //For the HCD-EThcD type of data, the EThcD spectrum also contains Oxonium Ion. However, it is likely we can ignore these Oxonium ions.
                        //The code here is not well designed. 
                        //if (GlycoPeptides.DissociationTypeContainHCD(CommonParameters.MS2ChildScanDissociationType))
                        //{
                        //    foreach (var c in theScan.ChildScans)
                        //    {
                        //        var _childoxo =  GlycoPeptides.ScanOxoniumIonFilter(c, ProductSearchMode, CommonParameters.MS2ChildScanDissociationType);
                        //        for (int i = 0; i < oxoniumIonIntensities.Length; i++)
                        //        {
                        //            oxoniumIonIntensities[i] += _childoxo[i];
                        //        }
                        //    }
                        //}

                        if (oxoniumIonIntensities[OxoniumIon204Index] == 0)
                        {
                            continue;
                        }
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
        private List<GlycoSpectralMatch> FindOGlycopeptideHashLocal(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
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

                    if (possibleGlycanMassHigh < GlycanBox.OGlycanBoxes.First().Mass || possibleGlycanMassLow > GlycanBox.OGlycanBoxes.Last().Mass)
                    {
                        continue;
                    }

                    //Filter by OxoniumIon
                    ////The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);

                    if (OxoniumIonFilter)
                    {
                        //Check FindNGlycopeptide for ChildScan Oxonium consideration.
                        if (oxoniumIonIntensities[OxoniumIon204Index] == 0)
                        {
                            continue;
                        }
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
                    ////The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);

                    if (OxoniumIonFilter)
                    {
                        //Check FindNGlycopeptide for ChildScan Oxonium consideration.
                        if (oxoniumIonIntensities[OxoniumIon204Index] == 0)
                        {
                            continue;
                        }
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

    }
}