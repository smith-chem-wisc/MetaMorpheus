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
        private readonly bool MixedGlycanAllowed;
        private readonly int TopN;
        private readonly int _maxOGlycanNum;
        private readonly int _maxNGlycanNum;
        private readonly bool OxoniumIonFilter; //To filt Oxonium Ion before searching a spectrum as glycopeptides. If we filter spectrum, it must contain oxonium ions such as 204 (HexNAc). 
        private readonly bool IndexChildScan; //Whether to index the child scan or not is a question.
        private readonly string _oglycanDatabase;
        private readonly string _nglycanDatabase;

        private readonly Tolerance PrecusorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;

        private readonly List<int>[] SecondFragmentIndex;

        public GlycoSearchEngine(List<GlycoSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
             string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, bool mixedGlycanAllowed, int glycoSearchTopNum, int maxOGlycanNum, int maxNGlycanNum, bool oxoniumIonFilter, bool indexChildScan, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, fileSpecificParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.GlycoSearchType = glycoSearchType;
            this.MixedGlycanAllowed = mixedGlycanAllowed;
            this.TopN = glycoSearchTopNum;
            this._maxOGlycanNum = maxOGlycanNum;
            this._maxNGlycanNum = maxNGlycanNum;
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
                GlycanBox.GlobalOGlycanMods = GlycanBox.BuildGlobalOGlycanMods(GlycanBox.GlobalOGlycans).ToArray();
                GlycanBox.OGlycanBoxes = GlycanBox.BuildGlycanBoxes(_maxOGlycanNum, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).OrderBy(p => p.Mass).ToArray();
            }
            else if (glycoSearchType == GlycoSearchType.NGlycanSearch)
            {
                GlycanBox.GlobalNGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false).OrderBy(p => p.Mass).ToArray();
                GlycanBox.GlobalNGlycanMods = GlycanBox.BuildGlobalNGlycanMods(GlycanBox.GlobalNGlycans).ToArray();
                GlycanBox.NGlycanBoxes = GlycanBox.BuildGlycanBoxes(_maxNGlycanNum, GlycanBox.GlobalNGlycans, GlycanBox.GlobalNGlycanMods).OrderBy(p => p.Mass).ToArray();
                //TO THINK: Glycan Decoy database.

            }
            else if (glycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.GlobalOGlycanMods = GlycanBox.BuildGlobalOGlycanMods(GlycanBox.GlobalOGlycans).ToArray();
                GlycanBox.OGlycanBoxes = GlycanBox.BuildGlycanBoxes(_maxOGlycanNum, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).OrderBy(p => p.Mass).ToArray();

                GlycanBox.GlobalNGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false).ToArray();
                GlycanBox.GlobalNGlycanMods = GlycanBox.BuildGlobalNGlycanMods(GlycanBox.GlobalNGlycans).ToArray();
                GlycanBox.NGlycanBoxes = GlycanBox.BuildGlycanBoxes(_maxNGlycanNum, GlycanBox.GlobalNGlycans, GlycanBox.GlobalNGlycanMods).OrderBy(p => p.Mass).ToArray();

                if (MixedGlycanAllowed)
                {
                    GlycanBox.GlobalMixedGlycans = GlycanBox.GlobalOGlycans.Concat(GlycanBox.GlobalNGlycans).ToArray();
                    GlycanBox.GlobalMixedGlycanMods = GlycanBox.GlobalOGlycanMods.Concat(GlycanBox.GlobalNGlycanMods).ToArray();
                    GlycanBox.MixedModBoxes = GlycanBox.BuildMixedGlycanBoxes(GlycanBox.GlobalMixedGlycans, GlycanBox.GlobalMixedGlycanMods, _maxOGlycanNum, _maxNGlycanNum, GlycanBox.GlobalOGlycans.Length, GlycanBox.GlobalNGlycans.Length).OrderBy(p => p.Mass).ToArray();
                }
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

                        gsms = FindGlycopeptide(scan, idsOfPeptidesTopN, scanIndex, (int)byteScoreCutoff);

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
            //TO DO: keep top 10 candidates in future for machine learning.
            double preScore = 0;
            int gsmsCount = 1;
            string preString = "";

            foreach (var gsm in gsms.Where(p => p != null).OrderByDescending(p => p.Score).ThenBy(c => c.FullSequence))
            {
                if (gsmsCount <= 1) //keep top 1 candidates currently.
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
        /// Plan to be Deprecated. (The function is in IndexingEngine, which is now used in Crosslink and Glyco Search)
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

        private GlycoSpectralMatch CreateGsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, 
            Route localization, double[] oxoniumIonIntensities, List<LocalizationGraph> localizationGraphs, 
            int[] n_modpos, GlycoType glycanType, Glycan[] globalglycans, Modification[] globalmods)
        {
            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(localization, peptide, globalmods);

            Glycan[] glycans = new Glycan[localizationGraphs.First().ModBox.ModCount];
            for (int i = 0; i < localizationGraphs.First().ModBox.ModCount; i++)
            {
                glycans[i] = globalglycans[((GlycanBox)localizationGraphs.First().ModBox).ModIds[i]];
            }
            List<int> npos = localization.Mods.Select(p => p.Item1).ToArray().Intersect(n_modpos).ToList();

            var fragmentsForEachGlycoPeptide = GlycoPeptides.GlyGetTheoreticalFragments(glycanType, CommonParameters.DissociationType, peptide, peptideWithMod, npos, glycans);

            var matchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);

            double score = CalculatePeptideScore(theScan.TheScan, matchedIons);

            var DiagnosticIonScore = CalculatePeptideScore(theScan.TheScan, matchedIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());

            var GlycanScore = CalculatePeptideScore(theScan.TheScan, matchedIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList());

            var PeptideScore = score - DiagnosticIonScore;

            var p = theScan.TheScan.MassSpectrum.Size * CommonParameters.ProductMassTolerance.GetRange(1000).Width / theScan.TheScan.MassSpectrum.Range.Width;

            int n = fragmentsForEachGlycoPeptide.Where(v => v.ProductType != ProductType.D && v.ProductType != ProductType.M).Count();

            var allMatchedChildIons = new Dictionary<int, List<MatchedFragmentIon>>();

            foreach (var childScan in theScan.ChildScans)
            {
                //People always use CID with low res. This special code works for Nic Scott's data. (High-HCD, High-EThcD, Low-CID, High-secHCD)
                if (childScan.TheScan.DissociationType == DissociationType.CID)
                {
                    continue;
                }
                var childFragments = GlycoPeptides.GlyGetTheoreticalFragments(glycanType, CommonParameters.MS2ChildScanDissociationType, peptide, peptideWithMod, npos, glycans);

                var matchedChildIons = MatchFragmentIons(childScan, childFragments, CommonParameters);

                n += childFragments.Where(v => v.ProductType != ProductType.D && v.ProductType != ProductType.M).Count();

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

            psmGlyco.GlycanType = glycanType;

            //TO THINK: Is this p and n properly calculated for different data type.
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

        private void GetModPosMotif(GlycoType glycoType, int[] n_modPos, int[] o_modPos, out int[] modPos, out string[] modMotifs)
        {
            int ip = 0;
            if (glycoType == GlycoType.OGlycoPep)
            {
                modPos = o_modPos;
                modMotifs = new string[o_modPos.Length];
               
                foreach (var o in o_modPos)
                {
                    modPos[ip] = o;
                    modMotifs[ip] = "S/T";
                    ip++;
                }
            }
            else if (glycoType == GlycoType.NGlycoPep)
            {
                modPos = n_modPos;
                modMotifs = new string[n_modPos.Length];

                foreach (var n in n_modPos)
                {
                    modPos[ip] = n;
                    modMotifs[ip] = "Nxs/t";
                    ip++;
                }
            }
            else
            {
                modPos = new int[n_modPos.Length + o_modPos.Length];
                modMotifs = new string[n_modPos.Length + o_modPos.Length];
                foreach (var n in n_modPos)
                {
                    modPos[ip] = n;
                    modMotifs[ip] = "Nxs/t";
                    ip++;
                }

                foreach (var o in o_modPos)
                {
                    modPos[ip] = o;
                    modMotifs[ip] = "S/T";
                    ip++;
                }
                Array.Sort(modPos, modMotifs);
            }
        }
        
        private void FindGlycan(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications theScanBestPeptide, 
            int ind, double possibleGlycanMassLow, double possibleGlycanMassHigh, double[] oxoniumIonIntensities, 
            GlycoType glycanType, Glycan[] globalglycans, Modification[] globalmods, GlycanBox[] globalBoxes, ref List<GlycoSpectralMatch> possibleMatches)
        {
            int[] n_modPos = GlycoPeptides.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxt", "Nxs" }).OrderBy(p => p).ToArray();

            if ((glycanType == GlycoType.NGlycoPep || glycanType == GlycoType.MixedGlycoPep) && n_modPos.Length < 1)
            {
                return;
            }

            int[] o_modPos = GlycoPeptides.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" }).OrderBy(p => p).ToArray();

            if ((glycanType == GlycoType.OGlycoPep || glycanType == GlycoType.MixedGlycoPep) && o_modPos.Length < 1)
            {
                return;
            }

            int[] modPos;
            string[] modMotifs;
            GetModPosMotif(glycanType, n_modPos, o_modPos, out modPos, out modMotifs);

            int iDLow = GlycoPeptides.BinarySearchGetIndex(globalBoxes.Select(p=>p.Mass).ToArray(), possibleGlycanMassLow);

            List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();

            double bestLocalizedScore = scoreCutOff;

            //while (iDLow < globalBoxes.Length && PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + globalBoxes[iDLow].Mass))
            while (iDLow < globalBoxes.Length && globalBoxes[iDLow].Mass <= possibleGlycanMassHigh)
            {
                if (glycanType == GlycoType.OGlycoPep || glycanType == GlycoType.NGlycoPep)
                {
                    if (globalBoxes[iDLow].ModCount > modPos.Length)
                    {
                        iDLow++;
                        continue;
                    }
                }
                else
                {
                    if (globalBoxes[iDLow].NGlycanCount == 0 || globalBoxes[iDLow].OGlycanCount == 0 ||
                        globalBoxes[iDLow].NGlycanCount > n_modPos.Length || globalBoxes[iDLow].OGlycanCount > o_modPos.Length)
                    {
                        iDLow++;
                        continue;
                    }
                }


                //TO DO: check the NGlyOxoniumIonsAnalysis for N-Glycan only or for both.
                if (!GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, globalBoxes[iDLow]))
                {
                    iDLow++;
                    continue;
                }

                //Peptide HCD fragments for graph localization.
                List<Product> hcdProducts = new List<Product>();
                theScanBestPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, hcdProducts);
                //Peptide ETD fragments for graph localization.
                List<Product> etdProducts = new List<Product>();
                theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, etdProducts);

                LocalizationGraph localizationGraph = new LocalizationGraph(modPos, modMotifs, globalBoxes[iDLow], globalBoxes[iDLow].ChildGlycanBoxes, iDLow);

                List<Product> mainProducts = new List<Product>();
                if (GlycoPeptides.DissociationTypeContainHCD(CommonParameters.DissociationType))
                {
                    mainProducts.AddRange(hcdProducts);
                }
                if (GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType))
                {
                    mainProducts.AddRange(etdProducts.Where(p=>p.ProductType != ProductType.y));
                }
                LocalizationGraph.LocalizeMod(localizationGraph, theScan, CommonParameters.ProductMassTolerance, mainProducts, GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
                
                if (theScan.ChildScans.Count > 0)
                {
                    foreach (var childScan in theScan.ChildScans)
                    {
                        List<Product> childProducts = new List<Product>();

                        if (GlycoPeptides.DissociationTypeContainHCD(CommonParameters.MS2ChildScanDissociationType))
                        {
                            childProducts.AddRange(hcdProducts);
                        }

                        if (GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType))
                        {
                            childProducts.AddRange(etdProducts.Where(p => p.ProductType != ProductType.y));
                        }

                        //run the same graph multiple times!
                        LocalizationGraph.LocalizeMod(localizationGraph, childScan, CommonParameters.ProductMassTolerance, childProducts, GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
                    }                   
                }

                double currentLocalizationScore = localizationGraph.TotalScore;
                if (currentLocalizationScore > bestLocalizedScore)
                {
                    bestLocalizedScore = currentLocalizationScore;
                    localizationGraphs.Clear();
                    localizationGraphs.Add(localizationGraph);
                }
                else if (currentLocalizationScore == bestLocalizedScore)
                {
                    localizationGraphs.Add(localizationGraph);
                }        

                iDLow++;
            }

            if (localizationGraphs.Count > 0)
            {
                var firstPath = LocalizationGraph.GetFirstPath(localizationGraphs[0].array, localizationGraphs[0].ChildModBoxes);
                var localizationCandidate = LocalizationGraph.GetLocalizedPath(localizationGraphs[0], firstPath);

                var psmGlyco = CreateGsm(theScan, scanIndex, ind, theScanBestPeptide, localizationCandidate, oxoniumIonIntensities, localizationGraphs, n_modPos, glycanType, globalglycans, globalmods);

                if (psmGlyco.Score > scoreCutOff)
                {
                    possibleMatches.Add(psmGlyco);
                }
            }

            return;
        }
        
        private List<GlycoSpectralMatch> FindGlycopeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
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

                    //Filter by glycanBoxes mass difference.
                    var possibleGlycanMassLow = PrecusorSearchMode.GetMinimumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    var possibleGlycanMassHigh = PrecusorSearchMode.GetMaximumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    if (GlycoSearchType == GlycoSearchType.OGlycanSearch || GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
                    {
                        if (possibleGlycanMassHigh < GlycanBox.OGlycanBoxes.First().Mass || possibleGlycanMassLow > GlycanBox.OGlycanBoxes.Last().Mass)
                        {
                            continue;
                        }

                        //Find O-Glycan use the original function.
                        FindOGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, possibleGlycanMassLow, oxoniumIonIntensities, ref possibleMatches);
                        //FindGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, possibleGlycanMassLow, possibleGlycanMassHigh, oxoniumIonIntensities,
                        //    GlycoType.OGlycoPep, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods, GlycanBox.OGlycanBoxes, ref possibleMatches);
                    }

                    if (GlycoSearchType == GlycoSearchType.NGlycanSearch || GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
                    {
                        if (possibleGlycanMassHigh < GlycanBox.NGlycanBoxes.First().Mass || possibleGlycanMassLow > GlycanBox.NGlycanBoxes.Last().Mass)
                        {
                            continue;
                        }

                        //Find N-Glycan 
                        FindGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, possibleGlycanMassLow, possibleGlycanMassHigh, oxoniumIonIntensities, 
                            GlycoType.NGlycoPep, GlycanBox.GlobalNGlycans, GlycanBox.GlobalNGlycanMods, GlycanBox.NGlycanBoxes, ref possibleMatches);
                    }

                    if (GlycoSearchType == GlycoSearchType.N_O_GlycanSearch && MixedGlycanAllowed)
                    {
                        if (possibleGlycanMassHigh < GlycanBox.MixedModBoxes.First().Mass || possibleGlycanMassLow > GlycanBox.MixedModBoxes.Last().Mass)
                        {
                            continue;
                        }

                        FindGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, ind, possibleGlycanMassLow, possibleGlycanMassHigh, oxoniumIonIntensities,
                            GlycoType.MixedGlycoPep, GlycanBox.GlobalMixedGlycans, GlycanBox.GlobalMixedGlycanMods, GlycanBox.MixedModBoxes, ref possibleMatches);                        
                    }
                }

                if (possibleMatches.Count != 0)
                {
                    possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
                }
            }

            return possibleMatches;
        }


        #region Original O-Pair functions, plan to deprecated.
        private GlycoSpectralMatch CreateOGsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, Route localization, double[] oxoniumIonIntensities, List<LocalizationGraph> localizationGraphs)
        {
            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(localization, peptide, GlycanBox.GlobalOGlycanMods);

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

            psmGlyco.GlycanType = GlycoType.OGlycoPep;
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

            int[] modPos = GlycoPeptides.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" }).OrderBy(p => p).ToArray();

            //Localization for O-glycopeptides only works on ETD related dissociationtype
            //No localization can be done with MS2-HCD spectrum
            bool is_HCD_only_data = !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType);
            if (is_HCD_only_data)
            {
                List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();

                while (iDLow < GlycanBox.OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBox.OGlycanBoxes[iDLow].Mass)))
                {
                    if (modPos.Length >= GlycanBox.OGlycanBoxes[iDLow].ModCount && GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, GlycanBox.OGlycanBoxes[iDLow]))
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
                    if (modPos.Length >= GlycanBox.OGlycanBoxes[iDLow].ModCount && GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, GlycanBox.OGlycanBoxes[iDLow]))
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

        #endregion

    }
}