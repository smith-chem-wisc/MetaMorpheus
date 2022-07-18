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

        protected readonly List<List<(double, int, double)>> Precursorss;
        protected readonly List<(int, int, int)>[] Candidates;

        private GlycoSearchType GlycoSearchType;
        private GlycoScoreType GlycoScoreType;
        private readonly bool MixedGlycanAllowed;
        private readonly int TopN;
        private readonly int _maxOGlycanNum;
        private readonly int _maxNGlycanNum;
        private readonly bool OxoniumIonFilter; //To filt Oxonium Ion before searching a spectrum as glycopeptides. If we filter spectrum, it must contain oxonium ions such as 204 (HexNAc). 
        private readonly bool FilterYcore = false;
        private readonly bool IndexChildScan; //Whether to index the child scan or not is a question.
        private readonly string _oglycanDatabase;
        private readonly string _nglycanDatabase;

        private readonly Tolerance PrecusorSearchMode;
        private readonly MassDiffAcceptor PrecusorSearchModeWithShift;
        private readonly MassDiffAcceptor ProductSearchMode;

        private readonly List<int>[] SecondFragmentIndex;

        public GlycoSearchEngine(List<GlycoSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters, 
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
             string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, GlycoScoreType glycoScoreType, bool mixedGlycanAllowed, int glycoSearchTopNum, 
             int maxOGlycanNum, int maxNGlycanNum, bool oxoniumIonFilter, bool indexChildScan, List<string> nestedIds,
            List<(int, int, int)>[] candidates, List<List<(double, int, double)>> precursorss)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, fileSpecificParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalCsms = globalCsms;
            this.GlycoSearchType = glycoSearchType;
            this.GlycoScoreType = glycoScoreType;
            this.MixedGlycanAllowed = mixedGlycanAllowed;
            this.TopN = glycoSearchTopNum;
            this._maxOGlycanNum = maxOGlycanNum;
            this._maxNGlycanNum = maxNGlycanNum;
            this.OxoniumIonFilter = oxoniumIonFilter;
            this.IndexChildScan = indexChildScan;
            this._oglycanDatabase = oglycanDatabase;
            this._nglycanDatabase = nglycanDatabase;

            SecondFragmentIndex = secondFragmentIndex;

            this.Candidates = candidates;
            this.Precursorss = precursorss;

            PrecusorSearchMode = commonParameters.PrecursorMassTolerance;
            PrecusorSearchModeWithShift = new DotMassDiffAcceptor(
                        "PlusOrMinus2Da",
                        new List<double>
                        {
                            -2 * Constants.C13MinusC12,
                            -1 * Constants.C13MinusC12,
                            0,
                            1 * Constants.C13MinusC12,
                            2 * Constants.C13MinusC12
                        },
                        commonParameters.PrecursorMassTolerance);
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


        public void FirstRoundSearch()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing glyco search ion index... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;
            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;

            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (scanIndex) =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();
                byte[] secondScoringTable = new byte[PeptideIndex.Count];
                List<int> childIdsOfPeptidesPossiblyObserved = new List<int>();

                byte scoreAtTopN = 0;
                int peptideCount = 0;

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

                    //Limit the high bound limitation, here assume it is possible to has max 3 Da shift. This allows for correcting precursor in the future.
                    var high_bound_limitation = scan.PrecursorMass + 5;

                    // first-pass scoring
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);

                    //child scan first - pass scoring
                    if (IndexChildScan && SecondFragmentIndex == null && scan.ChildScans != null && CommonParameters.MS2ChildScanDissociationType != DissociationType.Unknown && CommonParameters.MS2ChildScanDissociationType != DissociationType.LowCID)
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

                            if (Candidates[scanIndex] == null)
                            {
                                Candidates[scanIndex] = new List<(int, int, int)>();
                            }

                            Candidates[scanIndex].Add((CurrentPartition - 1, id, scoringTable[id]));
                        }

                        //Only keep TopN candidates.
                        if (CommonParameters.TotalPartitions > 1 && Candidates[scanIndex].Count() > TopN)
                        {
                            Candidates[scanIndex].Sort((x, y) => x.Item3.CompareTo(y.Item3));
                            Candidates[scanIndex].Reverse();
                            int minScore = Candidates[scanIndex][TopN - 1].Item3;
                            var id = Candidates[scanIndex].Count - 1;

                            var keepRemove = true;
                            while (id >= TopN && keepRemove)
                            {
                                if (Candidates[scanIndex][id].Item3 < minScore)
                                {
                                    Candidates[scanIndex].RemoveAt(id);
                                }
                                else
                                {
                                    keepRemove = false;
                                }
                                id--;
                            }
                        }

                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing glyco first round search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
                    }
                }
            });
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing glyco search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (scanIndex) =>
            {

                for (; scanIndex < ListOfSortedMs2Scans.Length; scanIndex += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    if (Candidates[scanIndex] == null)
                    {
                        continue;
                    }

                    var _candidates = Candidates[scanIndex].Where(p => p.Item1 == CurrentPartition - 1).Select(p => p.Item2).ToList();

                    var scan = ListOfSortedMs2Scans[scanIndex];

                    var precursors = Precursorss[scanIndex];

                    List<GlycoSpectralMatch> gsms;

                    gsms = _FindGlycopeptide(scan, precursors, _candidates, scanIndex, (int)byteScoreCutoff);

                    //gsms = FindOGlycopeptideHashLocal(scan, _candidates, scanIndex, (int)byteScoreCutoff);

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

        #region NOSearch. 

        private List<GlycoSpectralMatch> _FindGlycopeptide(Ms2ScanWithSpecificMass theScan, List<(double, int, double)> precursors, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
        {
            List<GlycoSpectralMatch> gsms = new List<GlycoSpectralMatch>();

            var oxoniumIonIntensities = new double[Glycan.AllOxoniumIons.Length];

            if (GlycoScoreType == GlycoScoreType.XcorrScore)
            {
                theScan.kojakSparseArray = GlycoXCorr.kojakXCorr(theScan, 0.01, 1);
                foreach (var childScan in theScan.ChildScans)
                {
                    if (childScan.TheScan.DissociationType == DissociationType.CID)
                    {
                        continue;
                    }
                    childScan.kojakSparseArray = GlycoXCorr.kojakXCorr(childScan, 0.01, 1);
                }
            }


            var gsmcds = _FindGsmCandidates(theScan, precursors, idsOfPeptidesPossiblyObserved, scanIndex, scoreCutOff, ref oxoniumIonIntensities);

            if (gsmcds.Count > 0)
            {
                //We only keep first gsm in the current version. The code could be expanded to consider keep more in the result for percolator in the future.

                var globalmods = GlycanBox.GlobalOGlycanMods;
                var globalBoxes = GlycanBox.OGlycanBoxes;
                var globalglycans = GlycanBox.GlobalOGlycans;

                if (gsmcds.First().GType == GlycoType.NGlycoPep)
                {
                    globalmods = GlycanBox.GlobalNGlycanMods;
                    globalBoxes = GlycanBox.NGlycanBoxes;
                    globalglycans = GlycanBox.GlobalNGlycans;
                }
                else if (gsmcds.First().GType == GlycoType.MixedGlycoPep)
                {
                    globalmods = GlycanBox.GlobalMixedGlycanMods;
                    globalBoxes = GlycanBox.MixedModBoxes;
                    globalglycans = GlycanBox.GlobalMixedGlycans;
                }

                if (gsmcds.First().GType == GlycoType.SinglePep)
                {
                    var gsm = _CreateSingleGsm(gsmcds.First(), theScan, scanIndex, scoreCutOff);
                    gsms.Add(gsm);
                }
                else
                {
                    var gsm = _CreateGsm(gsmcds.First(), theScan, scanIndex, oxoniumIonIntensities, scoreCutOff, globalBoxes, globalglycans, globalmods);
                    gsms.Add(gsm);
                }

            }

            return gsms;
        }

        private List<GlycoSpectraMatchCandidate> _FindGsmCandidates(Ms2ScanWithSpecificMass theScan, List<(double, int, double)> precursors, List<int> idsOfPeptidesPossiblyObserved, 
            int scanIndex, int scoreCutOff, ref double[] oxoniumIonIntensities)
        {
            List<GlycoSpectraMatchCandidate> gsmCandidates = new List<GlycoSpectraMatchCandidate>();
            

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                foreach (var pre in precursors)
                {
                    var pepId = idsOfPeptidesPossiblyObserved[ind];
                    var theScanBestPeptide = PeptideIndex[pepId];

                    if (PrecusorSearchModeWithShift.Accepts(pre.Item1, theScanBestPeptide.MonoisotopicMass)>=0)
                    {
                        _FindSingle(theScan, scanIndex, scoreCutOff, pepId, theScanBestPeptide, ind, ref gsmCandidates);

                    }

                    //The smallest glycan modification generally is HexNac 203. //Filter out unknow non-glycan modifications.
                    if (pre.Item1 - theScanBestPeptide.MonoisotopicMass < 150)
                    {
                        continue;
                    }

                    List<int> n_modPos = new List<int>();
                    List<int> o_modPos = new List<int>();
                    if (GlycoSearchType == GlycoSearchType.NGlycanSearch )
                    {
                        n_modPos = GlycoPeptides.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxt", "Nxs" }).OrderBy(p => p).ToList();
                        if (n_modPos.Count < 1)
                        {
                            continue;
                        }
                    }
                    else if (GlycoSearchType == GlycoSearchType.OGlycanSearch)
                    {
                        o_modPos = GlycoPeptides.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" }).OrderBy(p => p).ToList();
                        if (o_modPos.Count < 1)
                        {
                            continue;
                        }
                    }
                    else if (GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
                    {
                        n_modPos = GlycoPeptides.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxt", "Nxs" }).OrderBy(p => p).ToList();
                        o_modPos = GlycoPeptides.GetPossibleModSites(theScanBestPeptide, new string[] { "S", "T" }).OrderBy(p => p).ToList();
                        if (n_modPos.Count < 1 && o_modPos.Count < 1)
                        {
                            continue;
                        }
                    }

                    var oxo204check = false;

                    //Filter by glycanBoxes mass difference. Allow largest 2.052 mass deconvolution error.
                    var possibleGlycanMassLow = PrecusorSearchMode.GetMinimumValue(pre.Item1) - theScanBestPeptide.MonoisotopicMass - 2.052;

                    var possibleGlycanMassHigh = PrecusorSearchMode.GetMaximumValue(pre.Item1) - theScanBestPeptide.MonoisotopicMass + 2.052;

                    if (GlycoSearchType == GlycoSearchType.OGlycanSearch || GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
                    {
                        if (possibleGlycanMassHigh < GlycanBox.OGlycanBoxes.First().Mass || possibleGlycanMassLow > GlycanBox.OGlycanBoxes.Last().Mass)
                        {
                            continue;
                        }

                        //Filter by OxoniumIon //if oxo204check is ture. It will either continue of pass, the following OxoniumIonFilter must pass.
                        if (OxoniumIonFilter && !oxo204check)
                        {
                            ////The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                            oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);
                            oxo204check = true;
                        }

                        if (oxo204check && oxoniumIonIntensities[OxoniumIon204Index] == 0)
                        {
                            continue;
                        }

                        //Filter Ycore
                        if (FilterYcore)
                        {
                            //OGlycopeptide better to have at least one ycore ion.
                            var ycore_count = GlycoPeptides.YCoreIonsFilter(theScan, theScanBestPeptide, GlycoType.OGlycoPep, ProductSearchMode);
                            if (ycore_count < 1)
                            {
                                continue;
                            }
                        }

                        //Find O-Glycan use the original function.
                        //FindOGlycan(theScan, scanIndex, scoreCutOff, pepId, theScanBestPeptide, ind, possibleGlycanMassLow, oxoniumIonIntensities, ref gsmCandidates);

                        _FindGlycan(theScan, pre.Item1, scanIndex, scoreCutOff, pepId, theScanBestPeptide, n_modPos, o_modPos, ind, possibleGlycanMassLow, possibleGlycanMassHigh, oxoniumIonIntensities,
                            GlycoType.OGlycoPep, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods, GlycanBox.OGlycanBoxes, ref gsmCandidates);
                    }

                    if (GlycoSearchType == GlycoSearchType.NGlycanSearch || GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
                    {
                        if (possibleGlycanMassHigh < GlycanBox.NGlycanBoxes.First().Mass || possibleGlycanMassLow > GlycanBox.NGlycanBoxes.Last().Mass)
                        {
                            continue;
                        }

                        //Filter by OxoniumIon
                        if (OxoniumIonFilter && !oxo204check)
                        {
                            ////The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                            oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);
                            oxo204check = true;
                        }

                        if (oxo204check && oxoniumIonIntensities[OxoniumIon204Index] == 0)
                        {
                            continue;
                        }

                        //Filter Ycore
                        if (FilterYcore)
                        {
                            var ycore_count = GlycoPeptides.YCoreIonsFilter(theScan, theScanBestPeptide, GlycoType.NGlycoPep, ProductSearchMode);
                            if (ycore_count < 2)
                            {
                                continue;
                            }
                        }

                        //Find N-Glycan 
                        _FindGlycan(theScan, pre.Item1, scanIndex, scoreCutOff, pepId, theScanBestPeptide, n_modPos, o_modPos, ind, possibleGlycanMassLow, possibleGlycanMassHigh, oxoniumIonIntensities,
                            GlycoType.NGlycoPep, GlycanBox.GlobalNGlycans, GlycanBox.GlobalNGlycanMods, GlycanBox.NGlycanBoxes, ref gsmCandidates);
                    }

                    if (GlycoSearchType == GlycoSearchType.N_O_GlycanSearch && MixedGlycanAllowed)
                    {
                        if (possibleGlycanMassHigh < GlycanBox.MixedModBoxes.First().Mass || possibleGlycanMassLow > GlycanBox.MixedModBoxes.Last().Mass)
                        {
                            continue;
                        }

                        //Filter by OxoniumIon
                        if (OxoniumIonFilter && !oxo204check)
                        {
                            ////The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                            oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);
                            oxo204check = true;
                        }

                        if (oxo204check && oxoniumIonIntensities[OxoniumIon204Index] == 0)
                        {
                            continue;
                        }

                        //Filter Ycore
                        if (FilterYcore)
                        {
                            var ycore_count = GlycoPeptides.YCoreIonsFilter(theScan, theScanBestPeptide, GlycoType.MixedGlycoPep, ProductSearchMode);
                            if (ycore_count < 2)
                            {
                                continue;
                            }
                        }

                        _FindGlycan(theScan, pre.Item1, scanIndex, scoreCutOff, pepId, theScanBestPeptide, n_modPos, o_modPos, ind, possibleGlycanMassLow, possibleGlycanMassHigh, oxoniumIonIntensities,
                            GlycoType.MixedGlycoPep, GlycanBox.GlobalMixedGlycans, GlycanBox.GlobalMixedGlycanMods, GlycanBox.MixedModBoxes, ref gsmCandidates);

                    }
                }
            }

            if (gsmCandidates.Count != 0)
            {
                //Using which gsmcds for further search.
                if (GlycoScoreType == GlycoScoreType.MMScore)
                {
                    gsmCandidates = gsmCandidates.OrderByDescending(p => p.Score).ToList();

                }
                else if (GlycoScoreType == GlycoScoreType.pGlycoScore)
                {
                    gsmCandidates = gsmCandidates.OrderByDescending(p => p.PScore).ToList();
                }
                else if (GlycoScoreType == GlycoScoreType.XcorrScore)
                {
                    gsmCandidates = gsmCandidates.OrderByDescending(p => p.XcorrScore).ToList();
                }
            }

            return gsmCandidates;
        }

        private void _FindSingle(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, int pepId, PeptideWithSetModifications theScanBestPeptide, int ind, ref List<GlycoSpectraMatchCandidate> gsmCandidates)
        {
            List<Product> products = new List<Product>();
            theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, products);
            var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
            double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);
            double pScore = GlycoPeptides.Calc_pGlycoScore_SinglePep(matchedFragmentIons, products, CommonParameters);

            double xcorrScore = 0;
            if (GlycoScoreType == GlycoScoreType.XcorrScore)
            {
                xcorrScore = GlycoXCorr.KojakScoring(theScan, products, 0.01, 1, theScan.kojakSparseArray);
            }

            var allMatchedChildIons = new Dictionary<int, List<MatchedFragmentIon>>();

            foreach (var childScan in theScan.ChildScans)
            {
                //TO THINK: CID data may should not be excluded. 
                //People always use CID with low res. This special code works for Nick Scott's data. (High-HCD, High-EThcD, Low-CID, High-secHCD)
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
                double child_pScore = GlycoPeptides.Calc_pGlycoScore_SinglePep(matchedChildIons, childFragments, CommonParameters);

                double childXcorrScore = 0;
                if (GlycoScoreType == GlycoScoreType.XcorrScore)
                {
                    childXcorrScore = GlycoXCorr.KojakScoring(childScan, childFragments, 0.01, 1, childScan.kojakSparseArray);
                }
                //TO THINK:may think a different way to use childScore
                
                score += childScore;
                pScore += child_pScore;
                xcorrScore += childXcorrScore;
            }

            if (score > scoreCutOff)
            {
                var psmcd = new GlycoSpectraMatchCandidate(pepId, score, pScore, xcorrScore, GlycoType.SinglePep, ind, null, null, null, null);

                gsmCandidates.Add(psmcd);
            }
        }

        private void _FindGlycan(Ms2ScanWithSpecificMass theScan, double precursorMass, int scanIndex, int scoreCutOff, 
            int pepId, PeptideWithSetModifications theScanBestPeptide, List<int> n_modPos, List<int> o_modPos, int ind, 
            double possibleGlycanMassLow, double possibleGlycanMassHigh, double[] oxoniumIonIntensities,
            GlycoType glycanType, Glycan[] globalglycans, Modification[] globalmods, GlycanBox[] globalBoxes, ref List<GlycoSpectraMatchCandidate> gsmCandidates)
        {

            if (glycanType == GlycoType.NGlycoPep)
            {
                if (n_modPos.Count < 1)
                {
                    return;
                }
            }

            else if (glycanType == GlycoType.OGlycoPep)
            {
                if (o_modPos.Count < 1)
                {
                    return;
                }
            }
            else if (glycanType == GlycoType.MixedGlycoPep)
            {
                if (o_modPos.Count < 1 && n_modPos.Count < 1)
                {
                    return;
                }
            }

            int[] modPos;
            string[] modMotifs;
            GetModPosMotif(glycanType, n_modPos.ToArray(), o_modPos.ToArray(), out modPos, out modMotifs);

            int iDLow = GlycoPeptides.BinarySearchGetIndex(globalBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow);

            var glycanBoxes = new List<int>();
            //while (iDLow < globalBoxes.Length && PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + globalBoxes[iDLow].Mass))
            while (iDLow < globalBoxes.Length && globalBoxes[iDLow].Mass <= possibleGlycanMassHigh && PrecusorSearchModeWithShift.Accepts(precursorMass, theScanBestPeptide.MonoisotopicMass + globalBoxes[iDLow].Mass)>=0)
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
                        globalBoxes[iDLow].NGlycanCount > n_modPos.Count || globalBoxes[iDLow].OGlycanCount > o_modPos.Count)
                    {
                        iDLow++;
                        continue;
                    }
                }

                //TO DO: check the NGlyOxoniumIonsAnalysis for N-Glycan only or for both.
                if (OxoniumIonFilter && !GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, globalBoxes[iDLow]))
                {
                    iDLow++;
                    continue;
                }

                glycanBoxes.Add(iDLow);

                iDLow++;
            }

            if (glycanBoxes.Count() > 0)
            {
                //Here we do not know which glycanBoxes is the one that will generate most Y-ions. As we allow deconvolution error (-2.052, +2.052).
                //This may be OK as for N-glycopeptide, most share the same Ycore ions.
                Glycan[] glycans = new Glycan[globalBoxes[glycanBoxes.First()].ModCount];
                for (int i = 0; i < globalBoxes[glycanBoxes.First()].ModCount; i++)
                {
                    glycans[i] = globalglycans[globalBoxes[glycanBoxes.First()].ModIds[i]];
                }

                var fragments_hcd_peptide = GlycoPeptides.GlyGetPepHCDFragments(glycanType, CommonParameters.DissociationType, theScanBestPeptide, n_modPos, glycans);

                var matchedIons = MatchFragmentIons(theScan, fragments_hcd_peptide, CommonParameters);

                //Filter Y core ions.
                if(!GlycoPeptides.ScanTrimannosylCoreFilter(matchedIons, glycanType))
                {
                    return;
                }

                double score = CalculatePeptideScore(theScan.TheScan, matchedIons);
                double pScore = GlycoPeptides.Calc_pGlycoScore(matchedIons, fragments_hcd_peptide, CommonParameters);

                double xcorrScore = 0; 

                if (GlycoScoreType == GlycoScoreType.XcorrScore)
                {
                    xcorrScore = GlycoXCorr.KojakScoring(theScan, fragments_hcd_peptide, 0.01, 1, theScan.kojakSparseArray);
                }
 


                if (theScan.ChildScans.Count > 0)
                {
                    foreach (var childScan in theScan.ChildScans)
                    {
                        if (childScan.TheScan.DissociationType == DissociationType.CID)
                        {
                            continue;
                        }

                        var childFragments = GlycoPeptides.GlyGetPepHCDFragments(glycanType, CommonParameters.MS2ChildScanDissociationType, theScanBestPeptide, n_modPos, glycans);

                        var matchedChildIons = MatchFragmentIons(childScan, childFragments, CommonParameters);

                        double childScore = CalculatePeptideScore(theScan.TheScan, matchedChildIons);
                        double child_pScore = GlycoPeptides.Calc_pGlycoScore(matchedChildIons, childFragments, CommonParameters);

                        double childXcorrScore = 0;
                        if (GlycoScoreType == GlycoScoreType.XcorrScore)
                        {
                            childXcorrScore = GlycoXCorr.KojakScoring(childScan, childFragments, 0.01, 1, childScan.kojakSparseArray);

                        }

                        score += childScore;
                        pScore += child_pScore;
                        xcorrScore += childXcorrScore;
                    }
                }

                if (score > scoreCutOff)
                {
                    var psmcd = new GlycoSpectraMatchCandidate(pepId, score, pScore, xcorrScore, glycanType, ind, modPos, modMotifs, n_modPos.ToArray(), glycanBoxes);

                    gsmCandidates.Add(psmcd);
                }
            }

            return;
        }


        private GlycoSpectralMatch _CreateSingleGsm(GlycoSpectraMatchCandidate gsmcd, Ms2ScanWithSpecificMass theScan, int scanIndex, double scoreCutOff)
        {
            var theScanBestPeptide = PeptideIndex[gsmcd.PepId];
            List<Product> products = new List<Product>();
            theScanBestPeptide.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, products);
            var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);

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

            }

            double _score = 0;
            if (GlycoScoreType == GlycoScoreType.MMScore)
            {
                _score = gsmcd.Score;
            }
            else if (GlycoScoreType == GlycoScoreType.pGlycoScore)
            {
                _score = gsmcd.PScore;
            }
            else if (GlycoScoreType == GlycoScoreType.XcorrScore)
            {
                _score = gsmcd.XcorrScore;
            }

            var gsmSingle = new GlycoSpectralMatch(theScanBestPeptide, 0, _score, scanIndex, theScan, CommonParameters, matchedFragmentIons);
            gsmSingle.Rank = gsmcd.Rank;
            gsmSingle.MScore = gsmcd.Score;
            gsmSingle.PScore = gsmcd.PScore;
            gsmSingle.XcorrScore = gsmcd.XcorrScore;
            if (allMatchedChildIons.Count() > 0)
            {
                gsmSingle.ChildMatchedFragmentIons = allMatchedChildIons;
            }

            return gsmSingle;

        }


        private GlycoSpectralMatch _CreateGsm(GlycoSpectraMatchCandidate gsmcd, Ms2ScanWithSpecificMass theScan, int scanIndex, 
             double[] oxoniumIonIntensities, int scoreCutOff, GlycanBox[] globalBoxes, Glycan[] globalglycans, Modification[] globalmods)
        {
            var theScanBestPeptide = PeptideIndex[gsmcd.PepId];

            bool is_HCD_only_data = true;
            if (GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType))
            {
                is_HCD_only_data = false;
            }
            else if (theScan.ChildScans.Count > 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType))
            {
                is_HCD_only_data = false;
            }

            if (gsmcd.GType ==GlycoType.OGlycoPep && is_HCD_only_data)
            {
                List<LocalizationGraph> graphs = new List<LocalizationGraph>();

                foreach (var gb_id in gsmcd.GlycanBoxIds)
                {
                    LocalizationGraph localizationGraph = new LocalizationGraph( gsmcd.ModPos, gsmcd.ModMotifs, globalBoxes[gb_id], globalBoxes[gb_id].ChildGlycanBoxes, gb_id);
                    graphs.Add(localizationGraph);

                }

                var noSpecificRoute = LocalizationGraph.GetAnyOnePath(graphs[0]);

                var ogsm = _GenerateGsm(theScanBestPeptide, graphs, noSpecificRoute, gsmcd, theScan, scanIndex, oxoniumIonIntensities, globalBoxes, globalglycans, globalmods);

                return ogsm;

            }

            List<LocalizationGraph> localizationGraphs = _GenerateLocalizationGraph(theScanBestPeptide, gsmcd, theScan, globalBoxes);


            var firstPath = LocalizationGraph.GetFirstPath(localizationGraphs[0].array, localizationGraphs[0].ChildModBoxes);
            var route = LocalizationGraph.GetLocalizedPath(localizationGraphs[0], firstPath);

            var gsm = _GenerateGsm(theScanBestPeptide, localizationGraphs, route, gsmcd, theScan, scanIndex, oxoniumIonIntensities, globalBoxes, globalglycans, globalmods);

            return gsm;

        }

        private List<LocalizationGraph> _GenerateLocalizationGraph(PeptideWithSetModifications theScanBestPeptide, GlycoSpectraMatchCandidate gsmcd, 
            Ms2ScanWithSpecificMass theScan, GlycanBox[] globalBoxes)
        {
            List<Product> hcdProducts = new List<Product>();
            theScanBestPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, hcdProducts);
            //Peptide ETD fragments for graph localization.
            List<Product> etdProducts = new List<Product>();
            theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, etdProducts);

            List<Product> mainProducts = new List<Product>();
            if (gsmcd.GType!= GlycoType.OGlycoPep && GlycoPeptides.DissociationTypeContainHCD(CommonParameters.DissociationType))
            {
                mainProducts.AddRange(hcdProducts);
            }
            if (GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType))
            {
                //TO THINK: it is still unsure about what ions should be kept for etd or ethcd for glycopeptides.             
                mainProducts.AddRange(etdProducts.Where(p => p.ProductType != ProductType.y));
            }

            List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();

            double bestLocalizedScore = 0; 
            foreach (var iDLow in gsmcd.GlycanBoxIds)
            {
                LocalizationGraph localizationGraph = new LocalizationGraph(gsmcd.ModPos, gsmcd.ModMotifs, globalBoxes[iDLow], globalBoxes[iDLow].ChildGlycanBoxes, iDLow);

                if (mainProducts.Count >0)
                {
                    LocalizationGraph.LocalizeMod(localizationGraph, theScan, CommonParameters.ProductMassTolerance, mainProducts, GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
                }

                if (theScan.ChildScans.Count > 0)
                {
                    foreach (var childScan in theScan.ChildScans)
                    {
                        List<Product> childProducts = new List<Product>();

                        if (gsmcd.GType != GlycoType.OGlycoPep && GlycoPeptides.DissociationTypeContainHCD(CommonParameters.MS2ChildScanDissociationType))
                        {
                            childProducts.AddRange(hcdProducts);
                        }

                        if (GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType))
                        {
                            childProducts.AddRange(etdProducts.Where(p => p.ProductType != ProductType.y));
                        }

                        //run the same graph multiple times!
                        if (childProducts.Count >0)
                        {
                            LocalizationGraph.LocalizeMod(localizationGraph, childScan, CommonParameters.ProductMassTolerance, childProducts, GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
                        }
                    }
                }

                double currentLocalizationScore = localizationGraph.TotalScore;

                if (currentLocalizationScore > bestLocalizedScore)
                {
                    bestLocalizedScore = currentLocalizationScore;
                    localizationGraphs.Clear();
                    localizationGraphs.Add(localizationGraph);
                }
                else if ((currentLocalizationScore <= bestLocalizedScore + 0.00000001 && currentLocalizationScore >= bestLocalizedScore - 0.00000001))
                {
                    localizationGraphs.Add(localizationGraph);
                }
            }

            return localizationGraphs;          
        }

        private GlycoSpectralMatch _GenerateGsm(PeptideWithSetModifications theScanBestPeptide, List<LocalizationGraph> localizationGraphs, Route route,
            GlycoSpectraMatchCandidate gsmcd, Ms2ScanWithSpecificMass theScan, int scanIndex, double[] oxoniumIonIntensities, 
            GlycanBox[] globalBoxes, Glycan[] globalglycans, Modification[] globalmods)
        {
            
            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(route, theScanBestPeptide, globalmods);

            Glycan[] glycans = new Glycan[localizationGraphs.First().ModBox.ModCount];
            for (int i = 0; i < localizationGraphs.First().ModBox.ModCount; i++)
            {
                glycans[i] = globalglycans[((GlycanBox)localizationGraphs.First().ModBox).ModIds[i]];
            }

            List<int> npos = route.Mods.Select(p => p.ModSite).ToArray().Intersect(gsmcd.NPos).ToList();

            var fragmentsForEachGlycoPeptide = GlycoPeptides.GlyGetTheoreticalFragments(gsmcd.GType, CommonParameters.DissociationType, theScanBestPeptide, peptideWithMod, npos, glycans);
           
            //var matchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);
            var matchedIons = GlycoPeptides.GlyMatchOriginFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);

            double score = CalculatePeptideScore(theScan.TheScan, matchedIons.Where(p=>p.NeutralTheoreticalProduct.ProductType != ProductType.D).ToList());
            double pScore = GlycoPeptides.Calc_pGlycoScore(matchedIons, fragmentsForEachGlycoPeptide, CommonParameters);

            var DiagnosticIonScore = CalculatePeptideScore(theScan.TheScan, matchedIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());

            var GlycanScore = CalculatePeptideScore(theScan.TheScan, matchedIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.Ycore || v.NeutralTheoreticalProduct.ProductType == ProductType.Y).ToList());

            var p = theScan.TheScan.MassSpectrum.Size * CommonParameters.ProductMassTolerance.GetRange(1000).Width / theScan.TheScan.MassSpectrum.Range.Width;

            int n = fragmentsForEachGlycoPeptide.Where(v => v.ProductType != ProductType.D && v.ProductType != ProductType.Ycore && v.ProductType != ProductType.Y).Count();

            var allMatchedChildIons = new Dictionary<int, List<MatchedFragmentIon>>();

            foreach (var childScan in theScan.ChildScans)
            {
                //People always use CID with low res. This special code works for Nic Scott's data. (High-HCD, High-EThcD, Low-CID, High-secHCD)
                if (childScan.TheScan.DissociationType == DissociationType.CID)
                {
                    continue;
                }
                var childFragments = GlycoPeptides.GlyGetTheoreticalFragments(gsmcd.GType, CommonParameters.MS2ChildScanDissociationType, theScanBestPeptide, peptideWithMod, npos, glycans);

                var matchedChildIons = MatchFragmentIons(childScan, childFragments, CommonParameters);

                n += childFragments.Where(v => v.ProductType != ProductType.D && v.ProductType != ProductType.Ycore && v.ProductType != ProductType.Y).Count();

                if (matchedChildIons == null)
                {
                    continue;
                }

                allMatchedChildIons.Add(childScan.OneBasedScanNumber, matchedChildIons);
                double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons.Where(p => p.NeutralTheoreticalProduct.ProductType != ProductType.D).ToList());
                double child_pScore = GlycoPeptides.Calc_pGlycoScore(matchedChildIons, childFragments, CommonParameters);

                double childDiagnosticIonScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());
                double childGlycanScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.Ycore || v.NeutralTheoreticalProduct.ProductType == ProductType.Y).ToList());

                DiagnosticIonScore += childDiagnosticIonScore;
                GlycanScore += childGlycanScore;

                //TO THINK:may think a different way to use childScore
                score += childScore;
                pScore += child_pScore;

                p += childScan.TheScan.MassSpectrum.Size * CommonParameters.ProductMassTolerance.GetRange(1000).Width / childScan.TheScan.MassSpectrum.Range.Width;

            }

            double _score = 0;
            if (GlycoScoreType == GlycoScoreType.MMScore)
            {
                _score = gsmcd.Score;
            }
            else if (GlycoScoreType == GlycoScoreType.pGlycoScore)
            {
                _score = gsmcd.PScore;
            }
            else if (GlycoScoreType == GlycoScoreType.XcorrScore)
            {
                _score = gsmcd.XcorrScore;
            }

            var psmGlyco = new GlycoSpectralMatch(peptideWithMod, 0, _score, scanIndex, theScan, CommonParameters, matchedIons);

            psmGlyco.MScore = score;
            psmGlyco.PScore = pScore;
            psmGlyco.XcorrScore = gsmcd.XcorrScore;

            psmGlyco.GlycanType = gsmcd.GType;

            //TO THINK: Is this p and n properly calculated for different data type.
            psmGlyco.ScanInfo_p = p;

            psmGlyco.Thero_n = n;

            psmGlyco.Rank = gsmcd.Rank;

            psmGlyco.DiagnosticIonScore = DiagnosticIonScore;

            psmGlyco.GlycanScore = GlycanScore;

            psmGlyco.ChildMatchedFragmentIons = allMatchedChildIons;

            psmGlyco.LocalizationGraphs = localizationGraphs;

            psmGlyco.OxoniumIonIntensity = oxoniumIonIntensities;

            var product_nn = GlycoPeptides.GetIndicatorYIon(theScanBestPeptide.MonoisotopicMass, "NN");
            psmGlyco.PepNN = GlycoPeptides.MatchIndicatorYIon(theScan, product_nn, CommonParameters);

            var product_nh = GlycoPeptides.GetIndicatorYIon(theScanBestPeptide.MonoisotopicMass, "NH");
            psmGlyco.PepNH = GlycoPeptides.MatchIndicatorYIon(theScan, product_nh, CommonParameters);

            return psmGlyco;
        }

        #endregion

        #region Original O-Pair functions, plan to deprecated.
        private GlycoSpectralMatch CreateOGsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, Route localization, double[] oxoniumIonIntensities, List<LocalizationGraph> localizationGraphs)
        {
            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(localization, peptide, GlycanBox.GlobalOGlycanMods);

            var fragmentsForEachGlycoPeptide = GlycoPeptides.OGlyGetTheoreticalFragments(CommonParameters.DissociationType, peptide, peptideWithMod);

            var matchedIons = MatchFragmentIons(theScan, fragmentsForEachGlycoPeptide, CommonParameters);

            //double score = CalculatePeptideScore(theScan.TheScan, matchedIons);
            double score = GlycoPeptides.Calc_pGlycoScore(matchedIons, fragmentsForEachGlycoPeptide, CommonParameters);

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
                //double childScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons);
                double childScore = GlycoPeptides.Calc_pGlycoScore(matchedChildIons, childFragments, CommonParameters);

                double childDiagnosticIonScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.D).ToList());
                double childGlycanScore = CalculatePeptideScore(childScan.TheScan, matchedChildIons.Where(v => v.NeutralTheoreticalProduct.ProductType == ProductType.M).ToList());

                DiagnosticIonScore += childDiagnosticIonScore;
                GlycanScore += childGlycanScore;

                PeptideScore += childScore - childDiagnosticIonScore;
                //TO THINK:may think a different way to use childScore
                score += childScore;

                p += childScan.TheScan.MassSpectrum.Size * CommonParameters.ProductMassTolerance.GetRange(1000).Width / childScan.TheScan.MassSpectrum.Range.Width;

            }

            //var psmGlyco = new GlycoSpectralMatch(peptideWithMod, 0, PeptideScore, scanIndex, theScan, CommonParameters, matchedIons);
            var psmGlyco = new GlycoSpectralMatch(peptideWithMod, 0, score, scanIndex, theScan, CommonParameters, matchedIons);

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

            var localizationScan = theScan;
            List<Product> products = new List<Product>();

            //For HCD-pd-ETD or CD-pd-EThcD type of data
            if (theScan.ChildScans.Count > 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType))
            {
                localizationScan = theScan.ChildScans.First();
                theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
            }

            //For ETD type of data
            if (theScan.ChildScans.Count == 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType))
            {
                theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
            }

            //Localization for O-glycopeptides only works on ETD related dissociationtype
            //No localization can be done with MS2-HCD spectrum
            //TO THINK: there is a special situation. The HCD only scan from  HCD-pd-EThcD data can be a glycopeptide, but there is no ETD, so there is no localization. What to do with this?
            bool is_HCD_only_data = !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType);
            if (is_HCD_only_data)
            {
                theScanBestPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            }

            double bestLocalizedScore = 0;

            List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();

            string[] modMotifs = new string[1];
            //GetModPosMotif(glycanType, null, o_modPos, out modPos, out modMotifs);

            while (iDLow < GlycanBox.OGlycanBoxes.Count() && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBox.OGlycanBoxes[iDLow].Mass)))
            {
                if (modPos.Length >= GlycanBox.OGlycanBoxes[iDLow].ModCount && GlycoPeptides.OxoniumIonsAnalysis(oxoniumIonIntensities, GlycanBox.OGlycanBoxes[iDLow]))
                {
                    LocalizationGraph localizationGraph = new LocalizationGraph(modPos, modMotifs, GlycanBox.OGlycanBoxes[iDLow], GlycanBox.OGlycanBoxes[iDLow].ChildGlycanBoxes, iDLow);
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
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode, CommonParameters.DissociationType);

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

        #endregion
    
    }
}