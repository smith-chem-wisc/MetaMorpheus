using EngineLayer.GlycoSearch;
using EngineLayer.ModernSearch;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using MathNet.Numerics;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.ModSearch
{
    public class ModSearchEngine : ModernSearchEngine
    {
        protected readonly List<GlycoSpectralMatch>[] GlobalMsms; // Global mod spectral matches

        // For general modification settings
        private readonly int TopN; // DDA top Peak number
        private readonly int MaxModNumber;
        private readonly Tolerance PrecursorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;
        private readonly List<int>[] SecondFragmentIndex;
        private readonly ModBox[] ModBoxes; // The modBox for glyco search

        // For glyco settings
        private readonly bool OxoniumIonFilter;
        private readonly GlycoSearchType GlycoSearchType;
        private readonly string OGlycanDatabaseFile;
        private readonly string NGlycanDatabaseFile;

        public ModSearchEngine(List<GlycoSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans,
            List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition,
            CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
            string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, int modSearchTopNum,
            int maxModNum, bool oxoniumIonFilter, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters,
                fileSpecificParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalMsms = globalCsms;
            this.SecondFragmentIndex = secondFragmentIndex;
            this.TopN = modSearchTopNum;
            this.MaxModNumber = maxModNum;
            this.OxoniumIonFilter = oxoniumIonFilter;
            this.GlycoSearchType = glycoSearchType;
            this.OGlycanDatabaseFile = oglycanDatabase;
            this.NGlycanDatabaseFile = nglycanDatabase;

            this.SecondFragmentIndex = secondFragmentIndex;
            this.ProductSearchMode = new SinglePpmAroundZeroSearchMode(20); //For ScanOxoniumIonFilter only
            this.PrecursorSearchMode = commonParameters.PrecursorMassTolerance;

            //Load glycan databases and build the modBox 
            ModBox.GlobalModifications = GlycanDatabase
                .LoadGlycan(
                    GlobalVariables.OGlycanDatabasePaths.First(
                        p => System.IO.Path.GetFileName(p) == OGlycanDatabaseFile), true, true).ToArray();
            ModBoxes = ModBox.BuildModBoxes(MaxModNumber, false).OrderBy(p => p.Mass).ToArray();
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentComplete = 0;
            ReportProgress(new ProgressEventArgs(oldPercentComplete, "Performing modPair search..." + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;
            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, scanIndex =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();

                List<int> idsOfPeptidesTopN = new List<int>();
                byte scoreAtTopN = 0;
                int peptideCount = 0;

                for (; scanIndex < ListOfSortedMs2Scans.Length; scanIndex += maxThreadsPerFile)
                {
                    if (GlobalVariables.StopLoops) return; // Stop the loop if requested

                    Array.Clear(scoringTable, 0 , scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    idsOfPeptidesTopN.Clear();

                    var scan = ListOfSortedMs2Scans[scanIndex];

                    // Get fragments bins for this scan
                    List<int> allBinsToSearch = GetBinsToSearch(scan, FragmentIndex, CommonParameters.DissociationType);

                    //Limit the high bound limitation, here assume it is possible to has max 3 Da
                    var high_bound_limitation = scan.PrecursorMass + 1;

                    // Score the peptides candidates with binSearching. (Modern search)
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);
                    
                    // filtering the peptides candidate with the cufOff and limit the topN peptides.
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

                        gsms = FindModPeptide(scan, idsOfPeptidesPossiblyObserved, scanIndex, (int)byteScoreCutoff);


                        if (gsms.Count == 0)
                        {
                            progress++;
                            continue;
                        }

                        if (GlobalMsms[scanIndex] == null)
                        {
                            GlobalMsms[scanIndex] = new List<GlycoSpectralMatch>(); //the first one finished task, create teh new gsms list.
                        }
                        else
                        {
                            gsms.AddRange(GlobalMsms[scanIndex]);
                            GlobalMsms[scanIndex].Clear();
                        }

                        Add2GlobalMsms(ref gsms, scanIndex);

                    }
                }

                // report search progress
                progress++;
                var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                if (percentProgress > oldPercentComplete)
                {
                    oldPercentComplete = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Performing Mod search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
                }   //percentProgress = 100, "Performing glyco search...1/1", NestedIds = 3.
            });

            return new MetaMorpheusEngineResults(this); //Storage the result information into the result class.
        }

        /// <summary>
        /// Method to find the mod peptide. This method includes two searching steps:
        /// (1) Match the peptide candidates with the precursor mass to find the mass difference.
        /// (2) If mass difference is zero, then we match the peptide with unmodified precursor mass.
        /// If the mass difference is not zero, then search the possible modification.
        /// </summary>
        /// <returns></returns>
        private List<GlycoSpectralMatch> FindModPeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++) 
            {
                var thePeptideCandidate = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]]; // Get the peptide from the candidate list.

                if (PrecursorSearchMode.Within(theScan.PrecursorMass, thePeptideCandidate.MonoisotopicMass))
                {
                    MatchUnmodifiedPeptide(theScan, scanIndex, scoreCutOff, thePeptideCandidate, ind, ref possibleMatches);
                }
                else if (theScan.PrecursorMass - thePeptideCandidate.MonoisotopicMass >= 100) // 100 can be revised, it will be the minimum mass difference for the modification.
                {
                    double massDiffLow = PrecursorSearchMode.GetMinimumValue(theScan.PrecursorMass) - thePeptideCandidate.MonoisotopicMass;
                    double massDiffHigh = PrecursorSearchMode.GetMaximumValue(theScan.PrecursorMass) - thePeptideCandidate.MonoisotopicMass;

                    if (massDiffHigh < ModBoxes.First().Mass || massDiffLow > ModBoxes.Last().Mass)
                    {
                        continue; // if the mass difference is out of the range of the ModBox, we can skip this peptide.
                    }

                    // Search the possible modification
                    // (1) filter by Oxonium Ion
                    var oxoniumIonIntensities = GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode);
                    // (2) do the ModPair Search
                    MatchModifiedPeptide(theScan, scanIndex, scoreCutOff, thePeptideCandidate, ind, massDiffLow, oxoniumIonIntensities, ref possibleMatches);
                }

            }

            if (possibleMatches.Count != 0)
            {
                possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
            }

            return possibleMatches;
        }

        /// <summary>
        /// Match the mass of the peptide candidate with the precursor mass.
        /// </summary>
        private void MatchUnmodifiedPeptide(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications thePeptideCandidate, int ind, ref List<GlycoSpectralMatch> possibleMatches)
        {
            List<Product> products = new List<Product>();
            thePeptideCandidate.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both, products);
            var matchedFragmentIons = MatchFragmentIons(theScan, products, CommonParameters);
            double score = CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);

            if (score > scoreCutOff)
            {
                var psm = new GlycoSpectralMatch(thePeptideCandidate, 0, score, scanIndex, theScan, CommonParameters, matchedFragmentIons);
                psm.Rank = ind;

                possibleMatches.Add(psm);
            }
        }

        /// <summary>
        /// Match the mass difference with corresponding modification and localization.
        /// The generated Gsms will be stored in the possibleMatches.
        /// </summary>
        private void MatchModifiedPeptide(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications thePeptideCandidate, int ind, double possibleGlycanMassLow, double[] oxoniumIonIntensities, ref List<GlycoSpectralMatch> possibleMatches)
        {
            // try to find the index that closet match to the "possibleGlycanMassLow" within the glycanBox
            int boxIdForThisMass = GlycoPeptides.BinarySearchGetIndex(ModBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow);

            // We generate different products from different dissociation types.
            List<Product> products = new List<Product>(); // product list for the theoretical fragment ions
            var localizationScan = theScan; // The scan for localization, it can be the same as theScan or the child scan.

            // For HCD-pd-ETD or CD-pd-EThcD
            if (theScan.ChildScans.Count > 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType, CommonParameters.CustomIons))
            {
                localizationScan = theScan.ChildScans.First();
                thePeptideCandidate.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
            }

            //For ETD type of data
            if (theScan.ChildScans.Count == 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType, CommonParameters.CustomIons))
            {
                thePeptideCandidate.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
            }

            // For HCD, we don't need to localize the modification, because the HCD only generates b and c ions, which are not modified.
            bool is_HCD_only_data = !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType, CommonParameters.CustomIons) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType, CommonParameters.CustomIons);
            if (is_HCD_only_data)
            {
                thePeptideCandidate.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            }

            // Do the localization of the modification
            double bestLocalizedScore = 0;
            List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();

            // Iterate all modBoxes with the mass shift
            while (boxIdForThisMass < ModBoxes.Count() && (PrecursorSearchMode.Within(theScan.PrecursorMass, thePeptideCandidate.MonoisotopicMass + ModBoxes[boxIdForThisMass].Mass))) 
            {
                if (OxoniumIonFilter && !GlycoPeptides.DiagonsticFilter(oxoniumIonIntensities, ModBoxes[boxIdForThisMass])) // if the filter is turned on, we need to check does the oxoiums make sense.
                {
                    boxIdForThisMass++; // if the oxonium ions don't make sense (there is no 204, or without their diagnostic ion), we can skip this glycan.
                    continue;
                }

                // Generate the modPos, which is the possible modSite position on the peptide.
                HashSet<string> modMotifs = ModBoxes[boxIdForThisMass].ModIds
                    .Select(p => ModBox.GlobalModifications[p].Target.ToString())
                    .ToHashSet();
                SortedDictionary<int, string> modPos = GlycoSpectralMatch.GetPossibleModSites(thePeptideCandidate, modMotifs);

                if (GraphCheck(modPos, ModBoxes[boxIdForThisMass])) // the modSite number should be larger than the possible mod number. ex. we have two S in the peptide, then we can not localize the modBox with more than two S motifs
                {
                    // Create the localization graph with the mod mass and the possible modSite.
                    LocalizationGraph localizationGraph = new LocalizationGraph(modPos, ModBoxes[boxIdForThisMass], ModBoxes[boxIdForThisMass].ChildModBoxes, boxIdForThisMass);
                    LocalizationGraph.LocalizeOGlycan(localizationGraph, localizationScan, CommonParameters.ProductMassTolerance, products); // Temporary code, we will create the new LocalizedMod function in the future.

                    double currentLocalizationScore = localizationGraph.TotalScore;
                    if (currentLocalizationScore > bestLocalizedScore) //Try to find the best modBox with the highest score.
                    {
                        bestLocalizedScore = currentLocalizationScore;
                        localizationGraphs.Clear();
                        localizationGraphs.Add(localizationGraph); // we only keep the best modBox and its localizationGraph.
                    }
                    else if ((is_HCD_only_data || bestLocalizedScore > 0) && (currentLocalizationScore <= bestLocalizedScore + 0.00000001 && currentLocalizationScore >= bestLocalizedScore - 0.00000001))
                    {
                        localizationGraphs.Add(localizationGraph);
                    }
                }
                boxIdForThisMass++;
            }

            //In theory, the peptide_localization shouldn't be null, but it is possible that the real score is smaller than indexed score.
            if (localizationGraphs.Count > 0)
            {
                var firstPath = LocalizationGraph.GetFirstPath(localizationGraphs[0].array, localizationGraphs[0].ChildModBoxes); //Get the first path from the localization graph.
                var localizationCandidate = LocalizationGraph.GetLocalizedPath(localizationGraphs[0], firstPath); //Get the route of the localization from the first path inforation

                var psmMod = CreateMsm(theScan, scanIndex, ind, thePeptideCandidate, localizationCandidate, oxoniumIonIntensities, localizationGraphs); //Create the glycoSpectralMatch

                if (psmMod.Score > scoreCutOff)
                {
                    possibleMatches.Add(psmMod);
                }
            }
        }

        /// <summary>
        /// Valid the Graph created by this modPos and modBox.
        /// Check if the motif in peptide is sufficient to cover the motif in glycanBox.
        /// </summary>
        /// <param name="modPos"></param>
        /// <param name="glycanBox"></param>
        /// <returns></returns>
        private bool GraphCheck(SortedDictionary<int, string> modPos, ModBox modBox)
        {
            // If the motifs number is less than the glycanBox, we can skip this graph.
            if (modPos.Count < modBox.NumberOfMods)
                return false;

            // Calculate the motif in glycanBox.
            var motifInBox = new Dictionary<string, int>();
            foreach (var modId in modBox.ModIds)
            {
                var motif = ModBox.GlobalModifications[modId].Target.ToString();

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

        //For MatchModifiedPeptide, generate the msms (ModSpectralMatches) for Mod search
        private GlycoSpectralMatch CreateMsm(Ms2ScanWithSpecificMass theScan, int scanIndex, int rank, PeptideWithSetModifications peptide, Route localization, double[] oxoniumIonIntensities, List<LocalizationGraph> localizationGraphs)
        {
            return null;
        }

        private void Add2GlobalMsms(ref List<GlycoSpectralMatch> gsms, int scanIndex)
        {
        }

    }
}
