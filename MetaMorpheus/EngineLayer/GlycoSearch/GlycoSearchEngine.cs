using EngineLayer.ModernSearch;
using MzLibUtil;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MassSpectrometry;
using EngineLayer.SpectrumMatch;

namespace EngineLayer.GlycoSearch
{
    public class GlycoSearchEngine : ModernSearchEngine
    {
        public static readonly double ToleranceForMassDifferentiation = 1e-9;
        private readonly int OxoniumIon204Index = OxoniumIonReservedIndices.HexNAc204; // Check Glycan.AllOxoniumIons
        protected readonly List<GlycoSpectralMatch>[] GlobalGsms;  // Why don't we call it GlobalGsms?

        private GlycoSearchType GlycoSearchType;
        private readonly int TopN;              // DDA top Peak number.
        private readonly int _maxOGlycanNum;
        private readonly bool OxoniumIonFilter; // To filt Oxonium Ion before searching a spectrum as glycopeptides. If we filter spectrum, it must contain oxonium ions such as 204 (HexNAc). 
        private readonly string _oglycanDatabase;
        private readonly string _nglycanDatabase;
        private readonly GlycanBox[] GlycanBoxes; // GlycanBoxes for glycan search.
        private readonly double[] GlycanBoxMasses; // GlycanBoxes[i].Mass, built once so each candidate peptide does not copy it.

        private readonly Tolerance PrecusorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;

        private readonly Indexing.FragmentIndex SecondFragmentIndex;

        /// <summary>
        /// Per scan, the candidates that made the TopN cut across every partition searched so far, as
        /// (partition, peptide id within that partition, coarse score). Null for a single-partition search,
        /// which scores and matches in one pass. With more than one partition the cut has to be taken over the
        /// whole database rather than per partition -- otherwise N partitions send up to N x TopN candidates
        /// to glycan matching and the identifications depend on the partition count. Same two-round shape as
        /// CrosslinkSearchEngine: <see cref="FirstRoundSearch"/> over every partition, then Run() per partition.
        /// </summary>
        protected readonly List<(int Partition, int PeptideId, byte Score)>[] Candidates;

        /// <summary>
        /// Glyco search is proteomics-only, so it keeps a peptide-typed view of the index that
        /// ModernSearchEngine now holds as IBioPolymerWithSetMods. Same objects, narrower type.
        /// </summary>
        protected new readonly List<PeptideWithSetModifications> PeptideIndex;

        protected string[] Motifs
        {
            get
            {
                if (GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
                {
                    return new string[] { "S", "T", "Nxs", "Nxt" };
                }
                return new string[] { "S", "T" }; // motif for the OSearch
            }
        }

        // The constructor for GlycoSearchEngine, we can load the parameter for the searhcing like mode, topN, maxOGlycanNum, oxoniumIonFilter, datsbase, etc.
        public GlycoSearchEngine(List<GlycoSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, IEnumerable<IBioPolymerWithSetMods> peptideIndex,
            Indexing.FragmentIndex fragmentIndex, Indexing.FragmentIndex secondFragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
             string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, int glycoSearchTopNum, int maxOGlycanNum, bool oxoniumIonFilter, List<string> nestedIds,
             double maxGlycanBoxMass = GlycanBox.DefaultMaximumGlycanBoxMass, List<(int Partition, int PeptideId, byte Score)>[] candidates = null)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, fileSpecificParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.PeptideIndex = peptideIndex.Cast<PeptideWithSetModifications>().ToList();
            this.Candidates = candidates;
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
                GlycanBox.OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(_maxOGlycanNum, false, maxGlycanBoxMass).OrderBy(p => p.Mass).ToArray(); //generate glycan box for O-glycan search
                GlycanBoxes = GlycanBox.OGlycanBoxes;
                GlycoSpectralMatch.GlycanBoxes = GlycanBoxes;
            }
            else if (glycoSearchType == GlycoSearchType.NGlycanSearch) //because the there is only one glycan in N-glycanpeptide, so we don't need to build the n-glycanBox here.
            {
                // The single N-glycan is the whole box here, so the box mass cap applies to each glycan on its own.
                NGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false)
                    .Where(p => (double)p.Mass / 1E5 <= maxGlycanBoxMass).OrderBy(p => p.Mass).ToArray();
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
            }
            else if (glycoSearchType == GlycoSearchType.N_O_GlycanSearch) //search both N-glycan and O-glycan is still not tested and build completely yet.
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == _oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.GlobalNGlycans = new Dictionary<int, Glycan>();
                // For N-glycan, we use negative index to distinguish with O-glycan.
                var nGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanDatabasePaths.First(p => System.IO.Path.GetFileName(p) == _nglycanDatabase),
                        true, false).OrderBy(p => p.Mass);
                int indexForNGlycan = -1;
                foreach (var nGlycan in nGlycans)
                {
                    GlycanBox.GlobalNGlycans.Add(indexForNGlycan, nGlycan);
                    indexForNGlycan--;
                }

                GlycanBox.NOGlycanBoxes = GlycanBox.BuildNOGlycanBoxes(_maxOGlycanNum, false, maxGlycanBoxMass).OrderBy(p => p.Mass).ToArray();
                GlycanBoxes = GlycanBox.NOGlycanBoxes;
                GlycoSpectralMatch.GlycanBoxes = GlycanBoxes;
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
            }

            GlycanBoxMasses = GlycanBoxes?.Select(p => p.Mass).ToArray();
            NGlycanMasses = NGlycans?.Select(p => (double)p.Mass / 1E5).ToArray();
        }

        private Glycan[] NGlycans { get; }
        private double[] NGlycanMasses { get; } // NGlycans[i].Mass in Da, built once so each candidate peptide does not copy it.
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
            if (Candidates != null)
            {
                return SecondRoundSearch();
            }

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
                int[] candidateCountsByScore = new int[byte.MaxValue + 1];

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
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, base.PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);

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

                    //    IndexedScoring(SecondFragmentIndex, childBinsToSearch, secondScoringTable, byteScoreCutoff, childIdsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, base.PeptideIndex, MassDiffAcceptor, 0, CommonParameters.MS2ChildScanDissociationType);

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
                        SelectTopCandidates(idsOfPeptidesPossiblyObserved, scoringTable, byteScoreCutoff, TopN, candidateCountsByScore, idsOfPeptidesTopN);

                        List<GlycoSpectralMatch> gsms = MatchCandidates(scan, idsOfPeptidesTopN, scanIndex, (int)byteScoreCutoff, null);

                        if (gsms.Count == 0)
                        {
                            progress++;
                            continue;
                        }

                        MergeIntoGlobalGsms(gsms, scanIndex);

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

        /// <summary>
        /// Coarse scoring only, for a search split over more than one partition. Records each scan's
        /// candidates from this partition into <see cref="Candidates"/> and trims the running list back to the
        /// TopN cut over every partition seen so far. No glycan matching happens here; that is Run(), once every
        /// partition has been through this.
        /// </summary>
        public void FirstRoundSearch()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing glyco search first round... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;
            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (scanIndex) =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                List<int> idsOfPeptidesPossiblyObserved = new List<int>();
                List<int> idsOfPeptidesTopN = new List<int>();
                int[] candidateCountsByScore = new int[byte.MaxValue + 1];

                for (; scanIndex < ListOfSortedMs2Scans.Length; scanIndex += maxThreadsPerFile)
                {
                    if (GlobalVariables.StopLoops) { return; }

                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    idsOfPeptidesTopN.Clear();

                    var scan = ListOfSortedMs2Scans[scanIndex];
                    List<int> allBinsToSearch = GetBinsToSearch(scan, FragmentIndex, CommonParameters.DissociationType);

                    // same bound as the single-pass search in RunSpecific
                    var high_bound_limitation = scan.PrecursorMass + 1;
                    IndexedScoring(FragmentIndex, allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, Double.NegativeInfinity, high_bound_limitation, base.PeptideIndex, MassDiffAcceptor, 0, CommonParameters.DissociationType);

                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        // Cutting this partition to TopN first cannot lose a global survivor: a partition's TopN-th
                        // score is never above the database's, so anything at or above the database's is kept here.
                        SelectTopCandidates(idsOfPeptidesPossiblyObserved, scoringTable, byteScoreCutoff, TopN, candidateCountsByScore, idsOfPeptidesTopN);

                        if (idsOfPeptidesTopN.Count > 0)
                        {
                            var scanCandidates = Candidates[scanIndex] ?? new List<(int Partition, int PeptideId, byte Score)>();
                            foreach (int id in idsOfPeptidesTopN)
                            {
                                scanCandidates.Add((CurrentPartition - 1, id, scoringTable[id]));
                            }
                            Candidates[scanIndex] = KeepGlobalTopN(scanCandidates, TopN);
                        }
                    }

                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);
                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing glyco search first round... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));
                    }
                }
            });
        }

        /// <summary>
        /// Glycan matching for this partition's share of the candidates <see cref="FirstRoundSearch"/> kept.
        /// A candidate's Rank is its position in the scan's cut over the whole database, not within this partition.
        /// </summary>
        private MetaMorpheusEngineResults SecondRoundSearch()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing glyco search... " + CurrentPartition + "/" + CommonParameters.TotalPartitions, NestedIds));

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;
            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (scanIndex) =>
            {
                List<int> ids = new List<int>();
                List<int> ranks = new List<int>();

                for (; scanIndex < ListOfSortedMs2Scans.Length; scanIndex += maxThreadsPerFile)
                {
                    if (GlobalVariables.StopLoops) { return; }

                    var scanCandidates = Candidates[scanIndex];
                    if (scanCandidates != null)
                    {
                        ids.Clear();
                        ranks.Clear();
                        for (int rank = 0; rank < scanCandidates.Count; rank++)
                        {
                            if (scanCandidates[rank].Partition == CurrentPartition - 1)
                            {
                                ids.Add(scanCandidates[rank].PeptideId);
                                ranks.Add(rank);
                            }
                        }

                        if (ids.Count > 0)
                        {
                            List<GlycoSpectralMatch> gsms = MatchCandidates(ListOfSortedMs2Scans[scanIndex], ids, scanIndex, (int)byteScoreCutoff, ranks);
                            if (gsms.Count > 0)
                            {
                                MergeIntoGlobalGsms(gsms, scanIndex);
                            }
                        }
                    }

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

        /// <summary>
        /// The same cut as <see cref="SelectTopCandidates"/>, applied to candidates pooled from several partitions. The sort
        /// is stable, so equal scores stay in the order they were added: partition order, then the order
        /// SelectTopCandidates produced within a partition.
        /// </summary>
        internal static List<(int Partition, int PeptideId, byte Score)> KeepGlobalTopN(List<(int Partition, int PeptideId, byte Score)> candidates, int topN)
        {
            var ordered = candidates.OrderByDescending(c => c.Score).ToList();
            if (topN <= 0 || ordered.Count <= topN)
            {
                return ordered;
            }

            byte scoreAtTopN = ordered[topN - 1].Score;
            int keep = topN;
            while (keep < ordered.Count && ordered[keep].Score >= scoreAtTopN)
            {
                keep++;
            }
            ordered.RemoveRange(keep, ordered.Count - keep);
            return ordered;
        }

        private List<GlycoSpectralMatch> MatchCandidates(Ms2ScanWithSpecificMass scan, List<int> ids, int scanIndex, int scoreCutOff, List<int> ranks)
        {
            if (GlycoSearchType == GlycoSearchType.OGlycanSearch || GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                return MatchGlycopeptide(scan, ids, scanIndex, scoreCutOff, ranks); // Use the peptide candidate and the scan to generate the gsms.
            }
            return MatchNGlycopeptide(scan, ids, scanIndex, scoreCutOff, ranks);
        }

        private void MergeIntoGlobalGsms(List<GlycoSpectralMatch> gsms, int scanIndex)
        {
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

        /// <summary>
        /// Keeps the candidate peptides worth matching: every id scoring at least <paramref name="scoreCutoff"/>, cut after the
        /// <paramref name="topN"/>th best but keeping all ids tied with it, written to <paramref name="topCandidates"/> from highest score
        /// to lowest and, within a score, in the order the ids were observed. That is exactly what a stable descending sort by score,
        /// stopped below the topN-th score, produces. Scores are bytes, so counting them replaces sorting the whole candidate list, which
        /// can run to hundreds of thousands of ids per scan when only 50 are kept.
        /// </summary>
        /// <param name="countsByScore"> Scratch space of length 256, reused across scans. Left zeroed. </param>
        /// <param name="topCandidates"> Cleared and filled. </param>
        internal static void SelectTopCandidates(List<int> candidateIds, byte[] scores, int scoreCutoff, int topN, int[] countsByScore, List<int> topCandidates)
        {
            topCandidates.Clear();
            foreach (int id in candidateIds)
            {
                countsByScore[scores[id]]++;
            }

            // The lowest score kept: the score of the topN-th best candidate, or the cutoff when there are fewer than topN.
            int lowestKeptScore = Math.Max(scoreCutoff, 0);
            if (topN > 0)
            {
                int atOrAbove = 0;
                for (int score = byte.MaxValue; score >= lowestKeptScore; score--)
                {
                    atOrAbove += countsByScore[score];
                    if (atOrAbove >= topN)
                    {
                        lowestKeptScore = score;
                        break;
                    }
                }
            }

            // Where each kept score's ids start in the output, highest score first.
            int kept = 0;
            for (int score = byte.MaxValue; score >= lowestKeptScore; score--)
            {
                int count = countsByScore[score];
                countsByScore[score] = kept;
                kept += count;
            }

            System.Runtime.InteropServices.CollectionsMarshal.SetCount(topCandidates, kept);
            var output = System.Runtime.InteropServices.CollectionsMarshal.AsSpan(topCandidates);
            foreach (int id in candidateIds)
            {
                int score = scores[id];
                if (score >= lowestKeptScore)
                {
                    output[countsByScore[score]++] = id;
                }
            }

            Array.Clear(countsByScore);
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

            // Only Glycan-typed N-linked mods with generated Y ions contribute Y ions here.
            // Non-Glycan annotations (DB/GPTMD) and mod-loaded glycans without ions are intentionally
            // excluded; correct N-glycan ID is the NO-search's job.
            var nGlycan = peptideWithMod.AllModsOneIsNterminus.Values
                .OfType<Glycan>()
                .FirstOrDefault(g => g.Type == GlycanType.N_glycan && g.HasIons);
            if (nGlycan != null)
            {
                fragmentsForEachGlycoPeptide.AddRange(GlycoPeptides.GetGlycanYIons(theScan.PrecursorMass, nGlycan));
            }

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
                bool isIonTrapData = childScan.TheScan.MzAnalyzer == MZAnalyzerType.IonTrap2D || childScan.TheScan.MzAnalyzer == MZAnalyzerType.IonTrap3D;
                var matchedChildIons = MatchFragmentIons(childScan, childFragments, CommonParameters, isLowRes : isIonTrapData);

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

                var productTolerance = isIonTrapData ? CommonParameters.ProductMassTolerance_LowRes : CommonParameters.ProductMassTolerance;
                p += childScan.TheScan.MassSpectrum.Size * productTolerance.GetRange(1000).Width / childScan.TheScan.MassSpectrum.Range.Width;

            }

            var psmGlyco = new GlycoSpectralMatch(peptideWithMod, 0, PeptideScore, scanIndex, theScan, CommonParameters, matchedIons);

            //TO DO: This p is from childScan p, it works for HCD-pd-EThcD, which may not work for other type.
            psmGlyco.ScanInfo_p = p > 1? 1 : p;
            // p is the probability of randomly matching a single theoretical fragment ion.
            // With a wide mass tolerance (> 0.4 Da), the computed p may exceed 1,
            // which is not physically meaningful and can cause numerical issues.
            // In such cases, p is capped at 1.

            psmGlyco.Thero_n = n;

            psmGlyco.Rank = rank;

            psmGlyco.DiagnosticIonScore = DiagnosticIonScore;

            psmGlyco.GlycanScore = GlycanScore;

            psmGlyco.ChildMatchedFragmentIons = allMatchedChildIons;

            psmGlyco.LocalizationGraphs = localizationGraphs;

            if (oxoniumIonIntensities[OxoniumIonReservedIndices.R144] <= 0.00000001)
            {
                psmGlyco.R138vs144 = 100000000;
            }
            else
            {
                psmGlyco.R138vs144 = oxoniumIonIntensities[OxoniumIonReservedIndices.R138] / oxoniumIonIntensities[OxoniumIonReservedIndices.R144]; // if the ratio is high, that means the glycan is more likely to be N-glycan. Oppsitely, ration is small means close to O-glycan.
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


            int iDLow = GlycoPeptides.BinarySearchGetIndex(GlycanBoxMasses, possibleGlycanMassLow); // try to find the index that closet match to the "possibleGlycanMassLow" within the glycanBox

            // No glycan box fits the precursor: nothing below would run, so skip the site search and fragmentation altogether.
            if (iDLow >= GlycanBoxes.Length || !PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBoxes[iDLow].Mass))
            {
                return;
            }

            SortedDictionary<int, string> modPos = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, Motifs); //list all of the possible glycoslation site/postition

            // Where the protease says a glycan MUST be. Empty for every ordinary protease, and for a
            // glycoprotease it is derived from the peptide's own provenance rather than stored on it:
            // the cut that produced this peptide could not have happened unless a particular residue
            // carried a glycan. Derived, so it survives the index cache round trip and adds no field to
            // the serialized peptide.
            Dictionary<int, List<CleavageRequirement>> siteConditions = GetCleavageObligations(theScanBestPeptide);
            IReadOnlyCollection<int> obligatedSites = siteConditions.Keys;

            var localizationScan = theScan;
            var toleranceForLocalizationScan = CommonParameters.ProductMassTolerance;
            bool childScansCarryEtd = theScan.ChildScans.Count > 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType, CommonParameters.CustomIons);

            //For HCD-pd-ETD or CD-pd-EThcD type of data, we localize on the child scan.
            if (childScansCarryEtd)
            {
                localizationScan = GetLocalizationScan(theScan);
                // For the localization scan, if it is from ion trap, we will use a wider tolerance for the localization.
                toleranceForLocalizationScan = localizationScan.TheScan.MzAnalyzer == MZAnalyzerType.IonTrap2D ||
                    localizationScan.TheScan.MzAnalyzer == MZAnalyzerType.IonTrap3D ? CommonParameters.ProductMassTolerance_LowRes : CommonParameters.ProductMassTolerance;
            }

            //Localization for O-glycopeptides only works on ETD related dissociationtype
            //No localization can be done with MS2-HCD spectrum
            //TO THINK: there is a special situation. The HCD only scan from  HCD-pd-EThcD data can be a glycopeptide, but there is no ETD, so there is no localization. What to do with this?
            bool is_HCD_only_data = !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType, CommonParameters.CustomIons) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType, CommonParameters.CustomIons);

            // The theoretical fragments (product list) are only needed once a glycan box passes GraphCheck, so they are built then.
            List<Product> products = null;
            List<Product> GetProducts()
            {
                if (products != null)
                {
                    return products;
                }
                products = new List<Product>();

                //For HCD-pd-ETD or CD-pd-EThcD type of data, we generate the different rpoducts.
                if (childScansCarryEtd)
                {
                    theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
                }

                //For ETD type of data
                if (theScan.ChildScans.Count == 0 && GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType, CommonParameters.CustomIons))
                {
                    theScanBestPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
                }

                if (is_HCD_only_data) // In the HCD, there is no Y  ion, so we don't need to consider the modification here.
                {
                    theScanBestPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
                }
                return products;
            }

            double bestLocalizedScore = 0;

            List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>(); // if we also have ETD, then we will search the localization

            // These depend on the peptide alone, so every glycan box tried below shares them.
            string[] modPosMotifs = modPos.Values.ToArray();
            GlycoPeptides.SiteFragmentMasses[] siteFragments = null;

            while (iDLow < GlycanBoxes.Length && (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass + GlycanBoxes[iDLow].Mass))) // verify the glycan mass is invaild (within the range and match with mass shift)
            {
                if (OxoniumIonFilter && !GlycoPeptides.DiagonsticFilter(oxoniumIonIntensities, GlycanBoxes[iDLow])) // if the filter is turned on, we need to check does the oxoiums make sense.
                {
                    iDLow++; // if the oxonium ions don't make sense (there is no 204, or without their diagnostic ion), we can skip this glycan.
                    continue;
                }
                if (GraphCheck(modPosMotifs, GlycanBoxes[iDLow], obligatedSites.Count)) // the glycosite number should be larger than the possible glycan number.
                {
                    siteFragments ??= GlycoPeptides.GetSiteFragmentMasses(GetProducts(), modPos.Keys.ToArray());
                    LocalizationGraph localizationGraph = new LocalizationGraph(modPos, GlycanBoxes[iDLow], GlycanBoxes[iDLow].ChildGlycanBoxes, iDLow, obligatedSites, siteConditions);
                    LocalizationGraph.LocalizeOGlycan(localizationGraph, localizationScan, toleranceForLocalizationScan, GetProducts(), siteFragments); //create the localization graph with the glycan mass and the possible glycosite.

                    // No arrangement of this box satisfies the obligation, so it explains nothing about
                    // this spectrum. Skipping here rather than scoring it zero matters: a zero-scoring
                    // graph can still be the best one when every box fails, and GetFirstPath would then
                    // walk a graph that has no path.
                    if (!localizationGraph.HasRoute)
                    {
                        iDLow++;
                        continue;
                    }

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

        /// <summary>
        /// Chooses which child scan is used as the localization spectrum.
        /// Glycosite localization is scored with c/zDot ions, so it is only meaningful against an
        /// electron-based child scan. The child scans of a precursor are ordered by scan number, so
        /// simply taking the first one selects the lowest-numbered child whatever its activation. On an
        /// acquisition that places more than one activation on the same precursor (for example
        /// HCD-ETciD-CID) that can select the collision-based child and score c/zDot theoretical ions
        /// against a spectrum which cannot contain them.
        /// The child scan header is consulted first, and the previous first-child behavior is kept as the
        /// fallback when no child reports a usable activation. Files with a single child scan, which is
        /// every acquisition shape this method was originally written for, are unaffected either way.
        /// </summary>
        /// <param name="theScan">The parent scan whose children are candidates. Must have at least one child.</param>
        /// <returns>The child scan to localize against.</returns>
        private Ms2ScanWithSpecificMass GetLocalizationScan(Ms2ScanWithSpecificMass theScan)
        {
            foreach (var childScan in theScan.ChildScans)
            {
                DissociationType? childDissociationType = childScan.TheScan.DissociationType;

                // An absent or Autodetect header carries no activation information, so it cannot be trusted here.
                if (childDissociationType == null || childDissociationType == DissociationType.Autodetect)
                {
                    continue;
                }

                if (GlycoPeptides.DissociationTypeContainETD(childDissociationType.Value, CommonParameters.CustomIons))
                {
                    return childScan;
                }
            }

            // No child advertised an electron-based activation. Preserve the historical behavior rather than
            // dropping localization, because the declared MS2ChildScanDissociationType already asserted ETD.
            return theScan.ChildScans.First();
        }

        private void FindNGlycan(Ms2ScanWithSpecificMass theScan, int scanIndex, int scoreCutOff, PeptideWithSetModifications theScanBestPeptide, int ind, double possibleGlycanMassLow, double[] oxoniumIonIntensities, ref List<GlycoSpectralMatch> possibleMatches)
        {
            List<int> modPos_Nxs = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxs" }).Select(p => p.Key).ToList();
            List<int> modPos_Nxt = GlycoSpectralMatch.GetPossibleModSites(theScanBestPeptide, new string[] { "Nxt" }).Select(p => p.Key).ToList();
            if (modPos_Nxs.Count < 1 && modPos_Nxt.Count < 1) // if there is no possible glycosylation site, we can skip this peptide.
            {
                return;
            }

            int iDLow = GlycoPeptides.BinarySearchGetIndex(NGlycanMasses, possibleGlycanMassLow);
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

                    if (oxoniumIonIntensities[OxoniumIonReservedIndices.R144] <= 0.00000001)
                    {
                        psmGlyco.R138vs144 = 100000000;
                    }
                    else
                    {
                        psmGlyco.R138vs144 = oxoniumIonIntensities[OxoniumIonReservedIndices.R138] / oxoniumIonIntensities[OxoniumIonReservedIndices.R144];
                    }

                    possibleMatches.Add(psmGlyco);
                }

                iDLow++;
            }
        }

        // Conduct the search and generate the gsms for N-glycan search
        private List<GlycoSpectralMatch> MatchNGlycopeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff, List<int> ranks = null)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();
            double[] oxoniumIonIntensities = null;

            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]];
                int rank = ranks?[ind] ?? ind;

                //Considering coisolation, it doesn't mean it must from a glycopeptide even the scan contains oxonium ions.
                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass))
                {
                    FindSingle(theScan, scanIndex, scoreCutOff, theScanBestPeptide, rank, ref possibleMatches);
                }
                else
                {
                    //Filter by glycanBoxes mass difference.
                    var possibleGlycanMassLow = PrecusorSearchMode.GetMinimumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    var possibleGlycanMassHigh = PrecusorSearchMode.GetMaximumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    if (NGlycans.Length == 0 || possibleGlycanMassHigh < (double)NGlycans.First().Mass/1E5 || possibleGlycanMassLow > (double)NGlycans.Last().Mass/1E5)
                    {
                        continue;
                    }

                    //Filter by OxoniumIon. The intensities depend on the scan alone, so they are read once per scan.
                    oxoniumIonIntensities ??= GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode);

                    //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    if (OxoniumIonFilter && oxoniumIonIntensities[OxoniumIon204Index] == 0)
                    {
                        continue;
                    }

                    //Find N-Glycan 
                    FindNGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, rank, possibleGlycanMassLow, oxoniumIonIntensities, ref possibleMatches);

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
        private List<GlycoSpectralMatch> MatchGlycopeptide(Ms2ScanWithSpecificMass theScan, List<int> idsOfPeptidesPossiblyObserved, int scanIndex, int scoreCutOff, List<int> ranks = null)
        {
            List<GlycoSpectralMatch> possibleMatches = new List<GlycoSpectralMatch>();
            double[] oxoniumIonIntensities = null;


            for (int ind = 0; ind < idsOfPeptidesPossiblyObserved.Count; ind++)
            {
                var theScanBestPeptide = PeptideIndex[idsOfPeptidesPossiblyObserved[ind]]; // Get the peptide from the candidate list.
                int rank = ranks?[ind] ?? ind; // position in the TopN cut, which for a partitioned search spans every partition

                if (PrecusorSearchMode.Within(theScan.PrecursorMass, theScanBestPeptide.MonoisotopicMass)) // If the peptide mass is indentical to the precursor mass (or within the tolerance), we can directly search the glycopeptide.
                {
                    FindSingle(theScan, scanIndex, scoreCutOff, theScanBestPeptide, rank, ref possibleMatches);
                }
                else if (theScan.PrecursorMass - theScanBestPeptide.MonoisotopicMass >= 100) //If not, we need to consider the glycan mass difference.
                {
                    //Filter by glycanBoxes mass difference.
                    var possibleGlycanMassLow = PrecusorSearchMode.GetMinimumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    var possibleGlycanMassHigh = PrecusorSearchMode.GetMaximumValue(theScan.PrecursorMass) - theScanBestPeptide.MonoisotopicMass;

                    if (possibleGlycanMassHigh < GlycanBoxes.First().Mass || possibleGlycanMassLow > GlycanBoxes.Last().Mass)
                    {
                        continue; // if the glycan mass difference is out of the range of the glycan box, we can skip this peptide.
                    }

                    //Filter by OxoniumIon. The intensities depend on the scan alone, so they are read once per scan.
                    oxoniumIonIntensities ??= GlycoPeptides.ScanOxoniumIonFilter(theScan, ProductSearchMode);

                    //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons (the [9] is 204). A spectrum needs to have 204.0867 to be considered as a glycopeptide for now.
                    if (OxoniumIonFilter && oxoniumIonIntensities[OxoniumIon204Index] == 0)
                    {
                        continue;
                    }

                    //Find O-Glycan
                    FindOGlycan(theScan, scanIndex, scoreCutOff, theScanBestPeptide, rank, possibleGlycanMassLow, oxoniumIonIntensities, ref possibleMatches);
                }
            }

            // One stable sort after the loop orders the matches exactly as re-sorting after every candidate did.
            if (possibleMatches.Count != 0)
            {
                possibleMatches = possibleMatches.OrderByDescending(p => p.Score).ToList();
            }

            return possibleMatches;
        }

        /// <summary>
        /// Valid the Graph created by this modPos and glycanBox.
        /// Check if the motif in peptide is sufficient to cover the motif in glycanBox.
        /// </summary>
        /// <param name="modPosMotifs"> The motif at each candidate glycosite of the peptide. </param>
        /// <param name="glycanBox"></param>
        /// <returns></returns>
        /// <summary>
        /// The obligated sites this peptide's protease implies, in the same two-based key space as
        /// <see cref="GlycoSpectralMatch.GetPossibleModSites"/>. Empty whenever the peptide has no
        /// digestion provenance or its protease requires nothing.
        /// </summary>
        private static Dictionary<int, List<CleavageRequirement>> GetCleavageObligations(PeptideWithSetModifications peptide)
        {
            DigestionAgent agent = peptide?.DigestionParams?.DigestionAgent;
            if (agent == null || !agent.HasCleavageRequirement)
            {
                return new Dictionary<int, List<CleavageRequirement>>();
            }

            return peptide.GetCleavageObligations(agent);
        }

        private static bool GraphCheck(string[] modPosMotifs, GlycanBox glycanBox, int obligatedSiteCount)
        {
            // If the motifs number is less than the glycanBox, we can skip this graph.
            if (modPosMotifs.Length < glycanBox.NumberOfMods)
                return false;

            // A box with fewer glycans than the peptide has obligated sites cannot fill them all, so no
            // route through the graph would survive. Screening here matters rather than merely saving
            // work: FinishGraph throws when the terminal node is unreachable, and an unsatisfiable
            // obligation is exactly the case its own comment says never happens in a search.
            if (glycanBox.NumberOfMods < obligatedSiteCount)
                return false;

            // Check if the motif in peptide is sufficient to cover the motif in glycanBox.
            return glycanBox.GetMotifCount().CoveredBy(modPosMotifs, modPosMotifs.Length - 1);
        }

    }
}

