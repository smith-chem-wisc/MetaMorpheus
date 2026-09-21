using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using EngineLayer.Indexing;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using EngineLayer.Util;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class GlycoSearchTask : MetaMorpheusTask
    {
        public GlycoSearchTask() : base(MyTask.GlycoSearch)
        {
            //Default parameter setting which is different from SearchTask, can be overwriten
            var digestPara = new DigestionParams(
                minPeptideLength: 5,
                maxPeptideLength: 60

            );
            CommonParameters = new CommonParameters(
                precursorMassTolerance: new PpmTolerance(10),
                ms2childScanDissociationType: DissociationType.EThcD,
                scoreCutoff: 3,
                trimMsMsPeaks: false, 
                numberOfPeaksToKeepPerWindow: 1000,
                minimumAllowedIntensityRatioToBasePeak: 0.01,
                digestionParams: digestPara            
                
            );

            _glycoSearchParameters = new GlycoSearchParameters();
        }

        public GlycoSearchParameters _glycoSearchParameters { get; set; }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            List<List<GlycoSpectralMatch>> ListOfGsmsPerMS2Scan = new List<List<GlycoSpectralMatch>>();

            LoadModifications(taskId, out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);

            // load proteins
            var dbLoader = new DatabaseLoadingEngine(CommonParameters, this.FileSpecificParameters, [taskId], dbFilenameList, taskId, _glycoSearchParameters.DecoyType, true, localizeableModificationTypes);
            var loadingResults = dbLoader.Run() as DatabaseLoadingEngineResults;
            List<Protein> proteinList = loadingResults!.BioPolymers.Cast<Protein>().ToList();

            MyFileManager myFileManager = new (_glycoSearchParameters.DisposeOfFileWhenDone);
            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));

            int completedFiles = 0;

            Status("Searching files...", taskId);
            ProseCreatedWhileRunning.Append("\n");
            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.DigestionAgent + "; \n");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; \n");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinLength + "; \n");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxLength + "; \n");
            if (CommonParameters.DigestionParams is DigestionParams digestionParams)
                ProseCreatedWhileRunning.Append("initiator methionine behavior = " + digestionParams.InitiatorMethionineBehavior + "; \n");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; \n");
            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif)) + "; \n");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; \n");
            ProseCreatedWhileRunning.Append("parent mass tolerance(s) = " + CommonParameters.PrecursorMassTolerance + "; \n");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; \n");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. \n");
            if (_glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch)
            {
                ProseCreatedWhileRunning.Append("The O-glycan database: " + _glycoSearchParameters.OGlycanDatabasefile + "\n");
            }
            else if (_glycoSearchParameters.GlycoSearchType == GlycoSearchType.NGlycanSearch)
            {
                ProseCreatedWhileRunning.Append("The N-glycan database: " + _glycoSearchParameters.OGlycanDatabasefile + "\n");
            }
            else
            {
                ProseCreatedWhileRunning.Append("The O-glycan database: " + _glycoSearchParameters.OGlycanDatabasefile + "\n");
                ProseCreatedWhileRunning.Append("The N-glycan database: " + _glycoSearchParameters.NGlycanDatabasefile + "\n");
            }                
            
            ProseCreatedWhileRunning.Append("\n");

            FlashLfqResults flashLfqResults = null;
            int? decidedPartitions = null;
            int fileCount = currentRawFileList.Count;

            // Every file's parameters, built here on one thread: the CommonParameters constructor reads and clears a process-wide
            // list of custom product types, so building them from parallel searches would race.
            var combinedParamsPerFile = new CommonParameters[fileCount];
            var indexParamsPerFile = new CommonParameters[fileCount];
            for (int spectraFileIndex = 0; spectraFileIndex < fileCount; spectraFileIndex++)
            {
                combinedParamsPerFile[spectraFileIndex] = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                // scoped to indexing/searching only, so the settings the task reports stay as configured.
                // every TotalPartitions read below comes from indexParams, including the loop bound and the
                // protein-range slicing, which previously read from two different objects.
                indexParamsPerFile[spectraFileIndex] = RaisePartitionsToFitMemory(proteinList, combinedParamsPerFile[spectraFileIndex], fixedModifications,
                    variableModifications, null, null, null, 30000.0, ref decidedPartitions);
            }

            // The glycans and glycan boxes depend on the glycan databases and settings alone, so they are built once for every file.
            Status("Building glycan boxes...", taskId);
            GlycanSearchSpace glycanSearchSpace = GlycanSearchSpace.Build(_glycoSearchParameters.OGlycanDatabasefile, _glycoSearchParameters.NGlycanDatabasefile,
                _glycoSearchParameters.GlycoSearchType, _glycoSearchParameters.MaximumOGlycanAllowed, _glycoSearchParameters.MaximumGlycanBoxMass);

            var thisIds = currentRawFileList.Select(file => new List<string> { taskId, "Individual Spectra Files", file }).ToList();
            var databaseFileInfos = dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList();

            Ms2ScanWithSpecificMass[] LoadMs2Scans(int spectraFileIndex, CommonParameters loadParams)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                Status("Loading spectra file...", thisIds[spectraFileIndex]);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, loadParams);

                Status("Getting ms2 scans...", thisIds[spectraFileIndex]);
                Ms2ScanWithSpecificMass[] scans = GetMs2Scans(myMsDataFile, origDataFile, loadParams).OrderBy(b => b.PrecursorMass).ToArray();

                myFileManager.DoneWithFile(origDataFile);
                return scans;
            }

            // Files can be searched at once only when they share one unpartitioned index: every file's index-relevant settings must
            // agree, since a file-specific toml can change the protease, lengths, missed cleavages, mods or dissociation type.
            string sequentialReason = null;
            if (fileCount < 2)
            {
                sequentialReason = "there is one spectra file";
            }
            else if (indexParamsPerFile[0].TotalPartitions > 1)
            {
                sequentialReason = "the index was split into " + indexParamsPerFile[0].TotalPartitions + " partitions to fit in memory";
            }
            else if (indexParamsPerFile.Select(p => new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 0, _glycoSearchParameters.DecoyType, p, this.FileSpecificParameters, 30000.0, false, databaseFileInfos, TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId }).ToString()
                         + "\nCustom ions: " + string.Join(",", p.CustomIons ?? new List<Omics.Fragmentation.ProductType>())).Distinct().Count() > 1)
            {
                sequentialReason = "file-specific settings give the spectra files different indexes";
            }
            else if (_glycoSearchParameters.MaximumSpectraFilesInParallel == 1)
            {
                sequentialReason = "the task setting MaximumSpectraFilesInParallel = 1";
            }

            if (sequentialReason == null)
            {
                for (int spectraFileIndex = 0; spectraFileIndex < fileCount; spectraFileIndex++)
                {
                    NewCollection(Path.GetFileName(currentRawFileList[spectraFileIndex]), thisIds[spectraFileIndex]);
                }

                // One unpartitioned index for every file, so each file is searched in a single pass; the two rounds below are for partitions.
                Status("Getting fragment dictionary...", new List<string> { taskId });
                List<PeptideWithSetModifications> peptideIndex = null;
                FragmentIndex fragmentIndex = null;
                List<int>[] precursorIndex = null;
                var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 0, _glycoSearchParameters.DecoyType, indexParamsPerFile[0], this.FileSpecificParameters, 30000.0, false, databaseFileInfos, TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });
                GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);

                // Load the first file with the whole thread budget, and take what its scans hold in memory as the measure of a file.
                long heapBeforeFirstFile = GC.GetTotalMemory(forceFullCollection: true);
                var scansPerFile = new Ms2ScanWithSpecificMass[fileCount][];
                scansPerFile[0] = LoadMs2Scans(0, combinedParamsPerFile[0]);
                long scanBytes = Math.Max(0, GC.GetTotalMemory(forceFullCollection: true) - heapBeforeFirstFile);

                // A file being searched also holds its results and, per thread, two scoring tables the size of the index; doubling the
                // scans covers the results, and the tables are counted at the whole budget since the split is not known yet.
                int threadBudget = CommonParameters.MaxThreadsToUsePerFile;
                long bytesPerFile = 2 * scanBytes + 2L * Math.Max(1, threadBudget) * peptideIndex.Count;
                long freeBytes = IndexPartitioning.AvailableBytes();
                FileParallelismPlan plan = FileParallelism.Decide(fileCount, threadBudget, freeBytes, bytesPerFile, _glycoSearchParameters.MaximumSpectraFilesInParallel);

                // The split above sets how many files run at once and how many threads load each; the search itself draws on one budget
                // shared by every file, so threads move from files that finish to files still searching.
                var searchThreads = new SearchThreadBudget(threadBudget);
                ProseCreatedWhileRunning.Append("spectra files searched in parallel = " + plan.FilesInParallel + " of " + fileCount + " (limited by " + plan.LimitedBy + ")"
                    + ", each loaded with " + plan.ThreadsPerFile + " threads and searched with threads shared from a budget of " + searchThreads.TotalThreads
                    + " that move from files that finish to files still searching; \n");
                ProseCreatedWhileRunning.Append("memory free after building the index and loading the first spectra file = " + (freeBytes / 1e9).ToString("0.0") + " GB, estimated "
                    + (bytesPerFile / 1e9).ToString("0.0") + " GB for each further file; \n");

                // Each file's own copies of its parameters with the divided thread count, built here on one thread (see above).
                var loadParamsPerFile = combinedParamsPerFile.Select(p => p.CloneWithNewMaxThreadsToUsePerFile(plan.ThreadsPerFile)).ToArray();
                var searchParamsPerFile = indexParamsPerFile.Select(p => p.CloneWithNewMaxThreadsToUsePerFile(plan.ThreadsPerFile)).ToArray();

                var gsmsPerFile = new List<GlycoSpectralMatch>[fileCount][];
                Status("Searching files...", taskId);
                FileParallelism.ForEachFile(fileCount, plan.FilesInParallel, spectraFileIndex =>
                {
                    if (GlobalVariables.StopLoops) { return; }

                    Ms2ScanWithSpecificMass[] scans = scansPerFile[spectraFileIndex] ?? LoadMs2Scans(spectraFileIndex, loadParamsPerFile[spectraFileIndex]);
                    scansPerFile[spectraFileIndex] = null;

                    var gsms = new List<GlycoSpectralMatch>[scans.Length];
                    new GlycoSearchEngine(gsms, scans, peptideIndex, fragmentIndex, null, 0, searchParamsPerFile[spectraFileIndex], this.FileSpecificParameters, glycanSearchSpace,
                        _glycoSearchParameters.OGlycanDatabasefile, _glycoSearchParameters.NGlycanDatabasefile, _glycoSearchParameters.GlycoSearchTopNum, _glycoSearchParameters.MaximumOGlycanAllowed, _glycoSearchParameters.OxoniumIonFilt, thisIds[spectraFileIndex])
                    {
                        ThreadBudget = searchThreads
                    }.Run();
                    gsmsPerFile[spectraFileIndex] = gsms;

                    ReportProgress(new ProgressEventArgs(100, "Done with search 1/1!", thisIds[spectraFileIndex]));
                    int done = Interlocked.Increment(ref completedFiles);
                    ReportProgress(new ProgressEventArgs(100 * done / fileCount, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
                });

                // Results in file order, exactly as searching the files one after another adds them.
                foreach (var gsms in gsmsPerFile.Where(p => p != null))
                {
                    ListOfGsmsPerMS2Scan.AddRange(gsms.Where(p => p != null).ToList());
                }
            }
            else
            {
                ProseCreatedWhileRunning.Append("spectra files searched one after another because " + sequentialReason + "; \n");

                for (int spectraFileIndex = 0; spectraFileIndex < fileCount; spectraFileIndex++)
                {
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    CommonParameters combinedParams = combinedParamsPerFile[spectraFileIndex];

                    var thisId = thisIds[spectraFileIndex];
                    NewCollection(Path.GetFileName(origDataFile), thisId);

                    Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = LoadMs2Scans(spectraFileIndex, combinedParams);

                    List<GlycoSpectralMatch>[] newCsmsPerMS2ScanPerFile = new List<GlycoSpectralMatch>[arrayOfMs2ScansSortedByMass.Length];

                    CommonParameters indexParams = indexParamsPerFile[spectraFileIndex];

                    IndexingEngine MakeIndexEngine(int currentPartition)
                    {
                        //When partition, the proteinList will be split for each Thread.
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / indexParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / indexParams.TotalPartitions) - (currentPartition * proteinList.Count() / indexParams.TotalPartitions));

                        //Only reverse Decoy for glyco search has been tested and are set as fixed parameter.
                        return new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null, currentPartition, _glycoSearchParameters.DecoyType, indexParams, this.FileSpecificParameters, 30000.0, false, databaseFileInfos, TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });
                    }

                    GlycoSearchEngine MakeSearchEngine(List<PeptideWithSetModifications> peptideIndex, FragmentIndex fragmentIndex, int currentPartition, List<(int Partition, int PeptideId, byte Score)>[] candidates)
                    {
                        //The second Fragment index is for 'MS1-HCD_MS1-ETD_MS2s' type of data. If LowCID is used for MS1, ion-index is not allowed to use.
                        FragmentIndex secondFragmentIndex = null;

                        return new GlycoSearchEngine(newCsmsPerMS2ScanPerFile, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, secondFragmentIndex, currentPartition, indexParams, this.FileSpecificParameters, glycanSearchSpace,
                            _glycoSearchParameters.OGlycanDatabasefile, _glycoSearchParameters.NGlycanDatabasefile, _glycoSearchParameters.GlycoSearchTopNum, _glycoSearchParameters.MaximumOGlycanAllowed, _glycoSearchParameters.OxoniumIonFilt, thisId, candidates);
                    }

                    // The TopN candidate cut has to be taken over the whole database, or which glycopeptides are reported
                    // depends on the partition count -- and RaisePartitionsToFitMemory can raise that count on its own.
                    // So a partitioned search scores every partition first and glycan-matches afterwards, the same two
                    // rounds XLSearchTask uses. One partition keeps the single pass, which is the same cut already.
                    List<(int Partition, int PeptideId, byte Score)>[] candidates = indexParams.TotalPartitions > 1
                        ? new List<(int Partition, int PeptideId, byte Score)>[arrayOfMs2ScansSortedByMass.Length]
                        : null;

                    // Round 2 has only peptide ids from round 1, so it must see each partition's index in the same order.
                    // Recorded per partition in round 1 and checked in round 2.
                    int[] peptideOrderFingerprints = candidates == null ? null : new int[indexParams.TotalPartitions];

                    for (int currentPartition = 0; currentPartition < indexParams.TotalPartitions; currentPartition++)
                    {
                        List<PeptideWithSetModifications> peptideIndex = null;

                        Status("Getting fragment dictionary...", new List<string> { taskId });

                        var indexEngine = MakeIndexEngine(currentPartition);
                        FragmentIndex fragmentIndex = null;
                        List<int>[] precursorIndex = null;
                        GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);

                        Status("Searching files...", taskId);
                        var glycoSearchEngine = MakeSearchEngine(peptideIndex, fragmentIndex, currentPartition, candidates);
                        if (candidates == null)
                        {
                            glycoSearchEngine.Run();
                        }
                        else
                        {
                            glycoSearchEngine.FirstRoundSearch();
                            peptideOrderFingerprints[currentPartition] = IndexingEngine.PeptideOrderFingerprint(peptideIndex);
                        }

                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + indexParams.TotalPartitions + "!", thisId));
                        if (GlobalVariables.StopLoops) { break; }
                    }

                    // glycan matching, one partition at a time, for the candidates that survived the whole-database cut
                    for (int currentPartition = 0; candidates != null && currentPartition < indexParams.TotalPartitions; currentPartition++)
                    {
                        if (GlobalVariables.StopLoops) { break; }

                        Status("Getting peptide index...", new List<string> { taskId });

                        var indexEngine = MakeIndexEngine(currentPartition);
                        List<PeptideWithSetModifications> peptideIndex = null;
                        List<int>[] precursorIndex = null;
                        GenerateIndexes_PeptideOnly(indexEngine, dbFilenameList, ref peptideIndex, ref precursorIndex, proteinList, taskId);
                        if (peptideIndex == null)
                        {
                            // the first round's index could not be read back, so build it again
                            FragmentIndex unusedFragmentIndex = null;
                            GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref unusedFragmentIndex, ref precursorIndex, proteinList, taskId);
                        }

                        // A rebuild, or a different cache folder, can hold the same peptides in a different order: mass ties
                        // are broken by digestion order, which depends on the thread count, and the thread count is not part
                        // of the cache key. The candidates would then point at other peptides and report wrong
                        // glycopeptides without any error, so stop instead.
                        if (IndexingEngine.PeptideOrderFingerprint(peptideIndex) != peptideOrderFingerprints[currentPartition])
                        {
                            throw new MetaMorpheusException($"The peptide index for partition {currentPartition + 1} of {indexParams.TotalPartitions} " +
                                "changed between the two rounds of the glyco search, so the first round's candidates no longer point at the right " +
                                $"peptides. Delete the {IndexFolderName} folder next to the database and run the search again.");
                        }

                        Status("Searching files...", taskId);
                        MakeSearchEngine(peptideIndex, null, currentPartition, candidates).Run();

                        ReportProgress(new ProgressEventArgs(100, "Done with glycan matching " + (currentPartition + 1) + "/" + indexParams.TotalPartitions + "!", thisId));
                    }

                    ListOfGsmsPerMS2Scan.AddRange(newCsmsPerMS2ScanPerFile.Where(p => p != null).ToList());

                    completedFiles++;
                    ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
                }
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            //For every Ms2Scans, each have a list of candidates psms. The allPsms from GlycoSearchEngine is the list (all ms2scans) of list (each ms2scan) of psm (all candidate psm). 
            //Currently, only keep the first scan for consideration. 
            List<GlycoSpectralMatch> GsmPerScans = ListOfGsmsPerMS2Scan.Select(p => p.First()).ToList();

            var filteredAllPsms = new List<GlycoSpectralMatch>();

            //For each ms2scan, try to find the best candidate psm from the psms list. Do the localizaiton analysis. Add it into filteredAllPsms.
            foreach (var gsmsPerScan in GsmPerScans.GroupBy(p => (p.ScanNumber, p.FullFilePath)))
            {
                var glycos = RemoveSimilarSequenceDuplicates(gsmsPerScan.OrderByDescending(p=>p.Score).ToList());

                foreach (var glycoSpectralMatch in glycos)
                {
                    if (glycoSpectralMatch.LocalizationGraphs != null)
                    {

                        List<Route> localizationCandidates = new List<Route>();

                        for (int i = 0; i < glycoSpectralMatch.LocalizationGraphs.Count; i++)
                        {
                            var allPathWithMaxScore = LocalizationGraph.GetAllHighestScorePaths(glycoSpectralMatch.LocalizationGraphs[i].array, glycoSpectralMatch.LocalizationGraphs[i].ChildModBoxes);

                            foreach (var path in allPathWithMaxScore)
                            {
                                var local = LocalizationGraph.GetLocalizedPath(glycoSpectralMatch.LocalizationGraphs[i], path);
                                local.ModBoxId = glycoSpectralMatch.LocalizationGraphs[i].ModBoxId;
                                localizationCandidates.Add(local);
                            }
                        }

                        glycoSpectralMatch.Routes = localizationCandidates;

                    }

                    if (glycoSpectralMatch.Routes != null)
                    {
                        LocalizationLevel localLevel;
                        glycoSpectralMatch.LocalizedGlycan = GlycoSpectralMatch.GetLocalizedGlycan(glycoSpectralMatch.Routes, out localLevel);
                        glycoSpectralMatch.LocalizationLevel = localLevel;

                        //Localization PValue.
                        if (localLevel == LocalizationLevel.Level1 || localLevel == LocalizationLevel.Level2)
                        {
                            List<Route> allRoutes = new List<Route>();
                            foreach (var graph in glycoSpectralMatch.LocalizationGraphs)
                            {
                                allRoutes.AddRange(LocalizationGraph.GetAllPaths_CalP(graph, glycoSpectralMatch.ScanInfo_p, glycoSpectralMatch.Thero_n));
                            }
                            glycoSpectralMatch.ModSitePairProbDict = LocalizationGraph.CalProbabilityForModSitePair(allRoutes, glycoSpectralMatch.LocalizedGlycan);
                        }
                    }

                    filteredAllPsms.Add(glycoSpectralMatch);
                }
            }

            // Localization is done, so let go of what the search built on the glycan boxes before FDR and PEP (see ReleaseSearchCaches).
            GlycanBox.ReleaseSearchCaches(GlycanBox.OGlycanBoxes);
            GlycanBox.ReleaseSearchCaches(GlycanBox.NOGlycanBoxes);

            PostGlycoSearchAnalysisParameters pgsap = new()
            {
                GlycoSearchTaskResults = MyTaskResults,
                SearchTaskId = taskId,
                GlycoSearchParameters = _glycoSearchParameters,
                ProteinList = proteinList,
                VariableModifications = variableModifications,
                FixedModifications = fixedModifications,
                AllPsms = filteredAllPsms.OrderByDescending(p => p.Score).ToList(),
                OutputFolder = OutputFolder,
                IndividualResultsOutputFolder = Path.Combine(OutputFolder, "IndividualFileResults"),
                FlashLfqResults = flashLfqResults,
                FileSettingsList = fileSettingsList,
                DatabaseFilenameList = dbFilenameList,
                CurrentRawFileList = currentRawFileList,
                BioPolymerList = proteinList,
            };
            
            PostGlycoSearchAnalysisTask postGlycoSearchAnalysisTask = new PostGlycoSearchAnalysisTask()
            {
                Parameters = pgsap,
                FileSpecificParameters = this.FileSpecificParameters,
                CommonParameters = this.CommonParameters
            };

            return postGlycoSearchAnalysisTask.Run(OutputFolder, dbFilenameList, currentRawFileList, taskId, fileSettingsList, filteredAllPsms.OrderByDescending(p => p.Score).ToList(), CommonParameters, _glycoSearchParameters, proteinList, variableModifications, fixedModifications, localizeableModificationTypes, MyTaskResults);

        }

        //The coisolation works for general search doesn't work for glyco search workflow. Similar peptide with different glycan are identified because of poor precursor mass. 
        //glycoSpectralMatches must be OrderDecendingByScore.
        private static List<GlycoSpectralMatch> RemoveSimilarSequenceDuplicates(List<GlycoSpectralMatch> glycoSpectralMatches)
        {
            List<GlycoSpectralMatch> glycos = new List<GlycoSpectralMatch>();
            glycos.Add(glycoSpectralMatches.First());
            foreach (var g in glycoSpectralMatches)
            {
                if (glycoSpectralMatches.First().BaseSequence.Contains(g.BaseSequence) || g.BaseSequence.Contains(glycoSpectralMatches.First().BaseSequence))
                {
                    continue;
                }

                glycos.Add(g);
            }
            return glycos;
        }

    }
}
