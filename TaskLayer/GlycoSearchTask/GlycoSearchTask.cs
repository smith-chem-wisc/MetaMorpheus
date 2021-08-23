using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.GlycoSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;
using EngineLayer.FdrAnalysis;
using System;

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
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, true, _glycoSearchParameters.DecoyType, localizeableModificationTypes, CommonParameters);

            MyFileManager myFileManager = new MyFileManager(true);

            int completedFiles = 0;

            Status("Searching files...", taskId);
            ProseCreatedWhileRunning.Append("\n");
            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; \n");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; \n");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; \n");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; \n");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; \n");
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

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);

                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams);

                Status("Getting ms2 scans...", thisId);

                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy(b => b.PrecursorMass).ToArray();
                
                List<GlycoSpectralMatch>[] newCsmsPerMS2ScanPerFile = new List<GlycoSpectralMatch>[arrayOfMs2ScansSortedByMass.Length];
                
                myFileManager.DoneWithFile(origDataFile);

                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<PeptideWithSetModifications> peptideIndex = null;

                    //When partition, the proteinList will be split for each Thread.
                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                    Status("Getting fragment dictionary...", new List<string> { taskId });

                    //Only reverse Decoy for glyco search has been tested and are set as fixed parameter.
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null, currentPartition, _glycoSearchParameters.DecoyType, combinedParams, this.FileSpecificParameters, 30000.0, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });
                    List<int>[] fragmentIndex = null;
                    List<int>[] precursorIndex = null;
                    GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);

                    //The second Fragment index is for 'MS1-HCD_MS1-ETD_MS2s' type of data. If LowCID is used for MS1, ion-index is not allowed to use.
                    List<int>[] secondFragmentIndex = null;
                    //if (combinedParams.MS2ChildScanDissociationType != DissociationType.LowCID
                    //&& !CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(combinedParams.DissociationType, combinedParams.MS2ChildScanDissociationType))
                    //{
                    //    //Becuase two different type of dissociation methods are used, the parameters are changed with different dissociation type.
                    //    var secondCombinedParams = CommonParameters.CloneWithNewDissociationType(combinedParams.MS2ChildScanDissociationType);
                    //    var secondIndexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null, currentPartition, _glycoSearchParameters.DecoyType, secondCombinedParams, this.FileSpecificParameters, 30000.0, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), new List<string> { taskId });
                    //    GenerateSecondIndexes(indexEngine, secondIndexEngine, dbFilenameList, ref secondFragmentIndex, proteinList, taskId);
                    //}

                    Status("Searching files...", taskId);
                    new GlycoSearchEngine(newCsmsPerMS2ScanPerFile, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, secondFragmentIndex, currentPartition, combinedParams, this.FileSpecificParameters,
                        _glycoSearchParameters.OGlycanDatabasefile, _glycoSearchParameters.NGlycanDatabasefile, _glycoSearchParameters.GlycoSearchType, _glycoSearchParameters.GlycoSearchTopNum, _glycoSearchParameters.MaximumOGlycanAllowed, _glycoSearchParameters.OxoniumIonFilt, thisId).Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                    if (GlobalVariables.StopLoops) { break; }
                }

                ListOfGsmsPerMS2Scan.AddRange(newCsmsPerMS2ScanPerFile.Where(p => p != null).ToList());

                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            //For every Ms2Scans, each have a list of candidates psms. The allPsms from GlycoSearchEngine is the list (all ms2scans) of list (each ms2scan) of psm (all candidate psm). 
            //Currently, only keep the first scan for consideration. 
            List<GlycoSpectralMatch> GsmPerScans = ListOfGsmsPerMS2Scan.Select(p => p.First()).ToList();

            var filteredAllPsms = new List<GlycoSpectralMatch>();

            //For each ms2scan, try to find the best candidate psm from the psms list. Do the localizaiton analysis. Add it into filteredAllPsms.
            foreach (var gsmsPerScan in GsmPerScans.GroupBy(p => p.ScanNumber))
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
                            glycoSpectralMatch.SiteSpeciLocalProb = LocalizationGraph.CalSiteSpecificLocalizationProbability(allRoutes, glycoSpectralMatch.LocalizationGraphs.First().ModPos);
                        }
                    }

                    filteredAllPsms.Add(glycoSpectralMatch);
                }
            }

            PostGlycoSearchAnalysisTask postGlycoSearchAnalysisTask = new PostGlycoSearchAnalysisTask();

            postGlycoSearchAnalysisTask.FileSpecificParameters = this.FileSpecificParameters;
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
