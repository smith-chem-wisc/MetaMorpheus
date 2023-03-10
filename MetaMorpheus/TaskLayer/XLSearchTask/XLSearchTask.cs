using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public class XLSearchTask : MetaMorpheusTask
    {
        public XLSearchTask() : base(MyTask.XLSearch)
        {
            //Default parameter setting which is different from SearchTask, can be overwriten here.
            //maxMissedCleavage is set to 3 for crosslink generally induce long peptides.
            //maxPeptideLength is set to 60 in case generate peptides that are too long.
            //numberOfPeaksToKeepPerWindow is set to a large number to keep all peak with intensity > 0.01.
            var digestPara = new DigestionParams(
                minPeptideLength: 5,
                maxPeptideLength: 60,
                maxMissedCleavages: 3
            );
            CommonParameters = new CommonParameters(
                precursorMassTolerance: new PpmTolerance(10),
                scoreCutoff: 3,
                numberOfPeaksToKeepPerWindow: 200,
                //addCompIons: true,
                digestionParams: digestPara
            );

            XlSearchParameters = new XlSearchParameters();
        }

        public XlSearchParameters XlSearchParameters { get; set; }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            List<List<CrosslinkSpectralMatch>> ListOfCsmsPerMS2Scan = new List<List<CrosslinkSpectralMatch>>();

            LoadModifications(taskId, out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, true, XlSearchParameters.DecoyType, localizeableModificationTypes, CommonParameters);
            CommonParameters.TotalPartitions = proteinList.Count() / 250;
            if (CommonParameters.TotalPartitions == 0) { CommonParameters.TotalPartitions = 1; }

            MyFileManager myFileManager = new MyFileManager(true);

            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));
            HashSet<DigestionParams> ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);

            ProseCreatedWhileRunning.Append("The following crosslink discovery were used: ");
            ProseCreatedWhileRunning.Append("crosslinker name = " + XlSearchParameters.Crosslinker.CrosslinkerName + "; ");
            ProseCreatedWhileRunning.Append("crosslinker type = " + XlSearchParameters.Crosslinker.Cleavable + "; ");
            ProseCreatedWhileRunning.Append("crosslinker mass = " + XlSearchParameters.Crosslinker.TotalMass + "; ");
            ProseCreatedWhileRunning.Append("crosslinker modification site(s) = " + XlSearchParameters.Crosslinker.CrosslinkerModSites + "; ");

            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");

            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");

            ProseCreatedWhileRunning.Append("parent mass tolerance(s) = " + CommonParameters.PrecursorMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            List<CrosslinkSpectralMatch> AllCsms = new List<CrosslinkSpectralMatch>();

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);

                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams);

                Status("Getting ms2 scans...", thisId);

                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2ScansWrapByScanNum(myMsDataFile, origDataFile, combinedParams, out List<List<(double, int, double)>> precursorss).ToArray();

                List<CrosslinkSpectralMatch>[] newCsmsPerMS2ScanPerFile = new List<CrosslinkSpectralMatch>[arrayOfMs2ScansSortedByMass.Length];

                myFileManager.DoneWithFile(origDataFile);

                //Candidates: currentPartition, peptide Index, score 
                List<(int, int, int)>[] candidates = new List<(int, int, int)>[arrayOfMs2ScansSortedByMass.Length];

                //First round search.
                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<PeptideWithSetModifications> peptideIndex = null;

                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                    Status("Getting fragment dictionary...", new List<string> { taskId });

                    //Only reverse Decoy for crosslink search has been tested and are set as fixed parameter.
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null, currentPartition,
                        UsefulProteomicsDatabases.DecoyType.Reverse, combinedParams, this.FileSpecificParameters, 30000.0, false,
                        dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });

                    List<int>[] fragmentIndex = null;
                    List<int>[] precursorIndex = null;
                    GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);

                    //The second Fragment index is for 'MS1-HCD_MS1-ETD_MS2s' type of data. If LowCID is used for MS1, ion-index is not allowed to use.
                    List<int>[] secondFragmentIndex = null;
                    if (combinedParams.MS2ChildScanDissociationType != DissociationType.LowCID && combinedParams.MS2ChildScanDissociationType!= DissociationType.Unknown
                    && !CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(combinedParams.DissociationType, combinedParams.MS2ChildScanDissociationType))
                    {
                        //Becuase two different type of dissociation methods are used, the parameters are changed with different dissociation type.
                        var secondCombinedParams = CommonParameters.CloneWithNewDissociationType(combinedParams.MS2ChildScanDissociationType);
                        var secondIndexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null, currentPartition,
                            UsefulProteomicsDatabases.DecoyType.Reverse, secondCombinedParams, this.FileSpecificParameters, 30000.0, false,
                            dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });
                        GenerateSecondIndexes(indexEngine, secondIndexEngine, dbFilenameList, ref secondFragmentIndex, proteinList, taskId);
                    }

                    Status("Searching files...", taskId);
                    var crosslinkSearch1stRound = new CrosslinkSearchEngine(newCsmsPerMS2ScanPerFile, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, secondFragmentIndex, currentPartition,
                        combinedParams, this.FileSpecificParameters, XlSearchParameters.Crosslinker, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.CrosslinkAtCleavageSite, XlSearchParameters.XlQuench_H2O,
                        XlSearchParameters.XlQuench_NH2, XlSearchParameters.XlQuench_Tris, thisId, candidates, -1, null, precursorss);

                    crosslinkSearch1stRound.FirstRoundSearch();

                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                    if (GlobalVariables.StopLoops) { break; }
                }

                //Pair the crosslink between two partitions. 
                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<PeptideWithSetModifications> peptideIndex_a = null;

                    List<Protein> proteinListSubset_a = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                    Status("Getting fragment dictionary...", new List<string> { taskId });

                    var indexEngine_a = new IndexingEngine(proteinListSubset_a, variableModifications, fixedModifications, null, null, null, currentPartition,
                        UsefulProteomicsDatabases.DecoyType.Reverse, combinedParams, this.FileSpecificParameters, 30000.0, false,
                        dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });

                    List<int>[] precursorIndex_a = null;
                    GenerateIndexes_PeptideOnly(indexEngine_a, dbFilenameList, ref peptideIndex_a, ref precursorIndex_a, proteinList, taskId);

                    for (int nextPartition = 0; nextPartition < CommonParameters.TotalPartitions; nextPartition++)
                    {
                        List<PeptideWithSetModifications> peptideIndex_b = null;

                        if(currentPartition == nextPartition)
                        {
                            peptideIndex_b = peptideIndex_a;
                        }
                        else
                        {
                            List<Protein> proteinListSubset_b = proteinList.GetRange(nextPartition * proteinList.Count() / combinedParams.TotalPartitions, ((nextPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (nextPartition * proteinList.Count() / combinedParams.TotalPartitions));

                            Status("Getting fragment dictionary...", new List<string> { taskId });

                            var indexEngine_b = new IndexingEngine(proteinListSubset_b, variableModifications, fixedModifications, null, null, null, nextPartition,
                                UsefulProteomicsDatabases.DecoyType.Reverse, combinedParams, this.FileSpecificParameters, 30000.0, false,
                                dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });

                            List<int>[] precursorIndex_b = null;
                            GenerateIndexes_PeptideOnly(indexEngine_b, dbFilenameList, ref peptideIndex_b, ref precursorIndex_b, proteinList, taskId);

                        }

                        Status("Searching files...", taskId);
                        var crosslinkSearch2ndRound = new CrosslinkSearchEngine(newCsmsPerMS2ScanPerFile, arrayOfMs2ScansSortedByMass, peptideIndex_a, null, null, currentPartition,
                            combinedParams, this.FileSpecificParameters, XlSearchParameters.Crosslinker, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.CrosslinkAtCleavageSite, XlSearchParameters.XlQuench_H2O,
                            XlSearchParameters.XlQuench_NH2, XlSearchParameters.XlQuench_Tris, thisId, candidates, nextPartition, peptideIndex_b, precursorss);

                        crosslinkSearch2ndRound.Run();

                    }
                }


                // Add csms candidate for posttask analysis.
                List<List<CrosslinkSpectralMatch>> _ListOfCsmsPerMS2ScanParsimony = new List<List<CrosslinkSpectralMatch>>();

                //For every Ms2Scans, each have a list of candidates psms. The allPsms from CrosslinkSearchEngine is the list (all ms2scans) of list (each ms2scan) of psm (all candidate psm).
                //The allPsmsList is same as allPsms after ResolveAmbiguities.
                foreach (var csmsPerScan in newCsmsPerMS2ScanPerFile)
                {
                    if (csmsPerScan == null)
                    {
                        continue;
                    }
                    foreach (var csm in csmsPerScan)
                    {
                        csm.ResolveAllAmbiguities();
                        if (csm.BetaPeptide != null)
                        {
                            csm.BetaPeptide.ResolveAllAmbiguities();
                        }
                        csm.ResolveProteinPosAmbiguitiesForXl();
                    }

                    var orderedCsmsPerScan = RemoveDuplicateFromCsmsPerScan(csmsPerScan).ToList();

                    _ListOfCsmsPerMS2ScanParsimony.Add(orderedCsmsPerScan);
                }


                if (_ListOfCsmsPerMS2ScanParsimony.SelectMany(p=>p).Count() >0)
                {
                    _ListOfCsmsPerMS2ScanParsimony = SortListsOfCsms(_ListOfCsmsPerMS2ScanParsimony, CommonParameters);
                    AssignCrossType(_ListOfCsmsPerMS2ScanParsimony);


                    PostXLSearchAnalysisTask _postXLSearchAnalysisTask = new PostXLSearchAnalysisTask();
                    _postXLSearchAnalysisTask.FileSpecificParameters = this.FileSpecificParameters;
                    if (currentRawFileList.Count > 1)
                    {
                        var _filteredAllPsms = _postXLSearchAnalysisTask.RunSingleFile(OutputFolder + '/' + Path.GetFileNameWithoutExtension(origDataFile), taskId,
                            _ListOfCsmsPerMS2ScanParsimony.Select(p => p.First()).OrderByDescending(p => p.XLTotalScore).ToList(), XlSearchParameters);
                        AllCsms.AddRange(_filteredAllPsms);
                    }
                    else
                    {
                        AllCsms.AddRange(_ListOfCsmsPerMS2ScanParsimony.Select(p => p.First()).OrderByDescending(p => p.XLTotalScore));
                    }
                }


                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }     

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            PostXLSearchAnalysisTask postXLSearchAnalysisTask = new PostXLSearchAnalysisTask();
            postXLSearchAnalysisTask.FileSpecificParameters = this.FileSpecificParameters;

            return postXLSearchAnalysisTask.Run(OutputFolder, dbFilenameList, currentRawFileList, taskId, AllCsms.OrderByDescending(p => p.XLTotalScore).ToList(), CommonParameters, XlSearchParameters, proteinList, variableModifications, fixedModifications, localizeableModificationTypes, MyTaskResults);
        }


        public static List<List<CrosslinkSpectralMatch>> SortListsOfCsms(List<List<CrosslinkSpectralMatch>> ListOfCsmsPerMS2Scan, CommonParameters commonParameters)
        {
            List<List<CrosslinkSpectralMatch>> newLists = new List<List<CrosslinkSpectralMatch>>();
            foreach (var csmsPerScan in ListOfCsmsPerMS2Scan)
            {
                newLists.Add(SortOneListCsmsSetSecondBestScore(csmsPerScan, commonParameters));
            }

            return newLists.OrderByDescending(c => c.First().XLTotalScore).ThenByDescending(c => c.First().FullSequence + (c.First().BetaPeptide != null ? c.First().BetaPeptide.FullSequence : "")).ToList();
        }

        public static List<CrosslinkSpectralMatch> SortOneListCsmsSetSecondBestScore(List<CrosslinkSpectralMatch> csmsPerScan, CommonParameters commonParameters)
        {
            List<double> xlTotalScores = csmsPerScan.Select(s => s.XLTotalScore).OrderByDescending(s => s).ToList();
            xlTotalScores.RemoveAt(0);

            //This possibly needs to be doubled for xlinks. But, since each list can be a mix of xlinks and nonxlinks we just leave as is for now.
            double secondBestSore = commonParameters.ScoreCutoff;
            if (xlTotalScores.Count() > 0)
            {
                secondBestSore = xlTotalScores[0];
            }
            foreach (CrosslinkSpectralMatch csm in csmsPerScan)
            {
                csm.SecondBestXlScore = secondBestSore;
            }
            csmsPerScan = csmsPerScan.OrderByDescending(c => c.XLTotalScore).ThenBy(c => c.FullSequence + (c.BetaPeptide != null ? c.BetaPeptide.FullSequence : "")).ToList();

            return csmsPerScan;
        }

        public static void AssignCrossType(List<List<CrosslinkSpectralMatch>> ListOfCsmsPerMS2Scan)
        {
            foreach (var csmsPerScan in ListOfCsmsPerMS2Scan)
            {
                foreach (var csm in csmsPerScan)
                {
                    if (csm.CrossType == PsmCrossType.Cross)
                    {
                        if (CrosslinkSpectralMatch.IsIntraCsm(csm))
                        {
                            csm.CrossType = PsmCrossType.Intra;
                        }
                        else
                        {
                            csm.CrossType = PsmCrossType.Inter;
                        }
                    }
                }
            }
        }

        //Remove same peptide with from different protein. We search for every possible peptides from every protein for each scan, and do parsimony later.
        //Same peptides from different proteins are kept for csmsPerScan and parsimony is done here.
        public static List<CrosslinkSpectralMatch> RemoveDuplicateFromCsmsPerScan(List<CrosslinkSpectralMatch> crosslinkSpectralMatches)
        {
            //A dictionary is used for parsimony. the key of the dictionary is 'alphs_fullseq + beta_fullseq'.
            Dictionary<string, CrosslinkSpectralMatch> keyValuePairs = new Dictionary<string, CrosslinkSpectralMatch>();
            foreach (var csm in crosslinkSpectralMatches)
            {
                if (csm == null)
                {
                    continue;
                }
                string betaFullseq = "-";
                if (csm.BetaPeptide != null)
                {
                    betaFullseq += csm.BetaPeptide.FullSequence;
                }

                if (keyValuePairs.ContainsKey(csm.FullSequence + betaFullseq))
                {
                    keyValuePairs[csm.FullSequence + betaFullseq].AddProteinMatch(csm.BestMatchingPeptides.First(), csm.PeptidesToMatchingFragments[csm.BestMatchingPeptides.First().Peptide]);

                    if (csm.BetaPeptide != null)
                    {
                        keyValuePairs[csm.FullSequence + betaFullseq].BetaPeptide.AddProteinMatch(csm.BetaPeptide.BestMatchingPeptides.First(), csm.PeptidesToMatchingFragments[csm.BestMatchingPeptides.First().Peptide]);
                    }
                }
                else
                {
                    keyValuePairs.Add(csm.FullSequence + betaFullseq, csm);
                }
            }
            return keyValuePairs.Values.ToList();
        }

    }

}