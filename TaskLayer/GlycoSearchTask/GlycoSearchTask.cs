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
using Proteomics.RetentionTimePrediction;
using FlashLFQ;

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
                scoreCutoff: 2,
                trimMsMsPeaks: true, 
                numberOfPeaksToKeepPerWindow: 500,
                minimumAllowedIntensityRatioToBasePeak: 0.001,
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
            // Set TotalPartitions
            var partitionSize = 300;
            if (_glycoSearchParameters.GlycoSearchType == GlycoSearchType.NGlycanSearch)
            {
                partitionSize = partitionSize * 8;
            }
            else if (_glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch || _glycoSearchParameters.GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                partitionSize = partitionSize * 2;
            }
            CommonParameters.TotalPartitions = proteinList.Count() / partitionSize;
            if (CommonParameters.TotalPartitions == 0) { CommonParameters.TotalPartitions = 1; }

            MyFileManager myFileManager = new MyFileManager(true);

            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));

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
                ProseCreatedWhileRunning.Append("The N-glycan database: " + _glycoSearchParameters.NGlycanDatabasefile + "\n");
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

                // Ms2ScanWithSpecificMass[] _arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy(b => b.PrecursorMass).ToArray();
                //var _precursorss = new List<(double, int, double)>[_arrayOfMs2ScansSortedByMass.Length].ToList();
                Ms2ScanWithSpecificMass[] _arrayOfMs2ScansSortedByMass = GetMs2ScansWrapByScanNum(myMsDataFile, origDataFile, combinedParams, out List<List<(double, int, double)>> _precursorss).ToArray();
                myFileManager.DoneWithFile(origDataFile);
          
                
                var count = 0;
                while (count < _arrayOfMs2ScansSortedByMass.Length)
                {
                    var dataPartition = 30000;
                    if (_arrayOfMs2ScansSortedByMass.Length - count < dataPartition)
                    {
                        dataPartition = _arrayOfMs2ScansSortedByMass.Length - count;
                    }

                    var arrayOfMs2ScansSortedByMass = _arrayOfMs2ScansSortedByMass[count..(count + dataPartition)];

                    var precursorss = _precursorss.GetRange(count, dataPartition);

                    List<GlycoSpectralMatch>[] newCsmsPerMS2ScanPerFile = new List<GlycoSpectralMatch>[arrayOfMs2ScansSortedByMass.Length];

                    //Candidates: currentPartition, peptide Index, score 
                    List<(int, int, int)>[] candidates = new List<(int, int, int)>[arrayOfMs2ScansSortedByMass.Length];

                    //Ion-index
                    for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                    {

                        List<PeptideWithSetModifications> peptideIndex = null;

                        //When partition, the proteinList will be split for each Thread.
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                        Status("Getting fragment dictionary...", new List<string> { taskId });

                        //Indexing strategy for Glycopepitde.
                        //Generating a second fragment index slow down the program. 
                        //There are different type of fragmentation combination methods for N-glycopeptides and O-glycopeptides.

                        //We try to only index b/y fragment ions (HCD, CID, EThcD). Thus we will only generate b/y fragmentIndex. Even for HCD-trig-EThcD or EThcD-pd-HCD method.
                        //Due to the complexity of O-glycan, the ETD ions didn't generate a lot of peptide backbone ions.

                        //If there are paired spectrum, we will let the user to decide if to index the child spectrum if two indexes should be used. (HCD-trig-ETD)
                        //These data types require more testing: ETD-pd-HCD.  ETD-pd-EThcD.
                        var indexParams = combinedParams.Clone();
                        if (_glycoSearchParameters.Indexing_by_ion
                            && (combinedParams.DissociationType == DissociationType.HCD || combinedParams.DissociationType == DissociationType.CID || combinedParams.DissociationType == DissociationType.EThcD))
                        {
                            indexParams = CommonParameters.CloneWithNewDissociationType(DissociationType.HCD);
                        }

                        //Only reverse Decoy for glyco search has been tested and are set as fixed parameter.
                        var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null,
                            currentPartition, _glycoSearchParameters.DecoyType, indexParams, this.FileSpecificParameters, 30000.0,
                            false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });
                        List<int>[] fragmentIndex = null;
                        List<int>[] precursorIndex = null;
                        GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);

                        //If LowCID is used for MS1, ion-index is not allowed to use.
                        List<int>[] secondFragmentIndex = null;
                        if (combinedParams.MS2ChildScanDissociationType != DissociationType.LowCID
                            && _glycoSearchParameters.IndexingChildScan
                            && _glycoSearchParameters.IndexingChildScanDiffIndex
                            && !CrosslinkSearchEngine.DissociationTypeGenerateSameIons(combinedParams.DissociationType, combinedParams.MS2ChildScanDissociationType))
                        {
                            if (_glycoSearchParameters.Indexing_by_ion && CrosslinkSearchEngine.DissociationTypeGenerateSameIons(DissociationType.HCD, combinedParams.MS2ChildScanDissociationType))
                            {
                                continue;
                            }
                            //Becuase two different type of dissociation methods are used, the parameters are changed with different dissociation type.
                            var secondCombinedParams = CommonParameters.CloneWithNewDissociationType(combinedParams.MS2ChildScanDissociationType);
                            var secondIndexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null,
                                currentPartition, _glycoSearchParameters.DecoyType, secondCombinedParams, this.FileSpecificParameters, 30000.0,
                                false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });
                            GenerateSecondIndexes(indexEngine, secondIndexEngine, dbFilenameList, ref secondFragmentIndex, proteinList, taskId);
                        }

                        Status("Searching index files...", taskId);
                        var glycoSearch1stRound = new GlycoSearchEngine(newCsmsPerMS2ScanPerFile, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex,
                            secondFragmentIndex, currentPartition, combinedParams, this.FileSpecificParameters,
                            _glycoSearchParameters.OGlycanDatabasefile, _glycoSearchParameters.NGlycanDatabasefile,
                            _glycoSearchParameters.GlycoSearchType, _glycoSearchParameters.MixedGlycoAllowed,
                            _glycoSearchParameters.GlycoSearchTopNum, _glycoSearchParameters.MaximumOGlycanAllowed,
                            _glycoSearchParameters.MaximumNGlycanAllowed, _glycoSearchParameters.OxoniumIonFilt, _glycoSearchParameters.IndexingChildScan, thisId,
                            candidates, precursorss);

                        glycoSearch1stRound.FirstRoundSearch();

                        ReportProgress(new ProgressEventArgs(100, "Done with index search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                        if (GlobalVariables.StopLoops) { break; }
                    }

                    //Match peptides
                    for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                    {
                        List<PeptideWithSetModifications> peptideIndex = null;

                        List<Protein> proteinListSubset_a = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                        Status("Getting fragment dictionary...", new List<string> { taskId });

                        var indexEngine = new IndexingEngine(proteinListSubset_a, variableModifications, fixedModifications, null, null, null, currentPartition,
                            UsefulProteomicsDatabases.DecoyType.Reverse, combinedParams, this.FileSpecificParameters, 30000.0, false,
                            dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string> { taskId });

                        List<int>[] precursorIndex = null;
                        GenerateIndexes_PeptideOnly(indexEngine, dbFilenameList, ref peptideIndex, ref precursorIndex, proteinList, taskId);

                        Status("Searching files...", taskId);
                        var glycoSearch2ndRound = new GlycoSearchEngine(newCsmsPerMS2ScanPerFile, arrayOfMs2ScansSortedByMass, peptideIndex, null,
                            null, currentPartition, combinedParams, this.FileSpecificParameters,
                            _glycoSearchParameters.OGlycanDatabasefile, _glycoSearchParameters.NGlycanDatabasefile,
                            _glycoSearchParameters.GlycoSearchType, _glycoSearchParameters.MixedGlycoAllowed,
                            _glycoSearchParameters.GlycoSearchTopNum, _glycoSearchParameters.MaximumOGlycanAllowed,
                            _glycoSearchParameters.MaximumNGlycanAllowed, _glycoSearchParameters.OxoniumIonFilt, _glycoSearchParameters.IndexingChildScan, thisId,
                            candidates, precursorss);

                        glycoSearch2ndRound.Run();

                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                        if (GlobalVariables.StopLoops) { break; }
                    }

                    ListOfGsmsPerMS2Scan.AddRange(newCsmsPerMS2ScanPerFile.Where(p => p != null).ToList());

                    //// Does it release memory here?
                    //for (int i = count; i < count + dataPartition; i++)
                    //{
                    //    _arrayOfMs2ScansSortedByMass[i] = null;
                    //}

                    count += dataPartition;
                }


                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            //For every Ms2Scans, each have a list of candidates psms. The allPsms from GlycoSearchEngine is the list (all ms2scans) of list (each ms2scan) of psm (all candidate psm). 
            //Currently, only keep the first scan for consideration. 
            List<GlycoSpectralMatch> GsmPerScans = ListOfGsmsPerMS2Scan.Select(p => p.First()).ToList();

            var filteredAllPsms = new List<GlycoSpectralMatch>();

            //SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);

            //For each ms2scan, try to find the best candidate psm from the psms list. Do the localizaiton analysis. Add it into filteredAllPsms.
            foreach (var gsmsPerScan in GsmPerScans.GroupBy(p => p.ScanNumber))
            {
                var glycos = RemoveSimilarSequenceDuplicates(gsmsPerScan.OrderByDescending(p=>p.Score).ToList());

                foreach (var glycoSpectralMatch in glycos)
                {
                    //glycoSpectralMatch.PredictedHydrophobicity = calc.ScoreSequence(glycoSpectralMatch.BestMatchingPeptides.First().Peptide);

                    filteredAllPsms.Add(glycoSpectralMatch);
                }
            }

            FlashLfqResults flashLfqResults = null;

            PostSearchAnalysisParameters parameters = new PostSearchAnalysisParameters
            {
                SearchTaskResults = MyTaskResults,
                SearchTaskId = taskId,          
                ProteinList = proteinList,
                VariableModifications = variableModifications,
                FixedModifications = fixedModifications,
                ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams)),
                CurrentRawFileList = currentRawFileList,
                MyFileManager = myFileManager,
                OutputFolder = OutputFolder,
                IndividualResultsOutputFolder = Path.Combine(OutputFolder, "Individual File Results"),
                FlashLfqResults = flashLfqResults,
                FileSettingsList = fileSettingsList,
                DatabaseFilenameList = dbFilenameList
            };

            PostGlycoSearchAnalysisTask postGlycoSearchAnalysisTask = new PostGlycoSearchAnalysisTask 
            {
                Parameters = parameters,
                FileSpecificParameters = this.FileSpecificParameters,
                CommonParameters = CommonParameters,
                
            };

            return postGlycoSearchAnalysisTask.Run(filteredAllPsms.OrderByDescending(p => p.Score).ToList(),  _glycoSearchParameters, MyTaskResults);

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
