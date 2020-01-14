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
                childScanDissociationType: DissociationType.EThcD,
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

            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));
            HashSet<DigestionParams> ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);

            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");

            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif) + "; "));
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");

            ProseCreatedWhileRunning.Append("parent mass tolerance(s) = " + CommonParameters.PrecursorMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

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

                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<PeptideWithSetModifications> peptideIndex = null;

                    //When partition, the proteinList will be split for each Thread.
                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                    Status("Getting fragment dictionary...", new List<string> { taskId });

                    //Only reverse Decoy for glyco search has been tested and are set as fixed parameter.
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null, currentPartition, UsefulProteomicsDatabases.DecoyType.Reverse, combinedParams, 30000.0, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), new List<string> { taskId });
                    List<int>[] fragmentIndex = null;
                    List<int>[] precursorIndex = null;
                    GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);

                    //The second Fragment index is for 'MS1-HCD_MS1-ETD_MS2s' type of data. If LowCID is used for MS1, ion-index is not allowed to use.
                    List<int>[] secondFragmentIndex = null;
                    if (combinedParams.ChildScanDissociationType != DissociationType.LowCID
                    && !CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(combinedParams.DissociationType, combinedParams.ChildScanDissociationType))
                    {
                        //Becuase two different type of dissociation methods are used, the parameters are changed with different dissociation type.
                        var secondCombinedParams = CommonParameters.CloneWithNewDissociationType(combinedParams.ChildScanDissociationType);
                        var secondIndexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null, currentPartition, UsefulProteomicsDatabases.DecoyType.Reverse, secondCombinedParams, 30000.0, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), new List<string> { taskId });
                        GenerateSecondIndexes(indexEngine, secondIndexEngine, dbFilenameList, ref secondFragmentIndex, proteinList, taskId);
                    }

                    Status("Searching files...", taskId);
                    new GlycoSearchEngine(newCsmsPerMS2ScanPerFile, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, secondFragmentIndex, currentPartition, combinedParams,
                        _glycoSearchParameters.GlycanDatabasefileIndex, _glycoSearchParameters.IsOGlycoSearch, _glycoSearchParameters.GlycoSearchTopNum, _glycoSearchParameters.MaximumOGlycanAllowed, thisId).Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                    if (GlobalVariables.StopLoops) { break; }
                }

                ListOfGsmsPerMS2Scan.AddRange(newCsmsPerMS2ScanPerFile.Where(p => p != null).ToList());

                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            //For every Ms2Scans, each have a list of candidates psms. The allPsms from CrosslinkSearchEngine is the list (all ms2scans) of list (each ms2scan) of psm (all candidate psm). 
            //The allPsmsList is same as allPsms after ResolveAmbiguities. 
            foreach (var gsmsPerScan in ListOfGsmsPerMS2Scan)
            {
                foreach (var gsm in gsmsPerScan)
                {
                    gsm.ResolveAllAmbiguities();
                }
            }

            var filteredAllPsms = new List<GlycoSpectralMatch>();

            //For each ms2scan, try to find the best candidate psm from the psms list. Add it into filteredAllPsms
            //This function is for current usage, this can be replaced with PEP value. 
            foreach (var gsmsPerScan in ListOfGsmsPerMS2Scan)
            {
                GlycoSpectralMatch glycoSpectralMatch = gsmsPerScan[0];

                filteredAllPsms.Add(glycoSpectralMatch);
            }

            PostGlycoSearchAnalysisTask postGlycoSearchAnalysisTask = new PostGlycoSearchAnalysisTask();

            return postGlycoSearchAnalysisTask.Run(OutputFolder, dbFilenameList, currentRawFileList, taskId, fileSettingsList, filteredAllPsms.OrderByDescending(p => p.TotalScore).ToList(), CommonParameters, _glycoSearchParameters, proteinList, variableModifications, fixedModifications, localizeableModificationTypes, MyTaskResults);

        }


    }
}
