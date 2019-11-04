using EngineLayer;
using EngineLayer.CrosslinkSearch;
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
    public class XLSearchTask : MetaMorpheusTask
    {
        public XLSearchTask() : base(MyTask.XLSearch)
        {
            //Default parameter setting which is different from SearchTask, can be overwriten
            var digestPara = new DigestionParams(
                minPeptideLength: 5
            );
            CommonParameters = new CommonParameters(
                precursorMassTolerance: new PpmTolerance(10),
                scoreCutoff: 3,
                trimMsMsPeaks:false,
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

            var crosslinker = XlSearchParameters.Crosslinker;

            MyFileManager myFileManager = new MyFileManager(true);

            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));
            HashSet<DigestionParams> ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);

            ProseCreatedWhileRunning.Append("The following crosslink discovery were used: ");
            ProseCreatedWhileRunning.Append("crosslinker name = " + crosslinker.CrosslinkerName + "; ");
            ProseCreatedWhileRunning.Append("crosslinker type = " + crosslinker.Cleavable + "; ");
            ProseCreatedWhileRunning.Append("crosslinker mass = " + crosslinker.TotalMass + "; ");
            ProseCreatedWhileRunning.Append("crosslinker modification site(s) = " + crosslinker.CrosslinkerModSites + "; ");

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

                List<CrosslinkSpectralMatch>[] newCsmsPerMS2ScanPerFile = new List<CrosslinkSpectralMatch>[arrayOfMs2ScansSortedByMass.Length];
                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<PeptideWithSetModifications> peptideIndex = null;

                    //When partition, the proteinList will be split for each Thread.
                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                    Status("Getting fragment dictionary...", new List<string> { taskId });

                    //Only reverse Decoy for crosslink search has been tested and are set as fixed parameter.
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null, currentPartition, UsefulProteomicsDatabases.DecoyType.Reverse, combinedParams, 30000.0, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), new List<string> { taskId });

                    List<int>[] fragmentIndex = null;
                    List<int>[] precursorIndex = null;
                    GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);

                    //The second Fragment index is for 'MS1-HCD_MS1-ETD_MS2s' type of data. If LowCID is used for MS1, ion-index is not allowed to use.
                    List<int>[] secondFragmentIndex = null;
                    if (combinedParams.ChildScanDissociationType != DissociationType.Unknown)
                    {
                        if (combinedParams.ChildScanDissociationType != DissociationType.LowCID
                        && !CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(combinedParams.DissociationType, combinedParams.ChildScanDissociationType))
                        {
                            //Becuase two different type of dissociation methods are used, the parameters are changed with different dissociation type.
                            var secondCombinedParams = CommonParameters.CloneWithNewDissociationType(combinedParams.ChildScanDissociationType);
                            var secondIndexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, null, null, currentPartition, UsefulProteomicsDatabases.DecoyType.Reverse, secondCombinedParams, 30000.0, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), new List<string> { taskId });
                            GenerateSecondIndexes(indexEngine, secondIndexEngine, dbFilenameList, ref secondFragmentIndex, proteinList, taskId);
                        }
                        else
                        {
                            secondFragmentIndex = fragmentIndex;
                        }
                    }

                    Status("Searching files...", taskId);
                    new CrosslinkSearchEngine(newCsmsPerMS2ScanPerFile, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, secondFragmentIndex, currentPartition, combinedParams, crosslinker,
                        XlSearchParameters.RestrictToTopNHits, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.CrosslinkAtCleavageSite,
                        XlSearchParameters.XlQuench_H2O, XlSearchParameters.XlQuench_NH2, XlSearchParameters.XlQuench_Tris, thisId).Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                    if (GlobalVariables.StopLoops) { break; }
                }

                ListOfCsmsPerMS2Scan.AddRange(newCsmsPerMS2ScanPerFile.Where(p => p != null).ToList());

                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            List<List<CrosslinkSpectralMatch>> ListOfCsmsPerMS2ScanParsimony = new List<List<CrosslinkSpectralMatch>>();

            //For every Ms2Scans, each have a list of candidates psms. The allPsms from CrosslinkSearchEngine is the list (all ms2scans) of list (each ms2scan) of psm (all candidate psm). 
            //The allPsmsList is same as allPsms after ResolveAmbiguities. 
            foreach (var csmsPerScan in ListOfCsmsPerMS2Scan)
            {
                foreach (var csm in csmsPerScan)
                {
                    csm.ResolveAllAmbiguities();
                    if (csm.BetaPeptide != null)
                    {
                        csm.BetaPeptide.ResolveAllAmbiguities();
                    }
                }
                ListOfCsmsPerMS2ScanParsimony.Add(RemoveDuplicateFromCsmsPerScan(csmsPerScan));
            }

            var filteredAllPsms = new List<CrosslinkSpectralMatch>();

            //For each ms2scan, try to find the best candidate psm from the psms list. Add it into filteredAllPsms
            //This function is for current usage, this can be replaced with PEP value. 
            foreach (var csmsPerScan in ListOfCsmsPerMS2ScanParsimony)
            {
                filteredAllPsms.Add(csmsPerScan[0]);
            }

            PostXLSearchAnalysisTask postXLSearchAnalysisTask = new PostXLSearchAnalysisTask();

            return postXLSearchAnalysisTask.Run(OutputFolder, dbFilenameList, currentRawFileList, taskId, fileSettingsList, filteredAllPsms.OrderByDescending(p=>p.XLTotalScore).ToList(), CommonParameters, XlSearchParameters, proteinList, variableModifications, fixedModifications, localizeableModificationTypes, MyTaskResults);

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
                    keyValuePairs[csm.FullSequence + betaFullseq].AddProteinMatch(csm.BestMatchingPeptides.First());

                    if (csm.BetaPeptide != null)
                    {
                        keyValuePairs[csm.FullSequence + betaFullseq].BetaPeptide.AddProteinMatch(csm.BetaPeptide.BestMatchingPeptides.First());
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