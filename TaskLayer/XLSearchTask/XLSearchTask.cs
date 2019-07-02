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
            List<List<CrosslinkSpectralMatch>> allPsms = new List<List<CrosslinkSpectralMatch>>();

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

                List<CrosslinkSpectralMatch>[] newPsms = new List<CrosslinkSpectralMatch>[arrayOfMs2ScansSortedByMass.Length];
                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<PeptideWithSetModifications> peptideIndex = null;
                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                    Status("Getting fragment dictionary...", new List<string> { taskId });
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, currentPartition, UsefulProteomicsDatabases.DecoyType.Reverse, combinedParams, 30000.0, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), new List<string> { taskId });
                    List<int>[] fragmentIndex = null;
                    List<int>[] precursorIndex = null;
                    GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);

                    List<int>[] secondFragmentIndex = null;
                    if (combinedParams.ChildScanDissociationType != DissociationType.LowCID
                    && !CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(combinedParams.DissociationType, combinedParams.ChildScanDissociationType))
                    {
                        var secondCombinedParams = CommonParameters.CloneWithNewDissociationType(combinedParams.ChildScanDissociationType);
                        var secondIndexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, currentPartition, UsefulProteomicsDatabases.DecoyType.Reverse, secondCombinedParams, 30000.0, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), new List<string> { taskId });
                        GenerateSecondIndexes(indexEngine, secondIndexEngine, dbFilenameList, ref secondFragmentIndex, proteinList, taskId);
                    }

                    Status("Searching files...", taskId);
                    new CrosslinkSearchEngine(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, secondFragmentIndex, currentPartition, combinedParams, crosslinker,
                        XlSearchParameters.RestrictToTopNHits, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.XlQuench_H2O,
                        XlSearchParameters.XlQuench_NH2, XlSearchParameters.XlQuench_Tris, thisId).Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                    if (GlobalVariables.StopLoops) { break; }
                }

                allPsms.AddRange(newPsms.Where(p => p != null).ToList());

                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            List<List<CrosslinkSpectralMatch>> filteredAllPsms = new List<List<CrosslinkSpectralMatch>>();

            foreach (var psmsPerScan in allPsms)
            {
                foreach (var psm in psmsPerScan)
                {
                    psm.ResolveAllAmbiguities();
                    if (psm.BetaPeptide != null)
                    {
                        psm.BetaPeptide.ResolveAllAmbiguities();
                    }
                }
                filteredAllPsms.Add(RemoveDuplicateFromPsmsPerScan(psmsPerScan));
            }

            var allPsmsList = filteredAllPsms.Select(p => p.First()).OrderByDescending(p => p.XLTotalScore).ToList();

            PostXLSearchAnalysisTask postXLSearchAnalysisTask = new PostXLSearchAnalysisTask();

            return postXLSearchAnalysisTask.Run(OutputFolder, dbFilenameList, currentRawFileList, taskId, fileSettingsList, allPsmsList, CommonParameters, XlSearchParameters, proteinList, variableModifications, fixedModifications, localizeableModificationTypes, MyTaskResults);

        }

        //Remove same peptide with from different protein.
        public static List<CrosslinkSpectralMatch> RemoveDuplicateFromPsmsPerScan(List<CrosslinkSpectralMatch> crosslinkSpectralMatches)
        {
            Dictionary<string, CrosslinkSpectralMatch> keyValuePairs = new Dictionary<string, CrosslinkSpectralMatch>();
            foreach (var csm in crosslinkSpectralMatches)
            {
                if (csm == null)
                {
                    continue;
                }
                string betaFulseq = "-";
                if (csm.BetaPeptide != null)
                {
                    betaFulseq += csm.BetaPeptide.FullSequence;
                }

                if (keyValuePairs.ContainsKey(csm.FullSequence + betaFulseq))
                {
                    keyValuePairs[csm.FullSequence + betaFulseq].AddProteinMatch(csm.BestMatchingPeptides.First());

                    if (csm.BetaPeptide != null)
                    {
                        keyValuePairs[csm.FullSequence + betaFulseq].BetaPeptide.AddProteinMatch(csm.BetaPeptide.BestMatchingPeptides.First());
                    }
                }
                else
                {
                    keyValuePairs.Add(csm.FullSequence + betaFulseq, csm);
                }
            }
            return keyValuePairs.Values.ToList();
        }

    }
}