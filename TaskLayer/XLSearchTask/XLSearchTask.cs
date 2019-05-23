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
using System.Threading.Tasks;
using System.Collections.Concurrent;

namespace TaskLayer
{
    public partial class XLSearchTask : MetaMorpheusTask
    {
        public XLSearchTask() : base(MyTask.XLSearch)
        {
            CommonParameters = new CommonParameters(
                precursorMassTolerance: new PpmTolerance(10),
                scoreCutoff: 3
            );

            XlSearchParameters = new XlSearchParameters();
        }

        public XlSearchParameters XlSearchParameters { get; set; }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            List<CrosslinkSpectralMatch> allPsms = new List<CrosslinkSpectralMatch>();

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
                MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(combinedParams.PrecursorMassTolerance, XlSearchParameters.MassDiffAcceptorType, "");

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);

                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams);

                Status("Getting ms2 scans...", thisId);

                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy(b => b.PrecursorMass).ToArray();

                if (XlSearchParameters.AnalyzeOxiniumIons || XlSearchParameters.FilterScanOxiniumIons)
                {
                    MassDiffAcceptor massDiffAcceptor_oxiniumIons = new SinglePpmAroundZeroSearchMode(combinedParams.ProductMassTolerance.Value);
                    Tuple<int, double[]>[] tuples = new Tuple<int, double[]>[arrayOfMs2ScansSortedByMass.Length];
                    Parallel.ForEach(Partitioner.Create(0, arrayOfMs2ScansSortedByMass.Length), new ParallelOptions { MaxDegreeOfParallelism = combinedParams.MaxThreadsToUsePerFile }, (range, loopState) =>
                    {
                        for (int scanIndex = range.Item1; scanIndex < range.Item2; scanIndex++)
                        {
                            double[] oxoniumIonIntensities = GlycoPeptides.ScanGetOxoniumIons(arrayOfMs2ScansSortedByMass[scanIndex], massDiffAcceptor_oxiniumIons);
                            arrayOfMs2ScansSortedByMass[scanIndex].OxiniumIonNum = oxoniumIonIntensities.Where(p=>p>0).Count();
                            tuples[scanIndex] = new Tuple<int, double[]>(arrayOfMs2ScansSortedByMass[scanIndex].OneBasedScanNumber, oxoniumIonIntensities);
                        }
                    });
                    if (XlSearchParameters.AnalyzeOxiniumIons)
                    {
                        var writtenFile = Path.Combine(OutputFolder, "oxiniumIons" + ".tsv");
                        WriteOxoniumIons(tuples, writtenFile);
                    }

                    if (XlSearchParameters.FilterScanOxiniumIons)
                    {
                        arrayOfMs2ScansSortedByMass = arrayOfMs2ScansSortedByMass.Where(p => p.OxiniumIonNum >= 2).ToArray();
                    }
                }

                CrosslinkSpectralMatch[] newPsms = new CrosslinkSpectralMatch[arrayOfMs2ScansSortedByMass.Length];
                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<PeptideWithSetModifications> peptideIndex = null;
                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                    Status("Getting fragment dictionary...", new List<string> { taskId });
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, null, currentPartition, UsefulProteomicsDatabases.DecoyType.Reverse, combinedParams, 30000.0, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), new List<string> { taskId });
                    List<int>[] fragmentIndex = null;
                    List<int>[] precursorIndex = null;

                    GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, GlobalVariables.AllModsKnown.ToList(), taskId);                  

                    Status("Searching files...", taskId);
                    new CrosslinkSearchEngine(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition, combinedParams, massDiffAcceptor, XlSearchParameters.OpenSearchType, crosslinker,
                        XlSearchParameters.RestrictToTopNHits, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.XlQuench_H2O,
                        XlSearchParameters.XlQuench_NH2, XlSearchParameters.XlQuench_Tris, XlSearchParameters.SearchGlycan182, thisId).Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                    if (GlobalVariables.StopLoops) { break; }
                }

                allPsms.AddRange(newPsms.Where(p => p != null));

                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            allPsms = allPsms.OrderByDescending(p => p.XLTotalScore).ToList();

            if (XlSearchParameters.OpenSearchType == OpenSearchType.NGlyco)
            {
                //SingleFDRAnalysis(allPsms, new List<string> { taskId });             
                //var writtenFileInter = Path.Combine(OutputFolder, "all_fdr" + ".mytsv");
                //WritePsmCrossToTsv(allPsms, writtenFileInter, 3);

                var allPsmsSingle = allPsms.Where(p => p.Glycan == null && p.Score > 2).OrderByDescending(p => p.XLTotalScore).ToList();           
                SingleFDRAnalysis(allPsmsSingle, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_fdr" + ".mytsv");
                WritePsmCrossToTsv(allPsmsSingle, writtenFileInter1, 1);

                //TO DO: there may have a bug. I have to filter the following loopPsms, deadendPsms with a BestScore higher than 2, Or some of the Psms will have everything be 0!
                var allPsmsGly = allPsms.Where(p => p.Glycan != null && p.Score > 2).OrderByDescending(p => p.XLTotalScore).ToList();
                SingleFDRAnalysis(allPsmsGly, new List<string> { taskId });
                var writtenFileInter2 = Path.Combine(OutputFolder, "glyco_fdr" + ".mytsv");
                WritePsmCrossToTsv(allPsmsGly, writtenFileInter2, 3);

                //var allPsmsGlyForFDR = allPsms.SelectMany(p => p.crosslinkSpectralMatches).ToList();
                //var writtenFileInter3 = Path.Combine(OutputFolder, "glyco_for_fdr" + ".mytsv");
                //WritePsmCrossToTsv(allPsmsGlyForFDR, writtenFileInter3, 3);

                return MyTaskResults;
            }

            if (XlSearchParameters.OpenSearchType == OpenSearchType.OGlyco)
            {
                var allPsmsSingle = allPsms.Where(p => p.glycanBoxes == null && p.Score > 2).OrderByDescending(p => p.XLTotalScore).ToList();
                SingleFDRAnalysis(allPsmsSingle, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_fdr" + ".mytsv");
                WritePsmCrossToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsGly = allPsms.Where(p => p.glycanBoxes != null && p.Score > 2).OrderByDescending(p => p.XLTotalScore).ToList();
                SingleFDRAnalysis(allPsmsGly, new List<string> { taskId });
                var writtenFileInter2 = Path.Combine(OutputFolder, "glyco_fdr" + ".mytsv");
                WritePsmCrossToTsv(allPsmsGly, writtenFileInter2, 3);
                return MyTaskResults;
            }

            var allPsmsXL = allPsms.Where(p => p.CrossType == PsmCrossType.Cross).ToList();

            // inter-crosslinks; different proteins are linked
            var interCsms = allPsmsXL.Where(p => !p.ProteinAccession.Equals(p.BetaPeptide.ProteinAccession)).ToList();
            foreach (var item in interCsms)
            {
                item.CrossType = PsmCrossType.Inter;
            }

            // intra-crosslinks; crosslinks within a protein
            var intraCsms = allPsmsXL.Where(p => p.ProteinAccession.Equals(p.BetaPeptide.ProteinAccession)).ToList();
            foreach (var item in intraCsms)
            {
                item.CrossType = PsmCrossType.Intra;
            }

            // calculate FDR
            DoCrosslinkFdrAnalysis(interCsms);
            DoCrosslinkFdrAnalysis(intraCsms);
            SingleFDRAnalysis(allPsms, new List<string> { taskId });

            // calculate protein crosslink residue numbers
            foreach (var csm in allPsmsXL)
            {
                // alpha peptide crosslink residue in the protein
                csm.XlProteinPos = csm.OneBasedStartResidueInProtein.Value + csm.LinkPositions[0] - 1;

                // beta crosslink residue in protein
                csm.BetaPeptide.XlProteinPos = csm.BetaPeptide.OneBasedStartResidueInProtein.Value + csm.BetaPeptide.LinkPositions[0] - 1;
            }

            // write interlink CSMs
            if (interCsms.Any())
            {
                string file = Path.Combine(OutputFolder, "XL_Interlinks.tsv");
                WritePsmCrossToTsv(interCsms, file, 2);
                FinishedWritingFile(file, new List<string> { taskId });
            }
            MyTaskResults.AddNiceText("Target inter-crosslinks within 1% FDR: " + interCsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.BetaPeptide.IsDecoy));

            if (XlSearchParameters.WriteOutputForPercolator)
            {
                var interPsmsXLPercolator = interCsms.Where(p => p.Score >= 2 && p.BetaPeptide.Score >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteCrosslinkToTxtForPercolator(interPsmsXLPercolator, OutputFolder, "XL_Interlinks_Percolator", crosslinker, new List<string> { taskId });
            }

            // write intralink CSMs
            if (intraCsms.Any())
            {
                string file = Path.Combine(OutputFolder, "XL_Intralinks.tsv");
                WritePsmCrossToTsv(intraCsms, file, 2);
                FinishedWritingFile(file, new List<string> { taskId });
            }
            MyTaskResults.AddNiceText("Target intra-crosslinks within 1% FDR: " + intraCsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.BetaPeptide.IsDecoy));

            if (XlSearchParameters.WriteOutputForPercolator)
            {
                var intraPsmsXLPercolator = intraCsms.Where(p => p.Score >= 2 && p.BetaPeptide.Score >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteCrosslinkToTxtForPercolator(intraPsmsXLPercolator, OutputFolder, "XL_Intralinks_Percolator", crosslinker, new List<string> { taskId });
            }

            // write single peptides
            var singlePsms = allPsms.Where(p => p.CrossType == PsmCrossType.Single).ToList();
            if (singlePsms.Any())
            {
                string writtenFileSingle = Path.Combine(OutputFolder, "SinglePeptides" + ".tsv");
                WritePsmCrossToTsv(singlePsms, writtenFileSingle, 1);
                FinishedWritingFile(writtenFileSingle, new List<string> { taskId });
            }
            MyTaskResults.AddNiceText("Target single peptides within 1% FDR: " + singlePsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy));

            // write loops
            var loopPsms = allPsms.Where(p => p.CrossType == PsmCrossType.Loop).ToList();
            if (loopPsms.Any())
            {
                string writtenFileLoop = Path.Combine(OutputFolder, "Looplinks" + ".tsv");
                WritePsmCrossToTsv(loopPsms, writtenFileLoop, 1);
                FinishedWritingFile(writtenFileLoop, new List<string> { taskId });
            }
            MyTaskResults.AddNiceText("Target loop-linked peptides within 1% FDR: " + loopPsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy));

            // write deadends
            var deadendPsms = allPsms.Where(p => p.CrossType == PsmCrossType.DeadEnd 
                || p.CrossType == PsmCrossType.DeadEndH2O 
                || p.CrossType == PsmCrossType.DeadEndNH2 
                || p.CrossType == PsmCrossType.DeadEndTris).ToList();
            if (deadendPsms.Any())
            {
                string writtenFileDeadend = Path.Combine(OutputFolder, "Deadends" + ".tsv");
                WritePsmCrossToTsv(deadendPsms, writtenFileDeadend, 1);
                FinishedWritingFile(writtenFileDeadend, new List<string> { taskId });
            }
            MyTaskResults.AddNiceText("Target deadend peptides within 1% FDR: " + deadendPsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy));

            // write pepXML
            if (XlSearchParameters.WritePepXml)
            {
                List<CrosslinkSpectralMatch> writeToXml = new List<CrosslinkSpectralMatch>();
                writeToXml.AddRange(intraCsms.Where(p => !p.IsDecoy && !p.BetaPeptide.IsDecoy && p.FdrInfo.QValue <= 0.05));
                writeToXml.AddRange(interCsms.Where(p => !p.IsDecoy && !p.BetaPeptide.IsDecoy && p.FdrInfo.QValue <= 0.05));
                writeToXml.AddRange(singlePsms.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.05));
                writeToXml.AddRange(loopPsms.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.05));
                writeToXml.AddRange(deadendPsms.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.05));
                writeToXml = writeToXml.OrderBy(p => p.ScanNumber).ToList();

                foreach (var fullFilePath in currentRawFileList)
                {
                    string fileNameNoExtension = Path.GetFileNameWithoutExtension(fullFilePath);
                    WritePepXML_xl(writeToXml.Where(p => p.FullFilePath == fullFilePath).ToList(), proteinList, dbFilenameList[0].FilePath, variableModifications, fixedModifications, localizeableModificationTypes, OutputFolder, fileNameNoExtension, new List<string> { taskId });
                }
            }

            return MyTaskResults;
        }

        //Calculate the FDR of single peptide FP/TP
        private void SingleFDRAnalysis(List<CrosslinkSpectralMatch> items, List<string> taskIds)
        {
            // calculate single PSM FDR
            List<PeptideSpectralMatch> psms = items.Where(p => p.CrossType == PsmCrossType.Single).Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, CommonParameters, taskIds).Run();

            // calculate loop PSM FDR
            psms = items.Where(p => p.CrossType == PsmCrossType.Loop).Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, CommonParameters, taskIds).Run();

            // calculate deadend FDR
            psms = items.Where(p => p.CrossType == PsmCrossType.DeadEnd ||
                p.CrossType == PsmCrossType.DeadEndH2O ||
                p.CrossType == PsmCrossType.DeadEndNH2 ||
                p.CrossType == PsmCrossType.DeadEndTris).Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, CommonParameters, taskIds).Run();
        }

        //Calculate the FDR of crosslinked peptide FP/TP
        private void DoCrosslinkFdrAnalysis(List<CrosslinkSpectralMatch> csms)
        {
            int cumulativeTarget = 0;
            int cumulativeDecoy = 0;

            for (int i = 0; i < csms.Count; i++)
            {
                var csm = csms[i];
                if (csm.IsDecoy || csm.BetaPeptide.IsDecoy)
                {
                    cumulativeDecoy++;
                }
                else
                {
                    cumulativeTarget++;
                }

                double qValue = Math.Min(1, (double)cumulativeDecoy / cumulativeTarget);
                csm.SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, 0, 0, 0, 0, 0, 0, false);
            }

            double qValueThreshold = 1.0;
            for (int i = csms.Count - 1; i >= 0; i--)
            {
                CrosslinkSpectralMatch csm = csms[i];

                // threshold q-values
                if (csm.FdrInfo.QValue > qValueThreshold)
                {
                    csm.FdrInfo.QValue = qValueThreshold;
                }
                else if (csm.FdrInfo.QValue < qValueThreshold)
                {
                    qValueThreshold = csm.FdrInfo.QValue;
                }
            }
        }       
    }
}