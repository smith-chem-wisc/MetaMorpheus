using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.FdrAnalysis;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public class PostXLSearchAnalysisTask : MetaMorpheusTask
    {
        public PostXLSearchAnalysisTask()
    : base(MyTask.Search)
        {
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            return null;
        }

        public MyTaskResults Run(string outputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, List<CrosslinkSpectralMatch> allPsms, CommonParameters commonParameters, XlSearchParameters xlSearchParameters, List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications, List<string> localizeableModificationTypes, MyTaskResults MyTaskResults)
        {               
            // inter-crosslinks; different proteins are linked
            List<CrosslinkSpectralMatch> interCsms = allPsms.Where(p => p.CrossType == PsmCrossType.Inter).OrderByDescending(p => p.XLTotalScore).ToList();

            // intra-crosslinks; crosslinks within a protein
            List<CrosslinkSpectralMatch> intraCsms = allPsms.Where(p => p.CrossType == PsmCrossType.Intra).OrderByDescending(p => p.XLTotalScore).ToList();       

            var singlePsms = allPsms.Where(p => p.CrossType == PsmCrossType.Single).OrderByDescending(p => p.Score).ToList();

            var loopPsms = allPsms.Where(p => p.CrossType == PsmCrossType.Loop).OrderByDescending(p => p.Score).ToList();

            var deadendPsms = allPsms.Where(p => p.CrossType == PsmCrossType.DeadEnd
                    || p.CrossType == PsmCrossType.DeadEndH2O
                    || p.CrossType == PsmCrossType.DeadEndNH2
                    || p.CrossType == PsmCrossType.DeadEndTris).OrderByDescending(p => p.Score).ToList();

            ComputeXlinkQandPValues(allPsms, intraCsms, interCsms, commonParameters, taskId);

            WriteCsvFiles(outputFolder, interCsms, intraCsms, singlePsms, loopPsms, deadendPsms, taskId, xlSearchParameters);
            MyTaskResults.AddTaskSummaryText("Target inter-crosslinks within 1% FDR: " + interCsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.BetaPeptide.IsDecoy));
            MyTaskResults.AddTaskSummaryText("Target intra-crosslinks within 1% FDR: " + intraCsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.BetaPeptide.IsDecoy));
            MyTaskResults.AddTaskSummaryText("Target single peptides within 1% FDR: " + singlePsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy));
            MyTaskResults.AddTaskSummaryText("Target loop-linked peptides within 1% FDR: " + loopPsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy));
            MyTaskResults.AddTaskSummaryText("Target deadend peptides within 1% FDR: " + deadendPsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy));

            // write pepXML
            if (xlSearchParameters.WritePepXml)
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
                    WriteFile.WritePepXML_xl(writeToXml.Where(p => p.FullFilePath == fullFilePath).ToList(), proteinList, dbFilenameList[0].FilePath, variableModifications, fixedModifications, localizeableModificationTypes, outputFolder, fileNameNoExtension, commonParameters, xlSearchParameters);
                    FinishedWritingFile(Path.Combine(outputFolder, fileNameNoExtension + ".pep.XML"), new List<string> { taskId });
                }
            }

            return MyTaskResults;
        }

        public List<CrosslinkSpectralMatch> RunSingleFile(string outputFolder, string taskId, List<CrosslinkSpectralMatch> allPsms, XlSearchParameters xlSearchParameters)
        {
            List<CrosslinkSpectralMatch> allFilterPsms = new List<CrosslinkSpectralMatch>();

            if (!Directory.Exists(outputFolder)) { Directory.CreateDirectory(outputFolder); }


            // inter-crosslinks; different proteins are linked
            List<CrosslinkSpectralMatch> interCsms = allPsms.Where(p => p.CrossType == PsmCrossType.Inter).OrderByDescending(p => p.XLTotalScore).ToList();

            // intra-crosslinks; crosslinks within a protein
            List<CrosslinkSpectralMatch> intraCsms = allPsms.Where(p => p.CrossType == PsmCrossType.Intra).OrderByDescending(p => p.XLTotalScore).ToList();

            var singlePsms = allPsms.Where(p => p.CrossType == PsmCrossType.Single).OrderByDescending(p => p.Score).ToList();

            var loopPsms = allPsms.Where(p => p.CrossType == PsmCrossType.Loop).OrderByDescending(p => p.Score).ToList();

            var deadendPsms = allPsms.Where(p => p.CrossType == PsmCrossType.DeadEnd
                    || p.CrossType == PsmCrossType.DeadEndH2O
                    || p.CrossType == PsmCrossType.DeadEndNH2
                    || p.CrossType == PsmCrossType.DeadEndTris).OrderByDescending(p => p.Score).ToList();

            DoCrosslinkFdrAnalysis(interCsms);
            DoCrosslinkFdrAnalysis(intraCsms);
            DoCrosslinkFdrAnalysis(singlePsms);
            DoCrosslinkFdrAnalysis(loopPsms);
            DoCrosslinkFdrAnalysis(deadendPsms);

            WriteCsvFiles(outputFolder, interCsms, intraCsms, singlePsms, loopPsms, deadendPsms, taskId, xlSearchParameters);

            allFilterPsms.AddRange(interCsms.Where(p=>p.FdrInfo.QValue <= 0.2)); 
            allFilterPsms.AddRange(intraCsms.Where(p => p.FdrInfo.QValue <= 0.2));
            allFilterPsms.AddRange(singlePsms.Where(p => p.FdrInfo.QValue <= 0.2));
            allFilterPsms.AddRange(loopPsms.Where(p => p.FdrInfo.QValue <= 0.2));
            allFilterPsms.AddRange(deadendPsms.Where(p => p.FdrInfo.QValue <= 0.2));

            return allFilterPsms;
        }

        public void WriteCsvFiles(string outputFolder, List<CrosslinkSpectralMatch> interCsms, List<CrosslinkSpectralMatch> intraCsms,  List<CrosslinkSpectralMatch> singlePsms, 
            List<CrosslinkSpectralMatch> loopPsms, List<CrosslinkSpectralMatch> deadendPsms, string taskId, XlSearchParameters xlSearchParameters)
        {
            // write interlink CSMs
            if (interCsms.Any())
            {
                string file = Path.Combine(outputFolder, "XL_Interlinks.tsv");
                WriteFile.WritePsmCrossToTsv(interCsms, file, 2);
                FinishedWritingFile(file, new List<string> { taskId });
            }
            
            if (xlSearchParameters.WriteOutputForPercolator)
            {
                var interPsmsXLPercolator = interCsms.Where(p => p.Score >= 2 && p.BetaPeptide.Score >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteFile.WriteCrosslinkToTxtForPercolator(interPsmsXLPercolator, outputFolder, "XL_Interlinks_Percolator", xlSearchParameters.Crosslinker);
                FinishedWritingFile(Path.Combine(outputFolder, "XL_Interlinks_Percolator.txt"), new List<string> { taskId });
            }

            // write intralink CSMs
            if (intraCsms.Any())
            {
                string file = Path.Combine(outputFolder, "XL_Intralinks.tsv");
                WriteFile.WritePsmCrossToTsv(intraCsms, file, 2);
                FinishedWritingFile(file, new List<string> { taskId });
            }

            if (xlSearchParameters.WriteOutputForPercolator)
            {
                var intraPsmsXLPercolator = intraCsms.Where(p => p.Score >= 2 && p.BetaPeptide.Score >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteFile.WriteCrosslinkToTxtForPercolator(intraPsmsXLPercolator, outputFolder, "XL_Intralinks_Percolator", xlSearchParameters.Crosslinker);
                FinishedWritingFile(Path.Combine(outputFolder, "XL_Intralinks_Percolator.txt"), new List<string> { taskId });
            }

            // write single peptides
            
            if (singlePsms.Any())
            {
                string writtenFileSingle = Path.Combine(outputFolder, "SinglePeptides" + ".tsv");
                WriteFile.WritePsmCrossToTsv(singlePsms, writtenFileSingle, 1);
                FinishedWritingFile(writtenFileSingle, new List<string> { taskId });
            }

            // write loops        
            if (loopPsms.Any())
            {
                string writtenFileLoop = Path.Combine(outputFolder, "Looplinks" + ".tsv");
                WriteFile.WritePsmCrossToTsv(loopPsms, writtenFileLoop, 1);
                FinishedWritingFile(writtenFileLoop, new List<string> { taskId });
            }

            // write deadends
            if (deadendPsms.Any())
            {
                string writtenFileDeadend = Path.Combine(outputFolder, "Deadends" + ".tsv");
                WriteFile.WritePsmCrossToTsv(deadendPsms, writtenFileDeadend, 1);
                FinishedWritingFile(writtenFileDeadend, new List<string> { taskId });
            }
        }

        public void ComputeXlinkQandPValues(List<CrosslinkSpectralMatch> allPsms, List<CrosslinkSpectralMatch> intraCsms, List<CrosslinkSpectralMatch> interCsms, CommonParameters commonParameters, string taskId)
        {
            List<CrosslinkSpectralMatch> crossCsms = allPsms.Where(p => p.CrossType == PsmCrossType.Inter || p.CrossType == PsmCrossType.Intra).OrderByDescending(p => p.XLTotalScore).ToList();
            new FdrAnalysisEngine(crossCsms.ToList<PeptideSpectralMatch>(), 0, commonParameters, this.FileSpecificParameters, new List<string> { taskId }, "crosslink").Run();

            List<CrosslinkSpectralMatch> singles = allPsms.Where(p => p.CrossType != PsmCrossType.Inter).Where(p => p.CrossType != PsmCrossType.Intra).OrderByDescending(p => p.Score).ToList();
            new FdrAnalysisEngine(singles.ToList<PeptideSpectralMatch>(), 0, commonParameters, this.FileSpecificParameters, new List<string> { taskId }, "PSM").Run();
            SingleFDRAnalysis(singles, commonParameters, new List<string> { taskId });

            // calculate FDR
            DoCrosslinkFdrAnalysis(interCsms);
            DoCrosslinkFdrAnalysis(intraCsms);
        }

        //Calculate the FDR of single peptide FP/TP
        private void SingleFDRAnalysis(List<CrosslinkSpectralMatch> items, CommonParameters commonParameters, List<string> taskIds)
        {
            // calculate single PSM FDR
            List<PeptideSpectralMatch> psms = items.Where(p => p.CrossType == PsmCrossType.Single).Select(p => p as PeptideSpectralMatch).OrderByDescending(p => p.Score).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, this.FileSpecificParameters, taskIds, "skippep").Run();

            // calculate loop PSM FDR
            psms = items.Where(p => p.CrossType == PsmCrossType.Loop).Select(p => p as PeptideSpectralMatch).OrderByDescending(p => p.Score).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, this.FileSpecificParameters, taskIds, "skippep").Run();

            // calculate deadend FDR
            psms = items.Where(p => p.CrossType == PsmCrossType.DeadEnd ||
                p.CrossType == PsmCrossType.DeadEndH2O ||
                p.CrossType == PsmCrossType.DeadEndNH2 ||
                p.CrossType == PsmCrossType.DeadEndTris).Select(p => p as PeptideSpectralMatch).OrderByDescending(p => p.Score).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, this.FileSpecificParameters, taskIds, "skippep").Run();
        }

        //Calculate the FDR of crosslinked peptide FP/TP
        private void DoCrosslinkFdrAnalysis(List<CrosslinkSpectralMatch> csms)
        {
            int cumulativeTarget = 0;
            int cumulativeDecoy = 0;

            for (int i = 0; i < csms.Count; i++)
            {
                var csm = csms[i];
                if (csm.IsDecoy || (csm.BetaPeptide!=null && csm.BetaPeptide.IsDecoy))
                {
                    cumulativeDecoy++;
                }
                else
                {
                    cumulativeTarget++;
                }

                double qValue = Math.Min(1, (double)cumulativeDecoy / (cumulativeTarget <= 0 ? 1 : cumulativeTarget));
                double qValueNotch = 0; //maybe we should assign this some day?

                double pep = csm.FdrInfo == null ? double.NaN : csm.FdrInfo.PEP;
                double pepQValue = csm.FdrInfo == null ? double.NaN : csm.FdrInfo.PEP_QValue;

                csm.SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, 0, 0, qValueNotch, pep, pepQValue);
            }
        }
    }
}