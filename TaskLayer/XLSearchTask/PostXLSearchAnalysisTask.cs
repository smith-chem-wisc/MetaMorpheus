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

        public MyTaskResults Run(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList, List<CrosslinkSpectralMatch> allPsms, CommonParameters commonParameters, XlSearchParameters xlSearchParameters, List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications, List<string> localizeableModificationTypes, MyTaskResults MyTaskResults)
        {
            foreach (var csm in allPsms)
            {
                csm.ResolveProteinPosAmbiguitiesForXl();
            }

            var allPsmsXL = allPsms.Where(p => p.CrossType == PsmCrossType.Cross).ToList();

            // inter-crosslinks; different proteins are linked
            var interCsms = allPsmsXL.Where(p => !p.IsIntraCsm()).ToList();
            foreach (var item in interCsms)
            {
                item.CrossType = PsmCrossType.Inter;
            }

            // intra-crosslinks; crosslinks within a protein
            var intraCsms = allPsmsXL.Where(p => p.IsIntraCsm()).ToList();
            foreach (var item in intraCsms)
            {
                item.CrossType = PsmCrossType.Intra;
            }

            // calculate FDR
            DoCrosslinkFdrAnalysis(interCsms);
            DoCrosslinkFdrAnalysis(intraCsms);
            SingleFDRAnalysis(allPsms, commonParameters, new List<string> { taskId });

            // write interlink CSMs
            if (interCsms.Any())
            {
                string file = Path.Combine(OutputFolder, "XL_Interlinks.tsv");
                WriteFile.WritePsmCrossToTsv(interCsms, file, 2);
                FinishedWritingFile(file, new List<string> { taskId });
            }
            MyTaskResults.AddNiceText("Target inter-crosslinks within 1% FDR: " + interCsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.BetaPeptide.IsDecoy));

            if (xlSearchParameters.WriteOutputForPercolator)
            {
                var interPsmsXLPercolator = interCsms.Where(p => p.Score >= 2 && p.BetaPeptide.Score >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteFile.WriteCrosslinkToTxtForPercolator(interPsmsXLPercolator, OutputFolder, "XL_Interlinks_Percolator", xlSearchParameters.Crosslinker);
                FinishedWritingFile(Path.Combine(OutputFolder, "XL_Interlinks_Percolator.txt"), new List<string> { taskId });
            }

            // write intralink CSMs
            if (intraCsms.Any())
            {
                string file = Path.Combine(OutputFolder, "XL_Intralinks.tsv");
                WriteFile.WritePsmCrossToTsv(intraCsms, file, 2);
                FinishedWritingFile(file, new List<string> { taskId });
            }
            MyTaskResults.AddNiceText("Target intra-crosslinks within 1% FDR: " + intraCsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.BetaPeptide.IsDecoy));

            if (xlSearchParameters.WriteOutputForPercolator)
            {
                var intraPsmsXLPercolator = intraCsms.Where(p => p.Score >= 2 && p.BetaPeptide.Score >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteFile.WriteCrosslinkToTxtForPercolator(intraPsmsXLPercolator, OutputFolder, "XL_Intralinks_Percolator", xlSearchParameters.Crosslinker);
                FinishedWritingFile(Path.Combine(OutputFolder, "XL_Intralinks_Percolator.txt"), new List<string> { taskId });
            }

            // write single peptides
            var singlePsms = allPsms.Where(p => p.CrossType == PsmCrossType.Single).ToList();
            if (singlePsms.Any())
            {
                string writtenFileSingle = Path.Combine(OutputFolder, "SinglePeptides" + ".tsv");
                WriteFile.WritePsmCrossToTsv(singlePsms, writtenFileSingle, 1);
                FinishedWritingFile(writtenFileSingle, new List<string> { taskId });
            }
            MyTaskResults.AddNiceText("Target single peptides within 1% FDR: " + singlePsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy));

            // write loops
            var loopPsms = allPsms.Where(p => p.CrossType == PsmCrossType.Loop).ToList();
            if (loopPsms.Any())
            {
                string writtenFileLoop = Path.Combine(OutputFolder, "Looplinks" + ".tsv");
                WriteFile.WritePsmCrossToTsv(loopPsms, writtenFileLoop, 1);
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
                WriteFile.WritePsmCrossToTsv(deadendPsms, writtenFileDeadend, 1);
                FinishedWritingFile(writtenFileDeadend, new List<string> { taskId });
            }
            MyTaskResults.AddNiceText("Target deadend peptides within 1% FDR: " + deadendPsms.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy));

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
                    WriteFile.WritePepXML_xl(writeToXml.Where(p => p.FullFilePath == fullFilePath).ToList(), proteinList, dbFilenameList[0].FilePath, variableModifications, fixedModifications, localizeableModificationTypes, OutputFolder, fileNameNoExtension, commonParameters, xlSearchParameters);
                    FinishedWritingFile(Path.Combine(OutputFolder, fileNameNoExtension + ".pep.XML"), new List<string> { taskId });
                }
            }

            return MyTaskResults;
        }

        //Calculate the FDR of single peptide FP/TP
        private void SingleFDRAnalysis(List<CrosslinkSpectralMatch> items, CommonParameters commonParameters, List<string> taskIds)
        {
            // calculate single PSM FDR
            List<PeptideSpectralMatch> psms = items.Where(p => p.CrossType == PsmCrossType.Single).Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, taskIds).Run();

            // calculate loop PSM FDR
            psms = items.Where(p => p.CrossType == PsmCrossType.Loop).Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, taskIds).Run();

            // calculate deadend FDR
            psms = items.Where(p => p.CrossType == PsmCrossType.DeadEnd ||
                p.CrossType == PsmCrossType.DeadEndH2O ||
                p.CrossType == PsmCrossType.DeadEndNH2 ||
                p.CrossType == PsmCrossType.DeadEndTris).Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, taskIds).Run();
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
                csm.SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, 0, 0, 0, 0, 0);
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