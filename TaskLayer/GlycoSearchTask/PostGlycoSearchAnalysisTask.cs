using EngineLayer;
using EngineLayer.GlycoSearch;
using EngineLayer.FdrAnalysis;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public class PostGlycoSearchAnalysisTask : MetaMorpheusTask
    {
        public PostGlycoSearchAnalysisTask() : base(MyTask.Search)
        {
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            return null;
        }

        public MyTaskResults Run(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList, List<GlycoSpectralMatch> allPsms, CommonParameters commonParameters, GlycoSearchParameters glycoSearchParameters, List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications, List<string> localizeableModificationTypes, MyTaskResults MyTaskResults)
        {
            foreach (var csm in allPsms)
            {
                //TO DO: change this for glycopeptide
                //csm.ResolveProteinPosAmbiguitiesForXl();
            }

            if (glycoSearchParameters.OpenSearchType == OpenSearchType.NGlyco)
            {
                //SingleFDRAnalysis(allPsms, new List<string> { taskId });             
                //var writtenFileInter = Path.Combine(OutputFolder, "all_fdr" + ".mytsv");
                //WritePsmCrossToTsv(allPsms, writtenFileInter, 3);

                var allPsmsSingle = allPsms.Where(p => p.Glycan == null && p.Score > 2).OrderByDescending(p => p.XLTotalScore).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_fdr" + ".mytsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                //TO DO: there may have a bug. I have to filter the following loopPsms, deadendPsms with a BestScore higher than 2, Or some of the Psms will have everything be 0!
                var allPsmsGly = allPsms.Where(p => p.Glycan != null && p.Score > 2).OrderByDescending(p => p.XLTotalScore).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                var writtenFileInter2 = Path.Combine(OutputFolder, "glyco_fdr" + ".mytsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileInter2, 3);

                //var allPsmsGlyForFDR = allPsms.SelectMany(p => p.crosslinkSpectralMatches).ToList();
                //var writtenFileInter3 = Path.Combine(OutputFolder, "glyco_for_fdr" + ".mytsv");
                //WritePsmCrossToTsv(allPsmsGlyForFDR, writtenFileInter3, 3);

                return MyTaskResults;
            }

            if (glycoSearchParameters.OpenSearchType == OpenSearchType.OGlyco)
            {
                var allPsmsSingle = allPsms.Where(p => p.glycanBoxes == null && p.Score > 2).OrderByDescending(p => p.XLTotalScore).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_fdr" + ".mytsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsGly = allPsms.Where(p => p.glycanBoxes != null && p.Score > 2).OrderByDescending(p => p.XLTotalScore).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                var writtenFileInter2 = Path.Combine(OutputFolder, "glyco_fdr" + ".mytsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileInter2, 3);
                return MyTaskResults;
            }




            return MyTaskResults;
        }

        //Calculate the FDR of single peptide FP/TP
        private void SingleFDRAnalysis(List<GlycoSpectralMatch> items, CommonParameters commonParameters, List<string> taskIds)
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

    }
}

