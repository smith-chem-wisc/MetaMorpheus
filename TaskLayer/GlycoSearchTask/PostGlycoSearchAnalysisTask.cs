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

            if (!glycoSearchParameters.IsOGlycoSearch)
            {

                var allPsmsSingle = allPsms.Where(p => p.Glycan == null && p.Score > 2).OrderByDescending(p => p.TotalScore).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsGly = allPsms.Where(p => p.Glycan != null && p.Score > 2).OrderByDescending(p => p.TotalScore).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                var writtenFileInter2 = Path.Combine(OutputFolder, "glyco_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileInter2, 2);

                return MyTaskResults;
            }
            else
            {
                var allPsmsSingle = allPsms.Where(p => p.glycanBoxes == null && p.Score > 2).OrderByDescending(p => p.TotalScore).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsGly = allPsms.Where(p => p.glycanBoxes != null && p.Score > 2).OrderByDescending(p => p.TotalScore).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                var writtenFileInter2 = Path.Combine(OutputFolder, "glyco_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileInter2, 2);
                return MyTaskResults;
            }
        }

        //Calculate the FDR of single peptide FP/TP
        private void SingleFDRAnalysis(List<GlycoSpectralMatch> items, CommonParameters commonParameters, List<string> taskIds)
        {
            // calculate single PSM FDR
            List<PeptideSpectralMatch> psms = items.Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, taskIds).Run();

        }

    }
}

