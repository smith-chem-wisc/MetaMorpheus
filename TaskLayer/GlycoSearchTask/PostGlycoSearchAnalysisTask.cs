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
            if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.NGlycanSearch)
            {
                var allPsmsSingle = allPsms.Where(p => p.NGlycan == null && p.Score > 2).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsGly = allPsms.Where(p => p.NGlycan != null && p.Score > 2).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                var writtenFileInter2 = Path.Combine(OutputFolder, "glyco_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileInter2, 2);

                return MyTaskResults;
            }
            else if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch)
            {
                var allPsmsSingle = allPsms.Where(p => p.Routes == null && p.Score > 2).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_psm" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsGly = allPsms.Where(p => p.Routes != null && p.Score > 2).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });

                var writtenFileInter2 = Path.Combine(OutputFolder, "glyco_psm" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileInter2, 2);

                var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsGly);
                WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, Path.Combine(OutputFolder, "seen_glyco_localization" + ".tsv"));
                WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, Path.Combine(OutputFolder, "protein_glyco_localization" + ".tsv"));
                return MyTaskResults;
            }
            else
            {

                return MyTaskResults;
            }
        }

        //Calculate the FDR of single peptide FP/TP
        private void SingleFDRAnalysis(List<GlycoSpectralMatch> items, CommonParameters commonParameters, List<string> taskIds)
        {
            // calculate single PSM FDR
            List<PeptideSpectralMatch> psms = items.Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, this.FileSpecificParameters, taskIds).Run();

        }

    }
}

