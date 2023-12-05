using EngineLayer;
using EngineLayer.GlycoSearch;
using EngineLayer.FdrAnalysis;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MathNet.Numerics;
using System.Diagnostics;

namespace TaskLayer
{
    public class PostGlycoSearchAnalysisTask : MetaMorpheusTask
    {
        public PostGlycoSearchAnalysisParameters Parameters { get; set; }
        private List<ProteinGroup> ProteinGroups { get; set; }
        public PostGlycoSearchAnalysisTask() : base(MyTask.GlycoSearch)
        {
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            return null;
        }

        public MyTaskResults Run(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList, List<GlycoSpectralMatch> allPsms, CommonParameters commonParameters, GlycoSearchParameters glycoSearchParameters, List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications, List<string> localizeableModificationTypes, MyTaskResults MyTaskResults)
        {
            if (!Parameters.GlycoSearchParameters.WriteDecoys)
            {
                allPsms.RemoveAll(b => b.IsDecoy);
            }
            if (!Parameters.GlycoSearchParameters.WriteContaminants)
            {
                allPsms.RemoveAll(b => b.IsContaminant);
            }
            var allPSMs = allPsms.OrderByDescending(p => p.Score).ToList();

            SingleFDRAnalysis(allPSMs, commonParameters, new List<string> { taskId });

            List<GlycoSpectralMatch> filteredGsms = allPSMs.Where(p => p.FdrInfo.QValue < 0.01).ToList();

            if (Parameters.GlycoSearchParameters.DoParsimony)
            {
                if (Parameters.GlycoSearchParameters.WriteIndividualFiles)
                {
                    string individualFileResults = Path.Combine(OutputFolder, "IndividualFileResults");
                    if (!Directory.Exists(individualFileResults))
                    {
                        Directory.CreateDirectory(individualFileResults);
                    }

                    foreach (var fileSpecificGSMs in filteredGsms.GroupBy(p => p.FullFilePath))
                    {
                        string individualFileFolder = Path.GetFileNameWithoutExtension(fileSpecificGSMs.Key);
                        string individualFileFolderPath = Path.Combine(individualFileResults, individualFileFolder);
                        if (!Directory.Exists(individualFileFolderPath))
                        {
                            Directory.CreateDirectory(individualFileFolderPath);
                        }
                        var fsgList = fileSpecificGSMs.ToList();
                        GlycoProteinAnalysis(fsgList, individualFileFolderPath, individualFileFolder);
                        foreach (GlycoSpectralMatch gsm in fsgList) //maybe this needs to be the filterd list???
                        {
                            gsm.ResolveAllAmbiguities();
                        }
                        var individualFilePsmsPath = Path.Combine(individualFileFolderPath, individualFileFolder + "_AllPSMs.psmtsv");
                        WriteFile.WritePsmGlycoToTsv(fsgList, individualFilePsmsPath, 1);
                        DividePsmsIntoGroupsWriteToTsv(glycoSearchParameters.GlycoSearchType, fsgList, commonParameters, taskId, individualFileFolderPath, individualFileFolder);
                    }
                }

                GlycoProteinAnalysis(filteredGsms, OutputFolder);//Do the whole group last so inference is done on the whole group
                foreach (GlycoSpectralMatch gsm in allPSMs) //maybe this needs to be the filterd list???
                {
                    gsm.ResolveAllAmbiguities();
                }
                DividePsmsIntoGroupsWriteToTsv(glycoSearchParameters.GlycoSearchType, allPSMs, commonParameters, taskId, OutputFolder, null);
            }
            else
            {
                foreach (GlycoSpectralMatch gsm in allPSMs) //maybe this needs to be the filterd list???
                {
                    gsm.ResolveAllAmbiguities();
                }
                DividePsmsIntoGroupsWriteToTsv(glycoSearchParameters.GlycoSearchType, allPSMs, commonParameters, taskId, OutputFolder, null);
            }
            var writtenFileSingle = Path.Combine(OutputFolder, "AllPSMs.psmtsv");
            WriteFile.WritePsmGlycoToTsv(filteredGsms, writtenFileSingle, 1);
            FinishedWritingFile(writtenFileSingle, new List<string> { taskId });

            return MyTaskResults;
        }

        private void DividePsmsIntoGroupsWriteToTsv(GlycoSearchType glycoSearchType, List<GlycoSpectralMatch> gsms, CommonParameters commonParameters, string taskId, string individualFileFolderPath, string individualFileFolder)
        {
            switch (glycoSearchType)
            {
                case GlycoSearchType.OGlycanSearch:
                    var allPsmsOgly = gsms.Where(p => p.Routes != null).OrderByDescending(p => p.Score).ToList();
                    if (allPsmsOgly.Any())
                    {
                        SingleFDRAnalysis(allPsmsOgly, commonParameters, new List<string> { taskId });
                        var writtenFileOGlyco = Path.Combine(individualFileFolderPath, individualFileFolder + "oglyco" + ".psmtsv");
                        WriteFile.WritePsmGlycoToTsv(allPsmsOgly, writtenFileOGlyco, 2);

                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsOgly.Where(p => p.ProteinAccession != null && p.OneBasedStartResidueInProtein.HasValue).ToList());
                        var seen_oglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "seen_oglyco_localization" + ".tsv");
                        WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_oglyco_localization_file);

                        var protein_oglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "protein_oglyco_localization" + ".tsv");
                        WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_oglyco_localization_file);
                    }
                    break;
                case GlycoSearchType.NGlycanSearch:
                    var allPsmsNgly = gsms.Where(p => p.NGlycan != null).OrderByDescending(p => p.Score).ToList();
                    if (allPsmsNgly.Any())
                    {
                        SingleFDRAnalysis(allPsmsNgly, commonParameters, new List<string> { taskId });
                        var writtenFileOGlyco = Path.Combine(individualFileFolderPath, individualFileFolder + "nglyco" + ".psmtsv");
                        WriteFile.WritePsmGlycoToTsv(allPsmsNgly, writtenFileOGlyco, 2);

                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsNgly.Where(p => p.ProteinAccession != null && p.OneBasedStartResidueInProtein.HasValue).ToList());
                        var seen_nglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "seen_nglyco_localization" + ".tsv");
                        WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_nglyco_localization_file);

                        var protein_oglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "protein_nglyco_localization" + ".tsv");
                        WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_oglyco_localization_file);
                    }
                    break;
                default: //both N and O glycan search combined
                    var allPsmsGly = gsms.Where(p => p.NGlycan != null || p.Routes != null).OrderByDescending(p => p.Score).ToList();
                    if (allPsmsGly.Any())
                    {
                        SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                        var writtenFileOGlyco = Path.Combine(individualFileFolderPath, individualFileFolder + "nglyco" + ".psmtsv");
                        WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileOGlyco, 2);

                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsGly.Where(p => p.ProteinAccession != null && p.OneBasedStartResidueInProtein.HasValue).ToList());
                        var seen_oglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "seen_nglyco_localization" + ".tsv");
                        WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_oglyco_localization_file);

                        var protein_oglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "protein_nglyco_localization" + ".tsv");
                        WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_oglyco_localization_file);
                    }
                    break;
            }
        }

        //Calculate the FDR of single peptide FP/TP
        private void SingleFDRAnalysis(List<GlycoSpectralMatch> items, CommonParameters commonParameters, List<string> taskIds)
        {
            // calculate single PSM FDR
            List<PeptideSpectralMatch> psms = items.Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, this.FileSpecificParameters, taskIds).Run();

        }
        private void GlycoProteinAnalysis(List<GlycoSpectralMatch> gsms, string outputFolder, string individualFileFolder = null)
        {
            // convert gsms to psms
            List<PeptideSpectralMatch> psmsForProteinParsimony = gsms.Select(p => p as PeptideSpectralMatch).ToList();

            Status("Constructing protein groups...", Parameters.SearchTaskId);

            // run parsimony
            ProteinParsimonyResults proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(psmsForProteinParsimony, Parameters.GlycoSearchParameters.ModPeptidesAreDifferent, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run());

            // score protein groups and calculate FDR
            ProteinScoringAndFdrResults proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psmsForProteinParsimony,
                Parameters.GlycoSearchParameters.NoOneHitWonders, Parameters.GlycoSearchParameters.ModPeptidesAreDifferent, true, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run();

            ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

            Status("Done constructing protein groups!", Parameters.SearchTaskId);
            WriteProteinResults(outputFolder, individualFileFolder);
        }

        private void WriteProteinResults(string outputFolder, string individualFileFolder = null)
        {
            string fileName = "AllProteinGroups.tsv";
            string writtenFile = Path.Combine(outputFolder, individualFileFolder + "_"+ fileName);
            WriteProteinGroupsToTsv(ProteinGroups, writtenFile, new List<string> { Parameters.SearchTaskId });
        }
        private void WriteProteinGroupsToTsv(List<ProteinGroup> proteinGroups, string filePath, List<string> nestedIds)
        {
            if (proteinGroups != null && proteinGroups.Any())
            {
                using (StreamWriter output = new(filePath))
                {
                    output.WriteLine(proteinGroups.First().GetTabSeparatedHeader());
                    for (int i = 0; i < proteinGroups.Count; i++)
                    {
                        if ((!Parameters.GlycoSearchParameters.WriteDecoys && proteinGroups[i].IsDecoy) 
                            || (!Parameters.GlycoSearchParameters.WriteContaminants && proteinGroups[i].IsContaminant))
                        {
                            continue;
                        }
                        else
                        {
                            output.WriteLine(proteinGroups[i]);
                        }
                    }
                }
            }
        }
    }
}

