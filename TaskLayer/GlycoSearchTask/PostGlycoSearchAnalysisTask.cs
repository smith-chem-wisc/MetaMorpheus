using EngineLayer;
using EngineLayer.GlycoSearch;
using EngineLayer.FdrAnalysis;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MathNet.Numerics;

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
                    var gsmsGroupedByFile = filteredGsms.GroupBy(p => p.FullFilePath).ToList();
                    foreach (var fileSpecificGSMs in gsmsGroupedByFile)
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
            }

            var writtenFileSingle = Path.Combine(OutputFolder, "AllPSMs.psmtsv");
            WriteFile.WritePsmGlycoToTsv(filteredGsms, writtenFileSingle, 1);
            FinishedWritingFile(writtenFileSingle, new List<string> { taskId });

            //else if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch)
            //{
            //    var allPsmsGly = allPsms.Where(p => p.Routes != null ).OrderByDescending(p => p.Score).ToList();
            //    SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });

            //    var writtenFileOGlyco = Path.Combine(OutputFolder, "oglyco" + ".psmtsv");
            //    WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileOGlyco, 2);
            //    FinishedWritingFile(writtenFileOGlyco, new List<string> { taskId });

            //var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsGly.Where(p => p.ProteinAccession != null && p.OneBasedStartResidueInProtein.HasValue).ToList());

            //var seen_oglyco_localization_file = Path.Combine(OutputFolder, "seen_oglyco_localization" + ".tsv");
            //WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_oglyco_localization_file);
            //FinishedWritingFile(seen_oglyco_localization_file, new List<string> { taskId });

            //var protein_oglyco_localization_file = Path.Combine(OutputFolder, "protein_oglyco_localization" + ".tsv");
            //WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_oglyco_localization_file);
            //FinishedWritingFile(protein_oglyco_localization_file, new List<string> { taskId });

            //    return MyTaskResults;
            //}
            //else
            //{
            //    var allPsmsNGly = allPsms.Where(p => p.NGlycan != null).OrderByDescending(p => p.Score).ToList();
            //    SingleFDRAnalysis(allPsmsNGly, commonParameters, new List<string> { taskId });

            //    var writtenFileNGlyco = Path.Combine(OutputFolder, "nglyco" + ".psmtsv");
            //    WriteFile.WritePsmGlycoToTsv(allPsmsNGly, writtenFileNGlyco, 3);
            //    FinishedWritingFile(writtenFileNGlyco, new List<string> { taskId });

            //    var allPsmsOGly = allPsms.Where(p => p.Routes != null).OrderByDescending(p => p.Score).ToList();
            //    SingleFDRAnalysis(allPsmsOGly, commonParameters, new List<string> { taskId });

            //    var writtenFileOGlyco = Path.Combine(OutputFolder, "oglyco" + ".psmtsv");
            //    WriteFile.WritePsmGlycoToTsv(allPsmsOGly, writtenFileOGlyco, 2);
            //    FinishedWritingFile(writtenFileOGlyco, new List<string> { taskId });

            //    return MyTaskResults;
            //}
            return MyTaskResults;
        }

        private void DividePsmsIntoGroupsWriteToTsv(GlycoSearchType glycoSearchType, List<GlycoSpectralMatch> gsms, CommonParameters commonParameters, string taskId, string individualFileFolderPath, string individualFileFolder)
        {
            switch (glycoSearchType)
            {
                case GlycoSearchType.OGlycanSearch:
                    var allPsmsGly = gsms.Where(p => p.Routes != null).OrderByDescending(p => p.Score).ToList();
                    if (allPsmsGly.Any())
                    {
                        SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                        var writtenFileOGlyco = Path.Combine(individualFileFolderPath, individualFileFolder + "_oglyco" + ".psmtsv");
                        WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileOGlyco, 2);

                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsGly.Where(p => p.ProteinAccession != null && p.OneBasedStartResidueInProtein.HasValue).ToList());
                        var seen_oglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "_seen_oglyco_localization" + ".tsv");
                        WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_oglyco_localization_file);

                        var protein_oglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "_protein_oglyco_localization" + ".tsv");
                        WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_oglyco_localization_file);
                    }
                    break;
                case GlycoSearchType.NGlycanSearch:
                    break;
                default: //both N and O glycan search combined
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
            WriteProteinGroupsToTsv(ProteinGroups, writtenFile, new List<string> { Parameters.SearchTaskId }, Math.Min(CommonParameters.QValueOutputFilter, CommonParameters.PepQValueOutputFilter));
        }
        private void WriteProteinGroupsToTsv(List<ProteinGroup> proteinGroups, string filePath, List<string> nestedIds, double qValueCutoff)
        {
            if (proteinGroups != null && proteinGroups.Any())
            {
                using (StreamWriter output = new(filePath))
                {
                    output.WriteLine(proteinGroups.First().GetTabSeparatedHeader());
                    for (int i = 0; i < proteinGroups.Count; i++)
                    {
                        if(Parameters.GlycoSearchParameters.WriteDecoys && proteinGroups[i].IsDecoy && proteinGroups[i].QValue <= qValueCutoff)
                        {
                            output.WriteLine(proteinGroups[i]);
                            continue;
                        }
                        else if(Parameters.GlycoSearchParameters.WriteContaminants && proteinGroups[i].IsContaminant && proteinGroups[i].QValue <= qValueCutoff)
                        {
                            output.WriteLine(proteinGroups[i]);
                            continue;
                        }
                        else if(proteinGroups[i].QValue <= qValueCutoff)
                        {
                            output.WriteLine(proteinGroups[i]);
                        }
                    }
                }

                FinishedWritingFile(filePath, nestedIds);
            }
        }
    }
}

