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
        private List<EngineLayer.ProteinGroup> ProteinGroups { get; set; }
        public PostGlycoSearchAnalysisTask() : base(MyTask.GlycoSearch)
        {
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            return null;
        }

        public MyTaskResults Run(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList, List<GlycoSpectralMatch> allPsms, CommonParameters commonParameters, GlycoSearchParameters glycoSearchParameters, List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications, List<string> localizeableModificationTypes, MyTaskResults MyTaskResults)
        {
            var allPsmsSingle = allPsms.OrderByDescending(p => p.Score).ToList();

            if (!Parameters.GlycoSearchParameters.WriteDecoys)
            {
                allPsmsSingle.RemoveAll(b => b.IsDecoy);
            }
            if (!Parameters.GlycoSearchParameters.WriteContaminants)
            {
                allPsmsSingle.RemoveAll(b => b.IsContaminant);
            }
            bool qValueUsedAsFilter = !CommonParameters.QValueOutputFilter.IsLarger(CommonParameters.PepQValueOutputFilter, 3);
            string filterType = qValueUsedAsFilter ? "s with q-value = " : "s with pep q-value = ";
            double filterCutoffForResultsCounts = qValueUsedAsFilter ? Math.Min(0.01, CommonParameters.QValueOutputFilter) : Math.Min(0.01, CommonParameters.PepQValueOutputFilter);

            SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });

            var filteredPeptidesForOutput = allPsmsSingle.Where(p => p.FdrInfo.QValue <= CommonParameters.QValueOutputFilter
                && p.FdrInfo.QValueNotch <= CommonParameters.QValueOutputFilter
                && (p.FdrInfo.PEP_QValue <= CommonParameters.PepQValueOutputFilter || double.IsNaN(p.FdrInfo.PEP_QValue))).ToList();

            ProteinAnalysis(allPsms);
            WriteProteinResults(allPsms);
            var writtenFileSingle = Path.Combine(OutputFolder, "single" + ".psmtsv");
            WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileSingle, 1);
            FinishedWritingFile(writtenFileSingle, new List<string> { taskId });

            if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.NGlycanSearch)
            {
                var allPsmsGly = allPsms.Where(p => p.NGlycan != null ).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });

                var writtenFileNGlyco = Path.Combine(OutputFolder, "nglyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileNGlyco, 3);
                FinishedWritingFile(writtenFileNGlyco, new List<string> { taskId });

                return MyTaskResults;
            }
            else if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch)
            {
                var allPsmsGly = allPsms.Where(p => p.Routes != null ).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });

                var writtenFileOGlyco = Path.Combine(OutputFolder, "oglyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileOGlyco, 2);
                FinishedWritingFile(writtenFileOGlyco, new List<string> { taskId });

                var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsGly.Where(p=>p.ProteinAccession!=null && p.OneBasedStartResidueInProtein.HasValue).ToList());

                var seen_oglyco_localization_file = Path.Combine(OutputFolder, "seen_oglyco_localization" + ".tsv");
                WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_oglyco_localization_file);
                FinishedWritingFile(seen_oglyco_localization_file, new List<string> { taskId });

                var protein_oglyco_localization_file = Path.Combine(OutputFolder, "protein_oglyco_localization" + ".tsv");
                WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_oglyco_localization_file);
                FinishedWritingFile(protein_oglyco_localization_file, new List<string> { taskId });

                return MyTaskResults;
            }
            else
            {
                var allPsmsNGly = allPsms.Where(p => p.NGlycan != null).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsNGly, commonParameters, new List<string> { taskId });

                var writtenFileNGlyco = Path.Combine(OutputFolder, "nglyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsNGly, writtenFileNGlyco, 3);
                FinishedWritingFile(writtenFileNGlyco, new List<string> { taskId });

                var allPsmsOGly = allPsms.Where(p => p.Routes != null).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsOGly, commonParameters, new List<string> { taskId });

                var writtenFileOGlyco = Path.Combine(OutputFolder, "oglyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsOGly, writtenFileOGlyco, 2);
                FinishedWritingFile(writtenFileOGlyco, new List<string> { taskId });

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
        private void ProteinAnalysis(List<GlycoSpectralMatch> gsms)
        {
            // convert gsms to psms
            List<PeptideSpectralMatch> allPsms = gsms.Select(p => p as PeptideSpectralMatch).ToList();

            if (!Parameters.GlycoSearchParameters.DoParsimony)
            {
                return;
            }

            Status("Constructing protein groups...", Parameters.SearchTaskId);

            List<PeptideSpectralMatch> psmsForProteinParsimony = allPsms;

            // run parsimony
            ProteinParsimonyResults proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(psmsForProteinParsimony, Parameters.GlycoSearchParameters.ModPeptidesAreDifferent, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run());

            // score protein groups and calculate FDR
            ProteinScoringAndFdrResults proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psmsForProteinParsimony,
                Parameters.GlycoSearchParameters.NoOneHitWonders, Parameters.GlycoSearchParameters.ModPeptidesAreDifferent, true, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run();

            ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

            foreach (PeptideSpectralMatch psm in Parameters.AllPsms)
            {
                psm.ResolveAllAmbiguities();
            }

            Status("Done constructing protein groups!", Parameters.SearchTaskId);
        }

        private void WriteProteinResults(List<GlycoSpectralMatch> gsms)
        {
            if (Parameters.GlycoSearchParameters.DoParsimony)
            {
                string fileName = "AllProteinGroups.tsv";

                // convert gsms to psms
                List<PeptideSpectralMatch> allPsms = gsms.Select(p => p as PeptideSpectralMatch).ToList();

                // get PSMs to pass to FlashLFQ
                List<PeptideSpectralMatch> unambiguousPsmsBelowOnePercentFdr = new();
                if (allPsms.Count > 100)//PEP is not computed when there are fewer than 100 psms
                {
                    unambiguousPsmsBelowOnePercentFdr = allPsms.Where(p =>
                        p.FdrInfo.PEP_QValue <= 0.01
                        && !p.IsDecoy
                        && p.FullSequence != null).ToList(); //if ambiguous, there's no full sequence
                }
                else
                {
                    unambiguousPsmsBelowOnePercentFdr = allPsms.Where(p =>
                        p.FdrInfo.QValue <= 0.01
                        && !p.IsDecoy
                        && p.FullSequence != null).ToList(); //if ambiguous, there's no full sequence
                }

                //group psms by file
                var psmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy(p => p.FullFilePath);

                // write protein groups to tsv
                string writtenFile = Path.Combine(Parameters.OutputFolder, fileName);
                WriteProteinGroupsToTsv(ProteinGroups, writtenFile, new List<string> { Parameters.SearchTaskId }, Math.Min(CommonParameters.QValueOutputFilter, CommonParameters.PepQValueOutputFilter));

                // write all individual file results to subdirectory
                // local protein fdr, global parsimony, global psm fdr
                if (Parameters.CurrentRawFileList.Count > 1 && (Parameters.GlycoSearchParameters.WriteIndividualFiles))
                {
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);
                }

                //write the individual result files for each datafile
                foreach (var fullFilePath in psmsGroupedByFile.Select(v => v.Key))
                {
                    string strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                    List<PeptideSpectralMatch> psmsForThisFile = psmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList();
                    var subsetProteinGroupsForThisFile = ProteinGroups.Select(p => p.ConstructSubsetProteinGroup(fullFilePath)).ToList();

                    ProteinScoringAndFdrResults subsetProteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(subsetProteinGroupsForThisFile, psmsForThisFile,
                        Parameters.GlycoSearchParameters.NoOneHitWonders, Parameters.GlycoSearchParameters.ModPeptidesAreDifferent,
                        false, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }).Run();

                    subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults.SortedAndScoredProteinGroups;

                    Parameters.GlycoSearchTaskResults.AddTaskSummaryText("Target protein groups within 1 % FDR in " + strippedFileName + ": " + subsetProteinGroupsForThisFile.Count(b => b.QValue <= 0.01 && !b.IsDecoy));

                    // write individual spectra file protein groups results to tsv
                    if (Parameters.GlycoSearchParameters.WriteIndividualFiles && Parameters.CurrentRawFileList.Count > 1)
                    {
                        writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_ProteinGroups.tsv");
                        WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }, Math.Min(CommonParameters.QValueOutputFilter, CommonParameters.PepQValueOutputFilter));
                    }

                    ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }));
                }
            }
        }
        private void WriteProteinGroupsToTsv(List<EngineLayer.ProteinGroup> proteinGroups, string filePath, List<string> nestedIds, double qValueCutoff)
        {
            if (proteinGroups != null && proteinGroups.Any())
            {
                using (StreamWriter output = new StreamWriter(filePath))
                {
                    output.WriteLine(proteinGroups.First().GetTabSeparatedHeader());
                    for (int i = 0; i < proteinGroups.Count; i++)
                    {
                        if ((!Parameters.GlycoSearchParameters.WriteDecoys && proteinGroups[i].IsDecoy) ||
                            (!Parameters.GlycoSearchParameters.WriteContaminants && proteinGroups[i].IsContaminant))
                        {
                            continue;
                        }
                        else if (proteinGroups[i].QValue <= qValueCutoff)
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

