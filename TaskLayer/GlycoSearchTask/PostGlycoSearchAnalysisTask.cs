using EngineLayer;
using EngineLayer.GlycoSearch;
using EngineLayer.FdrAnalysis;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MathNet.Numerics.LinearRegression;
using FlashLFQ;
using Proteomics.Fragmentation;

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

        public PostSearchAnalysisParameters Parameters { get; set; }

        private IEnumerable<IGrouping<string, PeptideSpectralMatch>> PsmsGroupedByFile { get; set; }

        public MyTaskResults Run(List<GlycoSpectralMatch> allPsms, GlycoSearchParameters glycoSearchParameters, MyTaskResults MyTaskResults)
        {
            List<GlycoSpectralMatch> allgsms = new List<GlycoSpectralMatch>();

            var allPsmsSingle = allPsms.Where(p =>  p.LocalizationGraphs == null).OrderByDescending(p => p.Score).ToList();
            SingleFDRAnalysis(allPsmsSingle, CommonParameters, new List<string> { Parameters.SearchTaskId });
            //var allSinglePsmsFdr = allPsmsSingle.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01).ToList();
            var allSinglePsmsFdr = allPsmsSingle.Where(p => p.FdrInfo.QValue <= 0.05).ToList();
            allgsms.AddRange(allSinglePsmsFdr.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01));

            var writtenFileSingle = Path.Combine(Parameters.OutputFolder, "single" + ".psmtsv");
            WriteFile.WritePsmGlycoToTsv(allSinglePsmsFdr, writtenFileSingle, 1);
            FinishedWritingFile(writtenFileSingle, new List<string> { Parameters.SearchTaskId });

            if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch || glycoSearchParameters.GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                var allPsmsOGly = allPsms.Where(p => p.GlycanType == GlycoType.OGlycoPep).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsOGly, CommonParameters, new List<string> { Parameters.SearchTaskId });
                //var allOgsmsFdr = allPsmsOGly.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01).ToList();
                var allOgsmsFdr = allPsmsOGly.Where(p => p.FdrInfo.QValue <= 0.1).ToList();
                GlycoSite.GlycoLocalizationCalculation(allOgsmsFdr, CommonParameters);
                allgsms.AddRange(allOgsmsFdr.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01));

                var writtenFileOGlyco = Path.Combine(Parameters.OutputFolder, "oglyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allOgsmsFdr, writtenFileOGlyco, 2);
                FinishedWritingFile(writtenFileOGlyco, new List<string> { Parameters.SearchTaskId });

                bool is_HCD_only_data = !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType);

                if (!is_HCD_only_data)
                {
                    var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allOgsmsFdr.Where(p => p.ProteinAccession != null && p.OneBasedStartResidueInProtein.HasValue).ToList());

                    var seen_oglyco_localization_file = Path.Combine(Parameters.OutputFolder, "seen_oglyco_localization" + ".tsv");
                    WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_oglyco_localization_file);
                    FinishedWritingFile(seen_oglyco_localization_file, new List<string> { Parameters.SearchTaskId });

                    var protein_oglyco_localization_file = Path.Combine(Parameters.OutputFolder, "protein_oglyco_localization" + ".tsv");
                    WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_oglyco_localization_file);
                    FinishedWritingFile(protein_oglyco_localization_file, new List<string> { Parameters.SearchTaskId });
                }
            }
            
            if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.NGlycanSearch || glycoSearchParameters.GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                //TO THINK: a mixed glycopeptide has more properties similar to a NGlycopeptide.
                var allPsmsNGly = allPsms.Where(p => p.GlycanType == GlycoType.NGlycoPep || p.GlycanType == GlycoType.MixedGlycoPep).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsNGly, CommonParameters, new List<string> { Parameters.SearchTaskId });
                //var allNgsmsFdr = allPsmsNGly.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01).ToList();
                var allNgsmsFdr = allPsmsNGly.Where(p => p.FdrInfo.QValue <= 0.1).ToList();
                //NGlycoLocalizationCalculation(allNgsmsFdr, CommonParameters);
                GlycoSite.GlycoLocalizationCalculation(allNgsmsFdr, CommonParameters);
                allgsms.AddRange(allNgsmsFdr.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01));

                var writtenFileNGlyco = Path.Combine(Parameters.OutputFolder, "nglyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allNgsmsFdr, writtenFileNGlyco, 2);
                FinishedWritingFile(writtenFileNGlyco, new List<string> { Parameters.SearchTaskId });
            }

            {
                var writtenFileAllGlyco = Path.Combine(Parameters.OutputFolder, "all_glyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allgsms, writtenFileAllGlyco, 2);
                FinishedWritingFile(writtenFileAllGlyco, new List<string> { Parameters.SearchTaskId });

                var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allgsms.Where(p => p.GlycanType != GlycoType.SinglePep && p.ProteinAccession != null && p.OneBasedStartResidueInProtein.HasValue).ToList());

                var all_seen_glyco_localization_file = Path.Combine(Parameters.OutputFolder, "all_seen_glyco_localization" + ".tsv");
                WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, all_seen_glyco_localization_file);
                FinishedWritingFile(all_seen_glyco_localization_file, new List<string> { Parameters.SearchTaskId });

                var all_protein_glyco_localization = Path.Combine(Parameters.OutputFolder, "all_protein_glyco_localization" + ".tsv");
                WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, all_protein_glyco_localization);
                FinishedWritingFile(all_protein_glyco_localization, new List<string> { Parameters.SearchTaskId });

                var all_protein_info = Path.Combine(Parameters.OutputFolder, "all_protein_info" + ".tsv");
                WriteFile.WriteProteinInfo(allgsms, all_protein_info);
                FinishedWritingFile(all_protein_info, new List<string> { Parameters.SearchTaskId });
            }

            return MyTaskResults;
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

