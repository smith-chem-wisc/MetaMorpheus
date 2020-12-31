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
            var allPsmsSingle = allPsms.Where(p =>  p.LocalizationGraphs == null).OrderByDescending(p => p.Score).ToList();
            SingleFDRAnalysis(allPsmsSingle, CommonParameters, new List<string> { Parameters.SearchTaskId });
            //var allSinglePsmsFdr = allPsmsSingle.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01).ToList();
            var allSinglePsmsFdr = allPsmsSingle.Where(p => p.FdrInfo.QValue <= 0.05).ToList();

            var writtenFileSingle = Path.Combine(Parameters.OutputFolder, "single" + ".psmtsv");
            WriteFile.WritePsmGlycoToTsv(allSinglePsmsFdr, writtenFileSingle, 1);
            FinishedWritingFile(writtenFileSingle, new List<string> { Parameters.SearchTaskId });

            List<GlycoSpectralMatch> allgsms = new List<GlycoSpectralMatch>();

            if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch || glycoSearchParameters.GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                var allPsmsOGly = allPsms.Where(p => p.GlycanType == GlycoType.OGlycoPep).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsOGly, CommonParameters, new List<string> { Parameters.SearchTaskId });
                //var allOgsmsFdr = allPsmsOGly.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01).ToList();
                var allOgsmsFdr = allPsmsOGly.Where(p => p.FdrInfo.QValue <= 0.05).ToList();
                OGlycoLocalizationCalculation(allOgsmsFdr, CommonParameters);
                allgsms.AddRange(allOgsmsFdr);

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
                var allNgsmsFdr = allPsmsNGly.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.05).ToList();
                //NGlycoLocalizationCalculation(allNgsmsFdr, CommonParameters);
                OGlycoLocalizationCalculation(allNgsmsFdr, CommonParameters);
                allgsms.AddRange(allNgsmsFdr);

                var writtenFileNGlyco = Path.Combine(Parameters.OutputFolder, "nglyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allNgsmsFdr, writtenFileNGlyco, 3);
                FinishedWritingFile(writtenFileNGlyco, new List<string> { Parameters.SearchTaskId });
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

        //OGlyco Localization
        private static void OGlycoLocalizationCalculation(List<GlycoSpectralMatch> Ogsms, CommonParameters CommonParameters)
        {
            foreach (var glycoSpectralMatch in Ogsms)
            {
                if (glycoSpectralMatch.LocalizationGraphs == null)
                {
                    continue;
                }

                if (glycoSpectralMatch.LocalizationGraphs != null)
                {
                    bool is_HCD_only_data = !GlycoPeptides.DissociationTypeContainETD(CommonParameters.DissociationType) && !GlycoPeptides.DissociationTypeContainETD(CommonParameters.MS2ChildScanDissociationType);

                    if (is_HCD_only_data)
                    {
                        glycoSpectralMatch.LocalizationLevel = LocalizationLevel.Level3;
                        if (glycoSpectralMatch.LocalizationGraphs.Count == 1 && glycoSpectralMatch.LocalizationGraphs.First().ModPos.Length == 1)
                        {
                            glycoSpectralMatch.LocalizationLevel = LocalizationLevel.Level1b;
                        }

                    }
                    else
                    {
                        List<Route> localizationCandidates = new List<Route>();

                        for (int i = 0; i < glycoSpectralMatch.LocalizationGraphs.Count; i++)
                        {
                            var allPathWithMaxScore = LocalizationGraph.GetAllHighestScorePaths(glycoSpectralMatch.LocalizationGraphs[i].array, glycoSpectralMatch.LocalizationGraphs[i].ChildModBoxes);

                            foreach (var path in allPathWithMaxScore)
                            {
                                var local = LocalizationGraph.GetLocalizedPath(glycoSpectralMatch.LocalizationGraphs[i], path);
                                local.ModBoxId = glycoSpectralMatch.LocalizationGraphs[i].ModBoxId;
                                localizationCandidates.Add(local);
                            }
                        }

                        glycoSpectralMatch.Routes = localizationCandidates;
                    }
                }

                if (glycoSpectralMatch.Routes != null)
                {
                    LocalizationLevel localLevel;
                    glycoSpectralMatch.LocalizedGlycan = GlycoSpectralMatch.GetLocalizedGlycan(glycoSpectralMatch.Routes, out localLevel);
                    glycoSpectralMatch.LocalizationLevel = localLevel;

                    //Localization PValue.
                    if (localLevel == LocalizationLevel.Level1 || localLevel == LocalizationLevel.Level2)
                    {
                        List<Route> allRoutes = new List<Route>();
                        foreach (var graph in glycoSpectralMatch.LocalizationGraphs)
                        {
                            allRoutes.AddRange(LocalizationGraph.GetAllPaths_CalP(graph, glycoSpectralMatch.ScanInfo_p, glycoSpectralMatch.Thero_n));
                        }
                        glycoSpectralMatch.SiteSpeciLocalProb = LocalizationGraph.CalSiteSpecificLocalizationProbability(allRoutes, glycoSpectralMatch.LocalizationGraphs.First().ModPos);
                    }
                }

                CorrectLocalizationLevel(glycoSpectralMatch);

            }
        }

        //Correct Localization Level based on site specific probability. If LocalizationLevel = 1, and there are site probability lower than 0.75, Correct the level to 1b.
        private static void CorrectLocalizationLevel(GlycoSpectralMatch gsm)
        {
            if (gsm.SiteSpeciLocalProb == null || gsm.LocalizationLevel != LocalizationLevel.Level1)
            {
                return;
            }

            if (gsm.LocalizationGraphs.First().ModPos.Length == 1 && gsm.LocalizationGraphs.First().TotalScore == 0)
            {
                gsm.LocalizationLevel = LocalizationLevel.Level1b;
                return;
            }


            for (int i = 0; i < gsm.LocalizedGlycan.Count; i++)
            {
                var g = gsm.LocalizedGlycan[i];
                if (gsm.SiteSpeciLocalProb[g.Item1].Where(p => p.Item1 == g.Item2).First().Item2 < 0.75)
                {
                    gsm.LocalizationLevel = LocalizationLevel.Level1b;
                    return;
                }

                if (!gsm.Routes.First().Mods[i].Item3)
                {
                    gsm.LocalizationLevel = LocalizationLevel.Level1b;
                    return;
                }
            }
        }

        //Deprecated function. 
        //Glyco Localization
        private static void NGlycoLocalizationCalculation(List<GlycoSpectralMatch> gsms, CommonParameters CommonParameters)
        {
            foreach (var glycoSpectralMatch in gsms)
            {
                if (glycoSpectralMatch.LocalizationGraphs == null)
                {
                    continue;
                }

                if (glycoSpectralMatch.LocalizationGraphs != null)
                {
                    List<Route> localizationCandidates = new List<Route>();

                    for (int i = 0; i < glycoSpectralMatch.LocalizationGraphs.Count; i++)
                    {
                        var allPathWithMaxScore = LocalizationGraph.GetAllHighestScorePaths(glycoSpectralMatch.LocalizationGraphs[i].array, glycoSpectralMatch.LocalizationGraphs[i].ChildModBoxes);

                        foreach (var path in allPathWithMaxScore)
                        {
                            var local = LocalizationGraph.GetLocalizedPath(glycoSpectralMatch.LocalizationGraphs[i], path);
                            local.ModBoxId = glycoSpectralMatch.LocalizationGraphs[i].ModBoxId;
                            localizationCandidates.Add(local);
                        }
                    }

                    glycoSpectralMatch.Routes = localizationCandidates;

                }

                if (glycoSpectralMatch.Routes != null)
                {
                    LocalizationLevel localLevel;
                    glycoSpectralMatch.LocalizedGlycan = GlycoSpectralMatch.GetLocalizedGlycan(glycoSpectralMatch.Routes, out localLevel);
                    glycoSpectralMatch.LocalizationLevel = localLevel;

                    //Localization PValue.
                    if (localLevel == LocalizationLevel.Level1 || localLevel == LocalizationLevel.Level2)
                    {
                        List<Route> allRoutes = new List<Route>();
                        foreach (var graph in glycoSpectralMatch.LocalizationGraphs)
                        {
                            allRoutes.AddRange(LocalizationGraph.GetAllPaths_CalP(graph, glycoSpectralMatch.ScanInfo_p, glycoSpectralMatch.Thero_n));
                        }
                        glycoSpectralMatch.SiteSpeciLocalProb = LocalizationGraph.CalSiteSpecificLocalizationProbability(allRoutes, glycoSpectralMatch.LocalizationGraphs.First().ModPos);
                    }
                }

                CorrectLocalizationLevel(glycoSpectralMatch);

            }
        }

        //For N-Glycopeptide localzation Level.
        public static void NGlycoCorrectLocalizationLevel(List<GlycoSpectralMatch> ngsms)
        {
            foreach (var gsm in ngsms)
            {
                //Correct nglycoSpectrumMatch localizationLevel to Level1b based on matched fragment ions.
                if (gsm.ModPos.Count() == 1)
                {
                    gsm.LocalizationLevel = LocalizationLevel.Level1b;
                    //correct localization.
                    if (gsm.MatchedFragmentIons.Where(p => (p.NeutralTheoreticalProduct.ProductType == ProductType.b || p.NeutralTheoreticalProduct.ProductType == ProductType.y) && p.NeutralTheoreticalProduct.NeutralLoss > 0).Count() > 0
                        || gsm.MatchedFragmentIons.Where(p => (p.NeutralTheoreticalProduct.ProductType == ProductType.c && p.NeutralTheoreticalProduct.AminoAcidPosition >= gsm.ModPos[0] - 1) || (p.NeutralTheoreticalProduct.ProductType == ProductType.zDot && p.NeutralTheoreticalProduct.AminoAcidPosition <= gsm.ModPos[0] - 1)).Count() > 0)
                    {
                        gsm.LocalizationLevel = LocalizationLevel.Level1;
                    }
                    else
                    {
                        foreach (var childIons in gsm.ChildMatchedFragmentIons)
                        {
                            if (childIons.Value.Where(p => (p.NeutralTheoreticalProduct.ProductType == ProductType.b || p.NeutralTheoreticalProduct.ProductType == ProductType.y) && p.NeutralTheoreticalProduct.NeutralLoss > 0).Count() > 0
                                || childIons.Value.Where(p => (p.NeutralTheoreticalProduct.ProductType == ProductType.c && p.NeutralTheoreticalProduct.AminoAcidPosition >= gsm.ModPos[0] - 1) || (p.NeutralTheoreticalProduct.ProductType == ProductType.zDot && p.NeutralTheoreticalProduct.AminoAcidPosition <= gsm.ModPos[0] - 1)).Count() > 0)
                            {
                                gsm.LocalizationLevel = LocalizationLevel.Level1;
                            }
                        }
                    }
                }
            }
        }

    }
}

