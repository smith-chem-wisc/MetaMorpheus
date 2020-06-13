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
            var allPsmsSingle = allPsms.Where(p => p.NGlycan == null).OrderByDescending(p => p.Score).ToList();
            SingleFDRAnalysis(allPsmsSingle, CommonParameters, new List<string> { Parameters.SearchTaskId });
            var allSinglePsmsFdr = allPsmsSingle.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01).ToList();

            var writtenFileSingle = Path.Combine(Parameters.OutputFolder, "single" + ".psmtsv");
            WriteFile.WritePsmGlycoToTsv(allSinglePsmsFdr, writtenFileSingle, 1);
            FinishedWritingFile(writtenFileSingle, new List<string> { Parameters.SearchTaskId });

            List<GlycoSpectralMatch> allgsms = new List<GlycoSpectralMatch>();

            if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.NGlycanSearch || glycoSearchParameters.GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                var allPsmsNGly = allPsms.Where(p => p.NGlycan != null).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsNGly, CommonParameters, new List<string> { Parameters.SearchTaskId });
                var allNgsmsFdr = allPsmsNGly.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01).ToList();
                NGlycoCorrectLocalizationLevel(allNgsmsFdr);
                allgsms.AddRange(allNgsmsFdr);

                var writtenFileNGlyco = Path.Combine(Parameters.OutputFolder, "nglyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allNgsmsFdr, writtenFileNGlyco, 3);
                FinishedWritingFile(writtenFileNGlyco, new List<string> { Parameters.SearchTaskId });
            }

            if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch || glycoSearchParameters.GlycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                var allPsmsOGly = allPsms.Where(p => p.LocalizationGraphs != null).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsOGly, CommonParameters, new List<string> { Parameters.SearchTaskId });
                var allOgsmsFdr = allPsmsOGly.Where(p => !p.IsDecoy && p.FdrInfo.QValue <= 0.01).ToList();
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

            if (glycoSearchParameters.PerformQuantification || glycoSearchParameters.PerformRelativeRetentionTimePrediction)
            {
                QuantificationAnalysis(allgsms);
                if (glycoSearchParameters.PerformQuantification)
                {
                    WriteQuantificationResults();
                }
                if (glycoSearchParameters.PerformRelativeRetentionTimePrediction)
                {
                    var gps = GlycopeptideRRT.ConstructGlycopeptideRRT(allgsms, Parameters.FlashLfqResults);
                    string predictMessage = GlycopeptideRRT.PredictRRT(gps);
                    MyTaskResults.AddTaskSummaryText(predictMessage);

                    var writtenFileRRT = Path.Combine(Parameters.OutputFolder, "glycopeptideRT" + ".tsv");
                    GlycopeptideRRT.WriteGlycopeptideRRT(gps, writtenFileRRT);
                }
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

        //Quantification analysis
        private void QuantificationAnalysis(List<GlycoSpectralMatch> allPsms)
        {
            // pass quantification parameters to FlashLFQ
            Status("Quantifying...", Parameters.SearchTaskId);

            // get PSMs to pass to FlashLFQ
            var unambiguousPsmsBelowOnePercentFdr = allPsms.Where(p =>
                p.FdrInfo.QValue <= 0.01
                && p.FdrInfo.QValueNotch <= 0.01
                && !p.IsDecoy
                && p.FullSequence != null).ToList(); //if ambiguous, there's no full sequence

            //group psms by file
            var psmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy(p => p.FullFilePath);

            // construct file info for FlashLFQ
            var spectraFileInfo = new List<SpectraFileInfo>();
            foreach (var file in Parameters.CurrentRawFileList)
            {
                // experimental design info passed in here for each spectra file
                spectraFileInfo.Add(new SpectraFileInfo(fullFilePathWithExtension: file, condition: "", biorep: 0, fraction: 0, techrep: 0));
                Parameters.MyFileManager.DoneWithFile(file);
            }

            // pass PSM info to FlashLFQ
            var flashLFQIdentifications = new List<Identification>();


            foreach (var spectraFile in psmsGroupedByFile)
            {
                var rawfileinfo = spectraFileInfo.Where(p => p.FullFilePathWithExtension.Equals(spectraFile.Key)).First();

                foreach (var psm in spectraFile)
                {
                    flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence,
                        psm.PeptideMonisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, new List<FlashLFQ.ProteinGroup>(), useForProteinQuant:false));
                }
            }

            // run FlashLFQ
            var FlashLfqEngine = new FlashLfqEngine(
                allIdentifications: flashLFQIdentifications,
                normalize: false,
                ppmTolerance: 10,
                matchBetweenRuns: false,
                silent: true,
                maxThreads: CommonParameters.MaxThreadsToUsePerFile);

            if (flashLFQIdentifications.Any())
            {
                Parameters.FlashLfqResults = FlashLfqEngine.Run();
            }


        }

        //The three functions here are almost the same as in PostSearchAnalysisTask. It is possible to use the same functions.
        private void WriteQuantificationResults()
        {
            //if (Parameters.SearchParameters.DoQuantification && Parameters.FlashLfqResults != null)
            {
                // write peaks
                WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "AllQuantifiedPeaks", new List<string> { Parameters.SearchTaskId });

                // write peptide quant results
                string filename = "AllQuantified" + GlobalVariables.AnalyteType + "s";
                WritePeptideQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, filename, new List<string> { Parameters.SearchTaskId });

                //// write individual results
                //if (Parameters.CurrentRawFileList.Count > 1)
                //{
                //    foreach (var file in Parameters.FlashLfqResults.Peaks)
                //    {
                //        WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.IndividualResultsOutputFolder,
                //            file.Key.FilenameWithoutExtension + "_QuantifiedPeaks", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.Key.FullFilePathWithExtension });
                //    }
                //}
            }
        }

        private void WritePeptideQuantificationResultsToTsv(FlashLfqResults flashLFQResults, string outputFolder, string fileName, List<string> nestedIds)
        {
            var fullSeqPath = Path.Combine(outputFolder, fileName + ".tsv");

            flashLFQResults.WriteResults(null, fullSeqPath, null, null, true);

            FinishedWritingFile(fullSeqPath, nestedIds);
        }

        private void WritePeakQuantificationResultsToTsv(FlashLfqResults flashLFQResults, string outputFolder, string fileName, List<string> nestedIds)
        {
            var peaksPath = Path.Combine(outputFolder, fileName + ".tsv");

            flashLFQResults.WriteResults(peaksPath, null, null, null, true);

            FinishedWritingFile(peaksPath, nestedIds);
        }

    }
}

