using EngineLayer;
using EngineLayer.GlycoSearch;
using EngineLayer.FdrAnalysis;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MathNet.Numerics.LinearRegression;

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
                var allPsmsSingle = allPsms.Where(p => p.NGlycan == null ).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsGly = allPsms.Where(p => p.NGlycan != null ).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                var writtenFileInter2 = Path.Combine(OutputFolder, "nglyco_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileInter2, 3);

                return MyTaskResults;
            }
            else if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch)
            {
                var allPsmsSingle = allPsms.Where(p => p.Routes == null ).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_psm" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsGly = allPsms.Where(p => p.Routes != null ).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });

                //var regression = RetentionTimeRegression(allPsmsGly.Where(p=>p.FdrInfo.QValue <= 0.01 && !p.IsContaminant && !p.IsDecoy).ToList());
                //RetentionTimePrediction(allPsmsGly, regression);

                var writtenFileInter2 = Path.Combine(OutputFolder, "oglyco_psm" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileInter2, 2);

                var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsGly.Where(p=>p.ProteinAccession!=null && p.OneBasedStartResidueInProtein.HasValue).ToList());
                WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, Path.Combine(OutputFolder, "seen_oglyco_localization" + ".tsv"));
                WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, Path.Combine(OutputFolder, "protein_oglyco_localization" + ".tsv"));
                return MyTaskResults;
            }
            else
            {
                var allPsmsSingle = allPsms.Where(p => p.NGlycan == null && p.Routes == null).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_psm" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsNGly = allPsms.Where(p => p.NGlycan != null).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsNGly, commonParameters, new List<string> { taskId });
                var writtenFileInter2 = Path.Combine(OutputFolder, "nglyco_psm" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsNGly, writtenFileInter2, 3);

                var allPsmsOGly = allPsms.Where(p => p.Routes != null).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsOGly, commonParameters, new List<string> { taskId });
                var writtenFileInter3 = Path.Combine(OutputFolder, "oglyco_psm" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsOGly, writtenFileInter3, 2);

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

        public double[] RetentionTimeRegression(List<GlycoSpectralMatch> glycoSpectralMatches)
        {
            double[][] xs = new double[glycoSpectralMatches.Count][];

            double[] ys = glycoSpectralMatches.Select(p => p.ScanRetentionTime).ToArray();
       
            for (int i = 0; i < glycoSpectralMatches.Count; i++)
            {
                var glycanKind = GlycanBox.OGlycanBoxes[glycoSpectralMatches[i].Routes.First().ModBoxId].Kind;

                xs[i] = new double[1 + glycanKind.Length];

                xs[i][0] = glycoSpectralMatches[i].PredictedHydrophobicity;

                for (int j = 1; j <= glycanKind.Length; j++)
                {
                    xs[i][j] = glycanKind[j - 1];
                }   
            }

            using (StreamWriter output = new StreamWriter(@"E:\MassData\Glycan\Nick_2019_StcE\Rep1\_temp2\2020-05-19-17-12-20\Task1-GlycoSearchTask\RetentionTime.csv"))
            {
                for (int i = 0; i < xs.Length; i++)
                {
                    string line = "";
                    for (int j = 0; j < xs[0].Length; j++)
                    {
                        line += xs[i][j].ToString() + "\t";
                    }
                    line += ys[i].ToString();
                    output.WriteLine(line);
                }
            }

            double[] p = MultipleRegression.QR(xs, ys, intercept: false);

            //double[] p = Fit.MultiDim(xs, ys, intercept: true);

            return p;
        } 

        public void RetentionTimePrediction(List<GlycoSpectralMatch> glycoSpectralMatches, double[] regression)
        {
            foreach (var gsm in glycoSpectralMatches)
            {
                var glycanKind = GlycanBox.OGlycanBoxes[gsm.Routes.First().ModBoxId].Kind;

                var xs = new double[1 + glycanKind.Length];

                xs[0] = gsm.PredictedHydrophobicity;

                for (int j = 1; j <= glycanKind.Length; j++)
                {
                    xs[j] = glycanKind[j - 1];
                }

                double prt = regression[0];

                for (int j = 1; j <= xs.Length; j++)
                {
                    prt += regression[j] * xs[j - 1];
                }

                gsm.PredictedRT = prt;
            }
        }

    }
}

