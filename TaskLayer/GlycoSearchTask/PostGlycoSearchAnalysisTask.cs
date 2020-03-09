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
            if (!glycoSearchParameters.IsOGlycoSearch)
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
            else
            {
                var allPsmsSingle = allPsms.Where(p => p.OGlycanBoxLocalization == null && p.Score > 2).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                var writtenFileInter1 = Path.Combine(OutputFolder, "single_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileInter1, 1);

                var allPsmsGly = allPsms.Where(p => p.OGlycanBoxLocalization != null && p.Score > 2).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });

                var writtenFileInter2 = Path.Combine(OutputFolder, "glyco_fdr" + ".tsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileInter2, 2);

                var ProteinLevelLocalization = ProteinLevelGlycoParsimony(allPsmsGly);
                WriteFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, Path.Combine(OutputFolder, "seen_glyco_localization" + ".tsv"));
                WriteFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, Path.Combine(OutputFolder, "protein_glyco_localization" + ".tsv"));
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

        private Dictionary<string, Tuple<bool, double, double>> ProteinLevelGlycoParsimony(List<GlycoSpectralMatch> allPsmsGly)
        {
            //<id, <islocalized, minQValue, maxProb>>. id: ProteinAccession, ProtienPos, GlycanId. islocalized, minQValue, maxProb
            Dictionary<string, Tuple<bool, double, double>> localizedGlycan = new Dictionary<string, Tuple<bool, double, double>>();

            foreach (var gsm in allPsmsGly)
            {
                if (gsm.IsContaminant || gsm.IsDecoy)
                {
                    continue;
                }

                if (gsm.LocalizedGlycan.Count > 0)
                {
                    foreach (var local in gsm.LocalizedGlycan)
                    {
                        int proteinPos = local.Item1 + gsm.OneBasedStartResidueInProtein.Value;

                        string proPosId = gsm.ProteinAccession + "-" + proteinPos.ToString() + "-" + local.Item2;

                        double prob = -1;
                        if (gsm.SiteSpeciLocalProb!=null && gsm.SiteSpeciLocalProb.ContainsKey(local.Item1))
                        {
                            prob = gsm.SiteSpeciLocalProb[local.Item1].Where(p => p.Item1 == local.Item2).FirstOrDefault().Item2;
                        }


                        if (!localizedGlycan.ContainsKey(proPosId))
                        {
                            localizedGlycan.Add(proPosId, new Tuple<bool, double, double>(local.Item3, gsm.FdrInfo.QValue, prob));            
                        }
                        else
                        {
                            bool islocalized = (local.Item3 || localizedGlycan[proPosId].Item1);
                            double minQValue = localizedGlycan[proPosId].Item2 > gsm.FdrInfo.QValue ? gsm.FdrInfo.QValue : localizedGlycan[proPosId].Item2;
                            double maxProb = localizedGlycan[proPosId].Item3 > prob ? localizedGlycan[proPosId].Item3 : prob;

                            localizedGlycan[proPosId] = new Tuple<bool, double, double>(islocalized, minQValue, maxProb);                       
                        }
                    }
                }    
            }

            return localizedGlycan;
        }

    }
}

