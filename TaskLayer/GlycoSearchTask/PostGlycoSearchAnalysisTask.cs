using EngineLayer;
using EngineLayer.GlycoSearch;
using EngineLayer.FdrAnalysis;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using FlashLFQ;
using Proteomics.ProteolyticDigestion;
using Easy.Common.Extensions;

namespace TaskLayer
{
    public class PostGlycoSearchAnalysisTask : MetaMorpheusTask
    {
        public PostGlycoSearchAnalysisParameters Parameters { get; set; }
        private List<EngineLayer.ProteinGroup> ProteinGroups { get; set; }
        private IEnumerable<IGrouping<string, PeptideSpectralMatch>> PsmsGroupedByFile { get; set; }
        public PostGlycoSearchAnalysisTask() : base(MyTask.GlycoSearch)
        {
        }

        public MyTaskResults Run(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList, List<GlycoSpectralMatch> allPsms, CommonParameters commonParameters, GlycoSearchParameters glycoSearchParameters, List<Protein> proteinList, List<Modification> variableModifications, List<Modification> fixedModifications, List<string> localizeableModificationTypes, MyTaskResults MyTaskResults)
        {
            // Stop loop if canceled
            if (GlobalVariables.StopLoops) { return Parameters.GlycoSearchTaskResults; }

            if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.NGlycanSearch)
            {
                var allPsmsSingle = allPsms.Where(p => p.NGlycan == null).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                if (allPsmsSingle.Any())
                {
                    QuantificationAnalysis(allPsmsSingle);
                }
                WriteQuantificationResults("allSinglePsms");
                var writtenFileSingle = Path.Combine(OutputFolder, "single" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileSingle, 1);
                FinishedWritingFile(writtenFileSingle, new List<string> { taskId });

                var allPsmsGly = allPsms.Where(p => p.NGlycan != null ).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                if (allPsmsGly.Any())
                {
                    QuantificationAnalysis(allPsmsGly);
                }
                WriteQuantificationResults("allNGlycanPsms");
                var writtenFileNGlyco = Path.Combine(OutputFolder, "nglyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileNGlyco, 3);

                FinishedWritingFile(writtenFileNGlyco, new List<string> { taskId });
                return MyTaskResults;
            }
            else if (glycoSearchParameters.GlycoSearchType == GlycoSearchType.OGlycanSearch)
            {
                var allPsmsSingle = allPsms.Where(p => p.Routes == null ).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                if (allPsmsSingle.Any())
                {
                    QuantificationAnalysis(allPsmsSingle);
                }
                WriteQuantificationResults("allSinglePsms");
                var writtenFileSingle = Path.Combine(OutputFolder, "single" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileSingle, 1);
                FinishedWritingFile(writtenFileSingle, new List<string> { taskId });

                var allPsmsGly = allPsms.Where(p => p.Routes != null ).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                if (allPsmsGly.Any())
                {
                    QuantificationAnalysis(allPsmsGly);
                }
                WriteQuantificationResults("allOGlycoPsms");
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
                var allPsmsSingle = allPsms.Where(p => p.Routes == null).OrderByDescending(p => p.Score).ToList();
                SingleFDRAnalysis(allPsmsSingle, commonParameters, new List<string> { taskId });
                if (allPsmsSingle.Any())
                {
                    QuantificationAnalysis(allPsmsSingle);
                }
                WriteQuantificationResults("allSinglePsms");
                var writtenFileSingle = Path.Combine(OutputFolder, "single" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsSingle, writtenFileSingle, 1);
                FinishedWritingFile(writtenFileSingle, new List<string> { taskId });

                var allPsmsGly = allPsms.Where(p => p.NGlycan != null).ToList();
                SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                if (allPsmsGly.Any())
                {
                    QuantificationAnalysis(allPsmsGly);
                }
                WriteQuantificationResults("allPsmsGly");
                var writtenFileNGlyco = Path.Combine(OutputFolder, "n_o_glyco" + ".psmtsv");
                WriteFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileNGlyco, 3);
                FinishedWritingFile(writtenFileNGlyco, new List<string> { taskId });

                return MyTaskResults;
            }
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            return null;
        }
        //Calculate the FDR of single peptide FP/TP
        private void SingleFDRAnalysis(List<GlycoSpectralMatch> items, CommonParameters commonParameters, List<string> taskIds)
        {
            // calculate single PSM FDR
            List<PeptideSpectralMatch> psms = items.Select(p => p as PeptideSpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, this.FileSpecificParameters, taskIds).Run();

        }

        private void QuantificationAnalysis(List<GlycoSpectralMatch> gsms)
        {
            List<string> extensions = Parameters.CurrentRawFileList.Select(f => Path.GetExtension(f)).ToList();
            
            if (!Parameters.GlycoSearchParameters.DoQuantification || Parameters.CurrentRawFileList.Select(f => Path.GetExtension(f)).ToList().Any(f=>f == ".mgf"))//if there are mgf files then dont do quant
            {
                return;
            }

            // pass quantification parameters to FlashLFQ
            Status("Quantifying...", Parameters.SearchTaskId);

            foreach (var file in Parameters.CurrentRawFileList)
            {
                Parameters.MyFileManager.DoneWithFile(file);
            }

            // construct file info for FlashLFQ
            List<SpectraFileInfo> spectraFileInfo;

            // get experimental design info
            string pathToFirstSpectraFile = Directory.GetParent(Parameters.CurrentRawFileList.First()).FullName;
            string assumedExperimentalDesignPath = Path.Combine(pathToFirstSpectraFile, GlobalVariables.ExperimentalDesignFileName);

            if (File.Exists(assumedExperimentalDesignPath))
            {
                // copy experimental design file to output folder
                string writtenFile = Path.Combine(Parameters.OutputFolder, Path.GetFileName(assumedExperimentalDesignPath));
                try
                {
                    File.Copy(assumedExperimentalDesignPath, writtenFile, overwrite: true);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });
                }
                catch
                {
                    Warn("Could not copy Experimental Design file to search task output. That's ok, the search will continue");
                }

                spectraFileInfo = ExperimentalDesign.ReadExperimentalDesign(assumedExperimentalDesignPath, Parameters.CurrentRawFileList, out var errors);

                if (errors.Any())
                {
                    Warn("Error reading experimental design file: " + errors.First() + ". Skipping quantification");
                    return;
                }
            }
            else if (Parameters.GlycoSearchParameters.Normalize)
            {
                Warn("Could not find experimental design file at " + assumedExperimentalDesignPath + ", which is required for normalization. Skipping quantification");
                return;
            }
            else
            {
                spectraFileInfo = new List<SpectraFileInfo>();

                for (int i = 0; i < Parameters.CurrentRawFileList.Count; i++)
                {
                    var file = Parameters.CurrentRawFileList[i];

                    // experimental design info passed in here for each spectra file
                    spectraFileInfo.Add(new SpectraFileInfo(fullFilePathWithExtension: file, condition: "", biorep: i, fraction: 0, techrep: 0));
                }
            }

            // get PSMs to pass to FlashLFQ
            List<PeptideSpectralMatch> unambiguousPsmsBelowOnePercentFdr = new();
            if (gsms.Count > 100)//PEP is not computed when there are fewer than 100 psms
            {
                //someday when we compute PEP for glyco we can amend this
                unambiguousPsmsBelowOnePercentFdr = gsms.Where(p =>
                    p.FdrInfo.PEP_QValue <= 0.01
                    && !p.IsDecoy
                    && p.FullSequence != null).Select(p=>(PeptideSpectralMatch)p).ToList(); //if ambiguous, there's no full sequence
            }
            else
            {
                //someday when we compute PEP for glyco we can amend this
                unambiguousPsmsBelowOnePercentFdr = gsms.Where(p =>
                    p.FdrInfo.QValue <= 0.01
                    && !p.IsDecoy
                    && p.FullSequence != null).Select(p => (PeptideSpectralMatch)p).ToList(); //if ambiguous, there's no full sequence
            }

            // pass protein group info for each PSM
            var psmToProteinGroups = new Dictionary<PeptideSpectralMatch, List<FlashLFQ.ProteinGroup>>();
            if (ProteinGroups != null && ProteinGroups.Count != 0) //ProteinGroups can be null if parsimony wasn't done, and it can be empty if you're doing the two peptide rule
            {
                foreach (var proteinGroup in ProteinGroups)
                {
                    var proteinsOrderedByAccession = proteinGroup.Proteins.OrderBy(p => p.Accession);

                    var flashLfqProteinGroup = new FlashLFQ.ProteinGroup(proteinGroup.ProteinGroupName,
                        string.Join("|", proteinsOrderedByAccession.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault())),
                        string.Join("|", proteinsOrderedByAccession.Select(p => p.Organism).Distinct()));

                    foreach (var psm in proteinGroup.AllPsmsBelowOnePercentFDR.Where(v => v.FullSequence != null))
                    {
                        if (psmToProteinGroups.TryGetValue(psm, out var flashLfqProteinGroups))
                        {
                            flashLfqProteinGroups.Add(flashLfqProteinGroup);
                        }
                        else
                        {
                            psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { flashLfqProteinGroup });
                        }
                    }
                }
            }
            else
            {
                // if protein groups were not constructed, just use accession numbers
                var accessionToPg = new Dictionary<string, FlashLFQ.ProteinGroup>();
                foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
                {
                    var proteins = psm.BestMatchingPeptides.Select(b => b.Peptide.Protein).Distinct();

                    foreach (var protein in proteins)
                    {
                        if (!accessionToPg.ContainsKey(protein.Accession))
                        {
                            accessionToPg.Add(protein.Accession, new FlashLFQ.ProteinGroup(protein.Accession, string.Join("|", protein.GeneNames.Select(p => p.Item2).Distinct()), protein.Organism));
                        }

                        if (psmToProteinGroups.TryGetValue(psm, out var proteinGroups))
                        {
                            proteinGroups.Add(accessionToPg[protein.Accession]);
                        }
                        else
                        {
                            psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { accessionToPg[protein.Accession] });
                        }
                    }
                }
            }

            //group psms by file
            var psmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy(p => p.FullFilePath);

            // some PSMs may not have protein groups (if 2 peptides are required to construct a protein group, some PSMs will be left over)
            // the peptides should still be quantified but not considered for protein quantification
            var undefinedPg = new FlashLFQ.ProteinGroup("UNDEFINED", "", "");
            //sort the unambiguous psms by protease to make MBR compatible with multiple proteases
            Dictionary<Protease, List<PeptideSpectralMatch>> proteaseSortedPsms = new Dictionary<Protease, List<PeptideSpectralMatch>>();
            Dictionary<Protease, FlashLfqResults> proteaseSortedFlashLFQResults = new Dictionary<Protease, FlashLfqResults>();

            foreach (DigestionParams dp in Parameters.ListOfDigestionParams)
            {
                if (!proteaseSortedPsms.ContainsKey(dp.Protease))
                {
                    proteaseSortedPsms.Add(dp.Protease, new List<PeptideSpectralMatch>());
                }
            }
            foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
            {
                if (!psmToProteinGroups.ContainsKey(psm))
                {
                    psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { undefinedPg });
                }

                proteaseSortedPsms[psm.DigestionParams.Protease].Add(psm);
            }

            // pass PSM info to FlashLFQ
            var flashLFQIdentifications = new List<Identification>();
            foreach (var spectraFile in psmsGroupedByFile)
            {
                var rawfileinfo = spectraFileInfo.Where(p => p.FullFilePathWithExtension.Equals(spectraFile.Key)).First();

                foreach (var psm in spectraFile)
                {
                    flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence,
                        psm.PeptideMonisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroups[psm]));
                }
            }

            // run FlashLFQ
            var FlashLfqEngine = new FlashLfqEngine(
                allIdentifications: flashLFQIdentifications,
                normalize: Parameters.GlycoSearchParameters.Normalize,
                ppmTolerance: Parameters.GlycoSearchParameters.QuantifyPpmTol,
                matchBetweenRunsPpmTolerance: Parameters.GlycoSearchParameters.QuantifyPpmTol,  // If these tolerances are not equivalent, then MBR will falsely classify peptides found in the initial search as MBR peaks
                matchBetweenRuns: false, //we'll probably add this back later
                silent: true,
                maxThreads: CommonParameters.MaxThreadsToUsePerFile);

            if (flashLFQIdentifications.Any())
            {
                Parameters.FlashLfqResults = FlashLfqEngine.Run();
            }

            // get protein intensity back from FlashLFQ
            if (ProteinGroups != null && Parameters.FlashLfqResults != null)
            {
                foreach (var proteinGroup in ProteinGroups)
                {
                    proteinGroup.FilesForQuantification = spectraFileInfo;
                    proteinGroup.IntensitiesByFile = new Dictionary<SpectraFileInfo, double>();

                    foreach (var spectraFile in proteinGroup.FilesForQuantification)
                    {
                        if (Parameters.FlashLfqResults.ProteinGroups.TryGetValue(proteinGroup.ProteinGroupName, out var flashLfqProteinGroup))
                        {
                            proteinGroup.IntensitiesByFile.Add(spectraFile, flashLfqProteinGroup.GetIntensity(spectraFile));
                        }
                        else
                        {
                            proteinGroup.IntensitiesByFile.Add(spectraFile, 0);
                        }
                    }
                }
            }
        }

        private void WriteQuantificationResults(string apppendTextToFilename)
        {
            if (Parameters.GlycoSearchParameters.DoQuantification && Parameters.FlashLfqResults != null)
            {
                // write peaks
                WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "AllQuantifiedPeaks" + "_" + apppendTextToFilename, new List<string> { Parameters.SearchTaskId });

                // write peptide quant results
                string filename = "AllQuantified" + GlobalVariables.AnalyteType + "s";
                WritePeptideQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, filename + "_" + apppendTextToFilename, new List<string> { Parameters.SearchTaskId });
                
                // write individual results
                if (Parameters.CurrentRawFileList.Count > 1 && Parameters.GlycoSearchParameters.WriteIndividualFiles)
                {
                    foreach (var file in Parameters.FlashLfqResults.Peaks)
                    {
                        WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.IndividualResultsOutputFolder,
                            file.Key.FilenameWithoutExtension + "_QuantifiedPeaks", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.Key.FullFilePathWithExtension });
                    }
                }
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

