using EngineLayer;
using EngineLayer.GlycoSearch;
using EngineLayer.FdrAnalysis;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using FlashLFQ;
using Omics.Modifications;
using Omics.SpectrumMatch;
using pepXML.Generated;
using Proteomics.ProteolyticDigestion;
using Readers;
using ThermoFisher.CommonCore.Data;
using ProteinGroup = EngineLayer.ProteinGroup;

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
            MyTaskResults = new MyTaskResults(this);
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

            //This is all psms for all files including glyco- and non-glyco psms.
            SingleFDRAnalysis(allPSMs, commonParameters, new List<string> { taskId });

            List<GlycoSpectralMatch> filteredPsms = allPSMs.Where(p => p.FdrInfo.QValue <= 0.01).ToList();

            //write individual file results
            if (Parameters.GlycoSearchParameters.WriteIndividualFiles)
            {
                string individualFileResults = Path.Combine(OutputFolder, "IndividualFileResults");
                if (!Directory.Exists(individualFileResults))
                {
                    Directory.CreateDirectory(individualFileResults);
                }

                

                foreach (var fileSpecificPSMs in filteredPsms.GroupBy(p => p.FullFilePath)) //group by file path, and the path will be the key for the dictionary
                {
                    string individualFileFolder = Path.GetFileNameWithoutExtension(fileSpecificPSMs.Key); //folder name.
                    string individualFileFolderPath = Path.Combine(individualFileResults, individualFileFolder);
                    if (!Directory.Exists(individualFileFolderPath))
                    {
                        Directory.CreateDirectory(individualFileFolderPath);
                    }
                    var fspList = fileSpecificPSMs.ToList();
                    if (Parameters.GlycoSearchParameters.DoParsimony)
                    {
                        GlycoProteinAnalysis(fspList, individualFileFolderPath, individualFileFolder); //Creat the proteinGroups file 
                    }
                    
                    foreach (GlycoSpectralMatch gsm in fspList) //maybe this needs to be the filterd list???
                    {
                        gsm.ResolveAllAmbiguities();
                    }
                    var individualFilePsmsPath = Path.Combine(individualFileFolderPath, individualFileFolder + "_AllPSMs.psmtsv");
                    WriteGlycoFile.WritePsmGlycoToTsv(fspList, individualFilePsmsPath, false);//this is everything, glyco and non-glyco
                    //the individual file AllPSMs was just written. The next method writes only those PSMs that have a glyco mod
                    DivideGlycoPsmsIntoGroupsWriteToTsv(glycoSearchParameters.GlycoSearchType, fspList, commonParameters, taskId, individualFileFolderPath, individualFileFolder);
                }
            }

            //write combined results
            switch (glycoSearchParameters.GlycoSearchType)
            {
                case GlycoSearchType.OGlycanSearch:
                    var OglyInAllPsms = filteredPsms.Where(p => p.Routes != null).ToList(); //Try to filter out the non-glyco psms
                    if (OglyInAllPsms.Any()) // Is there any gsms in the allPsms?
                    {
                        SingleFDRAnalysis(OglyInAllPsms, commonParameters, new List<string> { taskId });
                        var writtenFileOGlyco = Path.Combine(OutputFolder + "\\oglyco" + ".psmtsv");

                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(OglyInAllPsms.Where(p => p.Accession != null && p.OneBasedStartResidue.HasValue).ToList());
                        var seen_oglyco_localization_file = Path.Combine(OutputFolder + "\\seen_oglyco_localization" + ".tsv"); //generate the localization file
                        WriteGlycoFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_oglyco_localization_file);

                        var protein_oglyco_localization_file = Path.Combine(OutputFolder + "\\protein_oglyco_localization" + ".tsv");
                        WriteGlycoFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_oglyco_localization_file);
                        // Writing the oglyco results to a file and summary text
                        WriteGlycoFile.WritePsmGlycoToTsv(OglyInAllPsms, writtenFileOGlyco, true); //we write this last so localization can be attempted 
                        

                    }
                    break;
                case GlycoSearchType.NGlycanSearch:
                    var NglyInAllPsms = filteredPsms.Where(p => p.GlycanScore > 0 && p.Routes == null).ToList();
                    if (NglyInAllPsms.Any())
                    {
                        SingleFDRAnalysis(NglyInAllPsms, commonParameters, new List<string> { taskId });

                        var writtenFileNGlyco = Path.Combine(OutputFolder + "\\nglyco" + ".psmtsv");

                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(NglyInAllPsms.Where(p => p.Accession != null && p.OneBasedStartResidue.HasValue).ToList());
                        var seen_nglyco_localization_file = Path.Combine(OutputFolder + "\\seen_nglyco_localization" + ".tsv");
                        WriteGlycoFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_nglyco_localization_file);

                        var protein_nglyco_localization_file = Path.Combine(OutputFolder + "\\protein_nglyco_localization" + ".tsv");
                        WriteGlycoFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_nglyco_localization_file);
                        WriteGlycoFile.WritePsmGlycoToTsv(NglyInAllPsms, writtenFileNGlyco, true); //we write this last so localization can be attempted                        
                    }
                    break;
                case GlycoSearchType.N_O_GlycanSearch:
                default:
                    var allPsmsgly = filteredPsms.Where(p => p.GlycanScore > 0).ToList();
                    if (allPsmsgly.Any())
                    {
                        SingleFDRAnalysis(allPsmsgly, commonParameters, new List<string> { taskId });

                        var writtenFileNOGlyco = Path.Combine(OutputFolder + "\\no_glyco" + ".psmtsv");

                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsgly.Where(p => p.Accession != null && p.OneBasedStartResidue.HasValue).ToList());
                        var seen_no_glyco_localization_file = Path.Combine(OutputFolder + "\\seen_no_glyco_localization" + ".tsv");
                        WriteGlycoFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_no_glyco_localization_file);

                        var protein_no_glyco_localization_file = Path.Combine(OutputFolder + "\\protein_no_glyco_localization" + ".tsv");
                        WriteGlycoFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_no_glyco_localization_file);
                        WriteGlycoFile.WritePsmGlycoToTsv(allPsmsgly, writtenFileNOGlyco, true); //we write this last so localization can be attempted
                    
                    }
                    break;
            }

            if (glycoSearchParameters.DoParsimony)
            {
                GlycoProteinAnalysis(filteredPsms, OutputFolder, null, MyTaskResults);//Do the whole group last so inference is done on the whole group
            }
            else
            {
                GlycoAccessionAnalysis(filteredPsms, OutputFolder);//Do the whole group last so inference is done on the whole group
            }
            
            QuantificationAnalysis();
            WriteQuantificationResults();

            var writtenFileSingle = Path.Combine(OutputFolder, "AllPSMs.psmtsv");
            WriteGlycoFile.WritePsmGlycoToTsv(filteredPsms, writtenFileSingle, true);
            
     

            if (Parameters.GlycoSearchParameters.WriteSpectrumLibrary)
            {
                List<LibrarySpectrum> spectrumLibrary = new List<LibrarySpectrum>();
                foreach (var gsm in filteredPsms)
                {
                    spectrumLibrary.Add(new LibrarySpectrum(gsm.FullSequence, gsm.ScanPrecursorMonoisotopicPeakMz, gsm.ScanPrecursorCharge, gsm.MatchedFragmentIons,gsm.ScanRetentionTime,gsm.IsDecoy));
                }
                WriteSpectrumLibrary(spectrumLibrary, OutputFolder);
            }

            FinishedWritingFile(writtenFileSingle, new List<string> { taskId });
            
            WriteSummary(filteredPsms, glycoSearchParameters, MyTaskResults);

            return MyTaskResults;
        }

        /// <summary>
        /// Wirte the summary of the glyco search results to the results txt file
        /// </summary>
        /// <param name="targetPsms"></param>
        /// <param name="proteinGroups"></param>
        /// <param name="glycoSearchParameters"></param>
        /// <param name="MyTaskResults"></param>
        private void WriteSummary(List<GlycoSpectralMatch> targetPsms, GlycoSearchParameters glycoSearchParameters, MyTaskResults MyTaskResults)
        {
            var gsms = targetPsms.Where(p => p.Routes != null).ToList();
            var Level1gsms = gsms.Where(p => p.LocalizationLevel == LocalizationLevel.Level1).ToList();
            MyTaskResults.AddTaskSummaryText("All target PSMs within 1% FDR: " + (targetPsms?.
                           Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.IsContaminant) ?? 0));
            MyTaskResults.AddTaskSummaryText("All target protein groups within 1% FDR: " + (ProteinGroups?.
                            Count(p => p.QValue <= 0.01 && !p.IsDecoy && !p.IsContaminant) ?? 0));

            switch (glycoSearchParameters.GlycoSearchType)
            {
                case GlycoSearchType.OGlycanSearch:
                    MyTaskResults.AddTaskSummaryText("All target O-Glyco PSMs within 1% FDR: " + (gsms?.
                            Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.IsContaminant) ?? 0));
                    MyTaskResults.AddTaskSummaryText("All target Level 1 O-Glyco PSMs within 1% FDR: " + (Level1gsms
                        ?.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.IsContaminant && p.LocalizationLevel == LocalizationLevel.Level1) ?? 0));
                    break;
                case GlycoSearchType.NGlycanSearch:
                    MyTaskResults.AddTaskSummaryText("All target N-Glyco PSMs within 1% FDR: " + (gsms?.
                            Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.IsContaminant) ?? 0));
                    MyTaskResults.AddTaskSummaryText("All target Level 1 N-Glyco PSMs within 1% FDR: " + (Level1gsms
                        ?.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.IsContaminant && p.LocalizationLevel == LocalizationLevel.Level1) ?? 0));
                    break;
                case GlycoSearchType.N_O_GlycanSearch:
                    MyTaskResults.AddTaskSummaryText("All target Glyco PSMs within 1% FDR: " + (gsms?.
                            Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.IsContaminant) ?? 0));
                    MyTaskResults.AddTaskSummaryText("All target Level 1 Glyco PSMs within 1% FDR: " + (Level1gsms
                        ?.Count(p => p.FdrInfo.QValue <= 0.01 && !p.IsDecoy && !p.IsContaminant && p.LocalizationLevel == LocalizationLevel.Level1) ?? 0));
                    break;
            }

        }



        private void DivideGlycoPsmsIntoGroupsWriteToTsv(GlycoSearchType glycoSearchType, List<GlycoSpectralMatch> gsms, CommonParameters commonParameters, string taskId, string individualFileFolderPath, string individualFileFolder)
        {
            switch (glycoSearchType)
            {
                case GlycoSearchType.OGlycanSearch:
                    var allPsmsOgly = gsms.Where(p => p.Routes != null).ToList(); 
                    if (allPsmsOgly.Any()) //In the all gsms. is there any gsms contain localization informaiton(route)
                    {
                        SingleFDRAnalysis(allPsmsOgly, commonParameters, new List<string> { taskId });
                        var writtenFileOGlyco = Path.Combine(individualFileFolderPath, individualFileFolder + "oglyco" + ".psmtsv");
                        
                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsOgly.Where(p => p.Accession != null && p.OneBasedStartResidue.HasValue).ToList());
                        var seen_oglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "seen_oglyco_localization" + ".tsv");
                        WriteGlycoFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_oglyco_localization_file);

                        var protein_oglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "protein_oglyco_localization" + ".tsv");
                        WriteGlycoFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_oglyco_localization_file);
                        WriteGlycoFile.WritePsmGlycoToTsv(allPsmsOgly, writtenFileOGlyco, true); //we write this last so localization can be attempted
                    }
                    break;
                case GlycoSearchType.NGlycanSearch:
                    var allPsmsNgly = gsms.Where(p => p.GlycanScore > 0 && p.Routes == null).OrderByDescending(p => p.Score).ToList();
                    if (allPsmsNgly.Any())
                    {
                        SingleFDRAnalysis(allPsmsNgly, commonParameters, new List<string> { taskId });
                        var writtenFileNGlyco = Path.Combine(individualFileFolderPath, individualFileFolder + "nglyco" + ".psmtsv");

                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsNgly.Where(p => p.Accession != null && p.OneBasedStartResidue.HasValue).ToList());
                        var seen_nglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "seen_nglyco_localization" + ".tsv");
                        WriteGlycoFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_nglyco_localization_file);

                        var protein_nglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "protein_nglyco_localization" + ".tsv");
                        WriteGlycoFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_nglyco_localization_file);
                        WriteGlycoFile.WritePsmGlycoToTsv(allPsmsNgly, writtenFileNGlyco, true); //we write this last so localization can be attempted
                    }
                    break;
                default: //both N and O glycan search combined
                    var allPsmsGly = gsms.Where(p => p.GlycanScore > 0).OrderByDescending(p => p.Score).ToList();
                    if (allPsmsGly.Any())
                    {
                        SingleFDRAnalysis(allPsmsGly, commonParameters, new List<string> { taskId });
                        var writtenFileNOGlyco = Path.Combine(individualFileFolderPath, individualFileFolder + "no_glyco" + ".psmtsv");
                        WriteGlycoFile.WritePsmGlycoToTsv(allPsmsGly, writtenFileNOGlyco, true);

                        var ProteinLevelLocalization = GlycoProteinParsimony.ProteinLevelGlycoParsimony(allPsmsGly.Where(p => p.Accession != null && p.OneBasedStartResidue.HasValue).ToList());
                        var seen_no_glyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "seen_no_glyco_localization" + ".tsv");
                        WriteGlycoFile.WriteSeenProteinGlycoLocalization(ProteinLevelLocalization, seen_no_glyco_localization_file);

                        var protein_noglyco_localization_file = Path.Combine(individualFileFolderPath, individualFileFolder + "protein_no_glyco_localization" + ".tsv");
                        WriteGlycoFile.WriteProteinGlycoLocalization(ProteinLevelLocalization, protein_noglyco_localization_file);
                    }
                    break;
            }
        }

        //Calculate the FDR of single peptide FP/TP
        private void SingleFDRAnalysis(List<GlycoSpectralMatch> items, CommonParameters commonParameters, List<string> taskIds)
        {
            // calculate single PSM FDR
            List<SpectralMatch> psms = items.Select(p => p as SpectralMatch).ToList();
            new FdrAnalysisEngine(psms, 0, commonParameters, this.FileSpecificParameters, taskIds).Run();

        }
        private void GlycoProteinAnalysis(List<GlycoSpectralMatch> gsms, string outputFolder, string individualFileFolder = null, MyTaskResults myTaskResults = null )
        {
            // convert gsms to psms
            List<SpectralMatch> psmsForProteinParsimony = gsms.Select(p => p as SpectralMatch).ToList();

            Status("Constructing protein groups...", Parameters.SearchTaskId);

            // run parsimony
            ProteinParsimonyResults proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(psmsForProteinParsimony, Parameters.GlycoSearchParameters.ModPeptidesAreDifferent, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run());

            // score protein groups and calculate FDR
            ProteinScoringAndFdrResults proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psmsForProteinParsimony,
                Parameters.GlycoSearchParameters.NoOneHitWonders, Parameters.GlycoSearchParameters.ModPeptidesAreDifferent, true, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run();

            ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

            Status("Done constructing protein groups!", Parameters.SearchTaskId);
            WriteProteinResults(outputFolder, individualFileFolder, myTaskResults);
         
        }
        private void GlycoAccessionAnalysis(List<GlycoSpectralMatch> gsms, string individualFileFolderPath, string individualFileFolder = null)
        {
            List<ProteinGroup> _proteinGroups = new List<ProteinGroup>();

            // convert gsms to psms
            List<SpectralMatch> psmsForProteinParsimony = gsms.Select(p => p as SpectralMatch).ToList();

            foreach (var psm in psmsForProteinParsimony)
            {
                psm.ResolveAllAmbiguities();
            }

            //get non-ambigous psms above cutoffs that are normally used for parsimony
            List<SpectralMatch>  _filteredPsms = psmsForProteinParsimony.Where(p => p.FullSequence != null && p.FdrInfo.QValue <= 0.01 && p.FdrInfo.QValueNotch <= 0.01).ToList();

            // if there are no peptides observed, there are no proteins; return an empty list of protein groups
            if (_filteredPsms.Count != 0)
            {
                foreach (var psm in _filteredPsms)
                {
                    List<Protein> proteinList = psm.BestMatchingBioPolymersWithSetMods.Select(p => ((PeptideWithSetModifications)p.SpecificBioPolymer).Protein).ToList();
                    ProteinGroup newProteinGroup = new ProteinGroup(new HashSet<Protein>(proteinList),
                        new HashSet<PeptideWithSetModifications>(new List<PeptideWithSetModifications>(psm.BestMatchingBioPolymersWithSetMods.Select(p=> (PeptideWithSetModifications)p.SpecificBioPolymer).ToList())), new HashSet<PeptideWithSetModifications>());

                    if (_proteinGroups.Any(p => p.Equals(newProteinGroup)))
                    {
                        _proteinGroups.Where(g => g.Equals(newProteinGroup)).First().MergeProteinGroupWith(newProteinGroup);
                    }
                    else
                    {
                        _proteinGroups.Add(newProteinGroup);
                    }   
                }
            }
        }
        private void WriteProteinResults(string outputFolder, string individualFileFolder = null, MyTaskResults myTaskResults = null)
        {
            double qValueCutoff_FORDEBUGONLY = 0.01;
            
            string fileName = "AllProteinGroups.tsv";
            string writtenFile = Path.Combine(outputFolder, individualFileFolder + "_"+ fileName);
            WriteProteinGroupsToTsv(ProteinGroups, writtenFile, new List<string> { Parameters.SearchTaskId }, qValueCutoff_FORDEBUGONLY);
           

        }
        private void WriteProteinGroupsToTsv(List<EngineLayer.ProteinGroup> proteinGroups, string filePath, List<string> nestedIds, double qValueCutoff)
        {
            if (proteinGroups != null && proteinGroups.Any())
            {
                using (StreamWriter output = new(filePath))
                {
                    output.WriteLine(proteinGroups.First().GetTabSeparatedHeader());
                    for (int i = 0; i < proteinGroups.Count; i++)
                    {
                        if ((Parameters.GlycoSearchParameters.WriteDecoys && proteinGroups[i].IsDecoy) // write decoys if they exist and are desired
                            || (Parameters.GlycoSearchParameters.WriteContaminants && proteinGroups[i].IsContaminant) // write contaminants if they exist and are desired
                            || (!proteinGroups[i].IsDecoy && !proteinGroups[i].IsContaminant)) // write targets regardless
                        {
                            output.WriteLine(proteinGroups[i]);
                        }
                    }
                }
            }
        }
        private void QuantificationAnalysis()
        {
            if (!Parameters.GlycoSearchParameters.DoQuantification)
            {
                return;
            }

            // pass quantification parameters to FlashLFQ
            Status("Quantifying...", Parameters.SearchTaskId);

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
            var unambiguousPsmsBelowOnePercentFdr = Parameters.AllPsms.Where(p =>
                p.FdrInfo.QValue <= 0.01
                && p.FdrInfo.QValueNotch <= 0.01
                && !p.IsDecoy
                && p.FullSequence != null).ToList();

            // pass protein group info for each PSM
            var psmToProteinGroups = new Dictionary<SpectralMatch, List<FlashLFQ.ProteinGroup>>();
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
            else //we will use protein accession numbers instead of protein group names
            {
                var accessionToPg = new Dictionary<string, FlashLFQ.ProteinGroup>();
                foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
                {
                    var proteins = psm.BestMatchingBioPolymersWithSetMods.Select(b => ((PeptideWithSetModifications)b.SpecificBioPolymer).Protein).Distinct();

                    foreach (var protein in proteins)
                    {
                        if (!accessionToPg.ContainsKey(protein.Accession))
                        {
                            accessionToPg.Add(protein.Accession, new FlashLFQ.ProteinGroup(protein.Accession, string.Join("|", protein.GeneNames.Select(p => p.Item2).Distinct()), protein.Organism));
                        }
                        
                    }
                    psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { accessionToPg[psm.Accession] });
                }
            }

            //group psms by file
            var psmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy(p => p.FullFilePath);

            // some PSMs may not have protein groups (if 2 peptides are required to construct a protein group, some PSMs will be left over)
            // the peptides should still be quantified but not considered for protein quantification
            var undefinedPg = new FlashLFQ.ProteinGroup("UNDEFINED", "", "");
            //sort the unambiguous psms by protease to make MBR compatible with multiple proteases
            Dictionary<Protease, List<SpectralMatch>> proteaseSortedPsms = new Dictionary<Protease, List<SpectralMatch>>();
            Dictionary<Protease, FlashLfqResults> proteaseSortedFlashLFQResults = new Dictionary<Protease, FlashLfqResults>();

            foreach (DigestionParams dp in Parameters.ListOfDigestionParams)
            {
                if (!proteaseSortedPsms.ContainsKey(dp.Protease))
                {
                    proteaseSortedPsms.Add(dp.Protease, new List<SpectralMatch>());
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
                        psm.BioPolymerWithSetModsMonoisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroups[psm]));
                }
            }

            // run FlashLFQ
            var FlashLfqEngine = new FlashLfqEngine(
                allIdentifications: flashLFQIdentifications,
                normalize: Parameters.GlycoSearchParameters.Normalize,
                ppmTolerance: Parameters.GlycoSearchParameters.QuantifyPpmTol,
                matchBetweenRunsPpmTolerance: Parameters.GlycoSearchParameters.QuantifyPpmTol,  // If these tolerances are not equivalent, then MBR will falsely classify peptides found in the initial search as MBR peaks
                matchBetweenRuns: Parameters.GlycoSearchParameters.DoMbrAnalysis,
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

        private void WriteQuantificationResults()
        {
            if (Parameters.GlycoSearchParameters.DoQuantification && Parameters.FlashLfqResults != null)
            {
                // write peaks
                // once SpectralRecoveryResults is implemented, this will need to be updated with an if statement

                WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "AllQuantifiedPeaks", new List<string> { Parameters.SearchTaskId });

                // write peptide quant results
                string filename = "AllQuantified" + GlobalVariables.AnalyteType + "s";

                // write peptides
                // once SpectralRecoveryResults is implemented, this will need to be updated with an if statement

                WritePeptideQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, filename, new List<string> { Parameters.SearchTaskId });

                WriteProteinQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder,
                    "AllQuantifiedProteins", new List<string> { Parameters.SearchTaskId });

                // write individual results
                if (Parameters.CurrentRawFileList.Count > 0 && Parameters.GlycoSearchParameters.WriteIndividualFiles)
                {
                    foreach (var file in Parameters.FlashLfqResults.Peaks)
                    {
                        WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Path.Combine(Parameters.IndividualResultsOutputFolder, file.Key.FilenameWithoutExtension),
                            file.Key.FilenameWithoutExtension + "_QuantifiedPeaks", new List<string> { Parameters.SearchTaskId, "IndividualFileResults", file.Key.FullFilePathWithExtension });
                        WritePeptideQuantificationResultsToTsv(Parameters.FlashLfqResults, Path.Combine(Parameters.IndividualResultsOutputFolder, file.Key.FilenameWithoutExtension),
                            file.Key.FilenameWithoutExtension + "_QuantifiedPeptides", new List<string> { Parameters.SearchTaskId, "IndividualFileResults", file.Key.FullFilePathWithExtension });
                        if (true)
                        {
                            WriteProteinQuantificationResultsToTsv(Parameters.FlashLfqResults, Path.Combine(Parameters.IndividualResultsOutputFolder, file.Key.FilenameWithoutExtension),
                                file.Key.FilenameWithoutExtension + "_QuantifiedProteins", new List<string> { Parameters.SearchTaskId, "IndividualFileResults", file.Key.FullFilePathWithExtension });
                        }
                    }
                }
            }
        }

        private void WriteProteinQuantificationResultsToTsv(FlashLfqResults flashLFQResults, string outputFolder, string fileName, List<string> nestedIds)
        {
            var fullSeqPath = Path.Combine(outputFolder, fileName + ".tsv");

            flashLFQResults.WriteResults(null, null, fullSeqPath, null, true);

            FinishedWritingFile(fullSeqPath, nestedIds);
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

