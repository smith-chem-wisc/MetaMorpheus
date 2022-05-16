using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using EngineLayer.Localization;
using EngineLayer.ModificationAnalysis;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class PostSearchAnalysisTask : MetaMorpheusTask
    {
        public PostSearchAnalysisParameters Parameters { get; set; }
        private List<EngineLayer.ProteinGroup> ProteinGroups { get; set; }
        private IEnumerable<IGrouping<string, PeptideSpectralMatch>> PsmsGroupedByFile { get; set; }

        public PostSearchAnalysisTask()
            : base(MyTask.Search)
        {
        }

        public MyTaskResults Run()
        {
            // Stop loop if canceled
            if (GlobalVariables.StopLoops) { return Parameters.SearchTaskResults; }

            if (Parameters.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.ModOpen
                || Parameters.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Open
                || Parameters.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Custom
                )
            {
                // This only makes sense if there is a mass difference that you want to localize. No use for exact and missed monoisotopic mass searches.
                Parameters.SearchParameters.DoLocalizationAnalysis = true;
            }
            else
            {
                Parameters.SearchParameters.DoLocalizationAnalysis = false;
            }

            //update all psms with peptide info
            if (Parameters.SearchParameters.SearchType != SearchType.NonSpecific) //if it hasn't been done already
            {
                Parameters.AllPsms = Parameters.AllPsms.Where(psm => psm != null).ToList();
                Parameters.AllPsms.ForEach(psm => psm.ResolveAllAmbiguities());
                Parameters.AllPsms = Parameters.AllPsms.OrderByDescending(b => b.Score)
                   .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                   .GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

                CalculatePsmFdr();
            }

            DoMassDifferenceLocalizationAnalysis();
            ProteinAnalysis();
            QuantificationAnalysis();

            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files" }));

            HistogramAnalysis();
            WritePsmResults();
            WriteProteinResults();
            WriteQuantificationResults();
            WritePrunedDatabase();
            if (Parameters.SearchParameters.WriteSpectralLibrary)
            {
                SpectralLibraryGeneration();
                PostQuantificationMbrAnalysis();
            }
            if (Parameters.ProteinList.Any((p => p.AppliedSequenceVariations.Count > 0)))
            {
                WriteVariantResults();
            }

            WritePeptideResults(); // modifies the FDR results for PSMs, so do this last
            CompressIndividualFileResults();
            return Parameters.SearchTaskResults;
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            return null;
        }

        /// <summary>
        /// Calculate estimated false-discovery rate (FDR) for peptide spectral matches (PSMs)
        /// </summary>
        private void CalculatePsmFdr()
        {
            // TODO: because FDR is done before parsimony, if a PSM matches to a target and a decoy protein, there may be conflicts between how it's handled in parsimony and the FDR engine here
            // for example, here it may be treated as a decoy PSM, where as in parsimony it will be determined by the parsimony algorithm which is agnostic of target/decoy assignments
            // this could cause weird PSM FDR issues

            Status("Estimating PSM FDR...", Parameters.SearchTaskId);
            new FdrAnalysisEngine(Parameters.AllPsms, Parameters.NumNotches, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }, outputFolder: Parameters.OutputFolder).Run();

            // sort by q-value because of group FDR stuff
            // e.g. multiprotease FDR, non/semi-specific protease, etc
            Parameters.AllPsms = Parameters.AllPsms
                .OrderBy(p => p.FdrInfo.QValue)
                .ThenByDescending(p => p.Score)
                .ThenBy(p => p.FdrInfo.CumulativeTarget)
                .ToList();

            Status("Done estimating PSM FDR!", Parameters.SearchTaskId);
        }

        private void ProteinAnalysis()
        {
            if (!Parameters.SearchParameters.DoParsimony)
            {
                return;
            }

            Status("Constructing protein groups...", Parameters.SearchTaskId);

            //if SILAC, modify the proteins to appear only light (we want a protein sequence to look like PROTEINK instead of PROTEINa)
            if (Parameters.SearchParameters.SilacLabels != null)
            {
                Parameters.AllPsms = SilacConversions.UpdateProteinSequencesToLight(Parameters.AllPsms, Parameters.SearchParameters.SilacLabels);
            }

            List<PeptideSpectralMatch> psmsForProteinParsimony = Parameters.AllPsms;

            // run parsimony
            ProteinParsimonyResults proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(psmsForProteinParsimony, Parameters.SearchParameters.ModPeptidesAreDifferent, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run());

            // score protein groups and calculate FDR
            ProteinScoringAndFdrResults proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psmsForProteinParsimony,
                Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent, true, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run();

            ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

            foreach (PeptideSpectralMatch psm in Parameters.AllPsms)
            {
                psm.ResolveAllAmbiguities();
            }

            Status("Done constructing protein groups!", Parameters.SearchTaskId);
        }

        private void DoMassDifferenceLocalizationAnalysis()
        {
            if (Parameters.SearchParameters.DoLocalizationAnalysis)
            {
                Status("Running mass-difference localization analysis...", Parameters.SearchTaskId);
                for (int spectraFileIndex = 0; spectraFileIndex < Parameters.CurrentRawFileList.Count; spectraFileIndex++)
                {
                    CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, Parameters.FileSettingsList[spectraFileIndex]);

                    var origDataFile = Parameters.CurrentRawFileList[spectraFileIndex];
                    Status("Running mass-difference localization analysis...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile });
                    MsDataFile myMsDataFile = Parameters.MyFileManager.LoadFile(origDataFile, combinedParams);
                    new LocalizationEngine(Parameters.AllPsms.Where(b => b.FullFilePath.Equals(origDataFile)).ToList(),
                        myMsDataFile, combinedParams, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile }).Run();
                    Parameters.MyFileManager.DoneWithFile(origDataFile);
                    ReportProgress(new ProgressEventArgs(100, "Done with localization analysis!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile }));
                }
            }

            // count different modifications observed
            new ModificationAnalysisEngine(Parameters.AllPsms, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run();
        }

        private void QuantificationAnalysis()
        {
            if (!Parameters.SearchParameters.DoQuantification)
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
            else if (Parameters.SearchParameters.Normalize)
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
            var unambiguousPsmsBelowOnePercentFdr = Parameters.AllPsms.Where(p =>
                p.FdrInfo.QValue <= 0.01
                && p.FdrInfo.QValueNotch <= 0.01
                && !p.IsDecoy
                && p.FullSequence != null).ToList(); //if ambiguous, there's no full sequence

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

            //If SILAC (Pre-Quantification) preparation
            //setup variables
            List<SilacLabel> allSilacLabels = Parameters.SearchParameters.SilacLabels;
            SilacLabel startLabel = Parameters.SearchParameters.StartTurnoverLabel;
            SilacLabel endLabel = Parameters.SearchParameters.EndTurnoverLabel;
            bool quantifyUnlabeledPeptides = Parameters.ListOfDigestionParams.Any(x => x.GeneratehUnlabeledProteinsForSilac);
            if (Parameters.SearchParameters.SilacLabels != null)
            {
                bool turnoverWithMultipleLabels = startLabel != null && endLabel != null; //used to check multiple labels
                //go through all the psms and duplicate them until a psm copy exists for the unlabeled and labeled proteins
                //The number of psms should roughly increase by a factor of N, where N is the number of labels.
                //It may not increase exactly by a factor of N if the amino acid(s) that gets labeled doesn't exist in the peptide

                List<PeptideSpectralMatch> silacPsms = new List<PeptideSpectralMatch>(); //populate with duplicate psms for heavy/light

                //multiply the psms by the number of labels
                foreach (PeptideSpectralMatch psm in unambiguousPsmsBelowOnePercentFdr)
                {
                    //get the original proteinGroup to give to the other psm clones
                    List<FlashLFQ.ProteinGroup> originalProteinGroups = psmToProteinGroups.ContainsKey(psm) ? psmToProteinGroups[psm] : new List<FlashLFQ.ProteinGroup>();

                    //see which label, if any, this peptide has
                    string peptideBaseSequence = psm.BaseSequence;
                    SilacLabel observedLabel = SilacConversions.GetRelevantLabelFromBaseSequence(peptideBaseSequence, allSilacLabels); //returns null if no label

                    //if it's not the light form, make a light form to start as a base.
                    PeptideSpectralMatch lightPsm = observedLabel == null ? psm : SilacConversions.GetSilacPsm(psm, observedLabel);

                    //get easy access to values we need for new psm generation
                    string unlabeledBaseSequence = lightPsm.BaseSequence;
                    int notch = psm.BestMatchingPeptides.First().Notch;
                    PeptideWithSetModifications pwsm = psm.BestMatchingPeptides.First().Peptide;

                    //check if turnover or multiplex experiment
                    if (startLabel == null && endLabel == null) //if multiplex
                    {
                        //If we need the light form, then add it
                        if (quantifyUnlabeledPeptides)
                        {
                            silacPsms.Add(lightPsm);
                            if (originalProteinGroups != null)
                            {
                                psmToProteinGroups[lightPsm] = originalProteinGroups; //add proteingroup info
                            }
                        }
                        //take the unlabeled sequence and modify everything for each label
                        foreach (SilacLabel label in allSilacLabels)
                        {
                            if (!label.Equals(observedLabel)) //if we need to change the psm
                            {
                                string labeledBaseSequence = SilacConversions.GetLabeledBaseSequence(unlabeledBaseSequence, label);

                                //check we're not adding a duplicate (happens if none of the heavy amino acids are present in the peptide)
                                if (!labeledBaseSequence.Equals(unlabeledBaseSequence))
                                {
                                    PeptideSpectralMatch labeledPsm = SilacConversions.GetLabeledPsm(psm, notch, pwsm, labeledBaseSequence);
                                    if (originalProteinGroups != null)
                                    {
                                        psmToProteinGroups[labeledPsm] = originalProteinGroups; //add proteingroup info
                                    }
                                    silacPsms.Add(labeledPsm);
                                }
                            }
                            else
                            {
                                silacPsms.Add(psm);
                            }
                        }
                    }
                    else //if turnover
                    {
                        //if it's possible that mixtures exist, then multiple labels must be removed
                        if (turnoverWithMultipleLabels)
                        {
                            observedLabel = SilacConversions.GetRelevantLabelFromBaseSequence(unlabeledBaseSequence, allSilacLabels); //returns null if no label
                            lightPsm = observedLabel == null ? lightPsm : SilacConversions.GetSilacPsm(lightPsm, observedLabel);
                            unlabeledBaseSequence = lightPsm.BaseSequence;
                        }

                        //Convert everything to the startLabel
                        string startLabeledBaseSequence = SilacConversions.GetLabeledBaseSequence(unlabeledBaseSequence, startLabel);

                        //Convert everything to the endLabel
                        string endLabeledBaseSequence = SilacConversions.GetLabeledBaseSequence(unlabeledBaseSequence, endLabel);

                        //figure out which residues can be changed
                        List<int> differentResidues = new List<int>();
                        List<char> startUniqueResidues = new List<char>();
                        List<char> endUniqueResidues = new List<char>();
                        for (int i = 0; i < startLabeledBaseSequence.Length; i++)
                        {
                            char startChar = startLabeledBaseSequence[i];
                            char endChar = endLabeledBaseSequence[i];
                            if (startChar != endChar)
                            {
                                differentResidues.Add(i);
                                startUniqueResidues.Add(startChar);
                                endUniqueResidues.Add(endChar);
                            }
                        }

                        //create missed cleavage combinations (only create HL; don't make HL and LH or FlashLFQ freaks at the ambiguity)
                        List<string> labeledBaseSequences = new List<string> { startLabeledBaseSequence };
                        string sequenceToModify = startLabeledBaseSequence;
                        for (int i = 0; i < differentResidues.Count; i++)
                        {
                            char[] charArray = sequenceToModify.ToCharArray();
                            charArray[differentResidues[i]] = endUniqueResidues[i];
                            sequenceToModify = string.Concat(charArray);
                            labeledBaseSequences.Add(sequenceToModify);
                        }
                        //add them
                        foreach (string sequence in labeledBaseSequences)
                        {
                            PeptideSpectralMatch labeledPsm = SilacConversions.GetLabeledPsm(psm, notch, pwsm, sequence);
                            silacPsms.Add(labeledPsm);
                            psmToProteinGroups[labeledPsm] = originalProteinGroups; //add proteingroup info
                        }
                    }
                }
                //update the list for FlashLFQ
                silacPsms.ForEach(x => x.ResolveAllAmbiguities()); //update the monoisotopic mass
                unambiguousPsmsBelowOnePercentFdr = silacPsms;
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
                normalize: Parameters.SearchParameters.Normalize,
                ppmTolerance: Parameters.SearchParameters.QuantifyPpmTol,
                matchBetweenRuns: Parameters.SearchParameters.MatchBetweenRuns,
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

            //Silac stuff for post-quantification
            if (Parameters.SearchParameters.SilacLabels != null) //if we're doing silac
            {
                SilacConversions.SilacConversionsPostQuantification(allSilacLabels, startLabel, endLabel, spectraFileInfo, ProteinGroups, Parameters.ListOfDigestionParams,
                    Parameters.FlashLfqResults, Parameters.AllPsms, Parameters.SearchParameters.ModsToWriteSelection, quantifyUnlabeledPeptides);
            }
        }

        private void HistogramAnalysis()
        {
            if (Parameters.SearchParameters.DoHistogramAnalysis)
            {
                var limitedpsms_with_fdr = Parameters.AllPsms.Where(b => (b.FdrInfo.QValue <= 0.01)).ToList();
                if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                {
                    Status("Running histogram analysis...", new List<string> { Parameters.SearchTaskId });
                    var myTreeStructure = new BinTreeStructure();
                    myTreeStructure.GenerateBins(limitedpsms_with_fdr, Parameters.SearchParameters.HistogramBinTolInDaltons);
                    var writtenFile = Path.Combine(Parameters.OutputFolder, "MassDifferenceHistogram.tsv");
                    WriteTree(myTreeStructure, writtenFile);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });
                }
            }
        }

        private void WritePsmResults()
        {
            Status("Writing PSM results...", Parameters.SearchTaskId);
            var FilteredPsmListForOutput = Parameters.AllPsms
                .Where(p => p.FdrInfo.QValue <= CommonParameters.QValueOutputFilter
                && p.FdrInfo.QValueNotch <= CommonParameters.QValueOutputFilter).ToList();

            if (!Parameters.SearchParameters.WriteDecoys)
            {
                FilteredPsmListForOutput.RemoveAll(b => b.IsDecoy);
            }
            if (!Parameters.SearchParameters.WriteContaminants)
            {
                FilteredPsmListForOutput.RemoveAll(b => b.IsContaminant);
            }

            // write PSMs
            string writtenFile = Path.Combine(Parameters.OutputFolder, "AllPSMs.psmtsv");
            WritePsmsToTsv(FilteredPsmListForOutput, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            // write PSMs for percolator
            // percolator native read format is .tab
            writtenFile = Path.Combine(Parameters.OutputFolder, "AllPSMs_FormattedForPercolator.tab");
            WritePsmsForPercolator(FilteredPsmListForOutput, writtenFile);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            // write summary text
            Parameters.SearchTaskResults.AddPsmPeptideProteinSummaryText("All target PSMS within 1% FDR: " + Parameters.AllPsms.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy) + Environment.NewLine);
            if (Parameters.SearchParameters.DoParsimony)
            {
                Parameters.SearchTaskResults.AddTaskSummaryText("All target protein groups within 1% FDR: " + ProteinGroups.Count(b => b.QValue <= 0.01 && !b.IsDecoy)
                    + Environment.NewLine);
            }

            PsmsGroupedByFile = FilteredPsmListForOutput.GroupBy(p => p.FullFilePath);

            foreach (var file in PsmsGroupedByFile)
            {
                // write summary text
                var psmsForThisFile = file.ToList();
                string strippedFileName = Path.GetFileNameWithoutExtension(file.First().FullFilePath);

                Parameters.SearchTaskResults.AddTaskSummaryText("MS2 spectra in " + strippedFileName + ": " + Parameters.NumMs2SpectraPerFile[strippedFileName][0]);
                Parameters.SearchTaskResults.AddTaskSummaryText("Precursors fragmented in " + strippedFileName + ": " + Parameters.NumMs2SpectraPerFile[strippedFileName][1]);
                Parameters.SearchTaskResults.AddTaskSummaryText("Target PSMs within 1% FDR in " + strippedFileName + ": " + psmsForThisFile.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy));

                // writes all individual spectra file search results to subdirectory
                if (Parameters.CurrentRawFileList.Count > 1 && Parameters.SearchParameters.WriteIndividualFiles)
                {
                    // create individual files subdirectory
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);

                    // write PSMs
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_PSMs.psmtsv");
                    WritePsmsToTsv(psmsForThisFile, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.First().FullFilePath });

                    // write PSMs for percolator
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_PSMsFormattedForPercolator.tab");
                    WritePsmsForPercolator(psmsForThisFile, writtenFile);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.First().FullFilePath });
                }
            }
        }

        //for those spectra matching the same peptide/protein with same charge, save the one with highest score
        private void SpectralLibraryGeneration()
        {
            var FilteredPsmList = Parameters.AllPsms
               .Where(p => p.FdrInfo.PEP_QValue <= 0.01 && p.FdrInfo.QValueNotch <= CommonParameters.QValueOutputFilter).ToList();
            FilteredPsmList.RemoveAll(b => b.IsDecoy);
            FilteredPsmList.RemoveAll(b => b.IsContaminant);
            //group psms by peptide and charge, the psms having same sequence and same charge will be in the same group
            Dictionary<(String, int), List<PeptideSpectralMatch>> PsmsGroupByPeptideAndCharge = new Dictionary<(String, int), List<PeptideSpectralMatch>>();
            foreach (var x in FilteredPsmList)
            {
                List<PeptideSpectralMatch> psmsWithSamePeptideAndSameCharge = FilteredPsmList.Where(b => b.FullSequence == x.FullSequence && b.ScanPrecursorCharge == x.ScanPrecursorCharge).OrderByDescending(p => p.Score).ToList();
                (String, int) peptideWithChargeState = (x.FullSequence, x.ScanPrecursorCharge);

                if (!PsmsGroupByPeptideAndCharge.ContainsKey(peptideWithChargeState))
                {
                    PsmsGroupByPeptideAndCharge.Add(peptideWithChargeState, psmsWithSamePeptideAndSameCharge);
                }
            }
            var spectraLibrary = new List<LibrarySpectrum>();
            foreach (var psm in PsmsGroupByPeptideAndCharge)
            {
                var standardSpectrum = new LibrarySpectrum(psm.Value[0].FullSequence, psm.Value[0].ScanPrecursorMonoisotopicPeakMz, psm.Value[0].ScanPrecursorCharge, psm.Value[0].MatchedFragmentIons, psm.Value[0].ScanRetentionTime);
                spectraLibrary.Add(standardSpectrum);
            }
            WriteSpectralLibrary(spectraLibrary, Parameters.OutputFolder);
        }

        private void WriteProteinResults()
        {
            if (Parameters.SearchParameters.DoParsimony)
            {
                string fileName = "AllProteinGroups.tsv";

                if (Parameters.SearchParameters.DoQuantification)
                {
                    fileName = "AllQuantifiedProteinGroups.tsv";
                }

                //set peptide output values
                ProteinGroups.ForEach(x => x.GetIdentifiedPeptidesOutput(Parameters.SearchParameters.SilacLabels));
                // write protein groups to tsv
                string writtenFile = Path.Combine(Parameters.OutputFolder, fileName);
                WriteProteinGroupsToTsv(ProteinGroups, writtenFile, new List<string> { Parameters.SearchTaskId }, CommonParameters.QValueOutputFilter);

                //If SILAC and no light, we need to update the psms (which were found in the "light" file) and say they were found in the "heavy" file)
                if (Parameters.SearchParameters.SilacLabels != null)
                {
                    //get the light filenames
                    List<string> fileNamesThatHadPsms = PsmsGroupedByFile.Select(v => v.Key).ToList();
                    EngineLayer.ProteinGroup firstProteinGroup = ProteinGroups.FirstOrDefault();
                    if (firstProteinGroup != null)
                    {
                        List<string> allFileNames = firstProteinGroup.FilesForQuantification.Select(x => x.FullFilePathWithExtension).ToList();
                        string heavyFileToSet = allFileNames.Where(x => !fileNamesThatHadPsms.Contains(x)).First();
                        List<PeptideSpectralMatch> psms = PsmsGroupedByFile.SelectMany(g => g).ToList();
                        PsmsGroupedByFile = psms.GroupBy(x => heavyFileToSet); //set them all to the same file
                    }
                }

                // write all individual file results to subdirectory
                // local protein fdr, global parsimony, global psm fdr

                if (Parameters.CurrentRawFileList.Count > 1 && (Parameters.SearchParameters.WriteIndividualFiles
                    || Parameters.SearchParameters.WriteMzId || Parameters.SearchParameters.WritePepXml))
                {
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);
                }

                foreach (var fullFilePath in PsmsGroupedByFile.Select(v => v.Key))
                {
                    string strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                    List<PeptideSpectralMatch> psmsForThisFile = PsmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList();
                    var subsetProteinGroupsForThisFile = ProteinGroups.Select(p => p.ConstructSubsetProteinGroup(fullFilePath)).ToList();

                    ProteinScoringAndFdrResults subsetProteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(subsetProteinGroupsForThisFile, psmsForThisFile,
                        Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent,
                        false, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }).Run();

                    subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults.SortedAndScoredProteinGroups;

                    Parameters.SearchTaskResults.AddTaskSummaryText("Target protein groups within 1 % FDR in " + strippedFileName + ": " + subsetProteinGroupsForThisFile.Count(b => b.QValue <= 0.01 && !b.IsDecoy));

                    // write individual spectra file protein groups results to tsv
                    if (Parameters.SearchParameters.WriteIndividualFiles && Parameters.CurrentRawFileList.Count > 1)
                    {
                        writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_ProteinGroups.tsv");
                        WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }, CommonParameters.QValueOutputFilter);
                    }

                    // write mzID
                    if (Parameters.SearchParameters.WriteMzId)
                    {
                        Status("Writing mzID...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });

                        string mzidFilePath = Path.Combine(Parameters.OutputFolder, strippedFileName + ".mzID");
                        if (Parameters.CurrentRawFileList.Count > 1)
                        {
                            mzidFilePath = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + ".mzID");
                        }
                        MzIdentMLWriter.WriteMzIdentMl(psmsForThisFile.Where(p => p.FdrInfo.QValue <= CommonParameters.QValueOutputFilter),
                            subsetProteinGroupsForThisFile, Parameters.VariableModifications, Parameters.FixedModifications, Parameters.SearchParameters.SilacLabels,
                            new List<Protease> { CommonParameters.DigestionParams.Protease }, CommonParameters.QValueOutputFilter, CommonParameters.ProductMassTolerance,
                            CommonParameters.PrecursorMassTolerance, CommonParameters.DigestionParams.MaxMissedCleavages, mzidFilePath,
                            Parameters.SearchParameters.IncludeModMotifInMzid);

                        FinishedWritingFile(mzidFilePath, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                    }

                    // write pepXML
                    if (Parameters.SearchParameters.WritePepXml)
                    {
                        Status("Writing pepXML...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });

                        string pepXMLFilePath = Path.Combine(Parameters.OutputFolder, strippedFileName + ".pep.XML");
                        if (Parameters.CurrentRawFileList.Count > 1)
                        {
                            pepXMLFilePath = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + ".pep.XML");
                        }

                        PepXMLWriter.WritePepXml(psmsForThisFile.Where(p => p.FdrInfo.QValue <= CommonParameters.QValueOutputFilter).ToList(), Parameters.DatabaseFilenameList, Parameters.VariableModifications, Parameters.FixedModifications,
                            CommonParameters, pepXMLFilePath, CommonParameters.QValueOutputFilter);

                        FinishedWritingFile(pepXMLFilePath, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                    }

                    ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }));
                }
            }
        }

        private void WriteQuantificationResults()
        {
            if (Parameters.SearchParameters.DoQuantification && Parameters.FlashLfqResults != null)
            {
                // write peaks
                WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "AllQuantifiedPeaks", new List<string> { Parameters.SearchTaskId });

                // write peptide quant results
                string filename = "AllQuantified" + GlobalVariables.AnalyteType + "s";
                WritePeptideQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, filename, new List<string> { Parameters.SearchTaskId });

                // write individual results
                if (Parameters.CurrentRawFileList.Count > 1 && Parameters.SearchParameters.WriteIndividualFiles)
                {
                    foreach (var file in Parameters.FlashLfqResults.Peaks)
                    {
                        WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.IndividualResultsOutputFolder,
                            file.Key.FilenameWithoutExtension + "_QuantifiedPeaks", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.Key.FullFilePathWithExtension });
                    }
                }
            }
        }

        private void PostQuantificationMbrAnalysis()
        {
            if (!Parameters.SearchParameters.DoMbrAnalysis)
            {
                return;
            }
            List<SpectraFileInfo> spectraFiles = Parameters.FlashLfqResults.Peaks.Select(p => p.Key).ToList();
            List<PeptideSpectralMatch> allPeptides = GetAllPeptides();
            List<(ChromatographicPeak, PeptideSpectralMatch)> bestPsmsForPeaks = new();

            foreach (SpectraFileInfo spectraFile in spectraFiles)
            {
                List<ChromatographicPeak> fileSpecificMbrPeaks = Parameters.FlashLfqResults.Peaks[spectraFile].Where(p => p.IsMbrPeak).ToList();
                if (fileSpecificMbrPeaks == null || (!fileSpecificMbrPeaks.Any())) break;

                MyFileManager myFileManager = new(true);
                MsDataFile myMsDataFile = myFileManager.LoadFile(spectraFile.FullFilePathWithExtension, CommonParameters);
                MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance,
                    Parameters.SearchParameters.MassDiffAcceptorType, Parameters.SearchParameters.CustomMdac);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByRT = GetMs2Scans(myMsDataFile, spectraFile.FullFilePathWithExtension, CommonParameters)
                    .OrderBy(b => b.RetentionTime).ToArray();
                double[] arrayOfRTs = arrayOfMs2ScansSortedByRT.Select(p => p.TheScan.RetentionTime).ToArray();

                

                foreach (ChromatographicPeak mbrPeak in fileSpecificMbrPeaks)
                {
                    //TODO: check if this really is the best donor PSM. Is there a way to get the peptide from the flashLFQ results? Maybe this is better anyway, IDK
                    PeptideSpectralMatch bestDonorPsm = allPeptides.Where(p => p.FullSequence == mbrPeak.Identifications.First().ModifiedSequence).First();
                    PeptideWithSetModifications bestDonorPwsm = bestDonorPsm.BestMatchingPeptides.First().Peptide;
                    double monoIsotopicMass = bestDonorPsm.PeptideMonisotopicMass.Value;

                    // Find MS2 scans falling within the relevant time window.
                    double apexRT = mbrPeak.Apex.IndexedPeak.RetentionTime;
                    double peakHalfWidth = 1.0; //Placeholder value to determine retention time window
                    int startIndex = Array.BinarySearch(arrayOfRTs, apexRT - peakHalfWidth);
                    if (startIndex < 0)
                        startIndex = ~startIndex;
                    int endIndex = Array.BinarySearch(arrayOfRTs, apexRT + peakHalfWidth);
                    if (endIndex < 0)
                        endIndex = ~endIndex;
                    Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = arrayOfMs2ScansSortedByRT[startIndex..endIndex].OrderBy(b => b.PrecursorMass).ToArray();

                    string spectralLibraryPath = Path.Combine(Parameters.OutputFolder, @"spectralLibrary.msp");
                    SpectralLibrary library = new(new List<string>() { spectralLibraryPath });
                    PeptideSpectralMatch[] peptideSpectralMatches = new PeptideSpectralMatch[arrayOfMs2ScansSortedByRT.Count()];
                    MiniClassicSearchEngine mcse = new(bestDonorPwsm, peptideSpectralMatches, arrayOfMs2ScansSortedByMass, Parameters.VariableModifications, Parameters.FixedModifications,
                        massDiffAcceptor, CommonParameters, FileSpecificParameters, library, new List<string> { Parameters.SearchTaskId });
                    mcse.Run();

                    if (peptideSpectralMatches.Any())
                    {
                        bestPsmsForPeaks.Add((mbrPeak, BestPsmForMbrPeak(peptideSpectralMatches)));
                    }

                    //write result somewhere. 
                }
            }
            WriteMbrPsmResults(bestPsmsForPeaks);
        }

        private void WriteMbrPsmResults(List<(ChromatographicPeak, PeptideSpectralMatch)> bestPsmsForPeaks)
        {
            string mbrOutputPath = Path.Combine(Parameters.OutputFolder, @"MbrAnalysis.psmtsv");
            using (StreamWriter output = new StreamWriter(mbrOutputPath))
            {
                output.WriteLine(TaskLayer.MbrWriter.TabSeparatedHeader);
                foreach (var peak in bestPsmsForPeaks)
                {
                    if (peak.Item2 != null)
                    {
                        output.WriteLine(peak.Item2.ToString() + "\t" + peak.Item1.ToString());
                    } 
                }
            }

            string unmatchedPeaksPath = Path.Combine(Parameters.OutputFolder, @"UnmatchedPeaks.tsv");
            using (StreamWriter output = new StreamWriter(unmatchedPeaksPath))
            {
                output.WriteLine(ChromatographicPeak.TabSeparatedHeader);
                foreach (var peak in bestPsmsForPeaks)
                {
                    if (peak.Item2 == null)
                    {
                        output.WriteLine(peak.Item1.ToString());
                    }
                }
            }
        }

        private PeptideSpectralMatch BestPsmForMbrPeak(PeptideSpectralMatch[] peptideSpectralMatches)
        {
            List<PeptideSpectralMatch> nonNullPsms = peptideSpectralMatches.Where(p=>p != null).ToList();

            if(nonNullPsms.Any())
            {
                // Setting Qvalue, QValueNotch, PEP, and PEP_Qvalue equal to 0 is necessary for MetaDraw to read the .psmtsv
                foreach (PeptideSpectralMatch psm in nonNullPsms)
                {
                    psm.SetFdrValues(0,0,0,0,0,0,0,0);
                }
                if (nonNullPsms.Select(p=>p.SpectralAngle).Any(g=>g != double.NaN))
                {
                    double maxSpectralAngle = nonNullPsms.Select(p => p.SpectralAngle).Max();
                    return nonNullPsms.Where(p => p.SpectralAngle == maxSpectralAngle).FirstOrDefault();
                }
                return nonNullPsms.FirstOrDefault();  
            }
            return null;
        }

        private List<PeptideSpectralMatch> GetAllPeptides()
        {
            List<PeptideSpectralMatch> peptides = Parameters.AllPsms.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();

            new FdrAnalysisEngine(peptides, Parameters.NumNotches, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }, "Peptide").Run();

            if (!Parameters.SearchParameters.WriteDecoys)
            {
                peptides.RemoveAll(b => b.IsDecoy);
            }
            if (!Parameters.SearchParameters.WriteContaminants)
            {
                peptides.RemoveAll(b => b.IsContaminant);
            }
            peptides.RemoveAll(p => p.FdrInfo.QValue > CommonParameters.QValueOutputFilter);
            return peptides;
        }

        private void WritePrunedDatabase()
        {
            if (Parameters.SearchParameters.WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { Parameters.SearchTaskId });
                HashSet<Modification> modificationsToWriteIfBoth = new HashSet<Modification>();
                HashSet<Modification> modificationsToWriteIfInDatabase = new HashSet<Modification>();
                HashSet<Modification> modificationsToWriteIfObserved = new HashSet<Modification>();

                var confidentPsms = Parameters.AllPsms.Where(b => b.FdrInfo.QValueNotch <= 0.01 && b.FdrInfo.QValue <= 0.01 && !b.IsDecoy && b.BaseSequence != null).ToList();
                var proteinToConfidentBaseSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                // associate all confident PSMs with all possible proteins they could be digest products of (before or after parsimony)
                foreach (PeptideSpectralMatch psm in confidentPsms)
                {
                    var myPepsWithSetMods = psm.BestMatchingPeptides.Select(p => p.Peptide);

                    foreach (var peptide in myPepsWithSetMods)
                    {
                        if (proteinToConfidentBaseSequences.TryGetValue(peptide.Protein.NonVariantProtein, out var myPepList))
                        {
                            myPepList.Add(peptide);
                        }
                        else
                        {
                            proteinToConfidentBaseSequences.Add(peptide.Protein.NonVariantProtein, new List<PeptideWithSetModifications> { peptide });
                        }
                    }
                }

                // Add user mod selection behavours to Pruned DB
                foreach (var modType in Parameters.SearchParameters.ModsToWriteSelection)
                {
                    foreach (Modification mod in GlobalVariables.AllModsKnown.Where(b => b.ModificationType.Equals(modType.Key)))
                    {
                        if (modType.Value == 1) // Write if observed and in database
                        {
                            modificationsToWriteIfBoth.Add(mod);
                        }
                        if (modType.Value == 2) // Write if in database
                        {
                            modificationsToWriteIfInDatabase.Add(mod);
                        }
                        if (modType.Value == 3) // Write if observed
                        {
                            modificationsToWriteIfObserved.Add(mod);
                        }
                    }
                }

                //generates dictionary of proteins with only localized modifications
                var ModPsms = Parameters.AllPsms.Where(b => b.FdrInfo.QValueNotch <= 0.01 && b.FdrInfo.QValue <= 0.01 && !b.IsDecoy && b.FullSequence != null).ToList();
                var proteinToConfidentModifiedSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                foreach (PeptideSpectralMatch psm in ModPsms)
                {
                    var myPepsWithSetMods = psm.BestMatchingPeptides.Select(p => p.Peptide);

                    foreach (var peptide in myPepsWithSetMods)
                    {
                        if (proteinToConfidentModifiedSequences.TryGetValue(peptide.Protein.NonVariantProtein, out var myPepList))
                        {
                            myPepList.Add(peptide);
                        }
                        else
                        {
                            proteinToConfidentModifiedSequences.Add(peptide.Protein.NonVariantProtein, new List<PeptideWithSetModifications> { peptide });
                        }
                    }
                }

                Dictionary<Protein, Dictionary<int, List<Modification>>> proteinsOriginalModifications = new Dictionary<Protein, Dictionary<int, List<Modification>>>();
                Dictionary<SequenceVariation, Dictionary<int, List<Modification>>> originalSequenceVariantModifications = new Dictionary<SequenceVariation, Dictionary<int, List<Modification>>>();

                // mods included in pruned database will only be confidently localized mods (peptide's FullSequence != null)
                foreach (var nonVariantProtein in Parameters.ProteinList.Select(p => p.NonVariantProtein).Distinct())
                {
                    if (!nonVariantProtein.IsDecoy)
                    {
                        proteinToConfidentModifiedSequences.TryGetValue(nonVariantProtein, out var psms);
                        HashSet<(int, Modification, SequenceVariation)> modsObservedOnThisProtein = new HashSet<(int, Modification, SequenceVariation)>(); // sequence variant is null if mod is not on a variant
                        foreach (PeptideWithSetModifications psm in psms ?? new List<PeptideWithSetModifications>())
                        {
                            foreach (var idxModKV in psm.AllModsOneIsNterminus)
                            {
                                int proteinIdx = GetOneBasedIndexInProtein(idxModKV.Key, psm);
                                SequenceVariation relevantVariant = psm.Protein.AppliedSequenceVariations.FirstOrDefault(sv => VariantApplication.IsSequenceVariantModification(sv, proteinIdx));
                                SequenceVariation unappliedVariant =
                                    relevantVariant == null ? null : // it's not a sequence variant mod
                                    psm.Protein.SequenceVariations.FirstOrDefault(sv => sv.Description != null && sv.Description.Equals(relevantVariant.Description));
                                modsObservedOnThisProtein.Add((VariantApplication.RestoreModificationIndex(psm.Protein, proteinIdx), idxModKV.Value, unappliedVariant));
                            }
                        }

                        IDictionary<(SequenceVariation, int), List<Modification>> modsToWrite = new Dictionary<(SequenceVariation, int), List<Modification>>();

                        //Add if observed (regardless if in database)
                        foreach (var observedMod in modsObservedOnThisProtein)
                        {
                            var tempMod = observedMod.Item2;

                            if (modificationsToWriteIfObserved.Contains(tempMod))
                            {
                                var svIdxKey = (observedMod.Item3, observedMod.Item1);
                                if (!modsToWrite.ContainsKey(svIdxKey))
                                {
                                    modsToWrite.Add(svIdxKey, new List<Modification> { observedMod.Item2 });
                                }
                                else
                                {
                                    modsToWrite[svIdxKey].Add(observedMod.Item2);
                                }
                            }
                        }

                        // Add modification if in database (two cases: always or if observed)
                        foreach (var modkv in nonVariantProtein.OneBasedPossibleLocalizedModifications)
                        {
                            foreach (var mod in modkv.Value)
                            {
                                //Add if always In Database or if was observed and in database and not set to not include
                                if (modificationsToWriteIfInDatabase.Contains(mod) ||
                                    (modificationsToWriteIfBoth.Contains(mod) && modsObservedOnThisProtein.Contains((modkv.Key, mod, null))))
                                {
                                    if (!modsToWrite.ContainsKey((null, modkv.Key)))
                                    {
                                        modsToWrite.Add((null, modkv.Key), new List<Modification> { mod });
                                    }
                                    else
                                    {
                                        modsToWrite[(null, modkv.Key)].Add(mod);
                                    }
                                }
                            }
                        }

                        //TODO add unit test here
                        // Add variant modification if in database (two cases: always or if observed)
                        foreach (SequenceVariation sv in nonVariantProtein.SequenceVariations)
                        {
                            foreach (var modkv in sv.OneBasedModifications)
                            {
                                foreach (var mod in modkv.Value)
                                {
                                    //Add if always In Database or if was observed and in database and not set to not include
                                    if (modificationsToWriteIfInDatabase.Contains(mod) ||
                                        (modificationsToWriteIfBoth.Contains(mod) && modsObservedOnThisProtein.Contains((modkv.Key, mod, sv))))
                                    {
                                        if (!modsToWrite.ContainsKey((sv, modkv.Key)))
                                        {
                                            modsToWrite.Add((sv, modkv.Key), new List<Modification> { mod });
                                        }
                                        else
                                        {
                                            modsToWrite[(sv, modkv.Key)].Add(mod);
                                        }
                                    }
                                }
                            }
                        }

                        var oldMods = nonVariantProtein.OneBasedPossibleLocalizedModifications.ToDictionary(p => p.Key, v => v.Value);
                        if (proteinsOriginalModifications.ContainsKey(nonVariantProtein.NonVariantProtein))
                        {
                            foreach (var entry in oldMods)
                            {
                                if (proteinsOriginalModifications[nonVariantProtein.NonVariantProtein].ContainsKey(entry.Key))
                                {
                                    proteinsOriginalModifications[nonVariantProtein.NonVariantProtein][entry.Key].AddRange(entry.Value);
                                }
                                else
                                {
                                    proteinsOriginalModifications[nonVariantProtein.NonVariantProtein].Add(entry.Key, entry.Value);
                                }
                            }
                        }
                        else
                        {
                            proteinsOriginalModifications.Add(nonVariantProtein.NonVariantProtein, oldMods);
                        }

                        // adds confidently localized and identified mods
                        nonVariantProtein.OneBasedPossibleLocalizedModifications.Clear();
                        foreach (var kvp in modsToWrite.Where(kv => kv.Key.Item1 == null))
                        {
                            nonVariantProtein.OneBasedPossibleLocalizedModifications.Add(kvp.Key.Item2, kvp.Value);
                        }
                        foreach (var sv in nonVariantProtein.SequenceVariations)
                        {
                            var oldVariantModifications = sv.OneBasedModifications.ToDictionary(p => p.Key, v => v.Value);
                            if (originalSequenceVariantModifications.ContainsKey(sv))
                            {
                                foreach (var entry in oldVariantModifications)
                                {
                                    if (originalSequenceVariantModifications[sv].ContainsKey(entry.Key))
                                    {
                                        originalSequenceVariantModifications[sv][entry.Key].AddRange(entry.Value);
                                    }
                                    else
                                    {
                                        originalSequenceVariantModifications[sv].Add(entry.Key, entry.Value);
                                    }
                                }
                            }
                            else
                            {
                                originalSequenceVariantModifications.Add(sv, oldVariantModifications);
                            }

                            sv.OneBasedModifications.Clear();
                            foreach (var kvp in modsToWrite.Where(kv => kv.Key.Item1 != null && kv.Key.Item1.Equals(sv)))
                            {
                                sv.OneBasedModifications.Add(kvp.Key.Item2, kvp.Value);
                            }
                        }
                    }
                }

                //writes all proteins
                if (Parameters.DatabaseFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), Parameters.ProteinList.Select(p => p.NonVariantProtein).Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);
                    FinishedWritingFile(outputXMLdbFullName, new List<string> { Parameters.SearchTaskId });
                }
                if (Parameters.DatabaseFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), Parameters.ProteinList.Select(p => p.NonVariantProtein).Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);
                    FinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { Parameters.SearchTaskId });
                }

                //writes only detected proteins
                if (Parameters.DatabaseFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);
                    FinishedWritingFile(outputXMLdbFullName, new List<string> { Parameters.SearchTaskId });
                }
                if (Parameters.DatabaseFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(Parameters.OutputFolder, string.Join("-", Parameters.DatabaseFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");
                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);
                    FinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { Parameters.SearchTaskId });
                }

                foreach (var nonVariantProtein in Parameters.ProteinList.Select(p => p.NonVariantProtein).Distinct())
                {
                    if (!nonVariantProtein.IsDecoy)
                    {
                        nonVariantProtein.OneBasedPossibleLocalizedModifications.Clear();
                        foreach (var originalMod in proteinsOriginalModifications[nonVariantProtein.NonVariantProtein])
                        {
                            nonVariantProtein.OneBasedPossibleLocalizedModifications.Add(originalMod.Key, originalMod.Value);
                        }
                        foreach (var sv in nonVariantProtein.SequenceVariations)
                        {
                            sv.OneBasedModifications.Clear();
                            foreach (var originalVariantMods in originalSequenceVariantModifications[sv])
                            {
                                sv.OneBasedModifications.Add(originalVariantMods.Key, originalVariantMods.Value);
                            }
                        }
                    }
                }
            }
        }

        private void WritePeptideResults()
        {
            Status("Writing peptide results...", Parameters.SearchTaskId);

            // write best (highest-scoring) PSM per peptide
            string filename = "All" + GlobalVariables.AnalyteType + "s.psmtsv";
            string writtenFile = Path.Combine(Parameters.OutputFolder, filename);
            List<PeptideSpectralMatch> peptides = Parameters.AllPsms.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();

            new FdrAnalysisEngine(peptides, Parameters.NumNotches, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }, "Peptide").Run();

            if (!Parameters.SearchParameters.WriteDecoys)
            {
                peptides.RemoveAll(b => b.IsDecoy);
            }
            if (!Parameters.SearchParameters.WriteContaminants)
            {
                peptides.RemoveAll(b => b.IsContaminant);
            }
            peptides.RemoveAll(p => p.FdrInfo.QValue > CommonParameters.QValueOutputFilter);

            WritePsmsToTsv(peptides, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            Parameters.SearchTaskResults.AddPsmPeptideProteinSummaryText("All target " + GlobalVariables.AnalyteType.ToLower() + "s within 1% FDR: " + peptides.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy));

            foreach (var file in PsmsGroupedByFile)
            {
                // write summary text
                var psmsForThisFile = file.ToList();
                string strippedFileName = Path.GetFileNameWithoutExtension(file.First().FullFilePath);
                var peptidesForFile = psmsForThisFile.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).OrderByDescending(b => b.Score).ToList();
                new FdrAnalysisEngine(peptidesForFile, Parameters.NumNotches, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }, "Peptide").Run();
                Parameters.SearchTaskResults.AddTaskSummaryText("Target " + GlobalVariables.AnalyteType.ToLower() + "s within 1% FDR in " + strippedFileName + ": " + peptidesForFile.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy) + Environment.NewLine);

                // writes all individual spectra file search results to subdirectory
                if (Parameters.CurrentRawFileList.Count > 1 && Parameters.SearchParameters.WriteIndividualFiles)
                {
                    // create individual files subdirectory
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);

                    // write best (highest-scoring) PSM per peptide
                    filename = "_" + GlobalVariables.AnalyteType + "s.psmtsv";
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + filename);
                    WritePsmsToTsv(peptidesForFile, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.First().FullFilePath });
                }
            }
        }

        private void CompressIndividualFileResults()
        {
            if (Parameters.SearchParameters.CompressIndividualFiles && Directory.Exists(Parameters.IndividualResultsOutputFolder))
            {
                ZipFile.CreateFromDirectory(Parameters.IndividualResultsOutputFolder, Parameters.IndividualResultsOutputFolder + ".zip");
                Directory.Delete(Parameters.IndividualResultsOutputFolder, true);
            }
        }

        private void WriteVariantResults()
        {
            Status("Writing variant peptide results...", Parameters.SearchTaskId);
            string variantPsmFile = Path.Combine(Parameters.OutputFolder, "VariantPSMs.psmtsv");

            string filename = "Variant" + GlobalVariables.AnalyteType + "s.psmtsv";
            string variantPeptideFile = Path.Combine(Parameters.OutputFolder, filename);
            List<PeptideSpectralMatch> FDRPsms = Parameters.AllPsms
                .Where(p => p.FdrInfo.QValue <= CommonParameters.QValueOutputFilter
                && p.FdrInfo.QValueNotch <= CommonParameters.QValueOutputFilter && p.BaseSequence != null).ToList();
            var possibleVariantPsms = FDRPsms.Where(p => p.BestMatchingPeptides.Any(pep => pep.Peptide.IsVariantPeptide())).OrderByDescending(pep => pep.Score).ToList();

            if (!Parameters.SearchParameters.WriteDecoys)
            {
                possibleVariantPsms.RemoveAll(b => b.IsDecoy);
            }
            if (!Parameters.SearchParameters.WriteContaminants)
            {
                possibleVariantPsms.RemoveAll(b => b.IsContaminant);
            }

            new FdrAnalysisEngine(possibleVariantPsms, Parameters.NumNotches, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }, "variant_PSMs").Run();

            possibleVariantPsms.OrderBy(p => p.FdrInfo.QValue).ThenByDescending(p => p.Score).ThenBy(p => p.FdrInfo.CumulativeTarget).ToList();

            WritePsmsToTsv(possibleVariantPsms, variantPsmFile, Parameters.SearchParameters.ModsToWriteSelection);

            List<PeptideSpectralMatch> variantPeptides = possibleVariantPsms.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).OrderByDescending(b => b.Score).ToList();
            List<PeptideSpectralMatch> confidentVariantPeps = new List<PeptideSpectralMatch>();

            new FdrAnalysisEngine(variantPeptides, Parameters.NumNotches, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }, "variant_Peptides").Run();

            WritePsmsToTsv(variantPeptides, variantPeptideFile, Parameters.SearchParameters.ModsToWriteSelection);

            // if a potential variant peptide can be explained by a canonical protein seqeunce then should not be counted as a confident variant peptide
            //because it is most probable that the peptide originated from the canonical protien.
            foreach (var entry in variantPeptides)
            {
                var pwsm = entry.BestMatchingPeptides;
                var nonVariantOption = pwsm.Any(p => p.Peptide.IsVariantPeptide() == false);
                if (nonVariantOption == false)
                {
                    confidentVariantPeps.Add(entry);
                }
            }
            // count of peptides that contain at least 1 of the given variant type
            int SNVmissenseCount = 0;
            int MNVmissenseCount = 0;
            int insertionCount = 0;
            int deletionCount = 0;
            int frameshiftCount = 0;
            int stopGainCount = 0;
            int stopLossCount = 0;

            // dictionaries facilitate the determination of unique variant sites
            Dictionary<Protein, HashSet<SequenceVariation>> MNVmissenseVariants = new Dictionary<Protein, HashSet<SequenceVariation>>();
            Dictionary<Protein, HashSet<SequenceVariation>> SNVmissenseVariants = new Dictionary<Protein, HashSet<SequenceVariation>>();
            Dictionary<Protein, HashSet<SequenceVariation>> insertionVariants = new Dictionary<Protein, HashSet<SequenceVariation>>();
            Dictionary<Protein, HashSet<SequenceVariation>> deletionVariants = new Dictionary<Protein, HashSet<SequenceVariation>>();
            Dictionary<Protein, HashSet<SequenceVariation>> frameshiftVariants = new Dictionary<Protein, HashSet<SequenceVariation>>();
            Dictionary<Protein, HashSet<SequenceVariation>> stopGainVariants = new Dictionary<Protein, HashSet<SequenceVariation>>();
            Dictionary<Protein, HashSet<SequenceVariation>> stopLossVariants = new Dictionary<Protein, HashSet<SequenceVariation>>();

            double FdrFilterValue = CommonParameters.QValueOutputFilter != 1.0 ? CommonParameters.QValueOutputFilter : 0.01;

            List<PeptideSpectralMatch> modifiedVariantPeptides = confidentVariantPeps.Where(p => p.ModsIdentified != null && p.ModsIdentified.Count > 0 && p.FdrInfo.QValue <= FdrFilterValue && p.FdrInfo.QValueNotch <= FdrFilterValue && !p.IsDecoy && !p.IsContaminant).ToList(); //modification can be on any AA in variant peptide

            List<PeptideSpectralMatch> modifiedVariantSitePeptides = new List<PeptideSpectralMatch>();// modification is speciifcally on the variant residue within the peptide
            foreach (var entry in modifiedVariantPeptides)
            {
                var variantPWSM = entry.BestMatchingPeptides.FirstOrDefault().Peptide;
                var peptideMods = variantPWSM.AllModsOneIsNterminus.Values.ToList();
                var variantProteinModifications = variantPWSM.Protein.OneBasedPossibleLocalizedModifications.Where(k => k.Key >= variantPWSM.OneBasedStartResidueInProtein && k.Key <= variantPWSM.OneBasedEndResidueInProtein).ToList();
                var variants = entry.BestMatchingPeptides.FirstOrDefault().Peptide.Protein.AppliedSequenceVariations.Where(v => entry.BestMatchingPeptides.FirstOrDefault().Peptide.IntersectsAndIdentifiesVariation(v).identifies).ToList();
                bool modifiedVariant = false;
                foreach (var mod in variantProteinModifications)
                {
                    var residue = mod.Key;
                    var peptideModsIdentified = mod.Value.Intersect(peptideMods).ToList().Count;
                    var modOnVariant = variants.Where(p => p.OneBasedBeginPosition >= residue && p.OneBasedEndPosition <= residue);
                    if (modOnVariant.Any() && peptideModsIdentified != 0)
                    {
                        modifiedVariant = true;
                    }
                }
                if (modifiedVariant == true)
                {
                    modifiedVariantSitePeptides.Add(entry);
                }
            }
            foreach (var peptide in confidentVariantPeps)
            {
                if (!peptide.IsDecoy && !peptide.IsContaminant && peptide.FdrInfo.QValue <= FdrFilterValue && peptide.FdrInfo.QValueNotch <= FdrFilterValue)
                {
                    var variantPWSM = peptide.BestMatchingPeptides.FirstOrDefault();//TODO: expand to all peptide options not just the first
                    var variants = variantPWSM.Peptide.Protein.AppliedSequenceVariations;
                    var culture = CultureInfo.CurrentCulture;
                    // these bools allow for us to accurrately count the number of peptides that have at least one variants of a given type.
                    // they will prevent double counting if a variant type is found more than once in a given peptide (most typically missense, but all will be covered)
                    bool SNVmissenseIdentified = false;
                    bool MNVmissenseIdentified = false;
                    bool insertionIdentified = false;
                    bool deletionIdentified = false;
                    bool frameshiftIdentified = false;
                    bool stopGainIdentified = false;
                    bool stopLossIdentifed = false;

                    foreach (var variant in variants)
                    {
                        if (variantPWSM.Peptide.IntersectsAndIdentifiesVariation(variant).identifies == true)
                        {
                            if (culture.CompareInfo.IndexOf(variant.Description.Description, "missense_variant", CompareOptions.IgnoreCase) >= 0)
                            {
                                if (variant.Description.ReferenceAlleleString.Length == 1 && variant.Description.AlternateAlleleString.Length == 1)
                                {
                                    if (SNVmissenseIdentified == false)
                                    {
                                        SNVmissenseCount++;
                                        SNVmissenseIdentified = true;
                                    }
                                    if (SNVmissenseVariants.ContainsKey(variantPWSM.Peptide.Protein))
                                    {
                                        SNVmissenseVariants[variantPWSM.Peptide.Protein].Add(variant);
                                    }
                                    else
                                    {
                                        SNVmissenseVariants.Add(variantPWSM.Peptide.Protein, new HashSet<SequenceVariation> { variant });
                                    }
                                }
                                else
                                {
                                    if (MNVmissenseIdentified == false)
                                    {
                                        MNVmissenseCount++;
                                        MNVmissenseIdentified = true;
                                    }
                                    if (MNVmissenseVariants.ContainsKey(variantPWSM.Peptide.Protein))
                                    {
                                        MNVmissenseVariants[variantPWSM.Peptide.Protein].Add(variant);
                                    }
                                    else
                                    {
                                        MNVmissenseVariants.Add(variantPWSM.Peptide.Protein, new HashSet<SequenceVariation> { variant });
                                    }
                                }
                            }
                            else if (culture.CompareInfo.IndexOf(variant.Description.Description, "frameshift_variant", CompareOptions.IgnoreCase) >= 0)
                            {
                                if (frameshiftIdentified == false)
                                {
                                    frameshiftCount++;
                                    frameshiftIdentified = true;
                                }
                                if (frameshiftVariants.ContainsKey(variantPWSM.Peptide.Protein))
                                {
                                    frameshiftVariants[variantPWSM.Peptide.Protein].Add(variant);
                                }
                                else
                                {
                                    frameshiftVariants.Add(variantPWSM.Peptide.Protein, new HashSet<SequenceVariation> { variant });
                                }
                            }
                            else if (culture.CompareInfo.IndexOf(variant.Description.Description, "stop_gained", CompareOptions.IgnoreCase) >= 0)
                            {
                                if (stopGainIdentified == false)
                                {
                                    stopGainCount++;
                                    stopGainIdentified = true;
                                }
                                if (stopGainVariants.ContainsKey(variantPWSM.Peptide.Protein))
                                {
                                    stopGainVariants[variantPWSM.Peptide.Protein].Add(variant);
                                }
                                else
                                {
                                    stopGainVariants.Add(variantPWSM.Peptide.Protein, new HashSet<SequenceVariation> { variant });
                                }
                            }
                            else if ((culture.CompareInfo.IndexOf(variant.Description.Description, "conservative_inframe_insertion", CompareOptions.IgnoreCase) >= 0) || (culture.CompareInfo.IndexOf(variant.Description.Description, "disruptive_inframe_insertion", CompareOptions.IgnoreCase) >= 0))
                            {
                                if (insertionIdentified == false)
                                {
                                    insertionCount++;
                                    insertionIdentified = true;
                                }
                                if (insertionVariants.ContainsKey(variantPWSM.Peptide.Protein))
                                {
                                    insertionVariants[variantPWSM.Peptide.Protein].Add(variant);
                                }
                                else
                                {
                                    insertionVariants.Add(variantPWSM.Peptide.Protein, new HashSet<SequenceVariation> { variant });
                                }
                            }
                            else if ((culture.CompareInfo.IndexOf(variant.Description.Description, "conservative_inframe_deletion", CompareOptions.IgnoreCase) >= 0) || (culture.CompareInfo.IndexOf(variant.Description.Description, "disruptive_inframe_deletion", CompareOptions.IgnoreCase) >= 0))
                            {
                                if (deletionIdentified == false)
                                {
                                    deletionCount++;
                                    deletionIdentified = true;
                                }
                                if (deletionVariants.ContainsKey(variantPWSM.Peptide.Protein))
                                {
                                    deletionVariants[variantPWSM.Peptide.Protein].Add(variant);
                                }
                                else
                                {
                                    deletionVariants.Add(variantPWSM.Peptide.Protein, new HashSet<SequenceVariation> { variant });
                                }
                            }
                            else if (culture.CompareInfo.IndexOf(variant.Description.Description, "stop_loss", CompareOptions.IgnoreCase) >= 0)
                            {
                                if (stopLossIdentifed == false)
                                {
                                    stopLossCount++;
                                    stopLossIdentifed = true;
                                }
                                if (stopLossVariants.ContainsKey(variantPWSM.Peptide.Protein))
                                {
                                    stopLossVariants[variantPWSM.Peptide.Protein].Add(variant);
                                }
                                else
                                {
                                    stopLossVariants.Add(variantPWSM.Peptide.Protein, new HashSet<SequenceVariation> { variant });
                                }
                            }
                        }
                    }
                }
            }
            int SNVmissenseSites = 0;
            foreach (var entry in SNVmissenseVariants)
            {
                SNVmissenseSites += entry.Value.Count;
            }

            int MNVmissenseSites = 0;
            foreach (var entry in MNVmissenseVariants)
            {
                MNVmissenseSites += entry.Value.Count;
            }

            int insertionSites = 0;
            foreach (var entry in insertionVariants)
            {
                insertionSites += entry.Value.Count;
            }

            int deletionSites = 0;
            foreach (var entry in deletionVariants)
            {
                deletionSites += entry.Value.Count;
            }

            int frameshiftSites = 0;
            foreach (var entry in frameshiftVariants)
            {
                frameshiftSites += entry.Value.Count;
            }

            int stopGainSites = 0;
            foreach (var entry in stopGainVariants)
            {
                stopGainSites += entry.Value.Count;
            }

            int stopLossSites = 0;
            foreach (var entry in stopLossVariants)
            {
                stopLossSites += entry.Value.Count;
            }

            int totalVariantSites = SNVmissenseSites + MNVmissenseSites + insertionSites + deletionSites + frameshiftSites + stopGainSites + stopLossSites;

            string[] variantResults = new string[25];
            variantResults[0] = "Variant Result Summary";
            variantResults[2] = "--------------------------------------------------";
            variantResults[4] = "Number of potential variant containing peptides identified at " + FdrFilterValue * 100 + "% group FDR: " + variantPeptides.Where(p => !p.IsDecoy && !p.IsContaminant && p.FdrInfo.QValue <= FdrFilterValue && p.FdrInfo.QValueNotch <= FdrFilterValue).ToList().Count;
            variantResults[5] = "Number of unqiuely identified variant peptides at " + FdrFilterValue * 100 + "% group FDR: " + confidentVariantPeps.Where(p => !p.IsDecoy && !p.IsContaminant && p.FdrInfo.QValue <= FdrFilterValue && p.FdrInfo.QValueNotch <= FdrFilterValue).ToList().Count;
            variantResults[6] = "Number of unique variants: " + totalVariantSites;
            variantResults[7] = "Number of SNV missense variant containing peptides at " + FdrFilterValue * 100 + "% group FDR: " + SNVmissenseCount;
            variantResults[8] = "Number of unique SNV missense variants: " + SNVmissenseSites;
            variantResults[9] = "Number of MNV missense variant containing peptides at " + FdrFilterValue * 100 + "% group FDR: " + MNVmissenseCount;
            variantResults[10] = "Number of unique MNV missense variants: " + MNVmissenseSites;
            variantResults[11] = "Number of frameshift variant containing peptides at " + FdrFilterValue * 100 + "% group FDR: " + frameshiftCount;
            variantResults[12] = "Number of unique frameshift variants: " + frameshiftSites;
            variantResults[13] = "Number of inframe insertion variant containing peptides at " + FdrFilterValue * 100 + "% group FDR: " + insertionCount;
            variantResults[14] = "Number of unique inframe insertion variants: " + insertionSites;
            variantResults[15] = "Number of inframe deletion variant containing peptides at " + FdrFilterValue * 100 + "% group FDR: " + deletionCount;
            variantResults[16] = "Number of unique inframe deletion variants: " + deletionSites;
            variantResults[17] = "Number of stop gain variant containing peptides at " + FdrFilterValue * 100 + "% group FDR: " + stopGainCount;
            variantResults[18] = "Number of unique stop gain variants: " + stopGainSites;
            variantResults[19] = "Number of stop loss variant containing peptides at " + FdrFilterValue * 100 + "% group FDR: " + stopLossCount;
            variantResults[20] = "Number of unique stop loss variants: " + stopLossSites;
            variantResults[21] = "Number of variant peptides at " + FdrFilterValue * 100 + "% group FDR with unambiguous localized modifications: " + modifiedVariantPeptides.Count;
            variantResults[22] = "Number of variant peptides at " + FdrFilterValue * 100 + "% group FDR with unambiguous localized modifications at the variant sites : " + modifiedVariantSitePeptides.Count;

            string filePath = Path.Combine(Parameters.OutputFolder, "VariantAnalysisResultSummary.txt");
            File.WriteAllLines(filePath, variantResults);
        }

        private static int GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications peptideWithSetModifications)
        {
            if (oneIsNterminus == 1)
            {
                return peptideWithSetModifications.OneBasedStartResidueInProtein;
            }
            if (oneIsNterminus == peptideWithSetModifications.Length + 2)
            {
                return peptideWithSetModifications.OneBasedEndResidueInProtein;
            }
            return peptideWithSetModifications.OneBasedStartResidueInProtein + oneIsNterminus - 2;
        }

        private static void WriteTree(BinTreeStructure myTreeStructure, string writtenFile)
        {
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tUnimodDiffs\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tprotNtermLocFrac\tpepNtermLocFrac\tpepCtermLocFrac\tprotCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot");
                foreach (Bin bin in myTreeStructure.FinalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString("F4", CultureInfo.InvariantCulture)
                        + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.Count == 0 ? double.NaN : (double)bin.CountDecoy / bin.Count).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : (double)bin.LocalizeableTarget / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.Mine
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.UnimodDiffs
                        + "\t" + bin.AA
                        + "\t" + bin.Combos
                        + "\t" + string.Join(",", bin.ModsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.ResidueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.FracWithSingle).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + ((double)bin.Overlapping / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.MedianLength).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.UniprotID);
                }
            }
        }

        private void WritePsmsForPercolator(List<PeptideSpectralMatch> psmList, string writtenFileForPercolator)
        {
            using (StreamWriter output = new StreamWriter(writtenFileForPercolator))
            {
                string searchType;
                if (psmList.Where(p => p != null).Any() && psmList[0].DigestionParams.Protease.Name != null && psmList[0].DigestionParams.Protease.Name == "top-down")
                {
                    searchType = "top-down";
                }
                else
                {
                    searchType = "standard";
                }

                string header = "SpecId\tLabel\tScanNr\t";
                header = header + String.Join("\t", PsmData.trainingInfos[searchType]);
                header = header + "\tPeptide\tProteins";

                output.WriteLine(header);

                StringBuilder directions = new StringBuilder();
                directions.Append("DefaultDirection\t-\t-");

                foreach (var headerVariable in PsmData.trainingInfos[searchType])
                {
                    directions.Append("\t");
                    directions.Append(PsmData.assumedAttributeDirection[headerVariable]);
                }

                output.WriteLine(directions.ToString());

                int idNumber = 0;
                psmList.OrderByDescending(p => p.Score);
                foreach (PeptideSpectralMatch psm in psmList.Where(p => p.PsmData_forPEPandPercolator != null))
                {
                    foreach (var peptide in psm.BestMatchingPeptides)
                    {
                        output.Write(idNumber.ToString());

                        output.Write('\t' + (peptide.Peptide.Protein.IsDecoy ? -1 : 1).ToString());
                        output.Write('\t' + psm.ScanNumber.ToString());
                        output.Write(psm.PsmData_forPEPandPercolator.ToString(searchType));

                        output.Write('\t' + (peptide.Peptide.PreviousAminoAcid + "." + peptide.Peptide.FullSequence + "." + peptide.Peptide.NextAminoAcid).ToString());
                        output.Write('\t' + (peptide.Peptide.Protein.Accession).ToString());
                        output.WriteLine();
                    }
                    idNumber++;
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
                        if ((!Parameters.SearchParameters.WriteDecoys && proteinGroups[i].IsDecoy) ||
                            (!Parameters.SearchParameters.WriteContaminants && proteinGroups[i].IsContaminant))
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