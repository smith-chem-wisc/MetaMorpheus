using EngineLayer;
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
using System.Linq;
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
            WritePeptideResults(); // modifies the FDR results for PSMs, so do this last

            return Parameters.SearchTaskResults;
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
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
            new FdrAnalysisEngine(Parameters.AllPsms, Parameters.NumNotches, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run();

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

            //if SILAC
            List<PeptideSpectralMatch> psmsForProteinParsimony = Parameters.AllPsms;
            if (Parameters.SearchParameters.SilacLabels != null)
            {
                psmsForProteinParsimony = SilacConversions.UpdatePsmsForParsimony(Parameters.SearchParameters.SilacLabels, psmsForProteinParsimony);
            }

            // run parsimony
            ProteinParsimonyResults proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(psmsForProteinParsimony, Parameters.SearchParameters.ModPeptidesAreDifferent, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run());

            // score protein groups and calculate FDR
            ProteinScoringAndFdrResults proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psmsForProteinParsimony,
                Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent, true, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run();

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
                        myMsDataFile, combinedParams, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile }).Run();
                    Parameters.MyFileManager.DoneWithFile(origDataFile);
                    ReportProgress(new ProgressEventArgs(100, "Done with localization analysis!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", origDataFile }));
                }
            }

            // count different modifications observed
            new ModificationAnalysisEngine(Parameters.AllPsms, CommonParameters, new List<string> { Parameters.SearchTaskId }).Run();
        }

        private void QuantificationAnalysis()
        {
            if (!Parameters.SearchParameters.DoQuantification)
            {
                return;
            }

            // pass quantification parameters to FlashLFQ
            Status("Quantifying...", Parameters.SearchTaskId);

            // construct file info for FlashLFQ
            var spectraFileInfo = new List<SpectraFileInfo>();

            // get experimental design info for normalization
            if (Parameters.SearchParameters.Normalize)
            {
                string assumedExperimentalDesignPath = Directory.GetParent(Parameters.CurrentRawFileList.First()).FullName;
                assumedExperimentalDesignPath = Path.Combine(assumedExperimentalDesignPath, GlobalVariables.ExperimentalDesignFileName);

                if (File.Exists(assumedExperimentalDesignPath))
                {
                    var experimentalDesign = File.ReadAllLines(assumedExperimentalDesignPath)
                        .ToDictionary(p => p.Split('\t')[0], p => p);

                    foreach (var file in Parameters.CurrentRawFileList)
                    {
                        string filename = Path.GetFileNameWithoutExtension(file);

                        var expDesignForThisFile = experimentalDesign[filename];
                        var split = expDesignForThisFile.Split('\t');

                        string condition = split[1];
                        int biorep = int.Parse(split[2]);
                        int fraction = int.Parse(split[3]);
                        int techrep = int.Parse(split[4]);

                        // experimental design info passed in here for each spectra file
                        spectraFileInfo.Add(new SpectraFileInfo(fullFilePathWithExtension: file,
                            condition: condition,
                            biorep: biorep - 1,
                            fraction: fraction - 1,
                            techrep: techrep - 1));

                        Parameters.MyFileManager.DoneWithFile(file);
                    }
                }
                else
                {
                    throw new MetaMorpheusException("Could not find experimental design file at location:\n" + assumedExperimentalDesignPath);
                }
            }
            else
            {
                foreach (var file in Parameters.CurrentRawFileList)
                {
                    // experimental design info passed in here for each spectra file
                    spectraFileInfo.Add(new SpectraFileInfo(fullFilePathWithExtension: file, condition: "", biorep: 0, fraction: 0, techrep: 0));
                    Parameters.MyFileManager.DoneWithFile(file);
                }
            }

            // get PSMs to pass to FlashLFQ
            var unambiguousPsmsBelowOnePercentFdr = Parameters.AllPsms.Where(p =>
                p.FdrInfo.QValue <= 0.01
                && p.FdrInfo.QValueNotch <= 0.01
                && !p.IsDecoy
                && p.FullSequence != null).ToList(); //if ambiguous, there's no full sequence

            //If SILAC (Pre-Quantification) preparation
            Dictionary<PeptideSpectralMatch, List<PeptideSpectralMatch>> silacPsmMatcher = new Dictionary<PeptideSpectralMatch, List<PeptideSpectralMatch>>(); //use for getting back info once we've separated light/heavy. It is the light psm for the key, and the list of heavy psm replicas for the value.
            Dictionary<string, List<string>> silacProteinGroupMatcher = new Dictionary<string, List<string>>(); //use for getting back info once we've separated light/heavy. It is the light protein group for the key, and a list of heavy protein group replicas for the value.
            if (Parameters.SearchParameters.SilacLabels != null)
            {
                //PEPTIDE LEVEL CHANGES
                //go through all the psms and duplicate them until a psm copy exists for the unlabeled and labeled proteins
                //The number of psms should roughly increase by a factor of N, where N is the number of labels.
                //It may not increase exactly by a factor of N if the amino acid that gets labeled doesn't exist in the peptide
                List<SilacLabel> allSilacLabels = Parameters.SearchParameters.SilacLabels;
                List<PeptideSpectralMatch> silacPsms = new List<PeptideSpectralMatch>();
                foreach (PeptideSpectralMatch psm in unambiguousPsmsBelowOnePercentFdr)
                {
                    //see which label, if any, this peptide has
                    string peptideBaseSequence = psm.BaseSequence;
                    SilacLabel observedLabel = SilacConversions.GetRelevantLabelFromBaseSequence(peptideBaseSequence, allSilacLabels); //returns null if no label

                    //if it's not the light form, make a light form and add it
                    PeptideSpectralMatch lightPsm = observedLabel == null ? psm : SilacConversions.GetSilacPsm(psm, observedLabel, true);
                    silacPsms.Add(lightPsm);

                    //see which labels, if any, this peptide can have
                    peptideBaseSequence = SilacConversions.GetSilacLightBaseSequence(peptideBaseSequence, observedLabel); //update to the light sequence, if applicable
                    List<SilacLabel> labelsToModify = allSilacLabels.Where(x => peptideBaseSequence.Contains(x.OriginalAminoAcid)).ToList(); //see which residues can be labeled
                    //foreach of the possible labels this peptide can have, add a new psm for it
                    List<PeptideSpectralMatch> addedSilacPsms = labelsToModify.Select(x => SilacConversions.GetSilacPsm(lightPsm, x, false)).ToList(); //if we had a heavy psm to start with, this will create an identical heavy
                    addedSilacPsms.ForEach(x => silacPsms.Add(x)); //add to psm list
                    silacPsmMatcher[lightPsm] = addedSilacPsms; //add to match list
                }

                //update the list for FlashLFQ
                silacPsms.ForEach(x => x.ResolveAllAmbiguities()); //update the base/full sequences
                unambiguousPsmsBelowOnePercentFdr = silacPsms;


                //PROTEIN LEVEL CHANGES
                //Need to create the heavy proteins in our parsimony list and assign the heavy peptides to these heavy proteins
                //During the search and protein parsimony, heavy and light peptides were considered to come from the same protein (light)
                //Now, we want to quantify the light and heavy proteins separately and we need to differentiate them.
                if (ProteinGroups != null) //if we did parsimony
                {
                    List<EngineLayer.ProteinGroup> silacProteinGroups = new List<EngineLayer.ProteinGroup>();
                    foreach (EngineLayer.ProteinGroup proteinGroup in ProteinGroups)
                    {
                        //add the unlabeled protein group
                        //this method removes the heavy psms from the protein groups and adds their light replica
                        EngineLayer.ProteinGroup unlabeledProteinGroup = SilacConversions.GetSilacProteinGroups(unambiguousPsmsBelowOnePercentFdr, proteinGroup);
                        silacProteinGroups.Add(unlabeledProteinGroup);
                        //add the labeled protein group(s)
                        List<EngineLayer.ProteinGroup> addedProteinGroups = allSilacLabels.Select(x => SilacConversions.GetSilacProteinGroups(unambiguousPsmsBelowOnePercentFdr, proteinGroup, x)).ToList(); //foreach label, create a new heavy protein group
                        addedProteinGroups.ForEach(x => silacProteinGroups.Add(x)); //add to psm list
                        silacProteinGroupMatcher[proteinGroup.ProteinGroupName] = addedProteinGroups.Select(x => x.ProteinGroupName).ToList(); //add to match list
                    }
                    ProteinGroups = silacProteinGroups;
                }
            } //END SILAC CODE

            // pass protein group info for each PSM
            var psmToProteinGroups = new Dictionary<PeptideSpectralMatch, List<FlashLFQ.ProteinGroup>>();
            if (ProteinGroups != null)
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
                normalize: Parameters.SearchParameters.Normalize,
                ppmTolerance: Parameters.SearchParameters.QuantifyPpmTol,
                matchBetweenRuns: Parameters.SearchParameters.MatchBetweenRuns,
                silent: true,
                optionalPeriodicTablePath: GlobalVariables.ElementsLocation,
                maxThreads: CommonParameters.MaxThreadsToUsePerFile);

            if (flashLFQIdentifications.Any())
            {
                Parameters.FlashLfqResults = FlashLfqEngine.Run();
            }

            //MultiProtease MBR capability code
            //Parameters.FlashLfqResults = null;

            //foreach (var proteasePsms in proteaseSortedPsms)
            //{
            //    var flashLFQIdentifications = new List<Identification>();
            //    var proteasePsmsGroupedByFile = proteasePsms.Value.GroupBy(p => p.FullFilePath);
            //    foreach (var spectraFile in proteasePsmsGroupedByFile)
            //    {
            //        var rawfileinfo = spectraFileInfo.Where(p => p.FullFilePathWithExtension.Equals(spectraFile.Key)).First();

            //        foreach (var psm in spectraFile)
            //        {
            //            flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence,
            //                psm.PeptideMonisotopicMass.Value, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroups[psm]));
            //        }
            //    }

            //    // run FlashLFQ
            //    var FlashLfqEngine = new FlashLFQEngine(
            //        allIdentifications: flashLFQIdentifications,
            //        normalize: Parameters.SearchParameters.Normalize,
            //        ppmTolerance: Parameters.SearchParameters.QuantifyPpmTol,
            //        matchBetweenRuns: Parameters.SearchParameters.MatchBetweenRuns,
            //        silent: true,
            //        optionalPeriodicTablePath: GlobalVariables.ElementsLocation);

            //    if (flashLFQIdentifications.Any())
            //    {
            //        //make specific to protease
            //        var results = FlashLfqEngine.Run();

            //        if (Parameters.FlashLfqResults == null)
            //        {
            //            Parameters.FlashLfqResults = results;
            //        }
            //        else
            //        {
            //            Parameters.FlashLfqResults.MergeResultsWith(results);
            //        }
            //    }
            //}

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

            if (Parameters.SearchParameters.SilacLabels != null)
            {
                SilacConversions.SilacConversionsPostQuantification(Parameters.SearchParameters.SilacLabels, spectraFileInfo, ProteinGroups,
                    Parameters.ListOfDigestionParams, silacProteinGroupMatcher, Parameters.FlashLfqResults, Parameters.AllPsms,
                    Parameters.SearchParameters.ModsToWriteSelection, FlashLfqEngine.Integrate);
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
            writtenFile = Path.Combine(Parameters.OutputFolder, "AllPSMs_FormattedForPercolator.tsv");
            WritePsmsForPercolator(FilteredPsmListForOutput, writtenFile, CommonParameters.QValueOutputFilter);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            // write summary text
            Parameters.SearchTaskResults.AddNiceText("All target PSMS within 1% FDR: " + Parameters.AllPsms.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy));
            if (Parameters.SearchParameters.DoParsimony)
            {
                Parameters.SearchTaskResults.AddNiceText("All target protein groups within 1% FDR: " + ProteinGroups.Count(b => b.QValue <= 0.01 && !b.IsDecoy)
                    + Environment.NewLine);
            }

            PsmsGroupedByFile = FilteredPsmListForOutput.GroupBy(p => p.FullFilePath);

            foreach (var file in PsmsGroupedByFile)
            {
                // write summary text
                var psmsForThisFile = file.ToList();
                string strippedFileName = Path.GetFileNameWithoutExtension(file.First().FullFilePath);

                Parameters.SearchTaskResults.AddNiceText("MS2 spectra in " + strippedFileName + ": " + Parameters.NumMs2SpectraPerFile[strippedFileName][0]);
                Parameters.SearchTaskResults.AddNiceText("Precursors fragmented in " + strippedFileName + ": " + Parameters.NumMs2SpectraPerFile[strippedFileName][1]);
                Parameters.SearchTaskResults.AddNiceText("Target PSMs within 1% FDR in " + strippedFileName + ": " + psmsForThisFile.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy));

                // writes all individual spectra file search results to subdirectory
                if (Parameters.CurrentRawFileList.Count > 1)
                {
                    // create individual files subdirectory
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);

                    // write PSMs
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_PSMs.psmtsv");
                    WritePsmsToTsv(psmsForThisFile, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.First().FullFilePath });

                    // write PSMs for percolator
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_PSMsFormattedForPercolator.tsv");
                    WritePsmsForPercolator(psmsForThisFile, writtenFile, CommonParameters.QValueOutputFilter);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.First().FullFilePath });
                }
            }
        }

        private void WriteProteinResults()
        {
            if (Parameters.SearchParameters.DoParsimony)
            {
                // write protein groups to tsv
                string writtenFile = Path.Combine(Parameters.OutputFolder, "AllProteinGroups.tsv");
                WriteProteinGroupsToTsv(ProteinGroups, writtenFile, new List<string> { Parameters.SearchTaskId }, CommonParameters.QValueOutputFilter);

                // write all individual file results to subdirectory
                // local protein fdr, global parsimony, global psm fdr
                if (Parameters.CurrentRawFileList.Count > 1 || Parameters.SearchParameters.WriteMzId || Parameters.SearchParameters.WritePepXml)
                {
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);

                    //If SILAC and no light, we need to update the psms (which were found in the "light" file) and say they were found in the "heavy" file)
                    if (!Parameters.ListOfDigestionParams.Any(x => x.GeneratehUnlabeledProteinsForSilac))
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

                    foreach (var fullFilePath in PsmsGroupedByFile.Select(v => v.Key))
                    {
                        string strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                        List<PeptideSpectralMatch> psmsForThisFile = PsmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList();
                        var subsetProteinGroupsForThisFile = ProteinGroups.Select(p => p.ConstructSubsetProteinGroup(fullFilePath)).ToList();

                        ProteinScoringAndFdrResults subsetProteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(subsetProteinGroupsForThisFile, psmsForThisFile,
                            Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent,
                            false, CommonParameters, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }).Run();

                        subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults.SortedAndScoredProteinGroups;

                        Parameters.SearchTaskResults.AddNiceText("Target protein groups within 1 % FDR in " + strippedFileName + ": " + subsetProteinGroupsForThisFile.Count(b => b.QValue <= 0.01 && !b.IsDecoy));

                        // write individual spectra file protein groups results to tsv
                        if (Parameters.CurrentRawFileList.Count > 1)
                        {
                            writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_ProteinGroups.tsv");
                            WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }, CommonParameters.QValueOutputFilter);
                        }

                        // write mzID
                        if (Parameters.SearchParameters.WriteMzId)
                        {
                            Status("Writing mzID...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });

                            var mzidFilePath = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + ".mzID");
                            MzIdentMLWriter.WriteMzIdentMl(psmsForThisFile, subsetProteinGroupsForThisFile, Parameters.VariableModifications, Parameters.FixedModifications, Parameters.SearchParameters.SilacLabels,
                                new List<Protease> { CommonParameters.DigestionParams.Protease }, CommonParameters.QValueOutputFilter, CommonParameters.ProductMassTolerance,
                                CommonParameters.PrecursorMassTolerance, CommonParameters.DigestionParams.MaxMissedCleavages, mzidFilePath);

                            FinishedWritingFile(mzidFilePath, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                        }

                        // write pepXML
                        if (Parameters.SearchParameters.WritePepXml)
                        {
                            Status("Writing pepXML...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });

                            var pepXMLFilePath = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + ".pep.XML");
                            PepXMLWriter.WritePepXml(psmsForThisFile, Parameters.DatabaseFilenameList, Parameters.VariableModifications, Parameters.FixedModifications,
                                CommonParameters, pepXMLFilePath, CommonParameters.QValueOutputFilter);

                            FinishedWritingFile(pepXMLFilePath, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                        }

                        ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }));
                    }
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
                WritePeptideQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "AllQuantifiedPeptides", new List<string> { Parameters.SearchTaskId });

                // write individual results
                if (Parameters.CurrentRawFileList.Count > 1)
                {
                    foreach (var file in Parameters.FlashLfqResults.Peaks)
                    {
                        WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.IndividualResultsOutputFolder,
                            file.Key.FilenameWithoutExtension + "_QuantifiedPeaks", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.Key.FullFilePathWithExtension });
                    }
                }
            }
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

                        if (proteinToConfidentBaseSequences.ContainsKey(nonVariantProtein.NonVariantProtein))
                        {
                            // adds confidently localized and identified mods
                            nonVariantProtein.OneBasedPossibleLocalizedModifications.Clear();
                            foreach (var kvp in modsToWrite.Where(kv => kv.Key.Item1 == null))
                            {
                                nonVariantProtein.OneBasedPossibleLocalizedModifications.Add(kvp.Key.Item2, kvp.Value);
                            }
                            foreach (var sv in nonVariantProtein.SequenceVariations)
                            {
                                sv.OneBasedModifications.Clear();
                                foreach (var kvp in modsToWrite.Where(kv => kv.Key.Item1 != null && kv.Key.Item1.Equals(sv)))
                                {
                                    sv.OneBasedModifications.Add(kvp.Key.Item2, kvp.Value);
                                }
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
            }
        }

        private void WritePeptideResults()
        {
            Status("Writing peptide results...", Parameters.SearchTaskId);

            // write best (highest-scoring) PSM per peptide
            string writtenFile = Path.Combine(Parameters.OutputFolder, "AllPeptides.psmtsv");
            List<PeptideSpectralMatch> peptides = Parameters.AllPsms.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
            
            new FdrAnalysisEngine(peptides, Parameters.NumNotches, CommonParameters, new List<string> { Parameters.SearchTaskId }, "Peptide").Run();

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

            Parameters.SearchTaskResults.AddNiceText("All target peptides within 1% FDR: " + peptides.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy));

            foreach (var file in PsmsGroupedByFile)
            {
                // write summary text
                var psmsForThisFile = file.ToList();
                string strippedFileName = Path.GetFileNameWithoutExtension(file.First().FullFilePath);
                var peptidesForFile = psmsForThisFile.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
                new FdrAnalysisEngine(peptidesForFile, Parameters.NumNotches, CommonParameters, new List<string> { Parameters.SearchTaskId }, "Peptide").Run();
                Parameters.SearchTaskResults.AddNiceText("Target peptides within 1% FDR in " + strippedFileName + ": " + peptidesForFile.Count(a => a.FdrInfo.QValue <= 0.01 && !a.IsDecoy) + Environment.NewLine);

                // writes all individual spectra file search results to subdirectory
                if (Parameters.CurrentRawFileList.Count > 1)
                {
                    // create individual files subdirectory
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);

                    // write best (highest-scoring) PSM per peptide
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_Peptides.psmtsv");
                    WritePsmsToTsv(peptidesForFile, writtenFile, Parameters.SearchParameters.ModsToWriteSelection);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", file.First().FullFilePath });
                }
            }
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

        private void WritePsmsForPercolator(List<PeptideSpectralMatch> psmList, string writtenFileForPercolator, double qValueCutoff)
        {
            using (StreamWriter output = new StreamWriter(writtenFileForPercolator))
            {
                output.WriteLine("SpecId\tLabel\tScanNr\tF1\tF2\tPeptide\tProteins");
                output.WriteLine("DefaultDirection\t-\t-\t1\t1\t\t");
                for (int i = 0; i < psmList.Count; i++)
                {
                    var psm = psmList[i];

                    if (psm.FdrInfo.QValue > qValueCutoff || psm.FdrInfo.QValueNotch > qValueCutoff)
                    {
                        continue;
                    }

                    output.Write(i.ToString());
                    output.Write('\t' + (psm.IsDecoy ? -1 : 1).ToString());
                    output.Write('\t' + psm.ScanNumber.ToString());

                    // Features
                    output.Write('\t' + string.Join("\t", psm.Features));

                    // HACKY: Ignores all ambiguity
                    var pwsm = psm.BestMatchingPeptides.First().Peptide;

                    output.Write('\t' + (pwsm.PreviousAminoAcid + "." + pwsm.FullSequence + "." + pwsm.NextAminoAcid).ToString());
                    output.Write('\t' + (pwsm.Protein.Accession).ToString());
                    output.WriteLine();
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

            flashLFQResults.WriteResults(null, fullSeqPath, null);

            FinishedWritingFile(fullSeqPath, nestedIds);
        }

        private void WritePeakQuantificationResultsToTsv(FlashLfqResults flashLFQResults, string outputFolder, string fileName, List<string> nestedIds)
        {
            var peaksPath = Path.Combine(outputFolder, fileName + ".tsv");

            flashLFQResults.WriteResults(peaksPath, null, null);

            FinishedWritingFile(peaksPath, nestedIds);
        }
    }
}