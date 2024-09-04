using Easy.Common.Extensions;
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
using System.IO.Compression;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;
using TaskLayer.MbrAnalysis;
using Chemistry;
using MzLibUtil;
using Omics.Modifications;
using Omics.SpectrumMatch;

namespace TaskLayer
{
    public class PostSearchAnalysisTask : MetaMorpheusTask
    {
        public PostSearchAnalysisParameters Parameters { get; set; }
        private List<EngineLayer.ProteinGroup> ProteinGroups { get; set; }
        private SpectralRecoveryResults SpectralRecoveryResults { get; set; }

        /// <summary>
        /// Used for storage of results for writing to Results.tsv. It is explained in the method ConstructResultsDictionary()
        /// </summary>
        private Dictionary<(string,string),string> ResultsDictionary { get; set; }

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
                // The GroupBy statement below gets rid of duplicate PSMs that can occur when the same peptides is matched to the same spectrum multiple times,
                // just with slightly different precursor masses.
                Parameters.AllPsms = Parameters.AllPsms.OrderByDescending(b => b)
                    .GroupBy(b => (b.FullFilePath, b.ScanNumber, b.BioPolymerWithSetModsMonoisotopicMass)).Select(b => b.First()).ToList();
                CalculatePsmAndPeptideFdr(Parameters.AllPsms);
            }
            ConstructResultsDictionary();
            DoMassDifferenceLocalizationAnalysis();
            ProteinAnalysis();
            QuantificationAnalysis();

            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files" }));

            HistogramAnalysis();
            WritePsmResults();
            WritePeptideResults();
            if (Parameters.CurrentRawFileList.Count > 1 && (Parameters.SearchParameters.WriteIndividualFiles
                                                            || Parameters.SearchParameters.WriteMzId ||
                                                            Parameters.SearchParameters.WritePepXml))
            {
                // create individual files subdirectory
                Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);
                if (Parameters.SearchParameters.WriteIndividualFiles)
                {
                    WriteIndividualPsmResults();
                    WriteIndividualPeptideResults();
                }
            }
            WriteProteinResults();
            AddResultsTotalsToAllResultsTsv();
            WritePrunedDatabase();
            var k = CommonParameters;
            if (Parameters.SearchParameters.WriteSpectralLibrary)
            {
                SpectralLibraryGeneration();
                if (Parameters.SearchParameters.DoLabelFreeQuantification && Parameters.FlashLfqResults != null)
                {
                    SpectralRecoveryResults = SpectralRecoveryRunner.RunSpectralRecoveryAlgorithm(Parameters, CommonParameters, FileSpecificParameters);
                }      
            }

            if(Parameters.SearchParameters.UpdateSpectralLibrary)
            {
                UpdateSpectralLibrary();
            }
          
            WriteFlashLFQResults();

            if (Parameters.ProteinList.Any((p => p.AppliedSequenceVariations.Count > 0)))
            {
                WriteVariantResults();
            }

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
        private void CalculatePsmAndPeptideFdr(List<SpectralMatch> psms, string analysisType = "PSM", bool doPep = true)
        {
            // TODO: because FDR is done before parsimony, if a PSM matches to a target and a decoy protein, there may be conflicts between how it's handled in parsimony and the FDR engine here
            // for example, here it may be treated as a decoy PSM, where as in parsimony it will be determined by the parsimony algorithm which is agnostic of target/decoy assignments
            // this could cause weird PSM FDR issues

            Status("Estimating PSM FDR...", Parameters.SearchTaskId);
            new FdrAnalysisEngine(psms, Parameters.NumNotches, CommonParameters, this.FileSpecificParameters, 
                new List<string> { Parameters.SearchTaskId }, analysisType: analysisType, doPEP: doPep, outputFolder: Parameters.OutputFolder).Run();

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
            if (Parameters.SearchParameters.SilacLabels != null && Parameters.AllPsms.First() is PeptideSpectralMatch)
            {
                Parameters.AllPsms = SilacConversions.UpdateProteinSequencesToLight(Parameters.AllPsms.Cast<PeptideSpectralMatch>().ToList(), Parameters.SearchParameters.SilacLabels);
                foreach (SpectralMatch psm in Parameters.AllPsms)
                {
                    psm.ResolveAllAmbiguities();
                }
            }

            var psmForParsimony = FilteredPsms.Filter(Parameters.AllPsms,
                commonParams: CommonParameters,
                includeDecoys: true,
                includeContaminants: true,
                includeAmbiguous: false,
                includeHighQValuePsms: false);

            // run parsimony
            ProteinParsimonyResults proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(psmForParsimony.FilteredPsmsList, Parameters.SearchParameters.ModPeptidesAreDifferent, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run());

            // score protein groups and calculate FDR
            ProteinScoringAndFdrResults proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, psmForParsimony.FilteredPsmsList,
                Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent, true, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId }).Run();

            ProteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

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
            if (Parameters.SearchParameters.DoMultiplexQuantification)
            {
                List<Modification> multiplexMods = Parameters.FixedModifications.Where(m => m.ModificationType == "Multiplex Label").ToList();
                if (multiplexMods.IsNotNullOrEmpty())
                {
                    Parameters.MultiplexModification = multiplexMods.MaxBy(m => m.DiagnosticIons.Count);
                }
                return;
            }

            if (!Parameters.SearchParameters.DoLabelFreeQuantification)
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
            var psmsForQuantification = FilteredPsms.Filter(Parameters.AllPsms,
                CommonParameters,
                includeDecoys: false,
                includeContaminants: true,
                includeAmbiguous: false,
                includeAmbiguousMods: false,
                includeHighQValuePsms: false);

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
            else
            {
                // if protein groups were not constructed, just use accession numbers
                var accessionToPg = new Dictionary<string, FlashLFQ.ProteinGroup>();
                foreach (var psm in psmsForQuantification)
                {
                    var proteins = psm.BestMatchingBioPolymersWithSetMods.Select(b => b.Peptide.Parent).Distinct();

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
            //Setup variables
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

                List<SpectralMatch> silacPsms = new(); //populate with duplicate psms for heavy/light

                //multiply the psms by the number of labels
                foreach (PeptideSpectralMatch psm in psmsForQuantification)
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
                    int notch = psm.BestMatchingBioPolymersWithSetMods.First().Notch;
                    PeptideWithSetModifications pwsm = psm.BestMatchingBioPolymersWithSetMods.First().Peptide as PeptideWithSetModifications;

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
                        List<int> differentResidues = new();
                        List<char> startUniqueResidues = new();
                        List<char> endUniqueResidues = new();
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
                        List<string> labeledBaseSequences = new() { startLabeledBaseSequence };
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
                psmsForQuantification.SetSilacFilteredPsms(silacPsms);
            }

            //group psms by file
            var psmsGroupedByFile = psmsForQuantification.GroupBy(p => p.FullFilePath);

            // some PSMs may not have protein groups (if 2 peptides are required to construct a protein group, some PSMs will be left over)
            // the peptides should still be quantified but not considered for protein quantification
            var undefinedPg = new FlashLFQ.ProteinGroup("UNDEFINED", "", "");
            //sort the unambiguous psms by protease to make MBR compatible with multiple proteases
            Dictionary<Protease, List<SpectralMatch>> proteaseSortedPsms = new Dictionary<Protease, List<SpectralMatch>>();

            foreach (DigestionParams dp in Parameters.ListOfDigestionParams)
            {
                if (!proteaseSortedPsms.ContainsKey(dp.Protease))
                {
                    proteaseSortedPsms.Add(dp.Protease, new List<SpectralMatch>());
                }
            }
            foreach (var psm in psmsForQuantification)
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
            var flashLfqEngine = new FlashLfqEngine(
                allIdentifications: flashLFQIdentifications,
                normalize: Parameters.SearchParameters.Normalize,
                ppmTolerance: Parameters.SearchParameters.QuantifyPpmTol,
                matchBetweenRunsPpmTolerance: Parameters.SearchParameters.QuantifyPpmTol,  // If these tolerances are not equivalent, then MBR will falsely classify peptides found in the initial search as MBR peaks
                matchBetweenRuns: Parameters.SearchParameters.MatchBetweenRuns,
                useSharedPeptidesForProteinQuant: Parameters.SearchParameters.UseSharedPeptidesForLFQ,
                silent: true,
                maxThreads: CommonParameters.MaxThreadsToUsePerFile);

            if (flashLFQIdentifications.Any())
            {
                Parameters.FlashLfqResults = flashLfqEngine.Run();
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
            if (Parameters.SearchParameters.SilacLabels != null && Parameters.AllPsms.First() is PeptideSpectralMatch) //if we're doing silac
            {
                SilacConversions.SilacConversionsPostQuantification(allSilacLabels, startLabel, endLabel, spectraFileInfo, ProteinGroups, Parameters.ListOfDigestionParams,
                    Parameters.FlashLfqResults, Parameters.AllPsms.Cast<PeptideSpectralMatch>().ToList(), Parameters.SearchParameters.ModsToWriteSelection, quantifyUnlabeledPeptides);
            }
        }

        private void HistogramAnalysis()
        {
            if (Parameters.SearchParameters.DoHistogramAnalysis)
            {
                var limitedpsms_with_fdr = FilteredPsms.Filter(Parameters.AllPsms,
                    commonParams: CommonParameters,
                    includeDecoys: false,
                    includeContaminants: true,
                    includeAmbiguous: false,
                    includeHighQValuePsms: false);

                if (limitedpsms_with_fdr.Any())
                {
                    Status("Running histogram analysis...", new List<string> { Parameters.SearchTaskId });
                    var myTreeStructure = new BinTreeStructure();
                    myTreeStructure.GenerateBins(limitedpsms_with_fdr.FilteredPsmsList, Parameters.SearchParameters.HistogramBinTolInDaltons);
                    var writtenFile = Path.Combine(Parameters.OutputFolder, "MassDifferenceHistogram.tsv");
                    WriteTree(myTreeStructure, writtenFile);
                    FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });
                }
            }
        }

        /// <summary>
        /// Writes PSMs to a .psmtsv file. If multiplex labeling was used (e.g., TMT), the intensities of the diagnostic ions are
        /// included, with each ion being reported in a separate column.
        /// </summary>
        /// <param name="psms">PSMs to be written</param>
        /// <param name="filePath">Full file path, up to and including the filename and extensioh. </param>
        protected void WritePsmsToTsv(IEnumerable<SpectralMatch> psms, string filePath, bool writePeptideLevelResults = false)
        {
            if (Parameters.SearchParameters.DoMultiplexQuantification &&
                Parameters.MultiplexModification != null &&
                psms.Any(p => p.BestMatchingBioPolymersWithSetMods
                    .SelectMany(pwsm => pwsm.Peptide.AllModsOneIsNterminus.Values)
                    .Any(mod => mod.OriginalId.Equals(Parameters.MultiplexModification.OriginalId))))
            {
                WritePsmPlusMultiplexIons(psms, filePath);
            }
            else
            {
                WritePsmsToTsv(psms, filePath, Parameters.SearchParameters.ModsToWriteSelection, writePeptideLevelResults);
            }
        }
        private void WritePsmResults()
        {
            Status("Writing PSM results...", Parameters.SearchTaskId);
            var psmsForPsmResults = FilteredPsms.Filter(Parameters.AllPsms,
                CommonParameters,
                includeDecoys: Parameters.SearchParameters.WriteDecoys,
                includeContaminants: Parameters.SearchParameters.WriteContaminants,
                includeAmbiguous: true,
                includeHighQValuePsms: Parameters.SearchParameters.WriteHighQValuePsms);

            // write PSMs
            string writtenFile = Path.Combine(Parameters.OutputFolder, "AllPSMs.psmtsv");
            WritePsmsToTsv(psmsForPsmResults.OrderByDescending(p=>p).ToList(), writtenFile, writePeptideLevelResults: false);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            // write PSMs for percolator
            // percolator native read format is .tab
            writtenFile = Path.Combine(Parameters.OutputFolder, "AllPSMs_FormattedForPercolator.tab");
            WritePsmsForPercolator(psmsForPsmResults.OrderByDescending(p=>p).ToList(), writtenFile);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            // write summary text
            if (psmsForPsmResults.FilteringNotPerformed)
            {
                
                Parameters.SearchTaskResults.AddPsmPeptideProteinSummaryText(
                    "PEP could not be calculated due to an insufficient number of PSMs. Results were filtered by q-value." +
                    Environment.NewLine);
            }
            string psmResultsText = "All target PSMs with " + psmsForPsmResults.GetFilterTypeString() + " <= " + Math.Round(psmsForPsmResults.FilterThreshold, 2) + ": " +
                psmsForPsmResults.TargetPsmsAboveThreshold;
            ResultsDictionary[("All", "PSMs")] = psmResultsText;
        }
        private void WritePeptideResults()
        {
            Status("Writing peptide results...", Parameters.SearchTaskId);

            var peptidesForPeptideResults = FilteredPsms.Filter(Parameters.AllPsms,
                CommonParameters,
                includeDecoys: Parameters.SearchParameters.WriteDecoys,
                includeContaminants: Parameters.SearchParameters.WriteContaminants,
                includeAmbiguous: true,
                includeHighQValuePsms: Parameters.SearchParameters.WriteHighQValuePsms,
                filterAtPeptideLevel: true);

            // write PSMs
            string writtenFile = Path.Combine(Parameters.OutputFolder, $"All{GlobalVariables.AnalyteType}s.psmtsv");
            WritePsmsToTsv(peptidesForPeptideResults.OrderByDescending(p => p).ToList(), writtenFile, writePeptideLevelResults: true);
            FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId });

            // write summary text
            if (peptidesForPeptideResults.FilteringNotPerformed)
            {
                Parameters.SearchTaskResults.AddPsmPeptideProteinSummaryText(
                    "PEP could not be calculated due to an insufficient number of PSMs. Results were filtered by q-value." + Environment.NewLine);
            }
            string peptideResultsText = $"All target {GlobalVariables.AnalyteType.ToLower()}s with " + peptidesForPeptideResults.GetFilterTypeString() + " <= " + Math.Round(peptidesForPeptideResults.FilterThreshold, 2) + ": " +
                peptidesForPeptideResults.TargetPsmsAboveThreshold;
            ResultsDictionary[("All", GlobalVariables.AnalyteType)] = peptideResultsText;
        }

        private void WriteIndividualPsmResults()
        {
            Status("Writing Individual PSM results...", Parameters.SearchTaskId);

            var psmsGroupedByFile = Parameters.AllPsms.GroupBy(p => p.FullFilePath);
            foreach (var psmFileGroup in psmsGroupedByFile)
            {
                // FDR Analysis is performed again for each file. File specific results show the results that would be 
                // generated by analyzing one file by itself. Therefore, the FDR info should change between AllResults and FileSpecific
                string strippedFileName = Path.GetFileNameWithoutExtension(psmFileGroup.Key);
                var psmsForThisFile = psmFileGroup.ToList();
                CalculatePsmAndPeptideFdr(psmsForThisFile,"PSM", false);
                var psmsToWrite = FilteredPsms.Filter(psmsForThisFile,
                    CommonParameters, 
                    includeDecoys: Parameters.SearchParameters.WriteDecoys,
                    includeContaminants: Parameters.SearchParameters.WriteContaminants,
                    includeAmbiguous: true,
                    includeHighQValuePsms: Parameters.SearchParameters.WriteHighQValuePsms);

                int count = psmsToWrite.Where(psm => psm.PsmFdrInfo.PEP <= 0.01).Count();

                // write PSMs
                string writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_PSMs.psmtsv");
                WritePsmsToTsv(psmsToWrite, writtenFile);
                FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", psmFileGroup.Key });

                // write PSMs for percolator
                writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_PSMsFormattedForPercolator.tab");
                WritePsmsForPercolator(psmsToWrite.FilteredPsmsList, writtenFile);
                FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", psmFileGroup.Key });

                // write summary text
                string psmResultsText = strippedFileName + " - Target PSMs with " + psmsToWrite.GetFilterTypeString() + " <= " + Math.Round(psmsToWrite.FilterThreshold, 2) + ": " +
                                        psmsToWrite.TargetPsmsAboveThreshold;
                ResultsDictionary[(strippedFileName, "PSMs")] = psmResultsText;
            }
        }
        private void WriteIndividualPeptideResults()
        {
            Status("Writing Individual Peptide results...", Parameters.SearchTaskId);

            var peptidesGroupedByFile = Parameters.AllPsms.GroupBy(p => p.FullFilePath);
            foreach (var psmFileGroup in peptidesGroupedByFile)
            {
                var peptideFileGroup = psmFileGroup
                    .OrderByDescending(p => p)
                    .GroupBy(p => p.FullSequence)
                    .Select(group => group.FirstOrDefault())
                    .ToList();

                // FDR Analysis is performed again for each file. File specific results show the results that would be 
                // generated by analyzing one file by itself. Therefore, the FDR info should change between AllResults and FileSpecific
                string strippedFileName = Path.GetFileNameWithoutExtension(psmFileGroup.Key);
                CalculatePsmAndPeptideFdr(peptideFileGroup, "peptide", false);
                var peptidesToWrite = FilteredPsms.Filter(peptideFileGroup,
                    CommonParameters,
                    includeDecoys: Parameters.SearchParameters.WriteDecoys,
                    includeContaminants: Parameters.SearchParameters.WriteContaminants,
                    includeAmbiguous: true,
                    includeHighQValuePsms: Parameters.SearchParameters.WriteHighQValuePsms,
                    filterAtPeptideLevel: true);

                // write PSMs
                string writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + $"_{GlobalVariables.AnalyteType}s.psmtsv");
                WritePsmsToTsv(peptidesToWrite, writtenFile, writePeptideLevelResults: true);
                FinishedWritingFile(writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", psmFileGroup.Key });

                // write summary text
                string peptideResultsText = strippedFileName + $" - Target {GlobalVariables.AnalyteType.ToLower()}s with " + peptidesToWrite.GetFilterTypeString() + " <= " + Math.Round(peptidesToWrite.FilterThreshold, 2) + ": " +
                                        peptidesToWrite.TargetPsmsAboveThreshold;
                ResultsDictionary[(strippedFileName, GlobalVariables.AnalyteType)] = peptideResultsText;
            }

        }
        private void UpdateSpectralLibrary()
        {
            var peptidesForSpectralLibrary = FilteredPsms.Filter(Parameters.AllPsms,
                CommonParameters,
                includeDecoys: false,
                includeContaminants: false,
                includeAmbiguous: false,
                includeHighQValuePsms: false
                );


            //group psms by peptide and charge, then write highest scoring PSM to dictionary
            Dictionary<(string, int), SpectralMatch> psmSeqChargeDictionary = peptidesForSpectralLibrary
                .GroupBy(p => (p.FullSequence, p.ScanPrecursorCharge))
                .ToDictionary(
                    // Key is a (FullSequence, Charge) tuple
                    keySelector: g => g.Key,
                    // Value is the highest scoring psm in the group
                    elementSelector: g => g.MaxBy(p => p.Score)); 

            //load the original library
            var originalLibrarySpectra = Parameters.SpectralLibrary.GetAllLibrarySpectra();
            List<LibrarySpectrum> updatedLibrarySpectra = new();

            // Add the better spectrum (library spectrum vs. current Psm) to updatedLibrarySpectra
            foreach (var ogLibrarySpectrum in originalLibrarySpectra)
            {
                //this statement is false of the old spectrum is not among the spectra observed in this search
                if (psmSeqChargeDictionary.TryGetValue((ogLibrarySpectrum.Sequence, ogLibrarySpectrum.ChargeState),
                        out var bestPsm))
                {
                    if (ogLibrarySpectrum.MatchedFragmentIons.Count > Math.Truncate(bestPsm.Score))
                    {
                        updatedLibrarySpectra.Add(ogLibrarySpectrum);
                    }
                    else
                    {
                        updatedLibrarySpectra.Add(new LibrarySpectrum(
                            bestPsm.FullSequence,
                            bestPsm.ScanPrecursorMonoisotopicPeakMz,
                            bestPsm.ScanPrecursorCharge,
                            bestPsm.MatchedFragmentIons,
                            bestPsm.ScanRetentionTime));
                    }
                    // once the spectrum is added, it is removed from the dictionary
                    // later we will add the remaining spectra to the updated library
                    psmSeqChargeDictionary.Remove((ogLibrarySpectrum.Sequence, ogLibrarySpectrum.ChargeState));
                }
                else // if the spectrum is not in the dictionary, we keep the original spectrum
                {
                    updatedLibrarySpectra.Add(ogLibrarySpectrum);
                }
            }

            // if we don't a spectrum in the original library, we add it to the updated library
            foreach (var bestPsm in psmSeqChargeDictionary.Values)
            {
                updatedLibrarySpectra.Add(new LibrarySpectrum(
                    bestPsm.FullSequence,
                    bestPsm.ScanPrecursorMonoisotopicPeakMz,
                    bestPsm.ScanPrecursorCharge,
                    bestPsm.MatchedFragmentIons,
                    bestPsm.ScanRetentionTime));
            }

            string updatedSpectralLibrary = UpdateSpectralLibrary(updatedLibrarySpectra, Parameters.OutputFolder);

            Parameters.SearchTaskResults.NewDatabases = new List<DbForTask> { new DbForTask(updatedSpectralLibrary, false) };

            DbForTask originalFastaDb = Parameters.DatabaseFilenameList.Where(p => p.IsSpectralLibrary == false && p.IsContaminant == false).First();
            Parameters.SearchTaskResults.NewDatabases.Add(originalFastaDb);
        }

        //for those spectra matching the same peptide/protein with same charge, save the one with highest score
        private void SpectralLibraryGeneration()
        {
            var peptidesForSpectralLibrary = FilteredPsms.Filter(Parameters.AllPsms,
                CommonParameters,
                includeDecoys: false,
                includeContaminants: false,
                includeAmbiguous: false,
                includeHighQValuePsms: false);

            //group psms by peptide and charge, the psms having same sequence and same charge will be in the same group
            var fullSeqChargeGrouping = peptidesForSpectralLibrary.GroupBy(p => (p.FullSequence, p.ScanPrecursorCharge));
            List<LibrarySpectrum> spectraLibrary = new();
            foreach (var matchGroup in fullSeqChargeGrouping)
            {
                SpectralMatch bestPsm = matchGroup.MaxBy(p => p.Score);
                if (bestPsm == null) continue;
                spectraLibrary.Add(new LibrarySpectrum(
                    bestPsm.FullSequence,
                    bestPsm.ScanPrecursorMonoisotopicPeakMz,
                    bestPsm.ScanPrecursorCharge,
                    bestPsm.MatchedFragmentIons,
                    bestPsm.ScanRetentionTime));
            }
            
            WriteSpectrumLibrary(spectraLibrary, Parameters.OutputFolder);
        }

        private void WriteProteinResults()
        {
            if (!Parameters.SearchParameters.DoParsimony)
            {
                return;
            }
            else
            {
                string proteinResultsText = "All target protein groups with q-value <= 0.01 (1% FDR): " + ProteinGroups.Count(b => b.QValue <= 0.01 && !b.IsDecoy);
                ResultsDictionary[("All", "Proteins")] = proteinResultsText;
            }
            
            string fileName = "AllProteinGroups.tsv";
            if (Parameters.SearchParameters.DoLabelFreeQuantification)
            {
                fileName = "AllQuantifiedProteinGroups.tsv";
            }

            //set peptide output values
            ProteinGroups.ForEach(x => x.GetIdentifiedPeptidesOutput(Parameters.SearchParameters.SilacLabels));
            // write protein groups to tsv
            string writtenFile = Path.Combine(Parameters.OutputFolder, fileName);
            WriteProteinGroupsToTsv(ProteinGroups, writtenFile, new List<string> { Parameters.SearchTaskId });

            var psmsGroupedByFile = FilteredPsms.Filter(Parameters.AllPsms,
                CommonParameters,
                includeDecoys: true,
                includeContaminants: true,
                includeAmbiguous: true,
                includeHighQValuePsms: false).FilteredPsmsList.GroupBy(f => f.FullFilePath);

            //if we're writing individual files, we need to reprocess the psms
            //If doing a SILAC search and no "unlabeled" labels were specified (i.e. multiple labels are used for multiplexing and no conditions are "unlabeled"),
            //then we need to update the psms (which were found in the data file that has the "unlabeled" named) and say they were found in the "heavy" file)
            if (Parameters.SearchParameters.SilacLabels != null) //if we have silac labels
            {
                //get the original filenames
                List<string> fileNamesThatHadPsms = psmsGroupedByFile.Select(v => v.Key).ToList();
                EngineLayer.ProteinGroup firstProteinGroup = ProteinGroups.FirstOrDefault(); //grab the first protein to extract the files used for quantification
                if (firstProteinGroup != null) //check that we even have a protein group to write
                {
                    var tempPsmsGroupedByFile = new List<IGrouping<string, SpectralMatch>>();
                    //foreach original file
                    foreach (string originalFile in fileNamesThatHadPsms)
                    {
                        //get all the "filenames" output by quantification. If no unlabeled condition was specified, the original datafile will not be present in the current grouping
                        //Example: the datafile "test.mzml" that was searched with +4 or +10 neutron mass difference on arginine would appear as "test(R+4).mzml" and "test(R+10).mzml".
                        //there would be no "test.mzml"
                        List<string> labeledFiles = new List<string> { originalFile };
                        foreach (SilacLabel label in Parameters.SearchParameters.SilacLabels)
                        {
                            //rediscover the previous naming conversion(s)
                            labeledFiles.Add(SilacConversions.GetHeavyFileInfo(new SpectraFileInfo(originalFile, "", 0, 0, 0), label).FullFilePathWithExtension);
                        }

                        //rename the file group for all of the relevant psms to their original file
                        List<SpectralMatch> psms = psmsGroupedByFile.Where(g => labeledFiles.Contains(g.Key)).SelectMany(x => x).ToList(); //grab all the psms
                        tempPsmsGroupedByFile.AddRange(psms.GroupBy(x => originalFile));
                    }
                    //overwrite the grouping for downstream processing
                    psmsGroupedByFile = tempPsmsGroupedByFile.ToList();
                }
            }

            //write the individual result files for each datafile
            foreach (var fullFilePath in psmsGroupedByFile.Select(v => v.Key))
            {
                string strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                List<SpectralMatch> psmsForThisFile = psmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList();
                var subsetProteinGroupsForThisFile = ProteinGroups.Select(p => p.ConstructSubsetProteinGroup(fullFilePath, Parameters.SearchParameters.SilacLabels)).ToList();

                ProteinScoringAndFdrResults subsetProteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(subsetProteinGroupsForThisFile, psmsForThisFile,
                    Parameters.SearchParameters.NoOneHitWonders, Parameters.SearchParameters.ModPeptidesAreDifferent,
                    false, CommonParameters, this.FileSpecificParameters, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }).Run();

                subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults.SortedAndScoredProteinGroups;

                if (Parameters.SearchParameters.WriteIndividualFiles && Parameters.CurrentRawFileList.Count > 1)
                {
                    // write summary text
                    string proteinResultsText = strippedFileName + " - Target protein groups within 1 % FDR: " + subsetProteinGroupsForThisFile.Count(b => b.QValue <= 0.01 && !b.IsDecoy);
                    ResultsDictionary[(strippedFileName, "Proteins")] = proteinResultsText;

                    // write result files
                    writtenFile = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + "_ProteinGroups.tsv");
                    WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, writtenFile, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                }

                psmsForThisFile = FilteredPsms.Filter(psmsForThisFile,
                    CommonParameters,
                    includeDecoys: Parameters.SearchParameters.WriteDecoys,
                    includeContaminants: Parameters.SearchParameters.WriteContaminants,
                    includeAmbiguous: true,
                    includeHighQValuePsms: true).FilteredPsmsList;

                // Filter psms in place before writing mzID
                if (Parameters.SearchParameters.WriteMzId)
                {
                    Status("Writing mzID...", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });

                    string mzidFilePath = Path.Combine(Parameters.OutputFolder, strippedFileName + ".mzID");
                    if (Parameters.CurrentRawFileList.Count > 1)
                    {
                        mzidFilePath = Path.Combine(Parameters.IndividualResultsOutputFolder, strippedFileName + ".mzID");
                    }

                    MzIdentMLWriter.WriteMzIdentMl(
                        psmsForThisFile,
                        subsetProteinGroupsForThisFile,
                        Parameters.VariableModifications,
                        Parameters.FixedModifications,
                        Parameters.SearchParameters.SilacLabels,
                        new List<Protease> { CommonParameters.DigestionParams.Protease },
                        CommonParameters.ProductMassTolerance,
                        CommonParameters.PrecursorMassTolerance,
                        CommonParameters.DigestionParams.MaxMissedCleavages,
                        mzidFilePath,
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

                    PepXMLWriter.WritePepXml(psmsForThisFile,
                        Parameters.DatabaseFilenameList,
                        Parameters.VariableModifications,
                        Parameters.FixedModifications,
                        CommonParameters, pepXMLFilePath);

                    FinishedWritingFile(pepXMLFilePath, new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath });
                }

                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { Parameters.SearchTaskId, "Individual Spectra Files", fullFilePath }));
            }
        }

        private void WriteFlashLFQResults()
        {
            if (Parameters.SearchParameters.DoLabelFreeQuantification && Parameters.FlashLfqResults != null)
            {
                // write peaks
                if (SpectralRecoveryResults != null)
                {
                    SpectralRecoveryResults.WritePeakQuantificationResultsToTsv(Parameters.OutputFolder, "AllQuantifiedPeaks");
                }
                else
                {
                    WritePeakQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, "AllQuantifiedPeaks", new List<string> { Parameters.SearchTaskId });
                }

                // write peptide quant results
                string filename = "AllQuantified" + GlobalVariables.AnalyteType + "s";
                if (SpectralRecoveryResults != null)
                {
                    SpectralRecoveryResults.WritePeptideQuantificationResultsToTsv(Parameters.OutputFolder, filename);
                }
                else
                {
                    WritePeptideQuantificationResultsToTsv(Parameters.FlashLfqResults, Parameters.OutputFolder, filename, new List<string> { Parameters.SearchTaskId });
                }

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

        private void WritePrunedDatabase()
        {
            if (Parameters.SearchParameters.WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { Parameters.SearchTaskId });
                HashSet<Modification> modificationsToWriteIfBoth = new HashSet<Modification>();
                HashSet<Modification> modificationsToWriteIfInDatabase = new HashSet<Modification>();
                HashSet<Modification> modificationsToWriteIfObserved = new HashSet<Modification>();

                var filteredPsms = FilteredPsms.Filter(Parameters.AllPsms,
                    CommonParameters,
                    includeDecoys: false,
                    includeContaminants: true,
                    includeAmbiguous: false,
                    includeHighQValuePsms: false);

                var proteinToConfidentBaseSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                // associate all confident PSMs with all possible proteins they could be digest products of (before or after parsimony)
                foreach (SpectralMatch psm in filteredPsms)
                {
                    var myPepsWithSetMods = psm.BestMatchingBioPolymersWithSetMods.Select(p => p.Peptide);

                    foreach (PeptideWithSetModifications peptide in myPepsWithSetMods)
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
                var originalModPsms = FilteredPsms.Filter(filteredPsms,
                        CommonParameters,
                        includeDecoys: false,
                        includeContaminants: true,
                        includeAmbiguous: false,
                        includeAmbiguousMods: false,
                        includeHighQValuePsms: false);


                var proteinToConfidentModifiedSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                HashSet<string> modPsmsFullSeq = originalModPsms.Select(p => p.FullSequence).ToHashSet();
                HashSet<string> originalModPsmsFullSeq = originalModPsms.Select(p => p.FullSequence).ToHashSet();
                modPsmsFullSeq.ExceptWith(originalModPsmsFullSeq);

                foreach (SpectralMatch psm in originalModPsms)
                {
                    var myPepsWithSetMods = psm.BestMatchingBioPolymersWithSetMods.Select(p => p.Peptide);

                    foreach (PeptideWithSetModifications peptide in myPepsWithSetMods)
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

        private void WritePsmPlusMultiplexIons(IEnumerable<SpectralMatch> psms, string filePath, bool writePeptideLevelResults = false)
        {
            PpmTolerance ionTolerance = new PpmTolerance(10);
            double[] reporterIonMzs = Parameters.MultiplexModification.DiagnosticIons.First().Value
                .Select(x => x.ToMz(1))
                .OrderBy(x => x)
                .ToArray();

            using (StreamWriter output = new StreamWriter(filePath))
            {
                string headerWithReporterIons = SpectralMatch.GetTabSeparatedHeader().Trim() + '\t' +
                                          GetMultiplexHeader();
                output.WriteLine(headerWithReporterIons);
                foreach (var psm in psms)
                {
                    IEnumerable<string> labelIonIntensities =
                        GetMultiplexIonIntensities(psm.MsDataScan.MassSpectrum, reporterIonMzs, ionTolerance)
                            .Select(d => d.ToString(CultureInfo.CurrentCulture));

                    output.Write(psm.ToString(Parameters.SearchParameters.ModsToWriteSelection, writePeptideLevelResults).Trim());
                    output.Write('\t');
                    output.WriteLine(String.Join('\t', labelIonIntensities));
                }
            }
        }

        private string GetMultiplexHeader()
        {
            List<string> ionLabels = new();
            var labelGroups = Parameters.MultiplexModification.DiagnosticIons.First().Value
                .Select(x => x.ToMz(1))
                .OrderBy(x => x)
                .GroupBy(x => (int)Math.Floor(x));

            if (Parameters.MultiplexModification.IdWithMotif.Contains("TMT"))
            {
                // TMT 126 contains no heavy isotopes. TMT 127N has one N15, TMT127C has one C15.
                // The "N" labels are slightly lighter than the "C" labels. 
                // Labels for the diagnostic ions are created accordingly
                foreach (var group in labelGroups)
                {
                    if (group.Count() == 1)
                    {
                        ionLabels.Add(group.Key.ToString());
                    }
                    else if (group.Count() == 2)
                    {
                        ionLabels.Add(group.Key + "N");
                        ionLabels.Add(group.Key + "C");
                    }
                }
            }
            else
            {
                foreach (var group in labelGroups)
                {
                    if (group.Count() == 1)
                    {
                        ionLabels.Add(group.Key.ToString());
                    }
                    else
                    {
                        ionLabels.AddRange(group.Select(mz => Math.Round(mz, 3).ToString(CultureInfo.CurrentCulture)));
                    }
                }
            }
            return String.Join('\t', ionLabels);
        }

        public static double[] GetMultiplexIonIntensities(MzSpectrum scan, double[] theoreticalIonMzs, Tolerance tolerance)
        {
            int peakIndex = scan.GetClosestPeakIndex(theoreticalIonMzs[0]);
            int lastPeakIndex = Math.Min(scan.GetClosestPeakIndex(theoreticalIonMzs.Last()) + 1, scan.XArray.Length - 1);
            double[] ionIntensities = new double[theoreticalIonMzs.Length];
            
            for (int ionIndex = 0; ionIndex < ionIntensities.Length; ionIndex++)
            {
                while (peakIndex <= lastPeakIndex && 
                       scan.XArray[peakIndex] < tolerance.GetMinimumValue(theoreticalIonMzs[ionIndex]))
                {
                    peakIndex++;
                }
                if (peakIndex > lastPeakIndex)
                {
                    break;
                }
                if (tolerance.Within(scan.XArray[peakIndex], theoreticalIonMzs[ionIndex]))
                {
                    ionIntensities[ionIndex] = scan.YArray[peakIndex];
                    peakIndex++;
                }
            }
            
            return ionIntensities;
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

            var fdrPsms = FilteredPsms.Filter(Parameters.AllPsms,
                        CommonParameters,
                        includeDecoys: true,
                        includeContaminants: true,
                        includeAmbiguous: true,
                        includeHighQValuePsms: false,
                        filterAtPeptideLevel: true);

            var possibleVariantPsms = fdrPsms.Where(p =>
                    p.BestMatchingBioPolymersWithSetMods.Any(pep => pep.Peptide is PeptideWithSetModifications pwsm && pwsm.IsVariantPeptide()))
                .OrderByDescending(pep => pep.Score)
                .ToList();

            new FdrAnalysisEngine(possibleVariantPsms, Parameters.NumNotches, CommonParameters, FileSpecificParameters,
                new List<string> { Parameters.SearchTaskId }, "variant_PSMs", doPEP: false).Run();

            possibleVariantPsms
                .OrderBy(p => p.FdrInfo.QValue)
                .ThenByDescending(p => p.Score)
                .ThenBy(p => p.FdrInfo.CumulativeTarget)
                .ToList();

            WritePsmsToTsv(possibleVariantPsms, variantPsmFile);

            List<SpectralMatch> variantPeptides = possibleVariantPsms
                .GroupBy(b => b.FullSequence)
                .Select(b => b.FirstOrDefault())
                .OrderByDescending(b => b.Score)
                .ToList();
            List<SpectralMatch> confidentVariantPeps = new List<SpectralMatch>();

            new FdrAnalysisEngine(variantPeptides, Parameters.NumNotches, CommonParameters, FileSpecificParameters,
                new List<string> { Parameters.SearchTaskId }, "variant_Peptides", doPEP: false).Run();

            WritePsmsToTsv(variantPeptides, variantPeptideFile);

            // if a potential variant peptide can be explained by a canonical protein seqeunce then should not be counted as a confident variant peptide
            //because it is most probable that the peptide originated from the canonical protien.
            foreach (var entry in variantPeptides)
            {
                var pwsm = entry.BestMatchingBioPolymersWithSetMods;
                var nonVariantOption = pwsm.Any(p => p.Peptide is PeptideWithSetModifications pwsm && pwsm.IsVariantPeptide() == false);
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
            Dictionary<Protein, HashSet<SequenceVariation>> MNVmissenseVariants = new();
            Dictionary<Protein, HashSet<SequenceVariation>> SNVmissenseVariants = new();
            Dictionary<Protein, HashSet<SequenceVariation>> insertionVariants = new();
            Dictionary<Protein, HashSet<SequenceVariation>> deletionVariants = new();
            Dictionary<Protein, HashSet<SequenceVariation>> frameshiftVariants = new();
            Dictionary<Protein, HashSet<SequenceVariation>> stopGainVariants = new();
            Dictionary<Protein, HashSet<SequenceVariation>> stopLossVariants = new();

            var filteredVariants = FilteredPsms.Filter(confidentVariantPeps,
                        CommonParameters,
                        includeDecoys: false,
                        includeContaminants: false,
                        includeAmbiguous: false,
                        includeHighQValuePsms: false);

            List<PeptideSpectralMatch> modifiedVariantPeptides = filteredVariants
                .Where(p => p.ModsIdentified != null && p.ModsIdentified.Count > 0 && p is PeptideSpectralMatch)
                .Select(p => (PeptideSpectralMatch)p)
                .ToList(); //modification can be on any AA in variant peptide

            List<PeptideSpectralMatch> modifiedVariantSitePeptides = new();// modification is speciifcally on the variant residue within the peptide
            foreach (PeptideSpectralMatch entry in modifiedVariantPeptides)
            {
                PeptideWithSetModifications firstOrDefault = entry.BestMatchingBioPolymersWithSetMods.FirstOrDefault().Peptide as PeptideWithSetModifications;

                var variantPWSM = firstOrDefault;
                var peptideMods = variantPWSM.AllModsOneIsNterminus.Values.ToList();
                var variantProteinModifications = variantPWSM.Parent.OneBasedPossibleLocalizedModifications.Where(k => k.Key >= variantPWSM.OneBasedStartResidue && k.Key <= variantPWSM.OneBasedEndResidue).ToList();
                var variants = firstOrDefault.Protein.AppliedSequenceVariations.Where(v => firstOrDefault.IntersectsAndIdentifiesVariation(v).identifies).ToList();
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
                var variantPWSM =
                    peptide.BestMatchingBioPolymersWithSetMods.FirstOrDefault() is (_, PeptideWithSetModifications)
                        ? ((int Notch, PeptideWithSetModifications Peptide))peptide.BestMatchingBioPolymersWithSetMods
                            .FirstOrDefault()
                        : (0, null);//TODO: expand to all peptide options not just the first
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
            variantResults[4] = "Number of potential variant containing peptides identified at " + fdrPsms.FilterThreshold * 100 + "% group FDR: " + fdrPsms.TargetPsmsAboveThreshold;
            variantResults[5] = "Number of unqiuely identified variant peptides at " + filteredVariants.FilterThreshold * 100 + "% group FDR: " + filteredVariants.TargetPsmsAboveThreshold;
            variantResults[6] = "Number of unique variants: " + totalVariantSites;
            variantResults[7] = "Number of SNV missense variant containing peptides at " + fdrPsms.FilterThreshold * 100 + "% group FDR: " + SNVmissenseCount;
            variantResults[8] = "Number of unique SNV missense variants: " + SNVmissenseSites;
            variantResults[9] = "Number of MNV missense variant containing peptides at " + fdrPsms.FilterThreshold * 100 + "% group FDR: " + MNVmissenseCount;
            variantResults[10] = "Number of unique MNV missense variants: " + MNVmissenseSites;
            variantResults[11] = "Number of frameshift variant containing peptides at " + fdrPsms.FilterThreshold * 100 + "% group FDR: " + frameshiftCount;
            variantResults[12] = "Number of unique frameshift variants: " + frameshiftSites;
            variantResults[13] = "Number of inframe insertion variant containing peptides at " + fdrPsms.FilterThreshold * 100 + "% group FDR: " + insertionCount;
            variantResults[14] = "Number of unique inframe insertion variants: " + insertionSites;
            variantResults[15] = "Number of inframe deletion variant containing peptides at " + fdrPsms.FilterThreshold * 100 + "% group FDR: " + deletionCount;
            variantResults[16] = "Number of unique inframe deletion variants: " + deletionSites;
            variantResults[17] = "Number of stop gain variant containing peptides at " + fdrPsms.FilterThreshold * 100 + "% group FDR: " + stopGainCount;
            variantResults[18] = "Number of unique stop gain variants: " + stopGainSites;
            variantResults[19] = "Number of stop loss variant containing peptides at " + fdrPsms.FilterThreshold * 100 + "% group FDR: " + stopLossCount;
            variantResults[20] = "Number of unique stop loss variants: " + stopLossSites;
            variantResults[21] = "Number of variant peptides at " + fdrPsms.FilterThreshold * 100 + "% group FDR with unambiguous localized modifications: " + modifiedVariantPeptides.Count;
            variantResults[22] = "Number of variant peptides at " + fdrPsms.FilterThreshold * 100 + "% group FDR with unambiguous localized modifications at the variant sites : " + modifiedVariantSitePeptides.Count;

            string filePath = Path.Combine(Parameters.OutputFolder, "VariantAnalysisResultSummary.txt");
            File.WriteAllLines(filePath, variantResults);
        }

        private static int GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications peptideWithSetModifications)
        {
            if (oneIsNterminus == 1)
            {
                return peptideWithSetModifications.OneBasedStartResidue;
            }
            if (oneIsNterminus == peptideWithSetModifications.Length + 2)
            {
                return peptideWithSetModifications.OneBasedEndResidue;
            }
            return peptideWithSetModifications.OneBasedStartResidue + oneIsNterminus - 2;
        }

        private static void WriteTree(BinTreeStructure myTreeStructure, string writtenFile)
        {
            using (StreamWriter output = new(writtenFile))
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

        private static void WritePsmsForPercolator(List<SpectralMatch> psmList, string writtenFileForPercolator)
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
                header += String.Join("\t", PsmData.trainingInfos[searchType]);
                header += "\tPeptide\tProteins";

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
                foreach (SpectralMatch psm in psmList.Where(p => p.PsmData_forPEPandPercolator != null))
                {
                    foreach (var peptide in psm.BestMatchingBioPolymersWithSetMods)
                    {
                        output.Write(idNumber.ToString());
                        output.Write('\t' + (peptide.Peptide.Parent.IsDecoy ? -1 : 1).ToString());
                        output.Write('\t' + psm.ScanNumber.ToString());
                        output.Write(psm.PsmData_forPEPandPercolator.ToString(searchType));
                        output.Write('\t' + (peptide.Peptide.PreviousResidue + "." + peptide.Peptide.FullSequence + "." + peptide.Peptide.NextResidue).ToString());
                        output.Write('\t' + (peptide.Peptide.Parent.Accession).ToString());
                        output.WriteLine();
                    }
                    idNumber++;
                }
            }
        }

        private void WriteProteinGroupsToTsv(List<EngineLayer.ProteinGroup> proteinGroups, string filePath, List<string> nestedIds)
        {
            if (proteinGroups != null && proteinGroups.Any())
            {
                double qValueThreshold = Math.Min(CommonParameters.QValueThreshold, CommonParameters.PepQValueThreshold);
                using (StreamWriter output = new StreamWriter(filePath))
                {
                    output.WriteLine(proteinGroups.First().GetTabSeparatedHeader());
                    for (int i = 0; i < proteinGroups.Count; i++)
                    {
                        if ((!Parameters.SearchParameters.WriteDecoys && proteinGroups[i].IsDecoy) ||
                            (!Parameters.SearchParameters.WriteContaminants && proteinGroups[i].IsContaminant) ||
                            (!Parameters.SearchParameters.WriteHighQValuePsms && proteinGroups[i].QValue > qValueThreshold))
                        {
                            continue;
                        }
                        else
                        {
                            output.WriteLine(proteinGroups[i]);
                        }
                    }
                }
                FinishedWritingFile(filePath, nestedIds);
            }
        }

        /// <summary>
        /// This is a handy dictionary to keep track of the PSM, peptide and protein count results at the
        ///  "All" level and at the individual raw file level.
        ///  The keys are a tuple such as ("All", "PSMs") or ("RawFileName", "Peptides")
        ///   The values are the results as a string
        /// </summary>
        private void ConstructResultsDictionary()
        {
            ResultsDictionary = new()
            {
                { ("All", "PSMs"), "" },
                { ("All", GlobalVariables.AnalyteType), "" }
            };

            if (Parameters.CurrentRawFileList.Count > 1 && Parameters.SearchParameters.WriteIndividualFiles)
            {
                foreach (var rawFile in Parameters.CurrentRawFileList)
                {
                    string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(rawFile);
                    ResultsDictionary.Add((fileNameWithoutExtension, "PSMs"), "");
                    ResultsDictionary.Add((fileNameWithoutExtension, GlobalVariables.AnalyteType), "");
                }
            }

            if (Parameters.SearchParameters.DoParsimony)
            {
                ResultsDictionary.Add(("All", "Proteins"), ""); 
                if (Parameters.CurrentRawFileList.Count > 1 && Parameters.SearchParameters.WriteIndividualFiles)
                {
                    foreach (var rawFile in Parameters.CurrentRawFileList)
                    {
                        string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(rawFile);
                        ResultsDictionary.Add((fileNameWithoutExtension, "Proteins"), "");
                    }
                }
            }
        }

        private string AllResultsTotals()
        {
            StringBuilder sb = new();
            foreach (var key in ResultsDictionary.Keys)
            {
                if (key.Item1 == "All")
                {
                    sb.AppendLine(ResultsDictionary[key]);
                }
            }
            
            var keys = ResultsDictionary.Keys.Where(k=>k.Item1 != "All").OrderBy(k=>k.Item1).ToList();
            if (keys.Any())
            {
                sb.AppendLine();
                var item1 = keys.First().Item1;
                foreach (var key in keys)
                {
                    if (key.Item1 != item1)
                    {
                        sb.AppendLine();
                        item1 = key.Item1;
                    }
                    sb.AppendLine(ResultsDictionary[key]);
                }
            }
            return sb.ToString();
        }

        private void AddResultsTotalsToAllResultsTsv()
        {
            Parameters.SearchTaskResults.AddPsmPeptideProteinSummaryText(AllResultsTotals());
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