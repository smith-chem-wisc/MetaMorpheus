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

            //if SILAC, modify the proteins to appear only light (we want a protein sequence to look like PROTEINK instead of PROTEINa)
            if (Parameters.SearchParameters.SilacLabels != null)
            {
                Parameters.AllPsms = SilacConversions.UpdateProteinSequencesToLight(Parameters.AllPsms, Parameters.SearchParameters.SilacLabels);
            }

            List<PeptideSpectralMatch> psmsForProteinParsimony = Parameters.AllPsms;

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
            (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = Parameters.SearchParameters.TurnoverLabels;
            bool quantifyUnlabeledPeptides = Parameters.ListOfDigestionParams.Any(x => x.GeneratehUnlabeledProteinsForSilac);
            if (Parameters.SearchParameters.SilacLabels != null)
            {
                //go through all the psms and duplicate them until a psm copy exists for the unlabeled and labeled proteins
                //The number of psms should roughly increase by a factor of N, where N is the number of labels.
                //It may not increase exactly by a factor of N if the amino acid(s) that gets labeled doesn't exist in the peptide

                List<PeptideSpectralMatch> silacPsms = new List<PeptideSpectralMatch>(); //populate with duplicate psms for heavy/light

                //multiply the psms by the number of labels
                foreach (PeptideSpectralMatch psm in unambiguousPsmsBelowOnePercentFdr)
                {
                    //get the original proteinGroup to give to the other psm clones
                    List<FlashLFQ.ProteinGroup> originalProteinGroups = psmToProteinGroups[psm];

                    //see which label, if any, this peptide has
                    string peptideBaseSequence = psm.BaseSequence;
                    SilacLabel observedLabel = SilacConversions.GetRelevantLabelFromBaseSequence(peptideBaseSequence, allSilacLabels); //returns null if no label

                    //if it's not the light form, make a light form to start as a base. 
                    PeptideSpectralMatch lightPsm = observedLabel == null ? psm : SilacConversions.GetSilacPsm(psm, observedLabel, true);

                    //get easy access to values we need for new psm generation
                    string unlabeledBaseSequence = lightPsm.BaseSequence;
                    int notch = psm.BestMatchingPeptides.First().Notch;
                    PeptideWithSetModifications pwsm = psm.BestMatchingPeptides.First().Peptide;

                    //check if turnover or multiplex experiment
                    if (turnoverLabels == null) //if multiplex
                    {
                        //If we need the light form, then add it
                        if (quantifyUnlabeledPeptides)
                        {
                            silacPsms.Add(lightPsm);
                            psmToProteinGroups[lightPsm] = originalProteinGroups; //add proteingroup info
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
                                    psmToProteinGroups[labeledPsm] = originalProteinGroups; //add proteingroup info
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
                        SilacLabel startLabel = turnoverLabels.Value.startLabel;
                        SilacLabel endLabel = turnoverLabels.Value.endLabel;

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
                            int index = differentResidues[i];
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
                //do protein quant if we had any results
                //if no results, we still may need to edit the psms
                if (Parameters.FlashLfqResults != null) //can be null if no unambiguous psms were found
                {
                    //after this point, we now have quantification values for the peptides, but they all belong to the same "unlabeled" protein and are in the same file
                    //We can remove "labeled" peptides from each file and put them in a new file as "unlabeled".

                    //MAKE NEW RAW FILES
                    //update number of spectra files to include a new file for each label*condition
                    List<SpectraFileInfo> silacSpectraFileInfo = new List<SpectraFileInfo>(); //new files
                    Dictionary<SpectraFileInfo, List<SpectraFileInfo>> originalToLabeledFileInfoDictionary = new Dictionary<SpectraFileInfo, List<SpectraFileInfo>>();
                    List<SpectraFileInfo> unlabeledFiles = new List<SpectraFileInfo>(); //null can't be a key in the dictionary, so this list is used for the unlabeled files.

                    Parameters.FlashLfqResults.SpectraFiles.Clear(); //clear existing files so we can replace them with labeled ones

                    //foreach existing file
                    //if multiplex
                    if (turnoverLabels == null)   
                    {
                        //populate dictionary
                        if (quantifyUnlabeledPeptides)
                        {
                            spectraFileInfo.ForEach(x => originalToLabeledFileInfoDictionary.Add(x, new List<SpectraFileInfo> { x }));
                            //get the unlabeled
                            unlabeledFiles = spectraFileInfo;
                            Parameters.FlashLfqResults.SpectraFiles.AddRange(spectraFileInfo);
                        }
                        else
                        {
                            spectraFileInfo.ForEach(x => originalToLabeledFileInfoDictionary.Add(x, new List<SpectraFileInfo>()));
                        }

                        //get the labeled
                        foreach (SilacLabel label in allSilacLabels)
                        {
                            List<SpectraFileInfo> labeledFiles = new List<SpectraFileInfo>();
                            foreach (SpectraFileInfo originalFile in spectraFileInfo)
                            {
                                //foreach label, add a new file with the label                      
                                //silacSpectraFileInfo.Add(silacFile);
                                SpectraFileInfo labeledInfo = SilacConversions.GetHeavyFileInfo(originalFile, label);
                                labeledFiles.Add(labeledInfo);
                                originalToLabeledFileInfoDictionary[originalFile].Add(labeledInfo);
                            }
                            Parameters.FlashLfqResults.SpectraFiles.AddRange(labeledFiles);
                        }
                    }

                    //we have the files, now let's reassign the psms.
                    //there are a few ways to do this, but we're going to generate the "base" peptide and assign to that

                    //Dictionary of protein accession to peptides
                    Dictionary<string, List<(SilacLabel label, Peptide peptide)>> unlabeledToPeptidesDictionary = new Dictionary<string, List<(SilacLabel label, Peptide peptide)>>();
                    Dictionary<string, Peptide> updatedPeptideModifiedSequences = new Dictionary<string, Peptide>();
                    IEnumerable<Peptide> lfqPwsms = Parameters.FlashLfqResults.PeptideModifiedSequences.Values;
                    if (turnoverLabels == null || turnoverLabels.Value.startLabel == null || turnoverLabels.Value.endLabel == null) //if multiplex, or if at least one of the turnover labels is unlabeled
                    {
                        foreach (Peptide peptide in lfqPwsms)
                        {
                            //convert to the unlabeled sequence
                            string labeledSequence = peptide.Sequence;
                            SilacLabel label = SilacConversions.GetRelevantLabelFromFullSequence(labeledSequence, allSilacLabels);
                            string unlabeledSequence = SilacConversions.GetSilacLightFullSequence(labeledSequence, label, false);
                            if (unlabeledToPeptidesDictionary.ContainsKey(unlabeledSequence))
                            {
                                unlabeledToPeptidesDictionary[unlabeledSequence].Add((label, peptide));
                            }
                            else
                            {
                                unlabeledToPeptidesDictionary.Add(unlabeledSequence, new List<(SilacLabel, Peptide)> { (label, peptide) });
                            }
                        }
                    }
                    else //if both start and end labels exist for turnover, then it can be trickier to find the start label
                    {
                        SilacLabel startLabel = turnoverLabels.Value.startLabel;
                        SilacLabel endLabel = turnoverLabels.Value.endLabel;
                        foreach (Peptide peptide in lfqPwsms)
                        {
                            string originalSequence = peptide.Sequence;
                            string partiallyCleanedSequence = originalSequence;
                            SilacLabel labelToReport = null;
                            //convert to the unlabeled sequence
                            if (endLabel != null)
                            {
                                partiallyCleanedSequence = SilacConversions.GetSilacLightFullSequence(originalSequence, endLabel, false);
                                if (!partiallyCleanedSequence.Equals(originalSequence)) //if any residue is new, then the whole peptide must be newly synthesized
                                {
                                    labelToReport = endLabel;
                                }
                            }
                            string fullyCleanedSequence = partiallyCleanedSequence;
                            if (startLabel != null) //there might also be some old residues in the new peptide, so this is an "if" and not an "else"
                            {
                                fullyCleanedSequence = SilacConversions.GetSilacLightFullSequence(partiallyCleanedSequence, startLabel, false);
                                if (!fullyCleanedSequence.Equals(partiallyCleanedSequence) && labelToReport == null)
                                {
                                    labelToReport = startLabel;
                                }
                            }

                            if (unlabeledToPeptidesDictionary.ContainsKey(fullyCleanedSequence))
                            {
                                unlabeledToPeptidesDictionary[fullyCleanedSequence].Add((labelToReport, peptide));
                            }
                            else
                            {
                                unlabeledToPeptidesDictionary.Add(fullyCleanedSequence, new List<(SilacLabel, Peptide)> { (labelToReport, peptide) });
                            }
                        }
                    }

                    //we now have a dictionary of unlabeledBaseSequence to the labeled peptides
                    //Better SILAC results can be obtained by using the summed intensities from ms1 scans where all peaks were found, rather than the apex
                    //foreach peptide, unlabeled peptide, get the isotopic envelope intensities for each labeled peptide in each file
                    //save the intensities from ms1s that are shared. If no ms1s contains all the peaks, then just use the apex intensity (default)
                    foreach (KeyValuePair<SpectraFileInfo, List<ChromatographicPeak>> kvp in Parameters.FlashLfqResults.Peaks)
                    {
                        //make a dictionary for easy peak lookup from peptide sequences
                        Dictionary<string, ChromatographicPeak> sequenceToPeakDictionary = new Dictionary<string, ChromatographicPeak>();
                        foreach(ChromatographicPeak peak in kvp.Value)
                        {
                            sequenceToPeakDictionary.Add(peak.Identifications.First().ModifiedSequence, peak);
                        }

                        foreach (List<(SilacLabel label, Peptide peptide)> peptideGroup in unlabeledToPeptidesDictionary.Select(x => x.Value))
                        {
                            List<string> sequences = peptideGroup.Select(x => x.peptide.Sequence).ToList();

                            //get peaks of interest for this peptide group
                            List<ChromatographicPeak> peaks = kvp.Value;
                            List<ChromatographicPeak> peaksOfInterest = new List<ChromatographicPeak>();
                            foreach(string sequence in sequences)
                            {
                                if(sequenceToPeakDictionary.TryGetValue(sequence, out ChromatographicPeak peakOfInterest))
                                {
                                    peaksOfInterest.Add(peakOfInterest);
                                }
                            }

                            //If there are fewer than 2 peaks, we can't do any comparisons
                            if (peaksOfInterest.Count > 1)
                            {
                                //get isotopic envelopes that are shared by all
                                List<int> scanIndex = peaksOfInterest.First().IsotopicEnvelopes.Select(x => x.IndexedPeak).Select(x => x.ZeroBasedMs1ScanIndex).ToList();
                                for (int i = 1; i < peaksOfInterest.Count; i++)
                                {
                                    List<int> currentScanIndexes = peaksOfInterest[i].IsotopicEnvelopes.Select(x => x.IndexedPeak).Select(x => x.ZeroBasedMs1ScanIndex).ToList();
                                    scanIndex = scanIndex.Intersect(currentScanIndexes).ToList();
                                    if (scanIndex.Count == 0) //if there's no overlap, then we're done!
                                    {
                                        break;
                                    }
                                }
                                //if we aren't sticking with the default values
                                if (scanIndex.Count != 0)
                                {
                                    //update peptides
                                    foreach (Peptide peptide in peptideGroup.Select(x => x.peptide))
                                    {
                                        ChromatographicPeak peakForThisPeptide = peaksOfInterest.Where(x => peptide.Sequence.Equals(x.Identifications.First().ModifiedSequence)).First();
                                        double summedIntensity = peakForThisPeptide.IsotopicEnvelopes.Where(x => scanIndex.Contains(x.IndexedPeak.ZeroBasedMs1ScanIndex)).Select(x => x.Intensity).Sum();
                                        peptide.SetIntensity(kvp.Key, summedIntensity);
                                    }
                                }
                            }
                        }
                    }


                    //SPLIT THE FILES
                    List<Peptide> updatedPeptides = new List<Peptide>();

                    //split the heavy/light peptides into separate raw files, remove the heavy peptide
                    if (turnoverLabels != null)
                    {
                        Dictionary<SpectraFileInfo, double[]> fileToRecycleDictionary = new Dictionary<SpectraFileInfo, double[]>();
                        //TODO just want the original files, no reason to look at the silac ones
                        foreach (SpectraFileInfo info in spectraFileInfo)
                        {
                            fileToRecycleDictionary[info] = new double[3];
                        }

                        //go through each peptide that has more than one label
                        foreach (KeyValuePair<string, List<(SilacLabel label, Peptide peptide)>> kvp in unlabeledToPeptidesDictionary.Where(x => x.Value.Count > 2))
                        {
                            //the order should always be from start to end
                            //TODO if that's not true, use the silac label
                            //TODO if that is true, remove the silaclabel from the dictionary
                            //TODO check that values are present for light and for double the label mass to prevent incorrectly id'd missed cleavages from swaying data?
                            Dictionary<SpectraFileInfo, List<double>> tempFileToRecycleDictionary = new Dictionary<SpectraFileInfo, List<double>>();
                            foreach (SpectraFileInfo info in spectraFileInfo)
                            {
                                tempFileToRecycleDictionary[info] = new List<double>();
                            }

                            List<Peptide> peptides = kvp.Value.Select(x => x.peptide).ToList();
                            foreach (Peptide peptide in peptides) //assumes sorted old->new
                            {
                                foreach (SpectraFileInfo info in spectraFileInfo)
                                {
                                    tempFileToRecycleDictionary[info].Add(peptide.GetIntensity(info));
                                }
                            }

                            //aggregate into LL, HL, and HH if they contain 2+ missed cleavages
                            foreach (SpectraFileInfo info in spectraFileInfo)
                            {
                                List<double> values = tempFileToRecycleDictionary[info];
                                if (values.Count > 6)
                                {
                                    //we can't really use data with this many missed cleavages (data isn't super accurate), so remove them
                                    values = new List<double> { 0, 0, 0 };
                                }
                                else if (values.Count > 3)
                                {
                                    //TODO: make sure this works with a ton of missed cleavages (5)
                                    if (values.Count == 6)
                                    {
                                        //LLLLL, LLLLH, LLLHH, LLHHH, LHHHH, HHHHHH have n combinations of 1, 5, 10, 10, 5, 1
                                        double LLLL = values[0] + values[1] / 5;
                                        double LLLH = values[1] * 4 / 5 + values[2] * 4 / 10;
                                        double LLHH = values[2] * 6 / 10 + values[3] * 6 / 10;
                                        double LHHH = values[3] * 4 / 10 + values[4] * 4 / 5;
                                        double HHHH = values[4] / 5 + values[5];
                                        values = new List<double> { LLLL, LLLH, LLHH, LHHH, HHHH };
                                    }
                                    if (values.Count == 5)
                                    {
                                        //LLLL, LLLH, LLHH, LHHH, HHHH have n combinations of 1, 4, 6, 4, 1
                                        double LLL = values[0] + values[1] / 4;
                                        double LLH = values[1] * 3 / 4 + values[2] / 2; // divide by two is multiply by three divide by six
                                        double LHH = values[2] / 2 + values[3] * 3 / 4;
                                        double HHH = values[3] / 4 + values[4];
                                        values = new List<double> { LLL, LLH, LHH, HHH };
                                    }
                                    if (values.Count == 4)
                                    {
                                        //LLL, LLH, LHH, HHH have probabilities of 1, 3, 3, 1
                                        double LL = values[0] + values[1] / 3; // divide by three for the probability difference (only one way to make LL from LLH)
                                        double LH = values[1] * 2 / 3 + values[2] * 2 / 3; //divide by three, multiply by two (two ways to make LH from LLH or from LHH)
                                        double HH = values[2] / 3 + values[3]; //same logic as LL
                                        values = new List<double> { LL, LH, HH };
                                    }

                                    //TODO: validate this
                                    //overwrite the current peptide quantification values with the total values shown here IF there was more than one missed cleavage
                                    for (int i = 0; i < 3; i++)
                                    {
                                        peptides[i].SetIntensity(info, values[i]);
                                        //TODO: does the detection type stay here, or do we need to set it?
                                    }
                                }
                                //else if three values, then we're fine
                                //add the values from this peptide to the file-specific-values of all peptides
                                double[] totalValues = fileToRecycleDictionary[info];
                                for (int i = 0; i < values.Count; i++)
                                {
                                    totalValues[i] += values[i];
                                }
                            }
                        }

                        //with the summed LL, LH, and HH values for each file, we can calculate the probability that a heavy amino acid is incorporated (Ph)
                        //the probability a light amino acid is incorporated (Pl) is simply Pl = 1 - Ph
                        //We can simplify the math if we assume that all protein started as completely light, @t=0: LL=100%, LH = 0%, and HH = 0%
                        //N = 1 - Qll + Qlh * Qlh / (4 * Qhh), (where N is the fraction of all protein that has turned over and Qll is the Quantity(fraction) of LL at the current timepoint
                        //Qhh = qhh -N*qhh+Ph*Ph*N, where Q is quantity(fraction) at the timepoint and q is the quantity(fraction) at the start point (time zero).
                        //Ph = sqrt((Qhh-qhh+N*qhh)/N), and qhh = 0, so we can simplify as
                        //Ph = sqrt(Qhh/N)
                        Dictionary<SpectraFileInfo, double> fileToHeavyProbabilityDictionary = new Dictionary<SpectraFileInfo, double>();
                        foreach (SpectraFileInfo info in spectraFileInfo)
                        {
                            //convert the absolute values of LL, LH, and HH to relative values
                            double[] values = fileToRecycleDictionary[info];
                            double sum = values.Sum();
                            for (int i = 0; i < values.Length; i++)
                            {
                                values[i] = values[i] / sum;
                            }
                            double Ph = Math.Sqrt(values[2] / (1 - values[0] + Math.Pow(values[1], 2) / (4 * values[2])));
                            fileToHeavyProbabilityDictionary[info] = Ph;
                        }


                        //foreach group, the labeled peptides should be split into their labeled files
                        //we're deleting the heavy results after we pull those results into a different file
                        //additionally, the Ph will be used to correct both values (L is originally too high, H too low)
                        foreach (SpectraFileInfo info in spectraFileInfo)
                        {
                            string fullPathWithExtension = info.FullFilePathWithExtension;
                            string[] pathArray = fullPathWithExtension.Split('.');
                            string extension = pathArray.Last();
                            string filePathWithoutExtension = fullPathWithExtension.Substring(0, fullPathWithExtension.Length - extension.Length - 1); //-1 removes the '.'
                            SpectraFileInfo lightInfo = new SpectraFileInfo(filePathWithoutExtension + "_Old." + extension, info.Condition, info.BiologicalReplicate, info.TechnicalReplicate, info.Fraction);
                            SpectraFileInfo heavyInfo = new SpectraFileInfo(filePathWithoutExtension + "_New." + extension, info.Condition, info.BiologicalReplicate, info.TechnicalReplicate, info.Fraction);
                            originalToLabeledFileInfoDictionary[info] = new List<SpectraFileInfo> { lightInfo, heavyInfo };
                            Parameters.FlashLfqResults.SpectraFiles.Add(lightInfo);
                            Parameters.FlashLfqResults.SpectraFiles.Add(heavyInfo);
                        }
                        foreach (KeyValuePair<string, List<(SilacLabel label, Peptide peptide)>> kvp in unlabeledToPeptidesDictionary)
                        {
                            string unlabeledSequence = kvp.Key; //this will be the key for the new quant entry
                            List<Peptide> peptides = kvp.Value.Select(x => x.peptide).ToList();
                            //if no missed cleavages
                            Peptide lightPeptide;
                            Peptide mixedPeptide = null; //null assignment needed to build
                            Peptide heavyPeptide;
                            bool mixedPeptideExists = peptides.Count != 2;
                            if (!mixedPeptideExists)
                            {
                                lightPeptide = peptides[0];
                                heavyPeptide = peptides[1];
                            }
                            else
                            {
                                lightPeptide = peptides[0];
                                mixedPeptide = peptides[1];
                                heavyPeptide = peptides[2];
                            }

                            foreach (SpectraFileInfo info in spectraFileInfo)
                            {
                                double Ph = fileToHeavyProbabilityDictionary[info];
                                double light = lightPeptide.GetIntensity(info);
                                double heavy = heavyPeptide.GetIntensity(info);
                                double mixed = mixedPeptideExists ? mixedPeptide.GetIntensity(info) : 0;

                                //all the heavy is new, but some of the light is also new protein
                                //Ph helps out here. The correction factor is Pl*Qh/Ph, or (1-Ph)*Qh/Ph.
                                //it's possible that we obtain a negative value for the start, which doesn't make sense.
                                //This occurs when Qh>Ph. In these cases, set the light to zero and the heavy to light+heavy.
                                double correction = (1 - Ph) * heavy / Ph;
                                var updatedInfo = originalToLabeledFileInfoDictionary[info];
                                SpectraFileInfo startInfo = updatedInfo[0];
                                SpectraFileInfo endInfo = updatedInfo[1];

                                Peptide updatedPeptide = new Peptide(lightPeptide.Sequence, lightPeptide.UseForProteinQuant);
                                updatedPeptide.ProteinGroups = SilacConversions.CleanPastProteinQuant(lightPeptide.ProteinGroups); //needed to keep protein info.

                                updatedPeptide.SetIntensity(startInfo, light - correction); //assign the corrected light intensity
                                updatedPeptide.SetDetectionType(startInfo, lightPeptide.GetDetectionType(info));
                                updatedPeptide.SetIntensity(endInfo, heavy + mixed + correction); //assign the corrected heavy intensity to the heavy file
                                updatedPeptide.SetDetectionType(endInfo, heavyPeptide.GetDetectionType(info)); //could include the mixed here if it really matters
                                updatedPeptides.Add(updatedPeptide);
                            }
                        }
                    }
                    else //multiplex
                    {
                        //foreach original file
                        foreach (SpectraFileInfo info in spectraFileInfo)
                        {
                            foreach (var kvp in unlabeledToPeptidesDictionary)
                            {
                                string unlabeledSequence = kvp.Key;
                                List<Peptide> peptides = kvp.Value.Select(x => x.peptide).ToList();
                                Peptide representativePeptide = peptides[0];
                                Peptide updatedPeptide = new Peptide(unlabeledSequence, representativePeptide.UseForProteinQuant);
                                updatedPeptide.ProteinGroups = SilacConversions.CleanPastProteinQuant(representativePeptide.ProteinGroups); //needed to keep protein info.
                                List<SpectraFileInfo> filesForThisFile = originalToLabeledFileInfoDictionary[info];
                                for (int i = 0; i < peptides.Count; i++) //the files and the peptides can use the same index, because there should be a distinct file for each label/peptide
                                {
                                    SpectraFileInfo currentInfo = filesForThisFile[i];
                                    Peptide currentPeptide = peptides[i];
                                    updatedPeptide.SetIntensity(currentInfo, currentPeptide.GetIntensity(info));
                                    updatedPeptide.SetDetectionType(currentInfo, currentPeptide.GetDetectionType(info));
                                }
                                updatedPeptides.Add(updatedPeptide);
                            }
                        }
                    }
                    //Update peptides
                    var peptideResults = Parameters.FlashLfqResults.PeptideModifiedSequences;
                    peptideResults.Clear();
                    foreach (Peptide peptide in updatedPeptides)
                    {
                        peptideResults.Add(peptide.Sequence, peptide);
                    }

                    //Do protein quant
                    Parameters.FlashLfqResults.ProteinGroups.Clear();
                    Parameters.FlashLfqResults.CalculateProteinResultsTop3();

                    //update proteingroups to have all files for quantification
                    if (ProteinGroups != null)
                    {
                        List<SpectraFileInfo> allInfo = originalToLabeledFileInfoDictionary.SelectMany(x => x.Value).ToList();
                        foreach (EngineLayer.ProteinGroup proteinGroup in ProteinGroups)
                        {
                            proteinGroup.FilesForQuantification = allInfo;
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

                    //Convert all lfqpeaks from heavy (a) to light (K+8.014) for output

                    if (Parameters.FlashLfqResults != null) //can be null if nothing was quantified (all peptides are ambiguous)
                    {
                        var lfqPeaks = Parameters.FlashLfqResults.Peaks;
                        List<SpectraFileInfo> peakKeys = lfqPeaks.Keys.ToList();

                        foreach (SpectraFileInfo key in peakKeys)
                        {
                            List<ChromatographicPeak> peaks = lfqPeaks[key];
                            for (int i = 0; i < peaks.Count; i++)
                            {
                                var peak = peaks[i];
                                List<Identification> identifications = new List<Identification>();
                                //check if we're removing light peaks and if it's a light peak
                                if (peak.Identifications.Any(x => SilacConversions.GetRelevantLabelFromBaseSequence(x.BaseSequence, allSilacLabels) != null)) //if no ids have any labels, remove them
                                {
                                    List<Identification> updatedIds = new List<Identification>();
                                    foreach (var id in peak.Identifications)
                                    {
                                        SilacLabel label = SilacConversions.GetRelevantLabelFromBaseSequence(id.BaseSequence, allSilacLabels);

                                        Identification updatedId = new Identification(
                                            id.FileInfo,
                                            SilacConversions.GetSilacLightBaseSequence(id.BaseSequence, label),
                                            SilacConversions.GetSilacLightFullSequence(id.ModifiedSequence, label),
                                            id.MonoisotopicMass,
                                            id.Ms2RetentionTimeInMinutes,
                                            id.PrecursorChargeState,
                                            id.ProteinGroups.ToList(),
                                            id.OptionalChemicalFormula,
                                            id.UseForProteinQuant
                                            );
                                        updatedIds.Add(updatedId);
                                    }
                                    peak.Identifications.Clear();
                                    peak.Identifications.AddRange(updatedIds);
                                    //TODO delete
                                    //ChromatographicPeak updatedPeak = new ChromatographicPeak(identifications.First(), peak.IsMbrPeak, peak.SpectraFileInfo);
                                    //for (int j = 1; j < identifications.Count; j++) //add all the original identification
                                    //{
                                    //    updatedPeak.MergeFeatureWith(new FlashLFQ.ChromatographicPeak(identifications[j], peak.IsMbrPeak, peak.SpectraFileInfo), FlashLfqEngine.Integrate);
                                    //}
                                    //updatedPeak.IsotopicEnvelopes = peak.IsotopicEnvelopes; //need to set isotopicEnevelopes, since the new identifications didn't have them.
                                    //updatedPeak.CalculateIntensityForThisFeature(FlashLfqEngine.Integrate); //needed to update info
                                    //peaks[i] = updatedPeak;
                                }
                            }
                        }
                    }
                }

                //convert all peaks into human readable format (PEPTIDE(K+8) and not PEPTIDEa)
                //convert all psms into human readable format
                //We can do this for all of the FlashLFQ peptides/peaks, because they use string sequences.
                //We are unable to do this for Parameters.AllPsms, because they store proteins and start/end residues instead
                //for Psms, we need to convert during the writing.
                List<PeptideSpectralMatch> allPsms = Parameters.AllPsms;
                var modsToWritePruned = Parameters.SearchParameters.ModsToWriteSelection;
                for (int i = 0; i < allPsms.Count; i++)
                {
                    allPsms[i].ResolveHeavySilacLabel(allSilacLabels, modsToWritePruned);
                }
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
                //set peptide output values
                ProteinGroups.ForEach(x => x.GetIdentifiedPeptidesOutput(Parameters.SearchParameters.SilacLabels));
                // write protein groups to tsv
                string writtenFile = Path.Combine(Parameters.OutputFolder, "AllProteinGroups.tsv");
                WriteProteinGroupsToTsv(ProteinGroups, writtenFile, new List<string> { Parameters.SearchTaskId }, CommonParameters.QValueOutputFilter);

                // write all individual file results to subdirectory
                // local protein fdr, global parsimony, global psm fdr
                if (Parameters.CurrentRawFileList.Count > 1 || Parameters.SearchParameters.WriteMzId || Parameters.SearchParameters.WritePepXml)
                {
                    Directory.CreateDirectory(Parameters.IndividualResultsOutputFolder);

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

            flashLFQResults.WriteResults(null, fullSeqPath, null,null);

            FinishedWritingFile(fullSeqPath, nestedIds);
        }

        private void WritePeakQuantificationResultsToTsv(FlashLfqResults flashLFQResults, string outputFolder, string fileName, List<string> nestedIds)
        {
            var peaksPath = Path.Combine(outputFolder, fileName + ".tsv");

            flashLFQResults.WriteResults(peaksPath, null, null,null);

            FinishedWritingFile(peaksPath, nestedIds);
        }
    }
}