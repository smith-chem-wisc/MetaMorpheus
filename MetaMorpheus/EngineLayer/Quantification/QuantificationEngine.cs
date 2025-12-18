using EngineLayer.SpectrumMatch;
using FlashLFQ;
using MzLibUtil;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.Quantification;

public class QuantificationEngine(
    QuantificationParameters parameters,
    List<SpectralMatch> allSpectralMatches,
    List<ProteinGroup> proteinGroups,
    CommonParameters commonParameters,
    List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
    List<string> nestedIds)
    : MetaMorpheusEngine(commonParameters, fileSpecificParameters, nestedIds)
{
    public QuantificationParameters Parameters { get; init; } = parameters;
    public List<SpectralMatch> AllSpectralMatches { get; init; } = allSpectralMatches;
    public List<ProteinGroup> ProteinGroups { get; init; } = proteinGroups;

    protected override MetaMorpheusEngineResults RunSpecific()
    {
        var results = new QuantificationEngineResults(this);

        var spectraFileInfos = GetSpectraFileInfos();

        // Issue loading experimental design file or normalization requested but no design file found
        if (spectraFileInfos.IsNullOrEmpty())
            return results;

        // get PSMs to pass to FlashLFQ
        var psmsForQuantification = FilteredPsms.Filter(AllSpectralMatches,
            CommonParameters,
            includeDecoys: Parameters.MatchBetweenRuns, // Decoys are required for PIP-ECHO, but are not written to the output file
            includeContaminants: true,
            includeAmbiguous: false,
            includeAmbiguousMods: false,
            includeHighQValuePsms: false);

        // Only these peptides will be written to the AllQuantifiedPeptides.tsv output file
        var peptideSequencesForQuantification = FilteredPsms.Filter(AllSpectralMatches,
            CommonParameters,
            includeDecoys: false,
            includeContaminants: true,
            includeAmbiguous: false,
            includeAmbiguousMods: false,
            includeHighQValuePsms: false,
            filterAtPeptideLevel: true).Select(p => p.FullSequence).ToList();

        var psmToProteinGroups = MapPsmsToProteinGroups(ProteinGroups, psmsForQuantification);

        // Silac adjustments
        List<SilacLabel> allSilacLabels = Parameters.SilacLabels;
        SilacLabel startLabel = Parameters.StartTurnoverLabel;
        SilacLabel endLabel = Parameters.EndTurnoverLabel;
        bool quantifyUnlabeledPeptides = Parameters.ListOfDigestionParams.Any(x => x is DigestionParams { GeneratehUnlabeledProteinsForSilac: true });
        if (Parameters.SilacLabels is { Count: > 0 })
        {
            MakePreQuantificationSilacAdjustments(quantifyUnlabeledPeptides, Parameters.SilacLabels, Parameters.StartTurnoverLabel, Parameters.EndTurnoverLabel, psmsForQuantification, ref psmToProteinGroups);
        }

        //group psms by file
        Dictionary<string, List<SpectralMatch>> psmsGroupedByFile = psmsForQuantification.GroupBy(p => p.FullFilePath)
            .ToDictionary(p => p.Key, p => p.ToList());
        Dictionary<DigestionAgent, List<SpectralMatch>> proteaseSortedPsms = GetPsmsGroupedByDigestionAgent(Parameters.ListOfDigestionParams.ToList(), psmsForQuantification, psmToProteinGroups);

        List<Identification> flashLFQIdentifications = GetFlashLFQIdentifications(spectraFileInfos, psmsForQuantification, psmsGroupedByFile, psmToProteinGroups);

        // run FlashLFQ
        var flashLfqEngine = new FlashLfqEngine(
            allIdentifications: flashLFQIdentifications,
            normalize: Parameters.Normalize,
            ppmTolerance: Parameters.QuantifyPpmTol,
            matchBetweenRunsPpmTolerance: Parameters.QuantifyPpmTol,  // If these tolerances are not equivalent, then MBR will falsely classify peptides found in the initial search as MBR peaks
            matchBetweenRuns: Parameters.MatchBetweenRuns,
            matchBetweenRunsFdrThreshold: Parameters.MbrFdrThreshold,
            useSharedPeptidesForProteinQuant: Parameters.UseSharedPeptides,
            peptideSequencesToQuantify: Parameters.SilacLabels == null ? peptideSequencesForQuantification : null, // Silac is doing it's own thing, no need to pass in peptide sequences
            silent: true,
            maxThreads: CommonParameters.MaxThreadsToUsePerFile);

        if (flashLFQIdentifications.Any())
        {
            results.FlashLfqResults = flashLfqEngine.Run();
        }

        // get protein intensity back from FlashLFQ
        if (ProteinGroups != null && results.FlashLfqResults != null)
        {
            foreach (var proteinGroup in ProteinGroups)
            {
                proteinGroup.FilesForQuantification = spectraFileInfos;
                proteinGroup.IntensitiesByFile = new Dictionary<SpectraFileInfo, double>();

                foreach (var spectraFile in proteinGroup.FilesForQuantification)
                {
                    if (results.FlashLfqResults.ProteinGroups.TryGetValue(proteinGroup.ProteinGroupName, out var flashLfqProteinGroup))
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
        if (Parameters.SilacLabels != null && AllSpectralMatches.First() is PeptideSpectralMatch) //if we're doing silac
        {
            SilacConversions.SilacConversionsPostQuantification(allSilacLabels, startLabel, endLabel, spectraFileInfos, ProteinGroups, Parameters.ListOfDigestionParams,
                results.FlashLfqResults, AllSpectralMatches.Cast<PeptideSpectralMatch>().ToList(), Parameters.ModsToWriteSelection, quantifyUnlabeledPeptides);
        }


        return results;
    }


    public List<SpectraFileInfo> GetSpectraFileInfos()
    {
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
                return null;
            }
        }
        else if (Parameters.Normalize)
        {
            Warn("Could not find experimental design file at " + assumedExperimentalDesignPath + ", which is required for normalization. Skipping quantification");
            return null;
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
        return spectraFileInfo;
    }

    public static List<Identification> GetFlashLFQIdentifications(List<SpectraFileInfo> spectraFileInfos, FilteredPsms psmsForQuantification, Dictionary<string, List<SpectralMatch>> psmsGroupedByFile, Dictionary<SpectralMatch, List<FlashLFQ.ProteinGroup>> psmToProteinGroups)
    {
        // pass PSM info to FlashLFQ
        var flashLFQIdentifications = new List<Identification>();
        foreach (var spectraFile in psmsGroupedByFile)
        {
            var rawfileinfo = spectraFileInfos.First(p => p.FullFilePathWithExtension.Equals(spectraFile.Key));

            foreach (var psm in spectraFile.Value)
            {
                flashLFQIdentifications.Add(
                    new Identification(
                        fileInfo: rawfileinfo,
                        psm.BaseSequence,
                        psm.FullSequence,
                        psm.BioPolymerWithSetModsMonoisotopicMass.Value,
                        psm.ScanRetentionTime,
                        psm.ScanPrecursorCharge,
                        psmToProteinGroups[psm],
                        psmScore: psm.Score,
                        qValue: psmsForQuantification.FilterType == FilterType.QValue ? psm.FdrInfo.QValue : psm.FdrInfo.PEP_QValue,
                        decoy: psm.IsDecoy));
            }
        }
        return flashLFQIdentifications;
    }

    public static Dictionary<DigestionAgent, List<SpectralMatch>> GetPsmsGroupedByDigestionAgent(List<IDigestionParams> listOfDigestionParameters, FilteredPsms psmsForQuantification, Dictionary<SpectralMatch, List<FlashLFQ.ProteinGroup>> psmToProteinGroups)
    {
        // some PSMs may not have protein groups (if 2 peptides are required to construct a protein group, some PSMs will be left over)
        // the peptides should still be quantified but not considered for protein quantification
        var undefinedPg = new FlashLFQ.ProteinGroup("UNDEFINED", "", "");
        //sort the unambiguous psms by protease to make MBR compatible with multiple proteases
        Dictionary<DigestionAgent, List<SpectralMatch>> proteaseSortedPsms = new Dictionary<DigestionAgent, List<SpectralMatch>>();

        foreach (IDigestionParams dp in listOfDigestionParameters)
        {
            if (!proteaseSortedPsms.ContainsKey(dp.DigestionAgent))
            {
                proteaseSortedPsms.Add(dp.DigestionAgent, new List<SpectralMatch>());
            }
        }
        foreach (var psm in psmsForQuantification)
        {
            if (!psmToProteinGroups.ContainsKey(psm))
            {
                psmToProteinGroups.Add(psm, new List<FlashLFQ.ProteinGroup> { undefinedPg });
            }

            proteaseSortedPsms[psm.DigestionParams.DigestionAgent].Add(psm);
        }

        return proteaseSortedPsms;
    }


    public static Dictionary<SpectralMatch, List<FlashLFQ.ProteinGroup>> MapPsmsToProteinGroups(List<EngineLayer.ProteinGroup> inputProteinGroups, FilteredPsms smsForQuantification)
    {
        var psmToProteinGroups = new Dictionary<SpectralMatch, List<FlashLFQ.ProteinGroup>>();
        if (inputProteinGroups != null && inputProteinGroups.Count != 0) //ProteinGroups can be null if parsimony wasn't done, and it can be empty if you're doing the two peptide rule
        {
            foreach (var proteinGroup in inputProteinGroups)
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
            foreach (var psm in smsForQuantification)
            {
                var proteins = psm.BestMatchingBioPolymersWithSetMods.Select(b => b.SpecificBioPolymer.Parent).Distinct();

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
        return psmToProteinGroups;
    }

    public void MakePreQuantificationSilacAdjustments(bool quantifyUnlabeledPeptides, List<SilacLabel> allSilacLabels, SilacLabel startLabel, SilacLabel endLabel, FilteredPsms psmsForQuantification, ref Dictionary<SpectralMatch, List<FlashLFQ.ProteinGroup>> psmToProteinGroups)
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
            PeptideWithSetModifications pwsm = psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer as PeptideWithSetModifications;

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
}
