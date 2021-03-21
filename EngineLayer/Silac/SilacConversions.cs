using FlashLFQ;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public static class SilacConversions
    {
        private static readonly string LABEL_DELIMITER = " & ";

        public static string GetLabeledBaseSequence(string unlabeledBaseSequence, SilacLabel label)
        {
            if (label != null)
            {
                string labeledBaseSequence = unlabeledBaseSequence.Replace(label.OriginalAminoAcid, label.AminoAcidLabel);
                if (label.AdditionalLabels != null)
                {
                    foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                    {
                        labeledBaseSequence = labeledBaseSequence.Replace(additionalLabel.OriginalAminoAcid, additionalLabel.AminoAcidLabel);
                    }
                }
                return labeledBaseSequence;
            }
            else
            {
                return unlabeledBaseSequence;
            }
        }

        public static PeptideSpectralMatch GetLabeledPsm(PeptideSpectralMatch psm, int notch, PeptideWithSetModifications pwsm, string labeledBaseSequence)
        {
            PeptideWithSetModifications labeledPwsm = new PeptideWithSetModifications(
                pwsm.Protein,
                pwsm.DigestionParams,
                pwsm.OneBasedStartResidueInProtein,
                pwsm.OneBasedEndResidueInProtein,
                pwsm.CleavageSpecificityForFdrCategory,
                pwsm.PeptideDescription,
                pwsm.MissedCleavages,
                pwsm.AllModsOneIsNterminus,
                pwsm.NumFixedMods,
                labeledBaseSequence);
            return psm.Clone(new List<(int Notch, PeptideWithSetModifications Peptide)> { (notch, labeledPwsm) });
        }

        public static PeptideSpectralMatch GetSilacPsm(PeptideSpectralMatch psm, SilacLabel silacLabel)
        {
            List<(int Notch, PeptideWithSetModifications Peptide)> updatedBestMatchingPeptides = new List<(int Notch, PeptideWithSetModifications Peptide)>();
            foreach ((int Notch, PeptideWithSetModifications Peptide) notchAndPwsm in psm.BestMatchingPeptides)
            {
                PeptideWithSetModifications modifiedPwsm = CreateSilacPwsm(silacLabel, notchAndPwsm.Peptide);
                updatedBestMatchingPeptides.Add((notchAndPwsm.Notch, modifiedPwsm));
            }
            return psm.Clone(updatedBestMatchingPeptides);
        }

        //modify the proteins to appear only light (we want a protein sequence to look like PROTEINK instead of PROTEINa)
        public static List<PeptideSpectralMatch> UpdateProteinSequencesToLight(List<PeptideSpectralMatch> originalPsms, List<SilacLabel> labels)
        {
            List<PeptideSpectralMatch> psmsToReturn = new List<PeptideSpectralMatch>();
            foreach (PeptideSpectralMatch psm in originalPsms)
            {
                List<(int Notch, PeptideWithSetModifications Peptide)> originalPeptides = psm.BestMatchingPeptides.ToList();
                List<(int Notch, PeptideWithSetModifications Peptide)> updatedPeptides = new List<(int Notch, PeptideWithSetModifications Peptide)>();
                foreach ((int Notch, PeptideWithSetModifications Peptide) notchPwsm in originalPeptides)
                {
                    PeptideWithSetModifications pwsm = notchPwsm.Peptide;
                    SilacLabel label = GetRelevantLabelFromBaseSequence(pwsm.BaseSequence, labels);
                    Protein updatedProtein = pwsm.Protein;
                    if (label != null)
                    {
                        string proteinLightSequence = updatedProtein.BaseSequence;
                        proteinLightSequence = proteinLightSequence.Replace(label.AminoAcidLabel, label.OriginalAminoAcid);
                        if (label.AdditionalLabels != null)
                        {
                            foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                            {
                                proteinLightSequence = proteinLightSequence.Replace(additionalLabel.AminoAcidLabel, additionalLabel.OriginalAminoAcid);
                            }
                        }
                        updatedProtein = new Protein(pwsm.Protein, proteinLightSequence);
                    }
                    PeptideWithSetModifications updatedPwsm = new PeptideWithSetModifications(
                        updatedProtein,
                        pwsm.DigestionParams,
                        pwsm.OneBasedStartResidueInProtein,
                        pwsm.OneBasedEndResidueInProtein,
                        pwsm.CleavageSpecificityForFdrCategory,
                        pwsm.PeptideDescription,
                        pwsm.MissedCleavages,
                        pwsm.AllModsOneIsNterminus,
                        pwsm.NumFixedMods,
                        pwsm.BaseSequence);
                    updatedPeptides.Add((notchPwsm.Notch, updatedPwsm));
                }

                psmsToReturn.Add(psm.Clone(updatedPeptides));
            }

            return psmsToReturn;
        }

        //Converts the heavy char label "a" into a human readable label "K+8.014"
        public static string GetSilacLightBaseSequence(string baseSequence, SilacLabel label)
        {
            if (label != null)
            {
                baseSequence = baseSequence.Replace(label.AminoAcidLabel.ToString(), HeavyStringForPeptides(label));
                if (label.AdditionalLabels != null)
                {
                    foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                    {
                        baseSequence = baseSequence.Replace(additionalLabel.AminoAcidLabel.ToString(), HeavyStringForPeptides(additionalLabel));
                    }
                }
            }
            return baseSequence;
        }

        //Converts the heavy char label "a" into a human readable label "K+8.014", but doesn't change mods containing the char like "Oxid'a'tion"
        public static string GetSilacLightFullSequence(string fullSequence, SilacLabel label, bool includeMassDifference = true)
        {
            //overwrite full sequence
            if (label != null)
            {
                List<SilacLabel> labels = new List<SilacLabel> { label };
                if (label.AdditionalLabels != null)
                {
                    labels.AddRange(label.AdditionalLabels);
                }

                foreach (SilacLabel additionalLabel in labels)
                {
                    string replacementSequence = includeMassDifference ? HeavyStringForPeptides(additionalLabel) : additionalLabel.OriginalAminoAcid.ToString();

                    bool inModification = false;
                    for (int i = 0; i < fullSequence.Length; i++)
                    {
                        if (inModification)
                        {
                            if (fullSequence[i] == ']')
                            {
                                inModification = false;
                            }
                        }
                        else
                        {
                            char currentChar = fullSequence[i];
                            if (currentChar == '[')
                            {
                                inModification = true;
                            }
                            else if (currentChar == additionalLabel.AminoAcidLabel)
                            {
                                fullSequence = fullSequence.Substring(0, i) + replacementSequence + fullSequence.Substring(i + 1, fullSequence.Length - i - 1);
                                i += replacementSequence.Length - 1; //-1 because we removed the label amino acid
                            }
                        }
                    }
                }
            }
            return fullSequence;
        }

        public static string GetAmbiguousLightSequence(string originalSequence, List<SilacLabel> allLabels, bool baseSequence)
        {
            string[] multipleSequences = originalSequence.Split('|').ToArray(); //if ambiguity
            string localSequence = "";
            foreach (string sequence in multipleSequences)
            {
                List<SilacLabel> labels = GetRelevantLabelsFromBaseSequenceForOutput(sequence, allLabels);
                string updatedSequence = sequence;
                if (labels != null)
                {
                    foreach (SilacLabel label in labels)
                    {
                        updatedSequence = baseSequence ? GetSilacLightBaseSequence(updatedSequence, label) : GetSilacLightFullSequence(updatedSequence, label);
                    }
                }
                localSequence += updatedSequence + "|";
            }
            if (localSequence.Length != 0)
            {
                localSequence = localSequence.Substring(0, localSequence.Length - 1); //remove last "|"
            }
            return localSequence;
        }

        public static string HeavyStringForPeptides(SilacLabel label)
        {
            return label.OriginalAminoAcid + "(" + label.MassDifference + ")";
        }

        public static SilacLabel GetRelevantLabelFromBaseSequence(string baseSequence, List<SilacLabel> labels)
        {
            return labels.Where(x => baseSequence.Contains(x.AminoAcidLabel) ||
                (x.AdditionalLabels != null && x.AdditionalLabels.Any(y => baseSequence.Contains(y.AminoAcidLabel)))).FirstOrDefault();
        }

        //This method is used specifically for turnover experiments where the start and the end conditions are labeled
        public static List<SilacLabel> GetRelevantLabelsFromBaseSequenceForOutput(string baseSequence, List<SilacLabel> labels)
        {
            return labels.Where(x => baseSequence.Contains(x.AminoAcidLabel) ||
                (x.AdditionalLabels != null && x.AdditionalLabels.Any(y => baseSequence.Contains(y.AminoAcidLabel)))).ToList();
        }

        public static SilacLabel GetRelevantLabelFromFullSequence(string fullSequence, List<SilacLabel> labels)
        {
            bool inModification = false;
            List<char> labelResidues = labels.Select(x => x.AminoAcidLabel).ToList();
            for (int i = 0; i < fullSequence.Length; i++)
            {
                if (inModification)
                {
                    if (fullSequence[i] == ']')
                    {
                        inModification = false;
                    }
                }
                else
                {
                    char currentChar = fullSequence[i];
                    if (currentChar == '[')
                    {
                        inModification = true;
                    }
                    else if (labelResidues.Contains(currentChar))
                    {
                        return labels.Where(x => currentChar == x.AminoAcidLabel).First();
                    }
                }
            }
            //if nothing
            return null;
        }

        public static PeptideWithSetModifications CreateSilacPwsm(SilacLabel silacLabel, PeptideWithSetModifications pwsm)
        {
            string baseSequence = pwsm.BaseSequence;

            baseSequence = baseSequence.Replace(silacLabel.AminoAcidLabel, silacLabel.OriginalAminoAcid); //create light sequence
            if (silacLabel.AdditionalLabels != null)
            {
                foreach (SilacLabel additionalLabel in silacLabel.AdditionalLabels)
                {
                    baseSequence = baseSequence.Replace(additionalLabel.AminoAcidLabel, additionalLabel.OriginalAminoAcid); //create light sequence
                }
            }

            return new PeptideWithSetModifications(
                pwsm.Protein,
                pwsm.DigestionParams,
                pwsm.OneBasedStartResidueInProtein,
                pwsm.OneBasedEndResidueInProtein,
                pwsm.CleavageSpecificityForFdrCategory,
                pwsm.PeptideDescription,
                pwsm.MissedCleavages,
                pwsm.AllModsOneIsNterminus,
                pwsm.NumFixedMods,
                baseSequence); //this is the only thing changing
        }

        public static SpectraFileInfo GetHeavyFileInfo(SpectraFileInfo originalFile, SilacLabel label)
        {
            string heavyFileName = originalFile.FilenameWithoutExtension + "(" + label.OriginalAminoAcid + label.MassDifference;
            string heavyCondition = originalFile.Condition + "(" + label.OriginalAminoAcid + label.MassDifference;
            if (label.AdditionalLabels != null)
            {
                foreach (SilacLabel additionaLabel in label.AdditionalLabels)
                {
                    heavyFileName += LABEL_DELIMITER + additionaLabel.OriginalAminoAcid + additionaLabel.MassDifference;
                    heavyCondition += LABEL_DELIMITER + additionaLabel.OriginalAminoAcid + additionaLabel.MassDifference;
                }
            }
            heavyFileName += ")." + originalFile.FullFilePathWithExtension.Split('.').Last(); //add extension
            heavyCondition += ")";

            return new SpectraFileInfo(heavyFileName, heavyCondition, originalFile.BiologicalReplicate, originalFile.TechnicalReplicate, originalFile.Fraction);
        }

        public static HashSet<FlashLFQ.ProteinGroup> CleanPastProteinQuant(HashSet<FlashLFQ.ProteinGroup> originalProteinGroups)
        {
            HashSet<FlashLFQ.ProteinGroup> cleanedProteinGroups = new HashSet<FlashLFQ.ProteinGroup>();
            foreach (FlashLFQ.ProteinGroup pg in originalProteinGroups)
            {
                cleanedProteinGroups.Add(new FlashLFQ.ProteinGroup(pg.ProteinGroupName, pg.GeneName, pg.Organism));
            }
            return cleanedProteinGroups;
        }

        //If SILAC (Post-Quantification), compress the light/heavy protein group pairs into the same light protein group but different files
        //Create new files for each silac label and file so that "file 1" now becomes "file 1 (light)" and "file 1 (heavy)"
        //Change heavy residue into the light residue plus a string label ("PEPTIDEa" -> "PEPTIDEK(+8.014)")
        //This light to heavy conversion needs to happen for the flashLFQ peptides here, but can't for the psm peptides, which are constrained to the protein
        //i.e. pwsms currently don't have sequences; they have start/end residues and a protein sequence. We have to change the output sequences when they're created.
        public static void SilacConversionsPostQuantification(List<SilacLabel> allSilacLabels, SilacLabel startLabel, SilacLabel endLabel,
            List<SpectraFileInfo> spectraFileInfo, List<ProteinGroup> proteinGroups, HashSet<DigestionParams> listOfDigestionParams, FlashLfqResults flashLfqResults,
            List<PeptideSpectralMatch> allPsms, Dictionary<string, int> modsToWriteSelection, bool quantifyUnlabeledPeptides)
        {
            //do protein quant if we had any results
            //if no results, we still may need to edit the psms
            if (flashLfqResults != null) //can be null if no unambiguous psms were found
            {
                //after this point, we now have quantification values for the peptides, but they all belong to the same "unlabeled" protein and are in the same file
                //We can remove "labeled" peptides from each file and put them in a new file as "unlabeled".

                //MAKE NEW RAW FILES

                //update number of spectra files to include a new file for each label/condition
                Dictionary<SpectraFileInfo, List<SpectraFileInfo>> originalToLabeledFileInfoDictionary = CreateSilacRawFiles(flashLfqResults, allSilacLabels, startLabel, endLabel, quantifyUnlabeledPeptides, spectraFileInfo);

                //we have the files, now let's reassign the psms.
                //there are a few ways to do this, but we're going to generate the "base" peptide and assign to that

                //Get Dictionary of protein accessions to peptides
                Dictionary<string, List<FlashLFQ.Peptide>> unlabeledToPeptidesDictionary = GetDictionaryOfProteinAccessionsToPeptides(flashLfqResults.PeptideModifiedSequences.Values, allSilacLabels, startLabel, endLabel);

                //we now have a dictionary of unlabeledBaseSequence to the labeled peptides
                //Better SILAC results can be obtained by using the summed intensities from ms1 scans where all peaks were found, rather than the apex
                //foreach peptide, unlabeled peptide, get the isotopic envelope intensities for each labeled peptide in each file
                //save the intensities from ms1s that are shared. If no ms1s contains all the peaks, then just use the apex intensity (default)
                CalculateSilacIntensities(flashLfqResults.Peaks, unlabeledToPeptidesDictionary);


                //SPLIT THE FILES
                List<FlashLFQ.Peptide> updatedPeptides = new List<FlashLFQ.Peptide>();

                //split the heavy/light peptides into separate raw files, remove the heavy peptide
                if (startLabel != null || endLabel != null) //if turnover
                {
                    //foreach group, the labeled peptides should be split into their labeled files
                    //we're deleting the heavy results after we pull those results into a different file
                    foreach (SpectraFileInfo info in spectraFileInfo)
                    {
                        string fullPathWithExtension = info.FullFilePathWithExtension;
                        string[] pathArray = fullPathWithExtension.Split('.');
                        string extension = pathArray.Last();
                        string filePathWithoutExtension = fullPathWithExtension.Substring(0, fullPathWithExtension.Length - extension.Length - 1); //-1 removes the '.'
                        SpectraFileInfo lightInfo = new SpectraFileInfo(filePathWithoutExtension + "_Original." + extension, info.Condition, info.BiologicalReplicate, info.TechnicalReplicate, info.Fraction);
                        SpectraFileInfo heavyInfo = new SpectraFileInfo(filePathWithoutExtension + "_NewlySynthesized." + extension, info.Condition + "_NewlySynthesized", info.BiologicalReplicate, info.TechnicalReplicate, info.Fraction);
                        originalToLabeledFileInfoDictionary[info] = new List<SpectraFileInfo> { lightInfo, heavyInfo };
                        flashLfqResults.SpectraFiles.Add(lightInfo);
                        flashLfqResults.SpectraFiles.Add(heavyInfo);
                    }

                    //This step converts the quantification intensities from light/heavy to original/newlySynthesized by splitting up the missed cleavage mixtures
                    foreach (KeyValuePair<string, List<FlashLFQ.Peptide>> kvp in unlabeledToPeptidesDictionary)
                    {
                        string unlabeledSequence = kvp.Key; //this will be the key for the new quant entry
                        List<FlashLFQ.Peptide> peptides = kvp.Value;
                        if (peptides.Count != 1) //sometimes it's one if there is no label site on the peptide (e.g. label K, peptide is PEPTIDER)
                        {
                            //Missed cleavages can yield multiple peptides (e.g. 1 missed = LL, LH, HH; 2 missed = LLL, LLH, LHH, HHH; etc)
                            //Compress into 2 values: Light and Heavy
                            FlashLFQ.Peptide updatedPeptide = new FlashLFQ.Peptide(unlabeledSequence, unlabeledSequence, peptides[0].UseForProteinQuant, CleanPastProteinQuant(peptides[0].ProteinGroups)); //needed to keep protein info.
                            foreach (SpectraFileInfo info in spectraFileInfo)
                            {
                                int maxNumberHeavyAminoAcids = peptides.Count - 1;
                                double lightIntensity = 0;
                                double heavyIntensity = 0;
                                int numUniquePeptidesQuantified = 0; 
                                for (int numHeavyAminoAcids = 0; numHeavyAminoAcids < peptides.Count; numHeavyAminoAcids++)
                                {
                                    double totalIntensity = peptides[numHeavyAminoAcids].GetIntensity(info);
                                    if (totalIntensity > 0)
                                    {
                                        //prevent confidence of a ratio if only the HL (and not the LL or HH) is observed. 
                                        //If LL or HH is observed (but not any other), the user knows the ratio is only from one peak.
                                        if (numHeavyAminoAcids == 0 || numHeavyAminoAcids == maxNumberHeavyAminoAcids)
                                        {
                                            numUniquePeptidesQuantified += 2;
                                        }
                                        else
                                        {
                                            numUniquePeptidesQuantified++;
                                        }
                                        double partHeavyIntensity = totalIntensity * numHeavyAminoAcids / maxNumberHeavyAminoAcids;
                                        lightIntensity += totalIntensity - partHeavyIntensity;
                                        heavyIntensity += partHeavyIntensity;
                                    }
                                }
                                //If only a mixed peptide with a missed cleavage was identified, reset the intensity values to zero so the user doesn't get a discreet, inaccurate measurement
                                if (numUniquePeptidesQuantified < 2)
                                {
                                    lightIntensity = 0;
                                    heavyIntensity = 0;
                                }

                                List<SpectraFileInfo> updatedInfo = originalToLabeledFileInfoDictionary[info];
                                SpectraFileInfo startInfo = updatedInfo[0];
                                SpectraFileInfo endInfo = updatedInfo[1];

                                updatedPeptide.SetIntensity(startInfo, lightIntensity); //assign the corrected light intensity
                                updatedPeptide.SetDetectionType(startInfo, peptides.First().GetDetectionType(info));
                                updatedPeptide.SetIntensity(endInfo, heavyIntensity); //assign the corrected heavy intensity to the heavy file
                                updatedPeptide.SetDetectionType(endInfo, peptides.Last().GetDetectionType(info)); //could include the mixed here if it really matters
                            }

                            //add the updated peptide to the list
                            updatedPeptides.Add(updatedPeptide);
                        }
                        else
                        {
                            updatedPeptides.Add(peptides[0]);
                        }
                    }
                }
                else //multiplex
                {
                    foreach (var kvp in unlabeledToPeptidesDictionary)
                    {
                        string unlabeledSequence = kvp.Key;
                        List<FlashLFQ.Peptide> peptides = kvp.Value;
                        FlashLFQ.Peptide representativePeptide = peptides[0];
                        FlashLFQ.Peptide updatedPeptide = new FlashLFQ.Peptide(unlabeledSequence, unlabeledSequence, representativePeptide.UseForProteinQuant, CleanPastProteinQuant(representativePeptide.ProteinGroups)); //needed to keep protein info.

                        //foreach original file
                        foreach (SpectraFileInfo info in spectraFileInfo)
                        {
                            List<SpectraFileInfo> filesForThisFile = originalToLabeledFileInfoDictionary[info];
                            for (int i = 0; i < peptides.Count; i++) //the files and the peptides can use the same index, because there should be a distinct file for each label/peptide
                            {
                                SpectraFileInfo currentInfo = filesForThisFile[i];
                                FlashLFQ.Peptide currentPeptide = peptides[i];
                                updatedPeptide.SetIntensity(currentInfo, currentPeptide.GetIntensity(info));
                                updatedPeptide.SetDetectionType(currentInfo, currentPeptide.GetDetectionType(info));
                            }
                        }
                        updatedPeptides.Add(updatedPeptide);
                    }
                }

                //Update peptides
                var peptideResults = flashLfqResults.PeptideModifiedSequences;
                peptideResults.Clear();
                foreach (FlashLFQ.Peptide peptide in updatedPeptides)
                {
                    peptideResults.Add(peptide.Sequence, peptide);
                }

                //Do protein quant
                flashLfqResults.CalculateProteinResultsMedianPolish(true);

                //update proteingroups to have all files for quantification
                if (proteinGroups != null)
                {
                    List<SpectraFileInfo> allInfo = originalToLabeledFileInfoDictionary.SelectMany(x => x.Value).ToList();
                    foreach (ProteinGroup proteinGroup in proteinGroups)
                    {
                        proteinGroup.FilesForQuantification = allInfo;
                        proteinGroup.IntensitiesByFile = new Dictionary<SpectraFileInfo, double>();

                        foreach (var spectraFile in allInfo)
                        {
                            if (flashLfqResults.ProteinGroups.TryGetValue(proteinGroup.ProteinGroupName, out var flashLfqProteinGroup))
                            {
                                proteinGroup.IntensitiesByFile.Add(spectraFile, flashLfqProteinGroup.GetIntensity(spectraFile));
                            }
                            else
                            {
                                //needed for decoys/contaminants/proteins that aren't quantified
                                proteinGroup.IntensitiesByFile.Add(spectraFile, 0);
                            }
                        }
                    }
                }

                //Convert all lfqpeaks from heavy (a) to light (K+8.014) for output
                if (flashLfqResults != null) //can be null if nothing was quantified (all peptides are ambiguous)
                {
                    var lfqPeaks = flashLfqResults.Peaks;
                    List<SpectraFileInfo> peakKeys = lfqPeaks.Keys.ToList();

                    foreach (SpectraFileInfo key in peakKeys)
                    {
                        List<ChromatographicPeak> peaks = lfqPeaks[key];
                        for (int i = 0; i < peaks.Count; i++)
                        {
                            var peak = peaks[i];
                            //check if we're removing light peaks and if it's a light peak
                            if (peak.Identifications.Any(x => GetRelevantLabelFromBaseSequence(x.BaseSequence, allSilacLabels) != null)) //if no ids have any labels, remove them
                            {
                                List<Identification> updatedIds = new List<Identification>();
                                foreach (var id in peak.Identifications)
                                {
                                    string baseSequence = id.BaseSequence;
                                    string fullSequence = id.ModifiedSequence;
                                    List<SilacLabel> labels = GetRelevantLabelsFromBaseSequenceForOutput(id.BaseSequence, allSilacLabels);
                                    if (labels != null)
                                    {
                                        foreach (SilacLabel label in labels)
                                        {
                                            baseSequence = GetSilacLightBaseSequence(baseSequence, label);
                                            fullSequence = GetSilacLightFullSequence(fullSequence, label);
                                        }
                                    }

                                    Identification updatedId = new Identification(
                                        id.FileInfo,
                                        baseSequence,
                                        fullSequence,
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
                            }
                        }
                    }
                }
            }

            //convert all psms into human readable format
            for (int i = 0; i < allPsms.Count; i++)
            {
                allPsms[i].ResolveHeavySilacLabel(allSilacLabels, modsToWriteSelection);
            }
        }

        private static Dictionary<SpectraFileInfo, List<SpectraFileInfo>> CreateSilacRawFiles(FlashLfqResults flashLfqResults, List<SilacLabel> allSilacLabels, SilacLabel startLabel, SilacLabel endLabel, bool quantifyUnlabeledPeptides, List<SpectraFileInfo> spectraFileInfo)
        {
            //update number of spectra files to include a new file for each label*condition
            Dictionary<SpectraFileInfo, List<SpectraFileInfo>> originalToLabeledFileInfoDictionary = new Dictionary<SpectraFileInfo, List<SpectraFileInfo>>();

            flashLfqResults.SpectraFiles.Clear(); //clear existing files so we can replace them with labeled ones

            //foreach existing file
            if (startLabel == null && endLabel == null) //if multiplex
            {
                //populate dictionary
                if (quantifyUnlabeledPeptides)
                {
                    spectraFileInfo.ForEach(x => originalToLabeledFileInfoDictionary.Add(x, new List<SpectraFileInfo> { x }));
                    flashLfqResults.SpectraFiles.AddRange(spectraFileInfo);
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
                        SpectraFileInfo labeledInfo = GetHeavyFileInfo(originalFile, label);
                        labeledFiles.Add(labeledInfo);
                        originalToLabeledFileInfoDictionary[originalFile].Add(labeledInfo);
                    }
                    flashLfqResults.SpectraFiles.AddRange(labeledFiles);
                }
            }
            return originalToLabeledFileInfoDictionary;
        }

        private static Dictionary<string, List<FlashLFQ.Peptide>> GetDictionaryOfProteinAccessionsToPeptides(IEnumerable<FlashLFQ.Peptide> lfqPwsms, List<SilacLabel> allSilacLabels, SilacLabel startLabel, SilacLabel endLabel)
        {
            Dictionary<string, List<FlashLFQ.Peptide>> unlabeledToPeptidesDictionary = new Dictionary<string, List<FlashLFQ.Peptide>>();
            if (startLabel == null || endLabel == null) //if multiplex, or if at least one of the turnover labels is unlabeled
            {
                foreach (FlashLFQ.Peptide peptide in lfqPwsms)
                {
                    //convert to the unlabeled sequence
                    string labeledSequence = peptide.Sequence;
                    SilacLabel label = GetRelevantLabelFromFullSequence(labeledSequence, allSilacLabels);
                    string unlabeledSequence = GetSilacLightFullSequence(labeledSequence, label, false);
                    if (unlabeledToPeptidesDictionary.ContainsKey(unlabeledSequence))
                    {
                        unlabeledToPeptidesDictionary[unlabeledSequence].Add(peptide);
                    }
                    else
                    {
                        unlabeledToPeptidesDictionary.Add(unlabeledSequence, new List<FlashLFQ.Peptide> { peptide });
                    }
                }
            }
            else //if both start and end labels exist for turnover, then it can be trickier to find the start label (if you're flooded enough to go from label A to label B)...
            {
                foreach (FlashLFQ.Peptide peptide in lfqPwsms)
                {
                    string originalSequence = peptide.Sequence;
                    //convert to the unlabeled sequence
                    string partiallyCleanedSequence = GetSilacLightFullSequence(originalSequence, endLabel, false);
                    string fullyCleanedSequence = GetSilacLightFullSequence(partiallyCleanedSequence, startLabel, false);

                    if (unlabeledToPeptidesDictionary.ContainsKey(fullyCleanedSequence))
                    {
                        unlabeledToPeptidesDictionary[fullyCleanedSequence].Add(peptide);
                    }
                    else
                    {
                        unlabeledToPeptidesDictionary.Add(fullyCleanedSequence, new List<FlashLFQ.Peptide> { peptide });
                    }
                }
            }
            return unlabeledToPeptidesDictionary;
        }

        //This function changes the LFQ (apex) intensities of the peaks to be SILAC intensities that only use intensities from MS1 scans where both heavy and light were quantified (only requires 2 for missed cleavages)
        //TODO: unclear how this affects data quality when deuterium shift is in play. A retention time calibration may be necessary.
        private static void CalculateSilacIntensities(Dictionary<SpectraFileInfo, List<ChromatographicPeak>> peakDictionary, Dictionary<string, List<FlashLFQ.Peptide>> unlabeledToPeptidesDictionary)
        {
            foreach (KeyValuePair<SpectraFileInfo, List<ChromatographicPeak>> kvp in peakDictionary)
            {
                //make a dictionary for easy peak lookup from peptide sequences
                Dictionary<string, ChromatographicPeak> sequenceToPeakDictionary = new Dictionary<string, ChromatographicPeak>();
                foreach (ChromatographicPeak peak in kvp.Value)
                {
                    //sometimes chromatographic peaks are separated (elute on both sides of the gradient, random blips, etc)
                    //in these situations, we want to use both chromatographic peaks for quantification, not just a single one.
                    string fullSequence = peak.Identifications.First().ModifiedSequence;
                    if (sequenceToPeakDictionary.ContainsKey(fullSequence))
                    {
                        ChromatographicPeak previousPeak = sequenceToPeakDictionary[fullSequence];
                        previousPeak.IsotopicEnvelopes.AddRange(peak.IsotopicEnvelopes);
                    }
                    else
                    {
                        sequenceToPeakDictionary.Add(fullSequence, peak);
                    }
                }

                foreach (List<FlashLFQ.Peptide> peptideGroup in unlabeledToPeptidesDictionary.Select(x => x.Value))
                {
                    List<string> sequences = peptideGroup.Select(x => x.Sequence).ToList();

                    //get peaks of interest for this peptide group
                    List<ChromatographicPeak> peaksOfInterest = new List<ChromatographicPeak>();
                    foreach (string sequence in sequences)
                    {
                        if (sequenceToPeakDictionary.TryGetValue(sequence, out ChromatographicPeak peakOfInterest))
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
                            foreach (FlashLFQ.Peptide peptide in peptideGroup)
                            {
                                ChromatographicPeak peakForThisPeptide = peaksOfInterest.Where(x => peptide.Sequence.Equals(x.Identifications.First().ModifiedSequence)).FirstOrDefault();
                                if (peakForThisPeptide != null)
                                {
                                    double summedIntensity = peakForThisPeptide.IsotopicEnvelopes.Where(x => scanIndex.Contains(x.IndexedPeak.ZeroBasedMs1ScanIndex)).Select(x => x.Intensity).Sum();
                                    peptide.SetIntensity(kvp.Key, summedIntensity);
                                }
                                else //rare instance, cause unknown. Crash identified using 180524_LMuscle_30d_bio3.raw, Mus_Canonical_180122.xml, 1 missed cleavage
                                {
                                    peptide.SetIntensity(kvp.Key, 0);
                                }
                            }
                        }
                    }
                }
            }
        }

        public static (SilacLabel updatedLabel, char nextHeavyLabel) UpdateAminoAcidLabel(SilacLabel currentLabel, char heavyLabel)
        {
            //make sure we're not overwriting something. , , and if it's a valid residue (not a motif/delimiter)
            while ((Residue.TryGetResidue(heavyLabel, out Residue residue) //Check if the amino acid exists. If it already exists, we don't want to overwrite it
                && !residue.ThisChemicalFormula.Formula.Equals(currentLabel.LabelChemicalFormula)) //if it exists but it's already the label (so we're not overwriting anything), then we're fine
                || GlobalVariables.InvalidAminoAcids.Contains(heavyLabel)) //If it didn't already exist, but it's invalid, we need to keep going
            {
                heavyLabel++;
            }
            SilacLabel updatedLabel = AssignValidHeavyCharacter(currentLabel, heavyLabel);
            heavyLabel++;
            if (currentLabel.AdditionalLabels != null)
            {
                foreach (SilacLabel additionalLabel in currentLabel.AdditionalLabels)
                {
                    updatedLabel.AddAdditionalSilacLabel(AssignValidHeavyCharacter(additionalLabel, heavyLabel));
                    heavyLabel++;
                }
            }
            return (updatedLabel, heavyLabel);
        }

        private static SilacLabel AssignValidHeavyCharacter(SilacLabel originalLabel, char heavyLabel)
        {
            double massDifference = Convert.ToDouble(originalLabel.MassDifference.Substring(1));
            if (originalLabel.MassDifference[0] == '-')
            {
                massDifference *= -1;
            }
            //Add the silac residues to the dictionary
            Residue.AddNewResiduesToDictionary(new List<Residue> { new Residue(originalLabel.MassDifference, heavyLabel, heavyLabel.ToString(), Chemistry.ChemicalFormula.ParseFormula(originalLabel.LabelChemicalFormula), ModificationSites.All) });

            return new SilacLabel(originalLabel.OriginalAminoAcid, heavyLabel, originalLabel.LabelChemicalFormula, massDifference);
        }
    }
}