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
        public static PeptideSpectralMatch GetSilacPsm(PeptideSpectralMatch psm, SilacLabel silacLabel, bool heavyToLight)
        {
            if (silacLabel == null)
            {
                return psm;
            }
            else
            {
                List<(int Notch, PeptideWithSetModifications Peptide)> updatedBestMatchingPeptides = new List<(int Notch, PeptideWithSetModifications Peptide)>();
                foreach ((int Notch, PeptideWithSetModifications Peptide) notchAndPwsm in psm.BestMatchingPeptides)
                {
                    PeptideWithSetModifications modifiedPwsm = CreateSilacPwsm(heavyToLight, silacLabel, notchAndPwsm.Peptide);
                    updatedBestMatchingPeptides.Add((notchAndPwsm.Notch, modifiedPwsm));
                }
                return psm.Clone(updatedBestMatchingPeptides);
            }
        }

        //Needed for parsimony, where there are ambiguous psms
        //Quantification ignores ambiguity
        public static PeptideSpectralMatch GetSilacPsmFromAmbiguousPsm(PeptideSpectralMatch psm, List<SilacLabel> silacLabels)
        {
            List<(int Notch, PeptideWithSetModifications Peptide)> updatedBestMatchingPeptides = new List<(int Notch, PeptideWithSetModifications Peptide)>();
            foreach ((int Notch, PeptideWithSetModifications Peptide) notchAndPwsm in psm.BestMatchingPeptides)
            {
                PeptideWithSetModifications pwsm = notchAndPwsm.Peptide;
                SilacLabel silacLabel = GetRelevantLabelFromBaseSequence(pwsm.Protein.BaseSequence, silacLabels);
                if (silacLabel == null)
                {
                    updatedBestMatchingPeptides.Add(notchAndPwsm);
                }
                else
                {
                    PeptideWithSetModifications modifiedPwsm = CreateSilacPwsm(true, silacLabel, pwsm); //create light pwsm
                    updatedBestMatchingPeptides.Add((notchAndPwsm.Notch, modifiedPwsm));
                }
            }
            return psm.Clone(updatedBestMatchingPeptides);
        }

        //This method creates a protein group based on the provided silac label. The input "proteinGroup" is expected to be light. If the label is null, then the light protein group will be output.
        //This method searches the provided psm list for which psms belong to the new/old protein group independent of the original proteinGroup's psms
        public static ProteinGroup GetSilacProteinGroups(List<PeptideSpectralMatch> unambiguousPsmsBelowOnePercentFdr, ProteinGroup proteinGroup, SilacLabel label = null)
        {
            //keep the proteins as light, or convert them all into the heavy versions
            HashSet<Protein> proteins = label == null ?
                proteinGroup.Proteins :
                new HashSet<Protein>(proteinGroup.Proteins.Select(x => CreateSilacProtein(false, label, x)));
            HashSet<PeptideWithSetModifications> allPeptides = new HashSet<PeptideWithSetModifications>();
            HashSet<PeptideWithSetModifications> uniquePeptides = new HashSet<PeptideWithSetModifications>();
            string firstAccession = proteins.First().Accession;
            HashSet<PeptideSpectralMatch> matchedPsms = new HashSet<PeptideSpectralMatch>();

            //go through all psms and find peptides that belong to this group
            foreach (PeptideSpectralMatch psm in unambiguousPsmsBelowOnePercentFdr)
            {
                var bestMatchingPeptides = psm.BestMatchingPeptides.ToList();
                if (bestMatchingPeptides.Count == 1)//if unique
                {
                    if (firstAccession.Equals(psm.ProteinAccession)) //since unique, we know there's only one protein for this sequence. Shared peptides are given multiple unique pwsms.
                    {
                        var peptide = bestMatchingPeptides.First().Peptide;
                        uniquePeptides.Add(peptide);
                        allPeptides.Add(peptide);
                        matchedPsms.Add(psm);
                    }
                }
                else //not unique
                {
                    foreach (var peptide in bestMatchingPeptides.Select(x => x.Peptide)) //go through all the peptides
                    {
                        if (firstAccession.Equals(peptide.Protein.Accession)) //if one of them matches, then we'll add it.
                        {
                            allPeptides.Add(peptide);
                            matchedPsms.Add(psm);
                            break;
                        }
                    }
                }
            }
            var updatedProtein = new ProteinGroup(proteins, allPeptides, uniquePeptides)
            {
                AllPsmsBelowOnePercentFDR = matchedPsms,
                CumulativeTarget = proteinGroup.CumulativeTarget,
                CumulativeDecoy = proteinGroup.CumulativeDecoy,
                QValue = proteinGroup.QValue,
                BestPeptideScore = proteinGroup.BestPeptideScore,
                BestPeptideQValue = proteinGroup.BestPeptideQValue
            };
            updatedProtein.CalculateSequenceCoverage();
            return updatedProtein;
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

        public static string GetAmbiguousLightSequence(string originalSequence, List<SilacLabel> labels, bool baseSequence)
        {
            string[] multipleSequences = originalSequence.Split('|').ToArray(); //if ambiguity
            string localSequence = "";
            foreach (string sequence in multipleSequences)
            {
                SilacLabel label = GetRelevantLabelFromBaseSequence(sequence, labels);
                localSequence += (baseSequence ? GetSilacLightBaseSequence(sequence, label) : GetSilacLightFullSequence(sequence, label)) + "|";
            }
            if (localSequence.Length != 0)
            {
                localSequence = localSequence.Substring(0, localSequence.Length - 1); //remove last "|"
            }
            return localSequence;
        }

        public static string GetProteinLightAccession(string proteinAccession, List<SilacLabel> labels)
        {
            foreach (SilacLabel label in labels)
            {
                proteinAccession = proteinAccession.Replace(label.MassDifference, "");
            }
            return proteinAccession;
        }

        public static List<PeptideSpectralMatch> UpdatePsmsForParsimony(List<SilacLabel> labels, List<PeptideSpectralMatch> psms)
        {
            //consider the heavy and light psms as being from the same protein
            //currently they have different accessions and sequences (PROTEIN and PROTEIN+8.014)
            List<PeptideSpectralMatch> psmsForProteinParsimony = new List<PeptideSpectralMatch>();
            foreach (PeptideSpectralMatch psm in psms)
            {
                if (psm.BaseSequence == null)
                {
                    psmsForProteinParsimony.Add(GetSilacPsmFromAmbiguousPsm(psm, labels));
                }
                else
                {
                    SilacLabel label = GetRelevantLabelFromBaseSequence(psm.BaseSequence, labels);
                    psmsForProteinParsimony.Add(GetSilacPsm(psm, label, true)); //if it's light, label will be null
                }
            }
            return psmsForProteinParsimony;
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

        public static Protein CreateSilacProtein(bool heavyToLight, SilacLabel silacLabel, Protein originalProtein)
        {
            string proteinSequence = originalProtein.BaseSequence;
            string proteinAccession = originalProtein.Accession;

            if (heavyToLight)
            {
                proteinSequence = proteinSequence.Replace(silacLabel.AminoAcidLabel, silacLabel.OriginalAminoAcid); //create light sequence
                int labelStart = proteinAccession.IndexOf(silacLabel.MassDifference);
                proteinAccession = proteinAccession.Substring(0, labelStart); //create light accession
                if (silacLabel.AdditionalLabels != null)
                {
                    foreach (SilacLabel additionalLabel in silacLabel.AdditionalLabels)
                    {
                        proteinSequence = proteinSequence.Replace(additionalLabel.AminoAcidLabel, additionalLabel.OriginalAminoAcid); //create light sequence
                    }
                }
            }
            else
            {
                proteinSequence = proteinSequence.Replace(silacLabel.OriginalAminoAcid, silacLabel.AminoAcidLabel); //create heavy sequence
                proteinAccession += "(" + silacLabel.OriginalAminoAcid + silacLabel.MassDifference; //add heavy accession
                if (silacLabel.AdditionalLabels != null)
                {
                    foreach (SilacLabel additionalLabel in silacLabel.AdditionalLabels)
                    {
                        proteinSequence = proteinSequence.Replace(additionalLabel.OriginalAminoAcid, additionalLabel.AminoAcidLabel); //create heavy sequence
                        proteinAccession += LABEL_DELIMITER + additionalLabel.OriginalAminoAcid + additionalLabel.MassDifference; //add heavy accession
                    }
                }
                proteinAccession += ")";
            }

            return new Protein(originalProtein, proteinSequence, proteinAccession);
        }

        public static PeptideWithSetModifications CreateSilacPwsm(bool heavyToLight, SilacLabel silacLabel, PeptideWithSetModifications pwsm)
        {
            Protein modifiedProtein = CreateSilacProtein(heavyToLight, silacLabel, pwsm.Protein);

            return new PeptideWithSetModifications(
                modifiedProtein,
                pwsm.DigestionParams,
                pwsm.OneBasedStartResidueInProtein,
                pwsm.OneBasedEndResidueInProtein,
                pwsm.CleavageSpecificityForFdrCategory,
                pwsm.PeptideDescription,
                pwsm.MissedCleavages,
                pwsm.AllModsOneIsNterminus,
                pwsm.NumFixedMods);
        }

        public static SilacLabel AssignValidHeavyCharacter(SilacLabel originalLabel, char heavyLabel)
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

        public static SpectraFileInfo GetHeavyFileInfo(SpectraFileInfo originalFile, SilacLabel label)
        {
            string heavyFileName = originalFile.FilenameWithoutExtension + "(" + label.OriginalAminoAcid + label.MassDifference;
            if (label.AdditionalLabels != null)
            {
                foreach (SilacLabel additionaLabel in label.AdditionalLabels)
                {
                    heavyFileName += LABEL_DELIMITER + additionaLabel.OriginalAminoAcid + additionaLabel.MassDifference;
                }
            }
            heavyFileName += ")." + originalFile.FullFilePathWithExtension.Split('.').Last(); //add extension

            return new SpectraFileInfo(heavyFileName, originalFile.Condition, originalFile.BiologicalReplicate, originalFile.TechnicalReplicate, originalFile.Fraction);
        }


        //If SILAC (Post-Quantification), compress the light/heavy protein group pairs into the same light protein group but different files
        //Create new files for each silac label and file so that "file 1" now becomes "file 1 (light)" and "file 1 (heavy)"
        //Change heavy residue into the light residue plus a string label ("PEPTIDEa" -> "PEPTIDEK(+8.014)")
        //This light to heavy conversion needs to happen for the flashLFQ peptides here, but can't for the psm peptides, which are constrained to the protein
        //i.e. pwsms currently don't have sequences; they have start/end residues and a protein sequence. We have to change the output sequences when they're created.
        public static void SilacConversionsPostQuantification(List<SilacLabel> silacLabels, List<SpectraFileInfo> spectraFileInfo, List<ProteinGroup> ProteinGroups,
            HashSet<DigestionParams> ListOfDigestionParams, Dictionary<string, List<string>> silacProteinGroupMatcher, FlashLfqResults FlashLfqResults,
            List<PeptideSpectralMatch> allPsms, Dictionary<string, int> ModsToWriteSelection, bool Integrate)
        {
            bool outputLightIntensities = ListOfDigestionParams.Any(x => x.GeneratehUnlabeledProteinsForSilac);


            //MAKE NEW RAW FILES
            //update number of spectra files to include a new file for each label*condition
            Dictionary<SpectraFileInfo, string> fileToLabelDictionary = new Dictionary<SpectraFileInfo, string>(); //figure out which file is which label, since some files will be only light and others only heavy. Key is file, value is the label string (label.MassDifference)
            Dictionary<SpectraFileInfo, SpectraFileInfo> labeledToUnlabeledFile = new Dictionary<SpectraFileInfo, SpectraFileInfo>(); //keep track of the heavy-to-light pairs. If multiple, looks like 3-1 and 2-1, but no 3-2 (only heavy to light, no heavy to heavy)
            List<SpectraFileInfo> silacSpectraFileInfo = new List<SpectraFileInfo>(); //new files

            //foreach existing file
            foreach (SpectraFileInfo originalFile in spectraFileInfo)
            {
                //add the existing file as the light
                silacSpectraFileInfo.Add(originalFile);
                //foreach label, add a new file with the label
                foreach (SilacLabel label in silacLabels)
                {
                    SpectraFileInfo silacFile = GetHeavyFileInfo(originalFile, label);
                    silacSpectraFileInfo.Add(silacFile);
                    fileToLabelDictionary[silacFile] = label.MassDifference;
                    labeledToUnlabeledFile[silacFile] = originalFile;
                }
            }


            //UPDATE PROTEIN GROUPS
            //remove the heavy protein groups so that there are only light ones
            //add the intensities of the heavy groups into the newly created heavy SpectraFileInfos
            HashSet<SpectraFileInfo> lightFilesToRemove = new HashSet<SpectraFileInfo>(); //this is only used when there user specified no unlabeled proteins
            if (ProteinGroups != null) //if we did parsimony
            {
                List<EngineLayer.ProteinGroup> silacProteinGroups = new List<EngineLayer.ProteinGroup>();
                //The light/unlabeled peptides/proteins were not searched if specified, but they were still quantified to keep track of the labels
                //we need to remove these unlabeled peptides/proteins before output
                //foreach protein group (which has its own quant for each file)
                foreach (EngineLayer.ProteinGroup proteinGroup in ProteinGroups)
                {
                    proteinGroup.FilesForQuantification = silacSpectraFileInfo; //update fileinfo for the group
                                                                                //grab the light groups. Using these light groups, find their heavy group pair(s), add them to the light group quant info, and then remove the heavy groups
                    if (silacProteinGroupMatcher.TryGetValue(proteinGroup.ProteinGroupName, out List<string> silacSubGroupNames)) //try to find the light protein groups. If it's not light, ignore it
                    {
                        //the out variable contains all the other heavy protein groups that were generated for this light protein group
                        //go through the files and see if any of them contain the same label. If not, put zeroes for those missing "files"
                        //If the user didn't specify to search light intensities, then don't output them
                        Dictionary<SpectraFileInfo, double> updatedIntensitiesByFile = proteinGroup.IntensitiesByFile; //light intensities
                        List<SpectraFileInfo> lightKeys = updatedIntensitiesByFile.Keys.ToList();

                        //go through all files (including "silac" files)
                        List<ProteinGroup> subGroup = ProteinGroups.Where(x => silacSubGroupNames.Contains(x.ProteinGroupName)).ToList(); //find the protein groups where the accession contains "light" accession of the current protein group
                        foreach (SpectraFileInfo fileInfo in silacSpectraFileInfo) //for every file (light and heavy)
                        {
                            //if it doesn't have a value, then it's a silac file (light missing values still have a value "0")
                            if (!updatedIntensitiesByFile.ContainsKey(fileInfo))
                            {
                                string labelSignature = fileToLabelDictionary[fileInfo]; //a string associated with a silac label
                                ProteinGroup foundGroup = subGroup.Where(x => x.Proteins.Any(y => y.Accession.Contains(labelSignature))).FirstOrDefault(); //get the protein groups containing this label
                                updatedIntensitiesByFile[fileInfo] = foundGroup == null ? 0 : foundGroup.IntensitiesByFile[labeledToUnlabeledFile[fileInfo]]; //update the intensity for that label in the light group
                            }
                            //else do nothing. The light version is already in the dictionary
                        }

                        //The light/unlabeled peptides/proteins were not searched if specified, but they were still quantified to keep track of the labels
                        //we need to remove these unlabeled peptides/proteins before output
                        if (!outputLightIntensities)
                        {
                            foreach (SpectraFileInfo info in lightKeys)
                            {
                                updatedIntensitiesByFile.Remove(info);
                                proteinGroup.FilesForQuantification.Remove(info);
                                lightFilesToRemove.Add(info);
                            }
                        }

                        silacProteinGroups.Add(proteinGroup);
                    }
                }

                //update
                ProteinGroups.Clear();
                ProteinGroups.AddRange(silacProteinGroups);
                //remove light files (if necessary)
                foreach (SpectraFileInfo info in lightFilesToRemove)
                {
                    FlashLfqResults.SpectraFiles.Remove(info);
                }

                //UPDATE FLASHLFQ PROTEINS
                if (FlashLfqResults != null) //can be null if nothing was quantified (all peptides are ambiguous)
                {
                    Dictionary<string, FlashLFQ.ProteinGroup> flashLfqProteins = FlashLfqResults.ProteinGroups; //dictionary of protein group names to protein groups
                                                                                                                //if the protein group is a heavy protein group, get rid of it. We already accounted for it above.
                    var keys = flashLfqProteins.Keys.ToList();
                    foreach (string key in keys)
                    {
                        if (silacLabels.Any(x => key.Contains(x.MassDifference)))
                        {
                            flashLfqProteins.Remove(key);
                        }
                    }
                }
            }

            ////UPDATE FLASHLFQ SPECTRA FILES
            if (FlashLfqResults != null) //can be null if nothing was quantified (all peptides are ambiguous)
            {
                List<SpectraFileInfo> originalFiles = FlashLfqResults.SpectraFiles; //pass reference
                foreach (SpectraFileInfo info in silacSpectraFileInfo)
                {
                    if (!originalFiles.Contains(info))
                    {
                        originalFiles.Add(info);
                    }
                }
            }

            //UPDATE PEPTIDE INFO
            //convert all psm/peptide/proteingroup sequences from the heavy label to the light label for output
            //We can do this for all of the FlashLFQ peptides/peaks, because they use string sequences.
            //We are unable to do this for Parameters.AllPsms, because they store proteins and start/end residues instead
            //for Psms, we need to convert during the writing.
            for (int i = 0; i < allPsms.Count; i++)
            {
                allPsms[i].ResolveHeavySilacLabel(silacLabels, ModsToWriteSelection);
            }

            //Convert all lfqpeaks from heavy (a) to light (K+8.014) for output
            if (FlashLfqResults != null) //can be null if nothing was quantified (all peptides are ambiguous)
            {
                var lfqPeaks = FlashLfqResults.Peaks;
                List<SpectraFileInfo> peakKeys = lfqPeaks.Keys.ToList();

                foreach (SpectraFileInfo key in peakKeys)
                {
                    List<FlashLFQ.ChromatographicPeak> peaks = lfqPeaks[key];
                    for (int i = 0; i < peaks.Count; i++)
                    {
                        var peak = peaks[i];
                        List<Identification> identifications = new List<Identification>();
                        //check if we're removing light peaks and if it's a light peak
                        if (!outputLightIntensities && !peak.Identifications.Any(x => GetRelevantLabelFromBaseSequence(x.BaseSequence, silacLabels) != null)) //if no ids have any labels, remove them
                        {
                            peaks.RemoveAt(i);
                            i--;
                        }
                        else
                        {
                            foreach (var id in peak.Identifications)
                            {
                                SilacLabel label = GetRelevantLabelFromBaseSequence(id.BaseSequence, silacLabels);
                                HashSet<FlashLFQ.ProteinGroup> originalGroups = id.proteinGroups;
                                List<FlashLFQ.ProteinGroup> updatedGroups = new List<FlashLFQ.ProteinGroup>();
                                foreach (FlashLFQ.ProteinGroup group in originalGroups)
                                {
                                    string groupName = group.ProteinGroupName;
                                    if (label == null) //if light
                                    {
                                        updatedGroups.Add(group);
                                    }
                                    else
                                    {
                                        string labelString = "(" + label.OriginalAminoAcid + label.MassDifference;
                                        int labelIndex = groupName.IndexOf(labelString);
                                        if (labelIndex != -1) //labelIndex == 1 if a) 2+ peptides are required per protein or b) somebody broke parsimony
                                        {
                                            groupName = groupName.Substring(0, labelIndex);
                                            updatedGroups.Add(new FlashLFQ.ProteinGroup(groupName, group.GeneName, group.Organism));
                                        }
                                    }
                                }

                                Identification updatedId = new Identification(
                                    id.fileInfo,
                                    GetSilacLightBaseSequence(id.BaseSequence, label),
                                    GetSilacLightFullSequence(id.ModifiedSequence, label),
                                    id.monoisotopicMass,
                                    id.ms2RetentionTimeInMinutes,
                                    id.precursorChargeState,
                                    updatedGroups,
                                    id.OptionalChemicalFormula,
                                    id.UseForProteinQuant
                                    );
                                identifications.Add(updatedId);
                            }
                            FlashLFQ.ChromatographicPeak updatedPeak = new FlashLFQ.ChromatographicPeak(identifications.First(), peak.IsMbrPeak, peak.SpectraFileInfo);
                            for (int j = 1; j < identifications.Count; j++) //add all the original identification
                            {
                                updatedPeak.MergeFeatureWith(new FlashLFQ.ChromatographicPeak(identifications[j], peak.IsMbrPeak, peak.SpectraFileInfo), Integrate);
                            }
                            updatedPeak.IsotopicEnvelopes = peak.IsotopicEnvelopes; //need to set isotopicEnevelopes, since the new identifications didn't have them.
                            updatedPeak.CalculateIntensityForThisFeature(Integrate); //needed to update info
                            peaks[i] = updatedPeak;
                        }
                    }
                }

                //convert all lfq peptides from heavy to light for output
                Dictionary<string, FlashLFQ.Peptide> lfqPwsms = FlashLfqResults.PeptideModifiedSequences;
                List<string> pwsmKeys = lfqPwsms.Keys.ToList();
                foreach (string key in pwsmKeys)
                {
                    FlashLFQ.Peptide currentPeptide = lfqPwsms[key];
                    SilacLabel label = GetRelevantLabelFromFullSequence(currentPeptide.Sequence, silacLabels);
                    if (label != null) //if it's a heavy peptide
                    {
                        lfqPwsms.Remove(key); //get rid of it
                                              //update the light version
                        string lightSequence = GetSilacLightFullSequence(currentPeptide.Sequence, label, false); //get the light sequence
                        List<SpectraFileInfo> heavyFiles = silacSpectraFileInfo.Where(x => x.FilenameWithoutExtension.Contains(label.MassDifference)).ToList(); //these are the heavy raw file names

                        //Find the light peptide (which has a value for the light datafile) and set the intensity for the heavy datafile from the current peptide
                        if (lfqPwsms.TryGetValue(lightSequence, out FlashLFQ.Peptide lightPeptide)) //this should always have a value, since we made replicas earlier, and yet it sometimes doesn't...
                        {
                            foreach (SpectraFileInfo heavyFile in heavyFiles)
                            {
                                SpectraFileInfo lightFile = labeledToUnlabeledFile[heavyFile];
                                lightPeptide.SetIntensity(heavyFile, currentPeptide.GetIntensity(lightFile));
                                lightPeptide.SetDetectionType(heavyFile, currentPeptide.GetDetectionType(lightFile));
                            }
                        }
                        else //if there's no light, create a new entry for the heavy
                        {
                            //new peptide
                            FlashLFQ.Peptide updatedPeptide = new FlashLFQ.Peptide(lightSequence, currentPeptide.UseForProteinQuant);
                            //update the heavy info, set the light values to zero
                            foreach (SpectraFileInfo info in heavyFiles)
                            {
                                updatedPeptide.SetIntensity(info, currentPeptide.GetIntensity(info));
                                updatedPeptide.SetDetectionType(info, currentPeptide.GetDetectionType(info));
                            }

                            //set the other values to zero
                            List<SpectraFileInfo> otherInfo = silacSpectraFileInfo.Where(x => !heavyFiles.Contains(x)).ToList();
                            foreach (SpectraFileInfo info in otherInfo)
                            {
                                updatedPeptide.SetIntensity(info, 0);
                                updatedPeptide.SetDetectionType(info, DetectionType.NotDetected);
                            }
                            HashSet<FlashLFQ.ProteinGroup> originalGroups = currentPeptide.proteinGroups;
                            HashSet<FlashLFQ.ProteinGroup> updatedGroups = new HashSet<FlashLFQ.ProteinGroup>();
                            foreach (FlashLFQ.ProteinGroup group in originalGroups)
                            {
                                string groupName = group.ProteinGroupName;
                                groupName = groupName.Replace(label.MassDifference, "");
                                updatedGroups.Add(new FlashLFQ.ProteinGroup(groupName, group.GeneName, group.Organism));
                            }
                            updatedPeptide.proteinGroups = updatedGroups;
                            lfqPwsms[updatedPeptide.Sequence] = updatedPeptide;
                        }
                    }
                }
            }
        }
    }
}