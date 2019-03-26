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
                proteinAccession += silacLabel.MassDifference; //add heavy accession
                if (silacLabel.AdditionalLabels != null)
                {
                    foreach (SilacLabel additionalLabel in silacLabel.AdditionalLabels)
                    {
                        proteinSequence = proteinSequence.Replace(additionalLabel.OriginalAminoAcid, additionalLabel.AminoAcidLabel); //create heavy sequence
                        proteinAccession += LABEL_DELIMITER + additionalLabel.MassDifference; //add heavy accession
                    }
                }
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
            if(label.AdditionalLabels!=null)
            {
                foreach(SilacLabel additionaLabel in label.AdditionalLabels)
                {
                    heavyFileName += LABEL_DELIMITER + additionaLabel.OriginalAminoAcid + additionaLabel.MassDifference;
                }
            }
            heavyFileName += ")." + originalFile.FullFilePathWithExtension.Split('.').Last(); //add extension

            return new SpectraFileInfo(heavyFileName, originalFile.Condition, originalFile.BiologicalReplicate, originalFile.TechnicalReplicate, originalFile.Fraction);
        }
    }
}