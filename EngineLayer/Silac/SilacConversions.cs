using EngineLayer;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{
    public static class SilacConversions
    {
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
                    PeptideWithSetModifications pwsm = notchAndPwsm.Peptide;
                    Protein originalProtein = pwsm.Protein;
                    Protein modifiedProtein = heavyToLight ?
                        new Protein(originalProtein, originalProtein.BaseSequence.Replace(silacLabel.AminoAcidLabel, silacLabel.OriginalAminoAcid), originalProtein.Accession.Replace(silacLabel.MassDifference, "")) : //create light protein
                        new Protein(originalProtein, originalProtein.BaseSequence.Replace(silacLabel.OriginalAminoAcid, silacLabel.AminoAcidLabel), originalProtein.Accession + silacLabel.MassDifference); //create heavy protein with different accessions
                    PeptideWithSetModifications modifiedPwsm = new PeptideWithSetModifications(
                        modifiedProtein,
                        pwsm.DigestionParams,
                        pwsm.OneBasedStartResidueInProtein,
                        pwsm.OneBasedEndResidueInProtein,
                        pwsm.CleavageSpecificityForFdrCategory,
                        pwsm.PeptideDescription,
                        pwsm.MissedCleavages,
                        pwsm.AllModsOneIsNterminus,
                        pwsm.NumFixedMods);
                    updatedBestMatchingPeptides.Add((notchAndPwsm.Notch, modifiedPwsm));
                }
                return psm.Clone(updatedBestMatchingPeptides);
            }
        }

        public static PeptideSpectralMatch GetSilacPsmFromAmbiguousPsm(PeptideSpectralMatch psm, List<SilacLabel> silacLabels)
        {
            List<(int Notch, PeptideWithSetModifications Peptide)> updatedBestMatchingPeptides = new List<(int Notch, PeptideWithSetModifications Peptide)>();
            foreach ((int Notch, PeptideWithSetModifications Peptide) notchAndPwsm in psm.BestMatchingPeptides)
            {
                PeptideWithSetModifications pwsm = notchAndPwsm.Peptide;
                Protein originalProtein = pwsm.Protein;
                SilacLabel silacLabel = silacLabels.Where(x => originalProtein.BaseSequence.Contains(x.AminoAcidLabel)).FirstOrDefault();
                if (silacLabel == null)
                {
                    updatedBestMatchingPeptides.Add(notchAndPwsm);
                }
                else
                {
                    Protein modifiedProtein = new Protein(originalProtein, originalProtein.BaseSequence.Replace(silacLabel.AminoAcidLabel, silacLabel.OriginalAminoAcid), originalProtein.Accession.Replace(silacLabel.MassDifference, "")); //create light protein
                    PeptideWithSetModifications modifiedPwsm = new PeptideWithSetModifications(
                        modifiedProtein,
                        pwsm.DigestionParams,
                        pwsm.OneBasedStartResidueInProtein,
                        pwsm.OneBasedEndResidueInProtein,
                        pwsm.CleavageSpecificityForFdrCategory,
                        pwsm.PeptideDescription,
                        pwsm.MissedCleavages,
                        pwsm.AllModsOneIsNterminus,
                        pwsm.NumFixedMods);
                    updatedBestMatchingPeptides.Add((notchAndPwsm.Notch, modifiedPwsm));
                }
            }
            return psm.Clone(updatedBestMatchingPeptides);
        }

        public static EngineLayer.ProteinGroup GetSilacProteinGroups(List<PeptideSpectralMatch> unambiguousPsmsBelowOnePercentFdr, EngineLayer.ProteinGroup proteinGroup, SilacLabel label = null)
        {
            HashSet<Protein> proteins = label == null ?
                proteinGroup.Proteins :
                new HashSet<Protein>(proteinGroup.Proteins.Select(x => new Protein(x, x.BaseSequence.Replace(label.OriginalAminoAcid, label.AminoAcidLabel), x.Accession + label.MassDifference)));
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
                    if (firstAccession.Equals(psm.ProteinAccession)) //since unique, we know there's only one protein for this sequence
                    {
                        var peptide = bestMatchingPeptides.First().Peptide;
                        uniquePeptides.Add(peptide);
                        allPeptides.Add(peptide);
                        matchedPsms.Add(psm);
                    }
                }
                else //not unique
                {
                    foreach (var peptide in bestMatchingPeptides.Select(x => x.Peptide))
                    {
                        if (firstAccession.Equals(peptide.Protein.Accession))
                        {
                            allPeptides.Add(peptide);
                            matchedPsms.Add(psm);
                        }
                    }
                }
            }
            return new EngineLayer.ProteinGroup(proteins, allPeptides, uniquePeptides)
            {
                AllPsmsBelowOnePercentFDR = matchedPsms
            };
        }

        public static string GetSilacLightBaseSequence(string baseSequence, SilacLabel label)
        {
            //overwrite base sequence
            return label == null ? baseSequence : baseSequence.Replace(label.AminoAcidLabel.ToString(), HeavyStringForPeptides(label));
        }

        public static string GetSilacLightFullSequence(string fullSequence, SilacLabel label)
        {
            //overwrite full sequence
            if (label == null)
            {
                return fullSequence;
            }
            else
            {
                string replacementSequence = HeavyStringForPeptides(label);

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
                        else if (currentChar == label.AminoAcidLabel)
                        {
                            fullSequence = fullSequence.Substring(0, i) + replacementSequence + fullSequence.Substring(i + 1, fullSequence.Length - i - 1);
                            i += replacementSequence.Length - 1; //-1 because we removed the label amino acid
                        }
                    }
                }
                return fullSequence;
            }
        }

        public static string HeavyStringForPeptides(SilacLabel label)
        {
            return label.OriginalAminoAcid + "(" + label.MassDifference + ")";
        }
    }
}
