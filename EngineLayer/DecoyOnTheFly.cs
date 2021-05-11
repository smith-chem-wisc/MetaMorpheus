using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class DecoyOnTheFly
    {
        //This function maintains the amino acids associated with the protease motif and reverses all other amino acids.
        //N-terminal modificatons are preserved. Other modifications travel with their respective amino acids. this results
        //in a decoy peptide composed the same amino acids and modifications as the original. 
        //Occasionally, this process results in peptide with exactly the same sequence. Therefore, there is a stop-gap measure
        //the returns the mirror image of the original. N-terminal mods are preserved, but other mods are also reversed. 
        //this should yield a unique decoy for each target sequence.
        public static PeptideWithSetModifications GetReverseDecoyFromTarget(PeptideWithSetModifications pwsm, int[] revisedAminoAcidOrder)
        {
            Dictionary<int, Modification> newModificationsDictionary = new Dictionary<int, Modification>();
            //Copy N-terminal modifications from target dictionary to decoy dictionary.
            if (pwsm.AllModsOneIsNterminus.ContainsKey(1))
            {
                newModificationsDictionary.Add(1, pwsm.AllModsOneIsNterminus[1]);
            }
            char[] newBase = new char[pwsm.BaseSequence.Length];
            ProteomicsExtenstionMethods.Fill(newBase, '0');
            char[] evaporatingBase = pwsm.BaseSequence.ToCharArray();
            List<DigestionMotif> motifs = pwsm.DigestionParams.Protease.DigestionMotifs;
            if (motifs != null && motifs.Count > 0)
            {
                foreach (var motif in motifs.Where(m => m.InducingCleavage != ""))//check the empty "" for topdown
                {
                    string cleavingMotif = motif.InducingCleavage;
                    List<int> cleavageMotifLocations = new List<int>();

                    for (int i = 0; i < pwsm.BaseSequence.Length; i++)
                    {
                        bool Fits;
                        bool Prevents;
                        (Fits, Prevents) = motif.Fits(pwsm.BaseSequence, i);

                        if (Fits && !Prevents)
                        {
                            cleavageMotifLocations.Add(i);
                        }
                    }

                    foreach (int location in cleavageMotifLocations)
                    {
                        char[] motifArray = pwsm.BaseSequence.Substring(location, cleavingMotif.Length).ToCharArray();
                        for (int i = 0; i < cleavingMotif.Length; i++)
                        {
                            newBase[location + i] = motifArray[i];
                            revisedAminoAcidOrder[location + i] = location + i;//
                            //directly copy mods that were on amino acids in the motif. Those amino acids don't change position.
                            if (pwsm.AllModsOneIsNterminus.ContainsKey(location + i + 2))
                            {
                                newModificationsDictionary.Add(location + i + 2, pwsm.AllModsOneIsNterminus[location + i + 2]);
                            }

                            evaporatingBase[location + i] = '0';//can null a char so i use a number which doesnt' appear in peptide string
                        }
                    }
                }
            }

            //We've kept amino acids in the digestion motif in the same position in the decoy peptide.
            //Now we will fill the remaining open positions in the decoy with the reverse of amino acids from the target.
            int fillPosition = 0;
            int extractPosition = pwsm.BaseSequence.Length - 1;
            while (fillPosition < pwsm.BaseSequence.Length && extractPosition >= 0)
            {
                if (evaporatingBase[extractPosition] != '0')
                {
                    while (newBase[fillPosition] != '0')
                    {
                        fillPosition++;
                    }
                    newBase[fillPosition] = evaporatingBase[extractPosition];
                    revisedAminoAcidOrder[fillPosition] = extractPosition;
                    if (pwsm.AllModsOneIsNterminus.ContainsKey(extractPosition + 2))
                    {
                        newModificationsDictionary.Add(fillPosition + 2, pwsm.AllModsOneIsNterminus[extractPosition + 2]);
                    }
                    fillPosition++;
                }
                extractPosition--;
            }

            string newBaseString = new string(newBase);

            var proteinSequence = pwsm.Protein.BaseSequence;
            var aStringBuilder = new StringBuilder(proteinSequence);
            aStringBuilder.Remove(pwsm.OneBasedStartResidueInProtein - 1, pwsm.BaseSequence.Length);
            aStringBuilder.Insert(pwsm.OneBasedStartResidueInProtein - 1, newBaseString);
            proteinSequence = aStringBuilder.ToString();

            Protein decoyProtein = new Protein(proteinSequence, "DECOY_" + pwsm.Protein.Accession, null, new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), null, null, null, true);

            DigestionParams d = pwsm.DigestionParams;

            if (newBaseString != pwsm.BaseSequence)
            {
                return new PeptideWithSetModifications(decoyProtein, d, pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedEndResidueInProtein, pwsm.CleavageSpecificityForFdrCategory, pwsm.PeptideDescription, pwsm.MissedCleavages, newModificationsDictionary, pwsm.NumFixedMods, newBaseString);
            }
            else
            {
                //The reverse decoy procedure failed to create a PeptideWithSetModificatons with a different sequence. Therefore,
                //we retrun the mirror image peptide.
                return GetPeptideMirror(pwsm, revisedAminoAcidOrder);
            }

        }

        //Returns a PeptideWithSetModifications mirror image. Used when reverse decoy sequence is same as target sequence
        public static PeptideWithSetModifications GetPeptideMirror(PeptideWithSetModifications pwsm, int[] revisedOrderNisOne)
        {
            Dictionary<int, Modification> newModificationsDictionary = new Dictionary<int, Modification>();
            //Copy N-terminal modifications from target dictionary to decoy dictionary.
            if (pwsm.AllModsOneIsNterminus.ContainsKey(1))
            {
                newModificationsDictionary.Add(1, pwsm.AllModsOneIsNterminus[1]);
            }

            //First step is to reverse the position of all modifications except the mod on the peptide N-terminus.
            if (pwsm.AllModsOneIsNterminus.Any())
            {
                int newPosition = pwsm.BaseSequence.Length + 1;
                for (int i = 2; i < pwsm.BaseSequence.Length + 2; i++)
                {
                    if (pwsm.AllModsOneIsNterminus.ContainsKey(i))
                    {
                        newModificationsDictionary.Add(newPosition, pwsm.AllModsOneIsNterminus[i]);
                        newPosition--;
                    }
                    else
                    {
                        newPosition--;
                    }
                }
            }

            //Second step is to reverse the sequence.
            string newBaseString = ProteomicsExtenstionMethods.Reverse(pwsm.BaseSequence);

            var proteinSequence = pwsm.Protein.BaseSequence;
            var aStringBuilder = new StringBuilder(proteinSequence);
            aStringBuilder.Remove(pwsm.OneBasedStartResidueInProtein - 1, pwsm.BaseSequence.Length);
            aStringBuilder.Insert(pwsm.OneBasedStartResidueInProtein - 1, newBaseString);
            proteinSequence = aStringBuilder.ToString();

            Protein decoyProtein = new Protein(proteinSequence, "DECOY_" + pwsm.Protein.Accession, null, new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), null, null, null, true);
            DigestionParams d = pwsm.DigestionParams;

            //now fill in the revised amino acid order
            int oldStringPosition = pwsm.BaseSequence.Length - 1;
            for (int i = 0; i < newBaseString.Length; i++)
            {
                revisedOrderNisOne[i] = oldStringPosition;
                oldStringPosition--;
            }
            return new PeptideWithSetModifications(decoyProtein, d, pwsm.OneBasedStartResidueInProtein, pwsm.OneBasedEndResidueInProtein, pwsm.CleavageSpecificityForFdrCategory, pwsm.PeptideDescription, pwsm.MissedCleavages, newModificationsDictionary, pwsm.NumFixedMods, newBaseString);
        }

        
    }
}
