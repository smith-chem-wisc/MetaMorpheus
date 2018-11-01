using FlashLFQ;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class ProteinGroup
    {
        public ProteinGroup(HashSet<Protein> proteins, HashSet<PeptideWithSetModifications> peptides, HashSet<PeptideWithSetModifications> uniquePeptides)
        {
            Proteins = proteins;
            ListOfProteinsOrderedByAccession = Proteins.OrderBy(p => p.Accession).ToList();
            ProteinGroupName = string.Join("|", ListOfProteinsOrderedByAccession.Select(p => p.Accession));
            AllPeptides = peptides;
            UniquePeptides = uniquePeptides;
            AllPsmsBelowOnePercentFDR = new HashSet<PeptideSpectralMatch>();
            SequenceCoveragePercent = new List<double>();
            SequenceCoverageDisplayList = new List<string>();
            SequenceCoverageDisplayListWithMods = new List<string>();
            ProteinGroupScore = 0;
            BestPeptideScore = 0;
            QValue = 0;
            IsDecoy = false;
            IsContaminant = false;
            ModsInfo = new List<string>();

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.IsDecoy)
                {
                    IsDecoy = true;
                    break;
                }
                if (protein.IsContaminant)
                {
                    IsContaminant = true;
                    break;
                }
            }
        }

        public bool IsDecoy { get; }

        public bool IsContaminant { get; }

        public List<SpectraFileInfo> FilesForQuantification { get; set; }

        public HashSet<Protein> Proteins { get; set; }

        public string ProteinGroupName { get; private set; }

        public double ProteinGroupScore { get; set; }

        public HashSet<PeptideWithSetModifications> AllPeptides { get; set; }

        public HashSet<PeptideWithSetModifications> UniquePeptides { get; set; }

        public HashSet<PeptideSpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        public List<double> SequenceCoveragePercent { get; private set; }

        public List<string> SequenceCoverageDisplayList { get; private set; }

        public List<string> SequenceCoverageDisplayListWithMods { get; private set; }

        public double QValue { get; set; }

        public double BestPeptideQValue { get; set; }

        public double BestPeptideScore { get; set; }

        public int CumulativeTarget { get; set; }

        public int CumulativeDecoy { get; set; }

        public bool DisplayModsOnPeptides { get; set; }

        public List<string> ModsInfo { get; private set; }

        public Dictionary<SpectraFileInfo, double> IntensitiesByFile { get; set; }

        private List<Protein> ListOfProteinsOrderedByAccession;

        public string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("Protein Accession" + '\t');
            sb.Append("Gene" + '\t');
            sb.Append("Organism" + '\t');
            sb.Append("Protein Full Name" + '\t');
            sb.Append("Protein Unmodified Mass" + '\t');
            sb.Append("Number of Proteins in Group" + '\t');
            sb.Append("Unique Peptides" + '\t');
            sb.Append("Shared Peptides" + '\t');
            sb.Append("Number of Peptides" + '\t');
            sb.Append("Number of Unique Peptides" + '\t');
            sb.Append("Sequence Coverage %" + '\t');
            sb.Append("Sequence Coverage" + '\t');
            sb.Append("Sequence Coverage with Mods" + '\t');
            sb.Append("Modification Info List" + "\t");
            if (FilesForQuantification != null)
            {
                for (int i = 0; i < FilesForQuantification.Count; i++)
                {
                    sb.Append("Intensity_" + FilesForQuantification[i].FilenameWithoutExtension + '\t');
                }
            }
            sb.Append("Number of PSMs" + '\t');
            sb.Append("Protein Decoy/Contaminant/Target" + '\t');
            sb.Append("Protein Cumulative Target" + '\t');
            sb.Append("Protein Cumulative Decoy" + '\t');
            sb.Append("Protein QValue" + '\t');
            sb.Append("Best Peptide Score" + '\t');
            sb.Append("Best Peptide Notch QValue");
            return sb.ToString();
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            // list of protein accession numbers
            sb.Append(ProteinGroupName);
            sb.Append("\t");

            // genes
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", ListOfProteinsOrderedByAccession.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault()))));
            sb.Append("\t");

            // organisms
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", ListOfProteinsOrderedByAccession.Select(p => p.Organism).Distinct())));
            sb.Append("\t");

            // list of protein names
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", ListOfProteinsOrderedByAccession.Select(p => p.FullName).Distinct())));
            sb.Append("\t");

            // list of masses
            var sequences = ListOfProteinsOrderedByAccession.Select(p => p.BaseSequence).Distinct();
            List<double> masses = new List<double>();
            foreach (var sequence in sequences)
            {
                try
                {
                    masses.Add(new Proteomics.AminoAcidPolymer.Peptide(sequence).MonoisotopicMass);
                }
                catch (System.Exception)
                {
                    masses.Add(double.NaN);
                }
            }
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", masses)));
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + Proteins.Count);
            sb.Append("\t");

            // list of unique peptides
            if (!DisplayModsOnPeptides)
                sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", UniquePeptides.Select(p => p.BaseSequence).Distinct())));
            else
                sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", UniquePeptides.Select(p => p.FullSequence).Distinct())));
            sb.Append("\t");

            // list of shared peptides
            var SharedPeptides = AllPeptides.Except(UniquePeptides);
            if (!DisplayModsOnPeptides)
                sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", SharedPeptides.Select(p => p.BaseSequence).Distinct())));
            else
                sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", SharedPeptides.Select(p => p.FullSequence).Distinct())));
            sb.Append("\t");

            // number of peptides
            if (!DisplayModsOnPeptides)
                sb.Append("" + AllPeptides.Select(p => p.BaseSequence).Distinct().Count());
            else
                sb.Append("" + AllPeptides.Select(p => p.FullSequence).Distinct().Count());
            sb.Append("\t");

            // number of unique peptides
            if (!DisplayModsOnPeptides)
                sb.Append("" + UniquePeptides.Select(p => p.BaseSequence).Distinct().Count());
            else
                sb.Append("" + UniquePeptides.Select(p => p.FullSequence).Distinct().Count());
            sb.Append("\t");

            // sequence coverage percent
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", SequenceCoveragePercent.Select(p => string.Format("{0:0}" + "%", (p * 100))))));
            sb.Append("\t");

            // sequence coverage
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", SequenceCoverageDisplayList)));
            sb.Append("\t");

            // sequence coverage with mods
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", SequenceCoverageDisplayListWithMods)));
            sb.Append("\t");

            //Detailed mods information list
            sb.Append(GlobalVariables.CheckLengthOfOutput(string.Join("|", ModsInfo)));
            sb.Append("\t");

            // MS1 intensity (retrieved from FlashLFQ in the SearchTask)
            if (IntensitiesByFile != null && FilesForQuantification != null)
            {
                foreach (var file in FilesForQuantification)
                {
                    if (IntensitiesByFile[file] > 0)
                    {
                        sb.Append(IntensitiesByFile[file]);
                    }
                    else
                    {
                        sb.Append("");
                    }
                    sb.Append("\t");
                }
            }

            // number of PSMs for listed peptides
            sb.Append("" + AllPsmsBelowOnePercentFDR.Count);
            sb.Append("\t");

            // isDecoy
            if (IsDecoy)
                sb.Append("D");
            else if (IsContaminant)
                sb.Append("C");
            else
                sb.Append("T");
            sb.Append("\t");

            // cumulative target
            sb.Append(CumulativeTarget);
            sb.Append("\t");

            // cumulative decoy
            sb.Append(CumulativeDecoy);
            sb.Append("\t");

            // q value
            sb.Append(QValue);
            sb.Append("\t");

            // best peptide score
            sb.Append(BestPeptideScore);
            sb.Append("\t");

            // best peptide q value
            sb.Append(BestPeptideQValue);
            sb.Append("\t");

            return sb.ToString();
        }

        // this method is only used internally, to make protein grouping faster
        // this is NOT an output and is NOT used for protein FDR calculations
        public void Score()
        {
            // sum the scores of the best PSM per base sequence
            ProteinGroupScore = AllPsmsBelowOnePercentFDR.GroupBy(p => p.BaseSequence).Select(p => p.Select(x => x.Score).Max()).Sum();
        }

        public void CalculateSequenceCoverage()
        {
            var proteinsWithUnambigSeqPsms = new Dictionary<Protein, List<PeptideWithSetModifications>>();
            var proteinsWithPsmsWithLocalizedMods = new Dictionary<Protein, List<PeptideWithSetModifications>>();

            foreach (var protein in Proteins)
            {
                proteinsWithUnambigSeqPsms.Add(protein, new List<PeptideWithSetModifications>());
                proteinsWithPsmsWithLocalizedMods.Add(protein, new List<PeptideWithSetModifications>());
            }

            foreach (var psm in AllPsmsBelowOnePercentFDR)
            {
                // null BaseSequence means that the amino acid sequence is ambiguous; do not use these to calculate sequence coverage
                if (psm.BaseSequence != null)
                {
                    var peptides = psm.BestMatchingPeptides.Select(p => p.Peptide);
                    foreach (var peptide in peptides)
                    {
                        // might be unambiguous but also shared; make sure this protein group contains this peptide+protein combo
                        if (Proteins.Contains(peptide.Protein))
                        {
                            proteinsWithUnambigSeqPsms[peptide.Protein].Add(peptide);

                            // null FullSequence means that mods were not successfully localized; do not display them on the sequence coverage mods info
                            if (psm.FullSequence != null)
                            {
                                proteinsWithPsmsWithLocalizedMods[peptide.Protein].Add(peptide);
                            }
                        }
                    }
                }
            }

            foreach (var protein in ListOfProteinsOrderedByAccession)
            {
                bool errorResult = false;
                var sequenceCoverageDisplay = protein.BaseSequence.ToLower(CultureInfo.InvariantCulture);
                HashSet<int> coveredOneBasedResidues = new HashSet<int>();

                // get residue numbers of each peptide in the protein and identify them as observed if the sequence is unambiguous
                foreach (var peptide in proteinsWithUnambigSeqPsms[protein])
                {
                    string sequenceExtractedFromProtein = "";
                    for (int i = peptide.OneBasedStartResidueInProtein; i <= peptide.OneBasedEndResidueInProtein; i++)
                    {
                        // check for bugs in sequence coverage; make sure we have the right amino acids!
                        sequenceExtractedFromProtein += sequenceCoverageDisplay[i - 1];
                        coveredOneBasedResidues.Add(i);
                    }

                    if (!sequenceExtractedFromProtein.ToUpper().Equals(peptide.BaseSequence))
                        errorResult = true;
                }

                // calculate sequence coverage percent
                double seqCoveragePercent = (double)coveredOneBasedResidues.Count / protein.Length;
                if (seqCoveragePercent > 1)
                    errorResult = true;

                // add the percent coverage or NaN if there was an error
                if (!errorResult)
                    SequenceCoveragePercent.Add(seqCoveragePercent);
                else
                    SequenceCoveragePercent.Add(double.NaN);

                // convert the observed amino acids to upper case if they are unambiguously observed
                var coverageArray = sequenceCoverageDisplay.ToCharArray();
                foreach (var obsResidueLocation in coveredOneBasedResidues)
                    coverageArray[obsResidueLocation - 1] = char.ToUpper(coverageArray[obsResidueLocation - 1]);
                sequenceCoverageDisplay = new string(coverageArray);

                // check to see if there was an errored result; if not, add the coverage display
                if (!errorResult)
                    SequenceCoverageDisplayList.Add(sequenceCoverageDisplay);
                else
                    SequenceCoverageDisplayList.Add("Error calculating sequence coverage");

                // put mods in the sequence coverage display
                if (!errorResult)
                {
                    // get mods to display in sequence (only unambiguously identified mods)
                    var modsOnThisProtein = new HashSet<KeyValuePair<int, Modification>>();
                    foreach (var pep in proteinsWithPsmsWithLocalizedMods[protein])
                    {
                        foreach (var mod in pep.AllModsOneIsNterminus)
                        {
                            if (!mod.Value.ModificationType.Contains("PeptideTermMod") && !mod.Value.ModificationType.Contains("Common Variable") && !mod.Value.ModificationType.Contains("Common Fixed"))
                                modsOnThisProtein.Add(new KeyValuePair<int, Modification>(pep.OneBasedStartResidueInProtein + mod.Key - 2, mod.Value));
                        }
                    }

                    var temp1 = modsOnThisProtein.OrderBy(p => p.Key).ToList();

                    foreach (var mod in temp1)
                    {
                        if (mod.Value.LocationRestriction.Equals("N-terminal."))
                            sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(0, "[" + mod.Value.IdWithMotif + "]-");
                        else if (mod.Value.LocationRestriction.Equals("Anywhere."))
                        {
                            int modStringIndex = sequenceCoverageDisplay.Length - (protein.Length - mod.Key);
                            sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(modStringIndex, "[" + mod.Value.IdWithMotif + "]");
                        }
                        else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                            sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(sequenceCoverageDisplay.Length, "-[" + mod.Value.IdWithMotif + "]");
                    }

                    SequenceCoverageDisplayListWithMods.Add(sequenceCoverageDisplay);

                    if (modsOnThisProtein.Any())
                    {
                        // calculate spectral count percentage of modified observation
                        string tempModStrings = ""; //The whole string
                        List<int> tempPepModTotals = new List<int>();  //The List of (For one mod, The Modified Pep Num)
                        List<int> tempPepTotals = new List<int>(); //The List of (For one mod, The total Pep Num)
                        List<string> tempPepModValues = new List<string>(); //The List of (For one mod, the Modified Name)
                        List<int> tempModIndex = new List<int>(); //The Index of the modified position.

                        foreach (var pep in proteinsWithPsmsWithLocalizedMods[protein])
                        {
                            foreach (var mod in pep.AllModsOneIsNterminus)
                            {
                                int tempPepNumTotal = 0; //For one mod, The total Pep Num
                                if (!mod.Value.ModificationType.Contains("Common Variable") && !mod.Value.ModificationType.Contains("Common Fixed") && !mod.Value.LocationRestriction.Equals(ModLocationOnPeptideOrProtein.PepC) && !mod.Value.LocationRestriction.Equals(ModLocationOnPeptideOrProtein.NPep))
                                {
                                    int tempIndexInProtein;
                                    if (mod.Value.LocationRestriction.Equals("N-terminal."))
                                        tempIndexInProtein = 1;
                                    else if (mod.Value.LocationRestriction.Equals("Anywhere."))
                                    {
                                        tempIndexInProtein = pep.OneBasedStartResidueInProtein + mod.Key - 2;
                                    }
                                    else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                                        tempIndexInProtein = protein.Length;
                                    else
                                        // In case it's a peptide terminal mod, skip!
                                        // we don't want this annotated in the protein's modifications
                                        continue;

                                    if (tempModIndex.Contains(tempIndexInProtein) && tempPepModValues[tempModIndex.IndexOf(tempIndexInProtein)] == mod.Value.IdWithMotif)
                                    {
                                        tempPepModTotals[tempModIndex.IndexOf(tempIndexInProtein)] += 1;
                                    }
                                    else
                                    {
                                        tempModIndex.Add(tempIndexInProtein);
                                        foreach (var pept in proteinsWithPsmsWithLocalizedMods[protein])
                                        {
                                            if (tempIndexInProtein >= pept.OneBasedStartResidueInProtein - (tempIndexInProtein == 1 ? 1 : 0) && tempIndexInProtein <= pept.OneBasedEndResidueInProtein)
                                            {
                                                tempPepNumTotal += 1;
                                            }
                                        }
                                        tempPepTotals.Add(tempPepNumTotal);
                                        tempPepModValues.Add(mod.Value.IdWithMotif);
                                        tempPepModTotals.Add(1);
                                    }
                                }
                            }
                        }
                        for (int i = 0; i < tempPepModTotals.Count; i++)
                        {
                            string tempString = ("#aa" + tempModIndex[i].ToString() + "[" + tempPepModValues[i].ToString() + ",info:occupancy=" + ((double)tempPepModTotals[i] / (double)tempPepTotals[i]).ToString("F2") + "(" + tempPepModTotals[i].ToString() + "/" + tempPepTotals[i].ToString() + ")" + "];");
                            tempModStrings += tempString;
                        }

                        if (!string.IsNullOrEmpty(tempModStrings))
                            ModsInfo.Add(tempModStrings);
                    }
                }
            }
        }

        public void MergeProteinGroupWith(ProteinGroup other)
        {
            this.Proteins.UnionWith(other.Proteins);
            this.AllPeptides.UnionWith(other.AllPeptides);
            this.UniquePeptides.UnionWith(other.UniquePeptides);
            this.AllPsmsBelowOnePercentFDR.UnionWith(other.AllPsmsBelowOnePercentFDR);
            other.ProteinGroupScore = 0;

            ListOfProteinsOrderedByAccession = Proteins.OrderBy(p => p.Accession).ToList();

            ProteinGroupName = string.Join("|", ListOfProteinsOrderedByAccession.Select(p => p.Accession));
        }

        public ProteinGroup ConstructSubsetProteinGroup(string fullFilePath)
        {
            var allPsmsForThisFile = new HashSet<PeptideSpectralMatch>(this.AllPsmsBelowOnePercentFDR.Where(p => p.FullFilePath.Equals(fullFilePath)));
            var allPeptidesForThisFile = new HashSet<PeptideWithSetModifications>(allPsmsForThisFile.SelectMany(p => p.BestMatchingPeptides.Select(v => v.Peptide)));
            var allUniquePeptidesForThisFile = new HashSet<PeptideWithSetModifications>(this.UniquePeptides.Intersect(allPeptidesForThisFile));

            ProteinGroup subsetPg = new ProteinGroup(this.Proteins, allPeptidesForThisFile, allUniquePeptidesForThisFile)
            {
                AllPsmsBelowOnePercentFDR = allPsmsForThisFile,
                DisplayModsOnPeptides = this.DisplayModsOnPeptides
            };

            SpectraFileInfo spectraFileInfo = null;
            if (FilesForQuantification != null)
            {
                spectraFileInfo = FilesForQuantification.Where(p => p.FullFilePathWithExtension == fullFilePath).First();
                subsetPg.FilesForQuantification = new List<SpectraFileInfo> { spectraFileInfo };
            }

            if (IntensitiesByFile == null)
            {
                subsetPg.IntensitiesByFile = null;
            }
            else
            {
                subsetPg.IntensitiesByFile = new Dictionary<SpectraFileInfo, double> { { spectraFileInfo, IntensitiesByFile[spectraFileInfo] } };
            }

            return subsetPg;
        }
    }
}