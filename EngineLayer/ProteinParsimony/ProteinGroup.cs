using Proteomics;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class ProteinGroup
    {
        #region Public Fields

        public readonly bool isDecoy;
        public readonly bool isContaminant;

        #endregion Public Fields

        #region Public Constructors

        public ProteinGroup(HashSet<Protein> proteins, HashSet<PeptideWithSetModifications> peptides, HashSet<PeptideWithSetModifications> uniquePeptides)
        {
            Proteins = proteins;
            ProteinGroupName = string.Join("|", new HashSet<string>(Proteins.Select(p => p.Accession)));
            AllPeptides = peptides;
            UniquePeptides = uniquePeptides;
            AllPsmsBelowOnePercentFDR = new HashSet<Psm>();
            SequenceCoveragePercent = new List<double>();
            SequenceCoverageDisplayList = new List<string>();
            SequenceCoverageDisplayListWithMods = new List<string>();
            ProteinGroupScore = 0;
            BestPeptideScore = 0;
            QValue = 0;
            isDecoy = false;
            isContaminant = false;
            ModsInfo = new List<string>();

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.IsDecoy)
                {
                    isDecoy = true;
                    break;
                }
                if (protein.IsContaminant)
                {
                    isContaminant = true;
                    break;
                }
            }
        }

        #endregion Public Constructors

        #region Public Properties

        public static string[] FilesForQuantification { get; set; }

        public HashSet<Protein> Proteins { get; set; }

        public string ProteinGroupName { get; private set; }

        public double ProteinGroupScore { get; set; }

        public HashSet<PeptideWithSetModifications> AllPeptides { get; set; }

        public HashSet<PeptideWithSetModifications> UniquePeptides { get; set; }

        public HashSet<Psm> AllPsmsBelowOnePercentFDR { get; set; }

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
        public double[] IntensitiesByFile { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static string GetTabSeparatedHeader(bool fileSpecificHeader)
        {
            var sb = new StringBuilder();
            sb.Append("Protein Accession" + '\t');
            sb.Append("Gene" + '\t');
            sb.Append("Protein Full Name" + '\t');
            sb.Append("Unmodified Mass" + '\t');
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
                if (!fileSpecificHeader)
                    for (int i = 0; i < FilesForQuantification.Length; i++)
                        sb.Append("Intensity_" + Path.GetFileNameWithoutExtension(FilesForQuantification[i]) + '\t');
                else
                    sb.Append("Intensity" + '\t');
            }
            sb.Append("Number of PSMs" + '\t');
            sb.Append("Summed Score" + '\t');

            sb.Append("Decoy/Contaminant/Target" + '\t');
            sb.Append("Cumulative Target" + '\t');
            sb.Append("Cumulative Decoy" + '\t');
            sb.Append("QValue" + '\t');

            sb.Append("Best Peptide Score" + '\t');
            sb.Append("Best Peptide QValue");
            return sb.ToString();
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            // list of protein accession numbers
            sb.Append(ProteinGroupName);
            sb.Append("\t");

            // genes
            var genes = new HashSet<string>(Proteins.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault()));
            sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", genes)));
            sb.Append("\t");

            // list of protein names
            sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", new HashSet<string>(Proteins.Select(p => p.FullName)))));
            sb.Append("\t");

            // list of masses
            IDigestionParams digestionParams = new TDdigest();
            sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", new HashSet<double>(Proteins.Select(p => p.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First().MonoisotopicMass)))));
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + Proteins.Count);
            sb.Append("\t");

            // list of unique peptides
            if (!DisplayModsOnPeptides)
                sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", new HashSet<string>(UniquePeptides.Select(p => p.BaseSequence)))));
            else
                sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", new HashSet<string>(UniquePeptides.Select(p => p.Sequence)))));
            sb.Append("\t");

            // list of shared peptides
            var SharedPeptides = AllPeptides.Except(UniquePeptides);
            if (!DisplayModsOnPeptides)
                sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", new HashSet<string>(SharedPeptides.Select(p => p.BaseSequence)))));
            else
                sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", new HashSet<string>(SharedPeptides.Select(p => p.Sequence)))));
            sb.Append("\t");

            // number of peptides
            if (!DisplayModsOnPeptides)
                sb.Append("" + new HashSet<string>(AllPeptides.Select(p => p.BaseSequence)).Count);
            else
                sb.Append("" + new HashSet<string>(AllPeptides.Select(p => p.Sequence)).Count);
            sb.Append("\t");

            // number of unique peptides
            if (!DisplayModsOnPeptides)
                sb.Append("" + new HashSet<string>(UniquePeptides.Select(p => p.BaseSequence)).Count);
            else
                sb.Append("" + new HashSet<string>(UniquePeptides.Select(p => p.Sequence)).Count);
            sb.Append("\t");

            // sequence coverage percent
            sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", SequenceCoveragePercent.Select(p => string.Format("{0:0}" + "%", (p * 100))))));
            sb.Append("\t");

            // sequence coverage
            sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", SequenceCoverageDisplayList)));
            sb.Append("\t");

            // sequence coverage with mods
            sb.Append(GlobalVariables.CheckLengthOfPSM(string.Join("|", SequenceCoverageDisplayListWithMods)));
            sb.Append("\t");

            //Detailed mods information list
            var modsInfoString = string.Join("|", ModsInfo);
            sb.Append(GlobalVariables.CheckLengthOfPSM(modsInfoString.Length < 32767 ? modsInfoString : "Too many mods to display"));
            sb.Append("\t");

            // MS1 intensity (retrieved from FlashLFQ in the SearchTask)
            if (IntensitiesByFile != null && FilesForQuantification != null)
            {
                for (int i = 0; i < IntensitiesByFile.Length; i++)
                {
                    if (IntensitiesByFile[i] > 0)
                        sb.Append(IntensitiesByFile[i]);
                    else
                        sb.Append("");
                    sb.Append("\t");
                }
            }

            // number of PSMs for listed peptides
            sb.Append("" + AllPsmsBelowOnePercentFDR.Count);
            sb.Append("\t");

            // summed metamorpheus score
            sb.Append(ProteinGroupScore);
            sb.Append("\t");

            // isDecoy
            if (isDecoy)
                sb.Append("D");
            else if (isContaminant)
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

        public void Score()
        {
            // sum the scores of the best PSM per base sequence
            ProteinGroupScore = AllPsmsBelowOnePercentFDR.GroupBy(p => p.BaseSequence).Select(p => p.Select(x => x.Score).Max()).Sum();
        }

        public void CalculateSequenceCoverage()
        {
            var peptidesGroupedByProtein = AllPeptides.GroupBy(p => p.Protein);

            var proteinsWithPsms = new Dictionary<Protein, List<PeptideWithSetModifications>>();

            foreach (var psm in AllPsmsBelowOnePercentFDR)
            {
                foreach (var pepWithSetMods in psm.CompactPeptides.SelectMany(b => b.Value.Item2))
                {
                    if (proteinsWithPsms.TryGetValue(pepWithSetMods.Protein, out List<PeptideWithSetModifications> temp))
                        temp.Add(pepWithSetMods);
                    else
                        proteinsWithPsms.Add(pepWithSetMods.Protein, new List<PeptideWithSetModifications> { pepWithSetMods });
                }
            }

            foreach (var protein in Proteins)
            {
                bool ModExist = false;
                foreach (var peptideGroup in peptidesGroupedByProtein)
                {
                    if (protein == peptideGroup.Key)
                    {
                        HashSet<int> coveredResidues = new HashSet<int>();

                        foreach (var peptide in peptideGroup)
                        {
                            if (peptide.Protein == peptideGroup.Key)
                            {
                                for (int i = peptide.OneBasedStartResidueInProtein; i <= peptide.OneBasedEndResidueInProtein; i++)
                                    coveredResidues.Add(i);
                            }
                        }

                        double sequenceCoverageHere = (double)coveredResidues.Count / peptideGroup.Key.Length;
                        SequenceCoveragePercent.Add(sequenceCoverageHere);

                        var sequenceCoverageDisplay = peptideGroup.Key.BaseSequence.ToLower(CultureInfo.InvariantCulture);
                        var coverageArray = sequenceCoverageDisplay.ToCharArray();

                        foreach (var residue in coveredResidues)
                        {
                            var temp = char.ToUpper(coverageArray[residue - 1]);
                            coverageArray[residue - 1] = temp;
                        }

                        var seq = new string(coverageArray);
                        SequenceCoverageDisplayList.Add(seq);

                        // get sequence coverage display with mods
                        var modsOnThisProtein = new HashSet<KeyValuePair<int, ModificationWithMass>>();
                        foreach (var pep in peptideGroup)
                        {
                            foreach (var mod in pep.allModsOneIsNterminus)
                            {
                                if (!mod.Value.modificationType.Contains("PeptideTermMod") && !mod.Value.modificationType.Contains("Common Variable") && !mod.Value.modificationType.Contains("Common Fixed"))
                                    modsOnThisProtein.Add(new KeyValuePair<int, ModificationWithMass>(pep.OneBasedStartResidueInProtein + mod.Key - 2, mod.Value));
                            }
                        }

                        var temp1 = modsOnThisProtein.OrderBy(p => p.Key).ToList();

                        foreach (var mod in temp1)
                        {
                            if (mod.Value.terminusLocalization.Equals(TerminusLocalization.NProt))
                                seq = seq.Insert(0, "[" + mod.Value.id + "]-");
                            else if (mod.Value.terminusLocalization.Equals(TerminusLocalization.Any))
                            {
                                int modStringIndex = seq.Length - (protein.Length - mod.Key);
                                seq = seq.Insert(modStringIndex, "[" + mod.Value.id + "]");
                            }
                            else if (mod.Value.terminusLocalization.Equals(TerminusLocalization.ProtC))
                                seq = seq.Insert(seq.Length, "-[" + mod.Value.id + "]");
                        }

                        SequenceCoverageDisplayListWithMods.Add(seq);

                        if (modsOnThisProtein.Any())
                        {
                            ModExist = true;
                        }
                    }
                }

                foreach (var aproteinWithPsms in proteinsWithPsms)
                {
                    if (aproteinWithPsms.Key == protein && ModExist == true)
                    {
                        string tempModStrings = ""; //The whole string
                        List<int> tempPepModTotals = new List<int>();  //The List of (For one mod, The Modified Pep Num)
                        List<int> tempPepTotals = new List<int>(); //The List of (For one mod, The total Pep Num)
                        List<string> tempPepModValues = new List<string>(); //The List of (For one mod, the Modified Name)
                        List<int> tempModIndex = new List<int>(); //The Index of the modified position.
                        foreach (var pep in aproteinWithPsms.Value)
                        {
                            foreach (var mod in pep.allModsOneIsNterminus)
                            {
                                int tempPepNumTotal = 0; //For one mod, The total Pep Num
                                if (!mod.Value.modificationType.Contains("Common Variable") && !mod.Value.modificationType.Contains("Common Fixed") && !mod.Value.terminusLocalization.Equals(TerminusLocalization.PepC) && !mod.Value.terminusLocalization.Equals(TerminusLocalization.NPep))
                                {
                                    int tempIndexInProtein;
                                    if (mod.Value.terminusLocalization.Equals(TerminusLocalization.NProt))
                                        tempIndexInProtein = 1;
                                    else if (mod.Value.terminusLocalization.Equals(TerminusLocalization.Any))
                                    {
                                        tempIndexInProtein = pep.OneBasedStartResidueInProtein + mod.Key - 2;
                                    }
                                    else if (mod.Value.terminusLocalization.Equals(TerminusLocalization.ProtC))
                                        tempIndexInProtein = protein.Length;
                                    else
                                        // In case it's a peptide mod, skip!
                                        break;

                                    if (tempModIndex.Contains(tempIndexInProtein) && tempPepModValues[tempModIndex.IndexOf(tempIndexInProtein)] == mod.Value.id)
                                    {
                                        tempPepModTotals[tempModIndex.IndexOf(tempIndexInProtein)] += 1;
                                    }
                                    else
                                    {
                                        tempModIndex.Add(tempIndexInProtein);
                                        foreach (var pept in aproteinWithPsms.Value)
                                        {
                                            if (tempIndexInProtein >= pept.OneBasedStartResidueInProtein - (tempIndexInProtein == 1 ? 1 : 0) && tempIndexInProtein <= pept.OneBasedEndResidueInProtein)
                                                tempPepNumTotal += 1;
                                        }
                                        tempPepTotals.Add(tempPepNumTotal);
                                        tempPepModValues.Add(mod.Value.id);
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
            ProteinGroupName = string.Join("|", new HashSet<string>(Proteins.Select(p => p.Accession)));
        }

        public ProteinGroup ConstructSubsetProteinGroup(string fullFilePath)
        {
            var allPsmsForThisFile = new HashSet<Psm>(this.AllPsmsBelowOnePercentFDR.Where(p => p.FullFilePath.Equals(fullFilePath)));
            var allPeptidesForThisFile = new HashSet<PeptideWithSetModifications>(allPsmsForThisFile.SelectMany(p => p.CompactPeptides.SelectMany(b => b.Value.Item2)));
            var allUniquePeptidesForThisFile = new HashSet<PeptideWithSetModifications>(this.UniquePeptides.Intersect(allPeptidesForThisFile));

            ProteinGroup subsetPg = new ProteinGroup(this.Proteins, allPeptidesForThisFile, allUniquePeptidesForThisFile)
            {
                AllPsmsBelowOnePercentFDR = allPsmsForThisFile,
                DisplayModsOnPeptides = this.DisplayModsOnPeptides
            };

            subsetPg.IntensitiesByFile = IntensitiesByFile != null ? new[] { IntensitiesByFile[System.Array.IndexOf(FilesForQuantification, fullFilePath)] } : new double[1];

            return subsetPg;
        }

        #endregion Public Methods

        #region Private Classes

        private class TDdigest : IDigestionParams
        {
            #region Public Properties

            public int MaxMissedCleavages => 0;

            public int? MinPeptideLength => 0;

            public int? MaxPeptideLength => null;

            public InitiatorMethionineBehavior InitiatorMethionineBehavior => InitiatorMethionineBehavior.Retain;

            public int MaxModificationIsoforms => 1;

            public int MaxModsForPeptide => 0;

            public Protease Protease => GlobalVariables.ProteaseDictionary["top-down"];

            public bool SemiProteaseDigestion => false;

            public TerminusType TerminusTypeSemiProtease => throw new System.NotImplementedException();

            #endregion Public Properties
        }

        #endregion Private Classes
    }
}