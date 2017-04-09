using Proteomics;
using System.Collections.Generic;
using System.Globalization;
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

        #region Internal Constructors

        internal ProteinGroup(HashSet<Protein> proteins, HashSet<PeptideWithSetModifications> peptides, HashSet<PeptideWithSetModifications> uniquePeptides)
        {
            Proteins = proteins;

            AllPeptides = peptides;
            UniquePeptides = uniquePeptides;
            RazorPeptides = new HashSet<PeptideWithSetModifications>();
            AllPsmsBelowOnePercentFDR = new HashSet<NewPsmWithFdr>();
            SequenceCoveragePercent = new List<double>();
            SequenceCoverageDisplayList = new List<string>();
            SequenceCoverageDisplayListWithMods = new List<string>();
            ProteinGroupScore = 0;
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

        #endregion Internal Constructors

        #region Public Properties

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("Protein Accession" + '\t');
                sb.Append("Gene" + '\t');
                sb.Append("Protein Full Name" + '\t');
                sb.Append("Number of proteins in group" + '\t');
                sb.Append("Unique peptides" + '\t');
                sb.Append("Shared peptides" + '\t');
                sb.Append("Razor peptides" + '\t');
                sb.Append("Number of peptides" + '\t');
                sb.Append("Number of unique peptides" + '\t');
                sb.Append("Sequence coverage %" + '\t');
                sb.Append("Sequence coverage" + '\t');
                sb.Append("Sequence coverage w Mods" + '\t');
                sb.Append("Modification Info List" + "\t");
                sb.Append("Intensity" + '\t');
                sb.Append("Number of PSMs" + '\t');
                sb.Append("Summed MetaMorpheus Score" + '\t');
                sb.Append("Decoy/Contaminant/Target" + '\t');
                sb.Append("Cumulative Target" + '\t');
                sb.Append("Cumulative Decoy" + '\t');
                sb.Append("Q-Value (%)");
                return sb.ToString();
            }
        }

        public double ProteinGroupScore { get; set; }
        public HashSet<Protein> Proteins { get; set; }
        public HashSet<PeptideWithSetModifications> AllPeptides { get; set; }
        public HashSet<PeptideWithSetModifications> UniquePeptides { get; set; }
        public HashSet<PeptideWithSetModifications> RazorPeptides { get; set; }
        public HashSet<NewPsmWithFdr> AllPsmsBelowOnePercentFDR { get; set; }
        public List<double> SequenceCoveragePercent { get; private set; }
        public List<string> SequenceCoverageDisplayList { get; private set; }
        public List<string> SequenceCoverageDisplayListWithMods { get; private set; }
        public double QValue { get; set; }
        public int CumulativeTarget { get; set; }
        public int CumulativeDecoy { get; set; }
        public double[] Intensity { get; set; }
        public bool DisplayModsOnPeptides { get; set; }
        public List<string> ModsInfo { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            // list of protein accession numbers
            sb.Append(string.Join("|", new HashSet<string>(Proteins.Select(p => p.Accession))));
            sb.Append("\t");

            // genes
            var genes = new HashSet<string>(Proteins.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault()));
            sb.Append(string.Join("|", genes));
            sb.Append("\t");

            // list of protein names
            sb.Append(string.Join("|", new HashSet<string>(Proteins.Select(p => p.FullName))));
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + Proteins.Count);
            sb.Append("\t");

            // list of unique peptides
            if (!DisplayModsOnPeptides)
                sb.Append(string.Join("|", new HashSet<string>(UniquePeptides.Select(p => p.BaseSequence))));
            else
                sb.Append(string.Join("|", new HashSet<string>(UniquePeptides.Select(p => p.Sequence))));
            sb.Append("\t");

            // list of shared peptides
            var SharedPeptides = AllPeptides.Except(UniquePeptides);
            if (!DisplayModsOnPeptides)
                sb.Append(string.Join("|", new HashSet<string>(SharedPeptides.Select(p => p.BaseSequence))));
            else
                sb.Append(string.Join("|", new HashSet<string>(SharedPeptides.Select(p => p.Sequence))));
            sb.Append("\t");

            // list of razor peptides
            if (!DisplayModsOnPeptides)
                sb.Append(string.Join("|", new HashSet<string>(RazorPeptides.Select(p => p.BaseSequence))));
            else
                sb.Append(string.Join("|", new HashSet<string>(RazorPeptides.Select(p => p.Sequence))));
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
            sb.Append(string.Join("|", SequenceCoveragePercent.Select(p => string.Format("{0:0}" + "%", (p * 100)))));
            sb.Append("\t");

            // sequence coverage
            sb.Append(string.Join("|", SequenceCoverageDisplayList));
            sb.Append("\t");

            // sequence coverage with mods
            sb.Append(string.Join("|", SequenceCoverageDisplayListWithMods));
            sb.Append("\t");

            //Detailed mods information list
            bool empty = false;
            if (ModsInfo != null)
            {foreach(var modsInfo in ModsInfo) { if (modsInfo != "") { empty = true; } }}
            if(empty == true) { sb.Append(string.Join("|", ModsInfo)); }
            sb.Append("\t");

            // summed MS1 intensity of razor and unique peptides
            sb.Append(string.Join("|", Intensity));
            sb.Append("\t");

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

            return sb.ToString();
        }

        public void Score()
        {
            // sum the scores of the best PSM per base sequence
            ProteinGroupScore = AllPsmsBelowOnePercentFDR.GroupBy(p => p.thisPSM.BaseSequence).Select(p => p.Select(x => x.thisPSM.Score).Max()).Sum();
        }

        public void CalculateSequenceCoverage()
        {
            var peptidesGroupedByProtein = AllPeptides.GroupBy(p => p.Protein);

            foreach (var protein in Proteins)
            {
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
                        string tempModStrings = ""; //The whole string 
                        List<int> tempPepModTotals = new List<int>();  //The List of (For one mod, The Modified Pep Num)
                        List<int> tempPepTotals = new List<int>(); //The List of (For one mod, The total Pep Num)
                        List<string> tempPepModValues = new List<string>(); //The List of (For one mod, the Modified Name)
                        List<int> tempModIndex = new List<int>(); //The Index of the modified position.

                        foreach (var pep in peptideGroup)
                        {            
                            foreach (var mod in pep.allModsOneIsNterminus)
                            {
                                int tempPepNumTotal = 1; //For one mod, The total Pep Num
                                int temp;
                                if (!mod.Value.modificationType.Contains("PeptideTermMod") && !mod.Value.modificationType.Contains("Common Variable") && !mod.Value.modificationType.Contains("Common Fixed"))
                                {     
                                    modsOnThisProtein.Add(new KeyValuePair<int, ModificationWithMass>(pep.OneBasedStartResidueInProtein + mod.Key - 2, mod.Value));

                                    if (tempModIndex.Contains(pep.OneBasedStartResidueInProtein + mod.Key - 2) && 
                                        tempPepModValues[tempModIndex.IndexOf(pep.OneBasedStartResidueInProtein + mod.Key - 2)]==mod.Value.id.ToString())
                                    { tempPepModTotals[tempModIndex.IndexOf(pep.OneBasedStartResidueInProtein + mod.Key - 2)] +=1; }
                                    else
                                    {
                                        temp = pep.OneBasedStartResidueInProtein + mod.Key - 2;
                                        tempModIndex.Add(temp);
                                        foreach (var pept in peptideGroup)
                                        {
                                            if (temp >= (pept.OneBasedStartResidueInProtein-1) && temp <= (pept.OneBasedEndResidueInProtein-1))
                                            {tempPepNumTotal += 1;}
                                        }
                                        tempPepTotals.Add(tempPepNumTotal);
                                        tempPepModValues.Add(mod.Value.id.ToString());
                                        tempPepModTotals.Add(1);
                                    }
                                }
                            }        
                        }
                        for(int i= 0; i < tempPepModTotals.Count(); i++)
                        {
                            string tempString = ("#aa" + tempModIndex[i].ToString() + "[" + tempPepModValues[i] + "|info:occupancy=" + ((double)tempPepModTotals[i] / (double)tempPepTotals[i]).ToString("F2") + "("+ tempPepModTotals[i].ToString()+"/"+ tempPepTotals[i].ToString() + ")"+ "]");
                            tempModStrings += tempString;
                        }
                        ModsInfo.Add(tempModStrings);
                        var temp1 = modsOnThisProtein.OrderBy(p => p.Key).ToList();

                        foreach (var mod in temp1)
                        {
                            int modStringIndex = seq.Length - (protein.Length - mod.Key);
                            seq = seq.Insert(modStringIndex, "[" + mod.Value.id + "]");
                        }

                        SequenceCoverageDisplayListWithMods.Add(seq);
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
        }

        public void Quantify()
        {
            Intensity = new double[AllPsmsBelowOnePercentFDR.First().thisPSM.newPsm.quantIntensity.Length];

            var psmsGroupedByBaseSequence = AllPsmsBelowOnePercentFDR.GroupBy(p => p.thisPSM.BaseSequence);
            var acceptedModTypesForProteinQuantification = new HashSet<string> { "Oxidation of M", "Carbamidomethyl of C", "TMT_tag_lysine", "TMT_tag_terminal" };

            foreach (var psmGroup in psmsGroupedByBaseSequence)
            {
                var psmsForThisBaseSeq = psmGroup.ToList();
                var psmsToIgnore = new List<NewPsmWithFdr>();

                // remove shared non-razor peptides
                foreach(var psm in psmGroup)
                {
                    var uniques = psm.thisPSM.PeptidesWithSetModifications.Intersect(UniquePeptides);
                    var razors = psm.thisPSM.PeptidesWithSetModifications.Intersect(RazorPeptides);

                    if (!uniques.Any() && !razors.Any())
                        psmsToIgnore.Add(psm);
                }

                psmsForThisBaseSeq = psmsForThisBaseSeq.Except(psmsToIgnore).ToList();

                // remove modified peptides that aren't used for quantification
                foreach (var psm in psmsForThisBaseSeq)
                {
                    var unacceptableModsForThisPsm = psm.thisPSM.PeptidesWithSetModifications.SelectMany(p => p.allModsOneIsNterminus.Values).Select(p => p.id).Except(acceptedModTypesForProteinQuantification);
                    if (unacceptableModsForThisPsm.Any())
                        psmsToIgnore.Add(psm);
                }

                psmsForThisBaseSeq = psmsForThisBaseSeq.Except(psmsToIgnore).ToList();

                if (psmsForThisBaseSeq.Any())
                    for(int i = 0; i < Intensity.Length; i++)
                        Intensity[i] += psmsForThisBaseSeq.Select(p => p.thisPSM.newPsm.quantIntensity[i]).Max();
            }
        }

        #endregion Public Methods

    }
}