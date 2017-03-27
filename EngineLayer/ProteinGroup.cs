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
            //sequenceCoverageDisplayListWithMods = new List<string>();
            ProteinGroupScore = 0;
            QValue = 0;
            isDecoy = false;
            isContaminant = false;

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
                //sb.Append("Razor peptides" + '\t');
                sb.Append("Number of peptides" + '\t');
                sb.Append("Number of unique peptides" + '\t');
                sb.Append("Sequence coverage %" + '\t');
                sb.Append("Sequence coverage" + '\t');
                //sb.Append("Sequence coverage w Localizable Mods" + '\t');
                //sb.Append("Intensity" + '\t');
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
        //public List<string> sequenceCoverageDisplayListWithMods { get; private set; }
        public double QValue { get; set; }
        public int CumulativeTarget { get; set; }
        public int CumulativeDecoy { get; set; }
        public double Intensity { get; private set; }

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
            sb.Append(string.Join("|", new HashSet<string>(UniquePeptides.Select(p => p.BaseSequence))));
            sb.Append("\t");

            // list of shared peptides
            var sharedPeptides = AllPeptides.Except(UniquePeptides);
            sb.Append(string.Join("|", new HashSet<string>(sharedPeptides.Select(p => p.BaseSequence))));
            sb.Append("\t");

            // list of razor peptides
            //sb.Append(string.Join("|", new HashSet<string>(RazorPeptides.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)))));
            //sb.Append("\t");

            // number of peptides
            sb.Append("" + new HashSet<string>(AllPeptides.Select(p => p.BaseSequence)).Count);
            sb.Append("\t");

            // number of unique peptides
            sb.Append("" + new HashSet<string>(UniquePeptides.Select(p => p.BaseSequence)).Count);
            sb.Append("\t");

            // sequence coverage percent
            sb.Append(string.Join("|", SequenceCoveragePercent.Select(p => string.Format("{0:0}" + "%", (p * 100)))));
            sb.Append("\t");

            // sequence coverage
            sb.Append(string.Join("|", SequenceCoverageDisplayList));
            sb.Append("\t");

            // summed MS1 intensity of razor and unique peptides
            //sb.Append(Intensity);
            //sb.Append("\t");

            // sequence coverage with mods
            //sb.Append(string.Join("|", sequenceCoverageDisplayListWithMods));
            //sb.Append("\t");

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
                        //var peptideModsBelow1FDR = new List<ModificationWithMass>();

                        foreach (var peptide in peptideGroup)
                        {
                            if (peptide.Protein == peptideGroup.Key)
                            {
                                for (int i = peptide.OneBasedStartResidueInProtein; i <= peptide.OneBasedEndResidueInProtein; i++)
                                    coveredResidues.Add(i);

                                //foreach(var mod in peptide.allModsOneIsNterminus)
                                //    peptideModsBelow1FDR.Add(mod.Value);
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

                        /*
                        foreach (var modsAtThisResidue in protein.OneBasedPossibleLocalizedModifications)
                        {
                            foreach (var mod in modsAtThisResidue.Value)
                            {
                                if (peptideModsBelow1FDR.Contains(mod))
                                {
                                    int modStringIndex = seq.Length - (protein.Length - modsAtThisResidue.Key);
                                    seq = seq.Insert(modStringIndex, "[" + mod.id + "]");
                                }
                            }
                        }
                        
                        sequenceCoverageDisplayListWithMods.Add(seq);
                        */
                    }
                }
            }
        }

        public void MergeProteinGroupWith(ProteinGroup other)
        {
            this.Proteins.UnionWith(other.Proteins);
            other.ProteinGroupScore = 0;
        }

        public void Quantify()
        {
            var groups = AllPsmsBelowOnePercentFDR.GroupBy(p => p.thisPSM.BaseSequence);
            var acceptedModTypesForProteinQuantification = new HashSet<string> { "Oxidation", "Carbamidomethyl", "Acetylation", "TMT_tag_lysine", "TMT_tag_terminal" };

            foreach (var group in groups)
            {
                var psmsForThisBaseSeq = group.ToList();
                var modsForThesePSMs = psmsForThisBaseSeq.Select(p => p.thisPSM.peptidesWithSetModifications.First().allModsOneIsNterminus.Values.ToList()).ToList();
                List<NewPsmWithFdr> psmsToIgnore = new List<NewPsmWithFdr>();

                for (int i = 0; i < modsForThesePSMs.Count; i++)
                {
                    var unacceptableMods = modsForThesePSMs[i].Select(p => p.id).Except(acceptedModTypesForProteinQuantification);
                    if (unacceptableMods.Any())
                        psmsToIgnore.Add(psmsForThisBaseSeq[i]);
                }

                psmsForThisBaseSeq = psmsForThisBaseSeq.Except(psmsToIgnore).ToList();

                if(psmsForThisBaseSeq.Any())
                    Intensity += psmsForThisBaseSeq.Select(p => p.thisPSM.newPsm.apexIntensity).Max();
            }
        }

        #endregion Public Methods

    }
}