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

        internal ProteinGroup(HashSet<Protein> proteins, HashSet<CompactPeptide> peptides, HashSet<CompactPeptide> uniquePeptides)
        {
            Proteins = proteins;

            AllPeptides = peptides;
            UniquePeptides = uniquePeptides;
            RazorPeptides = new HashSet<CompactPeptide>();
            AllPSMsBelow1PercentFDR = new HashSet<NewPsmWithFdr>();
            PeptidesWithSetMods = new HashSet<PeptideWithSetModifications>();
            sequenceCoveragePercent = new List<double>();
            sequenceCoverageDisplayList = new List<string>();
            //sequenceCoverageDisplayListWithMods = new List<string>();
            proteinGroupScore = 0;
            QValue = 0;
            isDecoy = false;
            isContaminant = false;

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.IsDecoy)
                    isDecoy = true;
                if (protein.IsContaminant)
                    isContaminant = true;
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
                sb.Append("Number of PSMs" + '\t');
                sb.Append("Summed MetaMorpheus Score" + '\t');
                sb.Append("Decoy/Contaminant/Target" + '\t');
                sb.Append("Cumulative Target" + '\t');
                sb.Append("Cumulative Decoy" + '\t');
                sb.Append("Q-Value (%)");
                return sb.ToString();
            }
        }

        public double proteinGroupScore { get; set; }
        public HashSet<Protein> Proteins { get; set; }
        public HashSet<CompactPeptide> AllPeptides { get; set; }
        public HashSet<NewPsmWithFdr> AllPSMsBelow1PercentFDR { get; set; } // all PSMs below 1% fdr
        public HashSet<CompactPeptide> UniquePeptides { get; set; }
        public HashSet<CompactPeptide> RazorPeptides { get; set; }
        public HashSet<PeptideWithSetModifications> PeptidesWithSetMods { get; set; }
        public List<double> sequenceCoveragePercent { get; private set; }
        public List<string> sequenceCoverageDisplayList { get; private set; }
        //public List<string> sequenceCoverageDisplayListWithMods { get; private set; }
        public double QValue { get; set; }
        public int cumulativeTarget { get; set; }
        public int cumulativeDecoy { get; set; }

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
            sb.Append(string.Join("|", new HashSet<string>(UniquePeptides.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)))));
            sb.Append("\t");

            // list of shared peptides
            var sharedPeptides = AllPeptides.Except(UniquePeptides);
            sb.Append(string.Join("|", new HashSet<string>(sharedPeptides.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)))));
            sb.Append("\t");

            // list of razor peptides
            //foreach (CompactPeptide razorPeptide in StrictRazorPeptideList)
            //{
            //    string peptideBaseSequence = string.Join("", razorPeptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));
            //
            //    sb.Append("" + peptideBaseSequence + " ;; ");
            //}
            //sb.Append("\t");

            // number of peptides
            sb.Append("" + new HashSet<string>(AllPeptides.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence))).Count);
            sb.Append("\t");

            // number of unique peptides
            sb.Append("" + new HashSet<string>(UniquePeptides.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence))).Count);
            sb.Append("\t");

            // sequence coverage percent
            sb.Append(string.Join("|", sequenceCoveragePercent.Select(p => string.Format("{0:0}" + "%", (p * 100)))));
            sb.Append("\t");

            // sequence coverage
            sb.Append(string.Join("|", sequenceCoverageDisplayList));
            sb.Append("\t");

            // sequence coverage with mods
            //sb.Append(string.Join("|", sequenceCoverageDisplayListWithMods));
            //sb.Append("\t");

            // number of PSMs for listed peptides
            sb.Append("" + AllPSMsBelow1PercentFDR.Count);
            sb.Append("\t");

            // summed metamorpheus score
            sb.Append(proteinGroupScore);
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
            sb.Append(cumulativeTarget);
            sb.Append("\t");

            // cumulative decoy
            sb.Append(cumulativeDecoy);
            sb.Append("\t");

            // q value
            sb.Append(QValue);
            sb.Append("\t");

            return sb.ToString();

        }

        public void Score()
        {
            // sum the scores of the best PSM per base sequence
            proteinGroupScore = AllPSMsBelow1PercentFDR.GroupBy(p => p.thisPSM.BaseSequence).Select(p => p.Select(x => x.thisPSM.Score).Max()).Sum();
        }

        public void CalculateSequenceCoverage()
        {
            var peptidesGroupedByProtein = PeptidesWithSetMods.GroupBy(p => p.Protein);

            foreach (var protein in Proteins)
            {
                foreach (var peptideGroup in peptidesGroupedByProtein)
                {
                    if (protein == peptideGroup.Key)
                    {
                        HashSet<int> coveredResidues = new HashSet<int>();
                        var peptideModsBelow1FDR = new List<ModificationWithMass>();

                        foreach (var peptide in peptideGroup)
                        {
                            if (peptide.Protein == peptideGroup.Key)
                            {
                                for (int i = peptide.OneBasedStartResidueInProtein; i <= peptide.OneBasedEndResidueInProtein; i++)
                                    coveredResidues.Add(i);

                                foreach(var mod in peptide.allModsOneIsNterminus)
                                    peptideModsBelow1FDR.Add(mod.Value);
                            }
                        }

                        double sequenceCoverageHere = (double)coveredResidues.Count / peptideGroup.Key.Length;
                        sequenceCoveragePercent.Add(sequenceCoverageHere);

                        var sequenceCoverageDisplay = peptideGroup.Key.BaseSequence.ToLower(CultureInfo.InvariantCulture);
                        var coverageArray = sequenceCoverageDisplay.ToCharArray();

                        foreach (var residue in coveredResidues)
                        {
                            var temp = char.ToUpper(coverageArray[residue - 1]);
                            coverageArray[residue - 1] = temp;
                        }

                        var seq = new string(coverageArray);
                        sequenceCoverageDisplayList.Add(seq);

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
            other.proteinGroupScore = 0;
        }

        #endregion Public Methods

    }
}