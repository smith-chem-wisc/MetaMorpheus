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

            TotalPeptideList = peptides;
            TotalUniquePeptideList = uniquePeptides;
            TotalPsmList = new HashSet<NewPsmWithFdr>();
            TotalPeptideWithSetModsList = new HashSet<PeptideWithSetModifications>();

            AllPsmsForStrictPeptideSequences = new HashSet<NewPsmWithFdr>();
            BestPsmPerBaseSeq = new HashSet<NewPsmWithFdr>();
            StrictPeptideList = new HashSet<CompactPeptide>();
            StrictRazorPeptideList = new HashSet<CompactPeptide>();
            StrictPeptideWithSetModsList = new HashSet<PeptideWithSetModifications>();
            StrictUniquePeptideList = new HashSet<CompactPeptide>();

            proteinGroupScore = 0;
            QValue = 0;
            isDecoy = false;

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
        public HashSet<CompactPeptide> TotalPeptideList { get; set; }
        public HashSet<NewPsmWithFdr> TotalPsmList { get; set; }
        public HashSet<CompactPeptide> TotalUniquePeptideList { get; set; }
        public HashSet<PeptideWithSetModifications> TotalPeptideWithSetModsList { get; set; }
        public HashSet<NewPsmWithFdr> BestPsmPerBaseSeq { get; private set; } // for scoring
        public HashSet<NewPsmWithFdr> AllPsmsForStrictPeptideSequences { get; private set; } // for PSMs per proteingroup output
        public HashSet<CompactPeptide> StrictPeptideList { get; private set; }
        public HashSet<CompactPeptide> StrictUniquePeptideList { get; private set; }
        public HashSet<CompactPeptide> StrictRazorPeptideList { get; private set; }
        public HashSet<PeptideWithSetModifications> StrictPeptideWithSetModsList { get; private set; }
        public List<double> sequenceCoveragePercent { get; private set; }
        public List<string> sequenceCoverageDisplayList { get; private set; }
        public double QValue { get; set; }
        public int cumulativeTarget { get; set; }
        public int cumulativeDecoy { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            // list of protein accession numbers
            sb.Append(string.Join(" | ", new List<string>(Proteins.Select(p => p.Accession))));
            sb.Append("\t");

            // genes
            var genes = new List<string>(Proteins.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault()));
            sb.Append(string.Join(" | ", genes));
            sb.Append("\t");

            // list of protein names
            sb.Append(string.Join(" | ", new List<string>(Proteins.Select(p => p.FullName))));
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + Proteins.Count);
            sb.Append("\t");

            // list of unique peptides
            sb.Append(string.Join(" | ", new HashSet<string>(StrictUniquePeptideList.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)))));
            sb.Append("\t");

            // list of shared peptides
            var sharedPeptides = StrictPeptideList.Except(TotalUniquePeptideList);
            sb.Append(string.Join(" | ", new HashSet<string>(sharedPeptides.Select(p => System.Text.Encoding.UTF8.GetString(p.BaseSequence)))));
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
            sb.Append("" + StrictPeptideList.Count);
            sb.Append("\t");

            // number of unique peptides
            sb.Append("" + StrictUniquePeptideList.Count);
            sb.Append("\t");

            // sequence coverage percent
            sb.Append(string.Join(" | ", sequenceCoveragePercent.Select(p => string.Format("{0:0}" + "%", (p * 100)))));
            sb.Append("\t");

            // sequence coverage
            sb.Append(string.Join(" | ", sequenceCoverageDisplayList));
            sb.Append("\t");

            // number of PSMs for listed peptides
            sb.Append("" + AllPsmsForStrictPeptideSequences.Count);
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

        public void ScoreThisProteinGroup(List<ModificationWithMass> variableModifications, List<ModificationWithMass> localizeableModifications, List<ModificationWithMass> fixedModifications)
        {
            // find the best psm per base sequence (peptide FDR must be <1%) for scoring
            Dictionary<string, NewPsmWithFdr> peptideBaseSeqToBestPsmMatching = new Dictionary<string, NewPsmWithFdr>();
            foreach (var psm in TotalPsmList)
            {
                if (psm.qValue < 0.01)
                {
                    CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications, fixedModifications);
                    string baseSeq = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    if (peptideBaseSeqToBestPsmMatching.ContainsKey(baseSeq))
                    {
                        NewPsmWithFdr psmHere;
                        peptideBaseSeqToBestPsmMatching.TryGetValue(baseSeq, out psmHere);
                        if (psm.thisPSM.Score > psmHere.thisPSM.Score)
                        {
                            peptideBaseSeqToBestPsmMatching[baseSeq] = psm;
                        }
                    }
                    else
                    {
                        peptideBaseSeqToBestPsmMatching.Add(baseSeq, psm);
                    }
                }
            }

            // create StrictPsmList (best psm per base seq)
            foreach (var kvp in peptideBaseSeqToBestPsmMatching)
            {
                BestPsmPerBaseSeq.Add(kvp.Value);
            }

            // create StrictPeptideList (only the CompactPeptide belonging to the best psm per base seq, and only if peptide FDR < 1%)
            foreach (var psm in BestPsmPerBaseSeq)
            {
                CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications, fixedModifications);
                StrictPeptideList.Add(peptide);
            }

            // create StrictUniquePeptideList
            foreach (var peptide in StrictPeptideList)
            {
                if (TotalUniquePeptideList.Contains(peptide))
                    StrictUniquePeptideList.Add(peptide);
            }

            // create AllPsmsForStrictPeptideSequences
            foreach (var peptide in StrictPeptideList)
            {
                string strictPepBaseSeq = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                // if the psm matches a base sequence of a peptide in the strict peptide list, add it to the psm list
                foreach (var psm in TotalPsmList)
                {
                    CompactPeptide psmPeptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications, fixedModifications);
                    string psmPeptideBaseSeq = string.Join("", psmPeptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    if (strictPepBaseSeq.Equals(psmPeptideBaseSeq))
                    {
                        AllPsmsForStrictPeptideSequences.Add(psm);
                    }
                }

                // TODO**: StrictPeptideWithSetModsList
                foreach (var pep in TotalPeptideWithSetModsList)
                {
                    if (pep.BaseSequence.Equals(strictPepBaseSeq))
                    {
                        StrictPeptideWithSetModsList.Add(pep);
                    }
                }
            }

            // score the protein group
            foreach (var psm in BestPsmPerBaseSeq)
            {
                proteinGroupScore += psm.thisPSM.Score;
            }
        }

        public void CalculateSequenceCoverage()
        {
            sequenceCoveragePercent = new List<double>();
            sequenceCoverageDisplayList = new List<string>();

            foreach (var protein in Proteins)
            {
                HashSet<int> coveredResidues = new HashSet<int>();

                foreach (var peptide in StrictPeptideWithSetModsList)
                {
                    if (peptide.Protein == protein)
                    {
                        for (int i = peptide.OneBasedStartResidueInProtein; i <= peptide.OneBasedEndResidueInProtein; i++)
                        {
                            coveredResidues.Add(i);
                        }
                    }
                }

                double sequenceCoverageHere = (double)coveredResidues.Count / protein.Length;
                sequenceCoveragePercent.Add(sequenceCoverageHere);

                var sequenceCoverageDisplay = protein.BaseSequence.ToLower(CultureInfo.InvariantCulture);
                var coverageArray = sequenceCoverageDisplay.ToCharArray();
                foreach (var residue in coveredResidues)
                {
                    var temp = char.ToUpper(coverageArray[residue - 1]);
                    coverageArray[residue - 1] = temp;
                }

                sequenceCoverageDisplay = new string(coverageArray);
                sequenceCoverageDisplayList.Add(sequenceCoverageDisplay);
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