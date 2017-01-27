using OldInternalLogic;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class ProteinGroup
    {

        #region Public Fields

        public readonly bool isDecoy;
        public double proteinGroupScore { get; private set; }

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
            }
        }

        #endregion Internal Constructors

        #region Public Properties

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("Proteins in group" + '\t');
                sb.Append("Number of proteins in group" + '\t');
                sb.Append("Unique peptides" + '\t');
                sb.Append("Shared peptides" + '\t');
                sb.Append("Razor peptides" + '\t');
                sb.Append("Number of peptides" + '\t');
                sb.Append("Number of unique peptides" + '\t');
                sb.Append("Sequence coverage" + '\t');
                sb.Append("Number of PSMs" + '\t');
                sb.Append("Summed MetaMorpheus Score" + '\t');
                sb.Append("Decoy?" + '\t');
                sb.Append("Cumulative Target" + '\t');
                sb.Append("Cumulative Decoy" + '\t');
                sb.Append("Q-Value (%)");
                return sb.ToString();
            }
        }

        public HashSet<Protein> Proteins { get; set; }
        public HashSet<CompactPeptide> TotalPeptideList { get; set; }
        public HashSet<NewPsmWithFdr> TotalPsmList { get; set; }
        public HashSet<CompactPeptide> TotalUniquePeptideList { get; set; }
        public HashSet<PeptideWithSetModifications> TotalPeptideWithSetModsList { get; set; }
        public HashSet<CompactPeptide> StrictPeptideList { get; private set; }
        public HashSet<NewPsmWithFdr> BestPsmPerBaseSeq { get; private set; } // for scoring
        public HashSet<NewPsmWithFdr> AllPsmsForStrictPeptideSequences { get; private set; } // for PSMs per protein group output
        public HashSet<CompactPeptide> StrictUniquePeptideList { get; private set; }
        public HashSet<CompactPeptide> StrictRazorPeptideList { get; private set; }
        public HashSet<PeptideWithSetModifications> StrictPeptideWithSetModsList { get; private set; }
        public List<double> sequenceCoverage { get; private set; }
        public double QValue { get; set; }
        public int cumulativeTarget { get; set; }
        public int cumulativeDecoy { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();

            // list of proteins in the group
            foreach (Protein protein in Proteins)
                sb.Append("" + protein.FullDescription + " ;; ");
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + Proteins.Count());
            sb.Append("\t");

            // list of unique peptides
            foreach (CompactPeptide uniquePeptide in StrictUniquePeptideList)
            {
                string peptideBaseSequence = string.Join("", uniquePeptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                sb.Append("" + peptideBaseSequence + " ;; ");
            }
            sb.Append("\t");

            // list of shared peptides
            foreach (CompactPeptide sharedPeptide in StrictPeptideList)
            {
                if (!TotalUniquePeptideList.Contains(sharedPeptide))
                {
                    string peptideBaseSequence = string.Join("", sharedPeptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    sb.Append("" + peptideBaseSequence + " ;; ");
                }
            }
            sb.Append("\t");

            // list of razor peptides
            foreach (CompactPeptide razorPeptide in StrictRazorPeptideList)
            {
                string peptideBaseSequence = string.Join("", razorPeptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                sb.Append("" + peptideBaseSequence + " ;; ");
            }
            sb.Append("\t");

            // number of peptides
            sb.Append("" + StrictPeptideList.Count());
            sb.Append("\t");

            // number of unique peptides
            sb.Append("" + StrictUniquePeptideList.Count());
            sb.Append("\t");

            // sequence coverage
            foreach (double coverage in sequenceCoverage)
            {
                double coverage1 = coverage * 100;
                string str = string.Format("{0:0}", coverage1);

                sb.Append("" + str + "% ;; ");
            }
            sb.Append("\t");

            // number of PSMs for listed peptides
            sb.Append("" + AllPsmsForStrictPeptideSequences.Count());
            sb.Append("\t");

            // number of PSMs for listed peptides
            //sb.Append("" + TotalPsmList.Count());
            //sb.Append("\t");

            // summed metamorpheus score
            sb.Append(proteinGroupScore);
            sb.Append("\t");

            // isDecoy
            sb.Append(isDecoy);
            sb.Append("\t");

            // cumulative target
            sb.Append(cumulativeTarget);
            sb.Append("\t");

            // cumulative decoy
            sb.Append(cumulativeDecoy);
            sb.Append("\t");

            // q value
            sb.Append(QValue * 100);
            sb.Append("\t");

            return sb.ToString();
        }

        public void ScoreThisProteinGroup(List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            // find the best psm per base sequence (peptide FDR must be <1%) for scoring
            Dictionary<string, NewPsmWithFdr> peptideBaseSeqToBestPsmMatching = new Dictionary<string, NewPsmWithFdr>();
            foreach (var psm in TotalPsmList)
            {
                if (psm.qValue < 0.01)
                {
                    CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications);
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
                CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications);
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
                    CompactPeptide psmPeptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications);
                    string psmPeptideBaseSeq = string.Join("", psmPeptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    if (strictPepBaseSeq.Equals(psmPeptideBaseSeq))
                    {
                        AllPsmsForStrictPeptideSequences.Add(psm);
                    }
                }
            }

            // score the protein group
            foreach (var psm in BestPsmPerBaseSeq)
            {
                proteinGroupScore += psm.thisPSM.Score;
            }

            //if (StrictPeptideList.Count == 1)
            //    proteinGroupScore = 0;
        }

        public void CalculateSequenceCoverage()
        {
            sequenceCoverage = new List<double>();

            foreach (var protein in Proteins)
            {
                HashSet<int> coveredResidues = new HashSet<int>();

                foreach (var peptide in TotalPeptideWithSetModsList)
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
                sequenceCoverage.Add(sequenceCoverageHere);
            }
        }

        #endregion Public Methods

    }
}
