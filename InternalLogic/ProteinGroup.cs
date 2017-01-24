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
        public double proteinGroupScore;

        #endregion Public Fields

        #region Internal Constructors

        internal ProteinGroup(HashSet<Protein> proteins, HashSet<CompactPeptide> peptides, HashSet<CompactPeptide> uniquePeptides)
        {
            Proteins = proteins;
            PeptideList = peptides;
            UniquePeptideList = uniquePeptides;
            BestPsmList = null;
            TotalPsmList = new List<NewPsmWithFdr>();
            RazorPeptideList = new HashSet<CompactPeptide>();
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
        public List<NewPsmWithFdr> BestPsmList { get; set; }
        public List<NewPsmWithFdr> TotalPsmList { get; set; }
        public HashSet<CompactPeptide> PeptideList { get; set; }
        public HashSet<CompactPeptide> UniquePeptideList { get; set; }
        public HashSet<CompactPeptide> RazorPeptideList { get; set; }
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
            foreach (CompactPeptide uniquePeptide in UniquePeptideList)
            {
                string peptideBaseSequence = string.Join("", uniquePeptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                sb.Append("" + peptideBaseSequence + " ;; ");
            }
            sb.Append("\t");

            // list of shared peptides
            foreach (CompactPeptide sharedPeptide in PeptideList)
            {
                if (!UniquePeptideList.Contains(sharedPeptide))
                {
                    string peptideBaseSequence = string.Join("", sharedPeptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    sb.Append("" + peptideBaseSequence + " ;; ");
                }
            }
            sb.Append("\t");

            // list of razor peptides
            foreach (CompactPeptide razorPeptide in RazorPeptideList)
            {
                string peptideBaseSequence = string.Join("", razorPeptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                sb.Append("" + peptideBaseSequence + " ;; ");
            }
            sb.Append("\t");

            // number of peptides
            sb.Append("" + PeptideList.Count());
            sb.Append("\t");

            // number of unique peptides
            sb.Append("" + UniquePeptideList.Count());
            sb.Append("\t");

            // number of PSMs for final base sequences
            sb.Append("" + TotalPsmList.Count());
            sb.Append("\t");

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

        public void scoreThisProteinGroup()
        {
            // score the protein group
            foreach (var psm in BestPsmList)
            {
                proteinGroupScore += psm.thisPSM.Score;
            }
        }

        #endregion Public Methods

    }
}
