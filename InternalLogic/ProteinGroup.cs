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
        public readonly double proteinGroupScore;

        #endregion Public Fields

        #region Private Fields

        private readonly double summedIntensity;
        private readonly double summedUniquePeptideIntensity;

        #endregion Private Fields

        #region Public Constructors

        public ProteinGroup(HashSet<Protein> proteins, List<NewPsmWithFDR> psmList, HashSet<CompactPeptide> allUniquePeptides, List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            this.proteins = proteins;
            this.psmList = psmList;
            peptideList = new List<CompactPeptide>();
            uniquePeptideList = new List<CompactPeptide>();
            proteinGroupScore = 0;

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.isDecoy)
                    isDecoy = true;
            }

            // build list of compact peptides associated with the protein group
            // all peptides in the group are associated with all proteins in the group
            // if encountering duplicate peptides, only use the best-scoring one
            foreach (var psm in psmList)
            {
                CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications);
                peptideList.Add(peptide);

                // calculate the protein group score
                if(allUniquePeptides.Contains(peptide))
                {
                    uniquePeptideList.Add(peptide);
                    proteinGroupScore += psm.thisPSM.Score;
                }
            }

            // calculate summed psm intensity
            summedIntensity = 0;

            // calculate intensity of only unique peptides
            summedUniquePeptideIntensity = 0;

            // q value
            QValue = 0;
        }

        #endregion Public Constructors

        #region Public Properties

        public HashSet<Protein> proteins { get; private set; }
        public List<NewPsmWithFDR> psmList { get; private set; }
        public List<CompactPeptide> peptideList { get; private set; }
        public List<CompactPeptide> uniquePeptideList { get; private set; }
        public double QValue { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("Protein Description" + '\t');
            sb.Append("Protein Sequence" + '\t');
            sb.Append("Protein Length" + '\t');
            sb.Append("Number of Proteins in Group" + '\t');
            sb.Append("Number of Peptide-Spectrum Matches" + '\t');
            sb.Append("Number of Unique Peptides" + '\t');
            sb.Append("Summed Peptide-Spectrum Match Precursor Intensity" + '\t');
            sb.Append("Summed Unique Peptide Precursor Intensity" + '\t');
            sb.Append("Summed Morpheus Score" + '\t');
            sb.Append("Decoy?" + '\t');
            sb.Append("Cumulative Target" + '\t');
            sb.Append("Cumulative Decoy" + '\t');
            sb.Append("Q-Value (%)");
            return sb.ToString();
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            // proteins in protein group
            foreach (Protein protein in proteins)
                sb.Append("" + protein.FullDescription + " ;; ");
            sb.Append("\t");

            // sequences of proteins in group
            foreach (Protein protein in proteins)
                sb.Append("" + protein.BaseSequence + " ;; ");
            sb.Append("\t");

            // length of each protein
            foreach (Protein protein in proteins)
                sb.Append("" + protein.BaseSequence.Length + " ;; ");
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + proteins.Count);
            sb.Append("\t");

            // number of psm's for the group
            sb.Append("" + psmList.Count);
            sb.Append("\t");

            // number of unique peptides
            sb.Append("" + uniquePeptideList.Count());
            sb.Append("\t");

            // summed psm precursor intensity
            sb.Append(summedIntensity);
            sb.Append("\t");

            // summed unique peptide precursor intensity
            sb.Append(summedUniquePeptideIntensity);
            sb.Append("\t");

            // summed metamorpheus score
            sb.Append(proteinGroupScore);
            sb.Append("\t");

            // isdecoy
            sb.Append(isDecoy);
            sb.Append("\t");

            return sb.ToString();
        }

        #endregion Public Methods
    }
}