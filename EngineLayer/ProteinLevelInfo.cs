using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class ProteinLinkedInfo
    {

        #region Public Constructors

        public ProteinLinkedInfo(HashSet<PeptideWithSetModifications> hashSet)
        {
            PeptidesWithSetModifications = hashSet;
            IsDecoy = PeptidesWithSetModifications.Any(bb => bb.Protein.IsDecoy);
            var representative = PeptidesWithSetModifications.First();

            PeptideMonoisotopicMass = representative.MonoisotopicMass;
            BaseSequence = representative.BaseSequence;
        }

        #endregion Public Constructors

        #region Public Properties

        public HashSet<PeptideWithSetModifications> PeptidesWithSetModifications { get; }
        public string BaseSequence { get; }
        public double PeptideMonoisotopicMass { get; }
        public bool IsDecoy { get; }

        #endregion Public Properties

        #region Internal Methods

        internal static string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();

            // Could have MANY options
            sb.Append("Same Sequence Ambiguity" + '\t');
            sb.Append("Protein Accession" + '\t');
            sb.Append("Protein Name" + '\t');
            sb.Append("Gene Name" + '\t');
            sb.Append("Peptide Description" + '\t');
            sb.Append("Start and End Residues In Protein" + '\t');
            sb.Append("Previous Amino Acid" + '\t');
            sb.Append("Next Amino Acid" + '\t');

            // Single info, common for all peptides/proteins
            sb.Append("Base Sequence" + '\t');
            sb.Append("Full Sequence" + '\t');
            sb.Append("Num Variable Mods" + '\t');
            sb.Append("Missed Cleavages" + '\t');
            sb.Append("Peptide Monoisotopic Mass" + '\t');
            sb.Append("Notch" + '\t');
            sb.Append("Decoy/Contaminant/Target");
            return sb.ToString();
        }

        #endregion Internal Methods

    }
}