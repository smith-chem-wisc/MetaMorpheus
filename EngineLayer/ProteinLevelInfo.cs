using System.Collections.Generic;
using System.Globalization;
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
            IsContaminant = PeptidesWithSetModifications.Any(bb => bb.Protein.IsContaminant);
            var representative = PeptidesWithSetModifications.First();

            PeptideMonoisotopicMass = representative.MonoisotopicMass;
            FullSequence = representative.Sequence;
            BaseSequence = representative.BaseSequence;
            MissedCleavages = representative.MissedCleavages;
            NumVariableMods = representative.NumMods - representative.numFixedMods;
            SequenceWithChemicalFormulas = representative.SequenceWithChemicalFormulas;
        }

        #endregion Public Constructors

        #region Public Properties

        public HashSet<PeptideWithSetModifications> PeptidesWithSetModifications { get; }
        public string FullSequence { get; }
        public string BaseSequence { get; }
        public int MissedCleavages { get; }
        public double PeptideMonoisotopicMass { get; }
        public int NumVariableMods { get; }
        public string SequenceWithChemicalFormulas { get; }
        public bool IsContaminant { get; }
        public bool IsDecoy { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            var s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.Protein.Accession));
            if (s.Length > 32000)
                s = "too many";
            sb.Append(s + "\t");

            s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.Protein.FullName));
            if (s.Length > 32000)
                s = "too many";
            sb.Append(s + "\t");

            s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.PeptideDescription));
            if (s.Length > 32000)
                s = "too many";
            sb.Append(s + "\t");

            s = string.Join(" or ", PeptidesWithSetModifications.Select(b => "[" + b.OneBasedStartResidueInProtein + " to " + b.OneBasedEndResidueInProtein + "]"));
            if (s.Length > 32000)
                s = "too many";
            sb.Append(s + "\t");

            s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.PreviousAminoAcid));
            if (s.Length > 32000)
                s = "too many";
            sb.Append(s + "\t");

            s = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.NextAminoAcid));
            if (s.Length > 32000)
                s = "too many";
            sb.Append(s + "\t");

            var representative = PeptidesWithSetModifications.First();

            sb.Append(representative.BaseSequence + "\t");
            sb.Append(representative.Sequence + "\t");
            sb.Append(NumVariableMods.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(MissedCleavages.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(PeptideMonoisotopicMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');

            if (IsDecoy)
                sb.Append("D");
            else if (IsContaminant)
                sb.Append("C");
            else
                sb.Append("T");

            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal static string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();

            // Could have MANY options
            sb.Append("Protein Accession" + '\t');
            sb.Append("Protein Name" + '\t');
            sb.Append("Peptide Description" + '\t');
            sb.Append("Start and End Residues In Protein" + '\t');
            sb.Append("Previous Amino Acid" + '\t');
            sb.Append("Next Amino Acid" + '\t');

            // Single info, common for all peptides/proteins
            sb.Append("Base Sequence" + '\t');
            sb.Append("Full Sequence" + '\t');
            sb.Append("Variable Mods" + '\t');
            sb.Append("Missed Cleavages" + '\t');
            sb.Append("Peptide Monoisotopic Mass" + '\t');
            sb.Append("Decoy/Contaminant/Target");
            return sb.ToString();
        }

        #endregion Internal Methods

    }
}