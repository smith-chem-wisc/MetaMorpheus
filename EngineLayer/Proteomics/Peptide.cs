using Proteomics;

namespace EngineLayer
{
    public abstract class Peptide
    {

        #region Private Fields

        private string baseSequence;
        private string baseLeucineSequence;

        #endregion Private Fields

        #region Protected Constructors

        protected Peptide(Protein protein, int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, bool addCompIons)
        {
            this.Protein = protein;
            this.OneBasedStartResidueInProtein = oneBasedStartResidueInProtein;
            this.OneBasedEndResidueInProtein = oneBasedEndResidueInProtein;
            this.addCompIons = addCompIons;
            Length = OneBasedEndResidueInProtein - OneBasedStartResidueInProtein + 1;
        }

        #endregion Protected Constructors

        #region Public Properties

        public Protein Protein { get; }
        public int OneBasedStartResidueInProtein { get; }
        public int OneBasedEndResidueInProtein { get; }
        public bool addCompIons { get; }

        public virtual string PeptideDescription { get; protected set; }

        public int Length { get; }

        public virtual char PreviousAminoAcid
        {
            get
            {
                return OneBasedStartResidueInProtein > 1 ? Protein[OneBasedStartResidueInProtein - 2] : '-';
            }
        }

        public virtual char NextAminoAcid
        {
            get
            {
                return OneBasedEndResidueInProtein < Protein.Length ? Protein[OneBasedEndResidueInProtein] : '-';
            }
        }

        public string BaseSequence
        {
            get
            {
                if (baseSequence == null)
                    baseSequence = Protein.BaseSequence.Substring(OneBasedStartResidueInProtein - 1, Length);
                return baseSequence;
            }
        }

        public string BaseLeucineSequence
        {
            get
            {
                if (baseLeucineSequence == null)
                    baseLeucineSequence = Protein.BaseSequence.Substring(OneBasedStartResidueInProtein - 1, Length).Replace('I', 'L');
                return baseLeucineSequence;
            }
        }

        #endregion Public Properties

        #region Public Indexers

        public char this[int zeroBasedIndex]
        {
            get
            {
                return Protein.BaseSequence[zeroBasedIndex + OneBasedStartResidueInProtein - 1];
            }
        }

        #endregion Public Indexers

    }
}