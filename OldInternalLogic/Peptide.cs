namespace OldInternalLogic
{
    public abstract class Peptide
    {
        #region Protected Fields

        protected string baseSequence;
        protected string baseLeucineSequence;

        #endregion Protected Fields

        #region Public Constructors

        protected Peptide(Protein protein, int OneBasedStartResidueInProtein, int OneBasedEndResidueInProtein)
        {
            this.protein = protein;
            this.OneBasedStartResidueInProtein = OneBasedStartResidueInProtein;
            this.OneBasedEndResidueInProtein = OneBasedEndResidueInProtein;
        }

        #endregion Public Constructors

        #region Public Properties

        public Protein protein { get; private set; }
        public int OneBasedStartResidueInProtein { get; private set; }
        public int OneBasedEndResidueInProtein { get; private set; }

        public virtual string PeptideDescription { get; protected set; }

        public int Length
        {
            get
            {
                return OneBasedEndResidueInProtein - OneBasedStartResidueInProtein + 1;
            }
        }

        public virtual char PreviousAminoAcid
        {
            get
            {
                return OneBasedStartResidueInProtein > 1 ? protein[OneBasedStartResidueInProtein - 2] : '-';
            }
        }

        public virtual char NextAminoAcid
        {
            get
            {
                return OneBasedEndResidueInProtein < protein.Length ? protein[OneBasedEndResidueInProtein] : '-';
            }
        }

        public string BaseSequence
        {
            get
            {
                if (baseSequence == null)
                    baseSequence = protein.BaseSequence.Substring(OneBasedStartResidueInProtein - 1, Length);
                return baseSequence;
            }
        }

        public string BaseLeucineSequence
        {
            get
            {
                if (baseLeucineSequence == null)
                    baseLeucineSequence = protein.BaseSequence.Substring(OneBasedStartResidueInProtein - 1, Length).Replace('I', 'L');
                return baseLeucineSequence;
            }
        }

        #endregion Public Properties

        #region Public Indexers

        public char this[int zeroBasedIndex]
        {
            get
            {
                return protein.BaseSequence[zeroBasedIndex + OneBasedStartResidueInProtein - 1];
            }
        }

        #endregion Public Indexers
    }
}