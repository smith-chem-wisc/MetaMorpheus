namespace MetaMorpheus
{
    public abstract class Peptide
    {
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

        public Peptide(Protein protein, int OneBasedStartResidueInProtein, int OneBasedEndResidueInProtein)
        {
            this.protein = protein;
            this.OneBasedStartResidueInProtein = OneBasedStartResidueInProtein;
            this.OneBasedEndResidueInProtein = OneBasedEndResidueInProtein;
        }

        public char this[int zeroBasedIndex]
        {
            get
            {
                return protein.BaseSequence[zeroBasedIndex + OneBasedStartResidueInProtein - 1];
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

        protected string baseSequence;

        public string BaseSequence
        {
            get
            {
                if (baseSequence == null)
                    baseSequence = protein.BaseSequence.Substring(OneBasedStartResidueInProtein - 1, Length);
                return baseSequence;
            }
        }

        protected string baseLeucineSequence;

        public string BaseLeucineSequence
        {
            get
            {
                if (baseLeucineSequence == null)
                    baseLeucineSequence = protein.BaseSequence.Substring(OneBasedStartResidueInProtein - 1, Length).Replace('I', 'L');
                return baseLeucineSequence;
            }
        }
    }
}