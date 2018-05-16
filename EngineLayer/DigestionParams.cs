using EngineLayer;

namespace EngineLayer
{
    public class DigestionParams 
    {
        public int MaxMissedCleavages { get; private set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; private set; }
        public int MinPeptideLength { get; private set; }
        public int MaxPeptideLength { get; private set; }
        public int MaxModificationIsoforms { get; private set; }
        public int MaxModsForPeptide { get; private set; }
        public Protease Protease { get; private set; }
        public bool SemiProteaseDigestion { get; private set; } //for nonspecific searching of proteases
        public TerminusType TerminusTypeSemiProtease { get; private set; }

        #region Public Constructors

        public DigestionParams(Protease Protease, int MaxMissedCleavages = 2, int MinPeptideLength = 7, int MaxPeptideLength=int.MaxValue, int MaxModificationIsoforms = 1024, 
            InitiatorMethionineBehavior InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable, int MaxModsForPeptides = 2, bool SemiProteaseDigestion= false, TerminusType TerminusTypeSemiProtease = TerminusType.N)
        {
            this.MaxMissedCleavages = MaxMissedCleavages;
            this.MinPeptideLength = MinPeptideLength;
            this.MaxPeptideLength = MaxPeptideLength;
            this.MaxModificationIsoforms = MaxModificationIsoforms;
            this.Protease = Protease;
            this.InitiatorMethionineBehavior = InitiatorMethionineBehavior;
            this.MaxModsForPeptide = MaxModsForPeptides;
            this.SemiProteaseDigestion = SemiProteaseDigestion;
            this.TerminusTypeSemiProtease = TerminusTypeSemiProtease;
        }

        #endregion Public Constructors

        

        

        

        #region Public Methods

        public DigestionParams Clone()
        {
            return (DigestionParams)this.MemberwiseClone();
        }
        public override bool Equals(object obj)
        {
            DigestionParams a = obj as DigestionParams;
            
            return a != null
                && this.MaxMissedCleavages.Equals(a.MaxMissedCleavages)
                && this.MinPeptideLength.Equals(a.MinPeptideLength)
                && this.MaxPeptideLength.Equals(a.MaxPeptideLength)
                && this.InitiatorMethionineBehavior.Equals(a.InitiatorMethionineBehavior)
                && this.MaxModificationIsoforms.Equals(a.MaxModificationIsoforms)
                && this.MaxModsForPeptide.Equals(a.MaxModsForPeptide)
                && this.Protease.Equals(a.Protease)
                && this.SemiProteaseDigestion.Equals(a.SemiProteaseDigestion)
                && this.TerminusTypeSemiProtease.Equals(a.TerminusTypeSemiProtease);
        }

        public override int GetHashCode()
        {
            return
                MaxMissedCleavages.GetHashCode()
                ^ InitiatorMethionineBehavior.GetHashCode()
                ^ MaxModificationIsoforms.GetHashCode()
                ^ MaxModsForPeptide.GetHashCode();
        }

        
        public void SetMaxMissedCleavages( int MaxMissedCleavages)
        {
            this.MaxMissedCleavages = MaxMissedCleavages;
        }

        public void SetInitiatorMethionineBehavior(InitiatorMethionineBehavior InitiatorMethionineBehavior)
        {
            this.InitiatorMethionineBehavior=InitiatorMethionineBehavior;
        }

        public void SetMinPeptideLength( int MinPeptideLength)
        {
            this.MinPeptideLength=MinPeptideLength;
        }

        public void SetMaxPeptideLength( int MaxPeptideLength)
        {
            this.MaxPeptideLength = MaxPeptideLength;
        }

        public void SetMaxModificationIsoforms(int MaxModificationIsoforms)
        {
            this.MaxModificationIsoforms = MaxModificationIsoforms;
        }

        public void SetMaxModsForPeptide( int MaxModsForPeptide)
        {
            this.MaxModsForPeptide = MaxModsForPeptide;
        }

        public void SetProtease( Protease Protease)
        {
            this.Protease = Protease;
        }

        public void SetSemiProteaseDigestion( bool SemiProteaseDigestion)
        {
            this.SemiProteaseDigestion = SemiProteaseDigestion;
        }

        public void SetTerminusTypeSemiProtease(TerminusType TerminusTypeSemiProtease)
        {
            this.TerminusTypeSemiProtease = TerminusTypeSemiProtease;
        }
        #endregion Public Methods
    }
}