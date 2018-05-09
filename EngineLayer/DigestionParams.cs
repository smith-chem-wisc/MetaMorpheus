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

        public DigestionParams(Protease Protease, int MaxMissedCleavages = 2, int MinPeptideLength = 7, int MaxPeptideLength=int.MaxValue, int MaxModificationIsoforms = 1024, InitiatorMethionineBehavior InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable, int MaxModsForPeptides = 2, bool SemiProteaseDigestion= false, TerminusType TerminusTypeSemiProtease = TerminusType.N)
        {
            this.MaxMissedCleavages = MaxMissedCleavages;
            this.MinPeptideLength = MinPeptideLength;
            this.MaxPeptideLength = MaxPeptideLength;
            this.MaxModificationIsoforms = MaxModificationIsoforms;
            this.Protease = Protease;//GlobalVariables.ProteaseDictionary["trypsin"];
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
            bool testMissed = MaxMissedCleavages.Equals(a.MaxMissedCleavages);
            bool testMin = MinPeptideLength.Equals(a.MinPeptideLength);
            bool testMax = MaxPeptideLength.Equals(a.MaxPeptideLength);
            bool testInit = InitiatorMethionineBehavior.Equals(a.InitiatorMethionineBehavior);
            bool testMaxModsIso = MaxModificationIsoforms.Equals(a.MaxModificationIsoforms);
            bool testMaxMods = MaxModsForPeptide.Equals(a.MaxModsForPeptide);
            bool testProtease = Protease.Equals(a.Protease);
            bool testSemiProtease = SemiProteaseDigestion.Equals(a.SemiProteaseDigestion);
            bool testTerminus = TerminusTypeSemiProtease.Equals(a.TerminusTypeSemiProtease);
            return a != null
                && MaxMissedCleavages.Equals(a.MaxMissedCleavages)
                && MinPeptideLength.Equals(a.MinPeptideLength)
                && MaxPeptideLength.Equals(a.MaxPeptideLength)
                && InitiatorMethionineBehavior.Equals(a.InitiatorMethionineBehavior)
                && MaxModificationIsoforms.Equals(a.MaxModificationIsoforms)
                && MaxModsForPeptide.Equals(a.MaxModsForPeptide)
                && Protease.Equals(a.Protease)
                && SemiProteaseDigestion.Equals(a.SemiProteaseDigestion)
                && TerminusTypeSemiProtease.Equals(a.TerminusTypeSemiProtease);
        }

        public override int GetHashCode()
        {
            return
                MaxMissedCleavages.GetHashCode()
                ^ InitiatorMethionineBehavior.GetHashCode()
                ^ MaxModificationIsoforms.GetHashCode()
                ^ MaxModsForPeptide.GetHashCode();
        }

        
        public void setMaxMissedCleavages( int MaxMissedCleavages)
        {
            this.MaxMissedCleavages = MaxMissedCleavages;
        }

        public void setInitiatorMethionineBehavior(InitiatorMethionineBehavior InitiatorMethionineBehavior)
        {
            this.InitiatorMethionineBehavior=InitiatorMethionineBehavior;
        }

        public void setMinPeptideLength( int MinPeptideLength)
        {
            this.MinPeptideLength=MinPeptideLength;
        }

        public void setMaxPeptideLength( int MaxPeptideLength)
        {
            this.MaxPeptideLength = MaxPeptideLength;
        }

        public void setMaxModificationIsoforms(int MaxModificationIsoforms)
        {
            this.MaxModificationIsoforms = MaxModificationIsoforms;
        }

        public void setMaxModsForPeptide( int MaxModsForPeptide)
        {
            this.MaxModsForPeptide = MaxModsForPeptide;
        }

        public void setProtease( Protease Protease)
        {
            this.Protease = Protease;
        }

        public void setSemiProteaseDigestion( bool SemiProteaseDigestion)
        {
            this.SemiProteaseDigestion = SemiProteaseDigestion;
        }

        public void setTerminusTypeSemiProtease(TerminusType TerminusTypeSemiProtease)
        {
            this.TerminusTypeSemiProtease = TerminusTypeSemiProtease;
        }
        #endregion Public Methods
    }
}