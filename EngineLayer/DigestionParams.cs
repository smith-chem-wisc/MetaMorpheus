using EngineLayer;

namespace EngineLayer
{
    public class DigestionParams
    {
        private int MaxMissedCleavages { get; set; }
        private InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }
        private int MinPeptideLength { get; set; }
        private int MaxPeptideLength { get; set; }
        private int MaxModificationIsoforms { get; set; }
        private int MaxModsForPeptide { get; set; }
        private Protease Protease { get; set; }
        private bool SemiProteaseDigestion { get; set; } //for nonspecific searching of proteases
        private TerminusType TerminusTypeSemiProtease { get; set; }

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

        public int getMaxMissedCleavages()
        {
            return MaxMissedCleavages;
        }

        public InitiatorMethionineBehavior getInitiatorMethionineBehavior()
        {
            return InitiatorMethionineBehavior;
        }
        
        public int getMinPeptideLength()
        {
            return MinPeptideLength;
        }

        public int getMaxPeptideLength()
        {
            return MaxPeptideLength;
        }

        public int getMaxModificationIsoforms()
        {
            return MaxModificationIsoforms;
        }

        public int getMaxModsForPeptide()
        {
            return MaxModsForPeptide;
        }

        public Protease getProtease()
        {
            return Protease;
        }

        public bool getSemiProteaseDigestion()
        {
            return SemiProteaseDigestion;
        } 

        public TerminusType getTerminusTypeSemiProtease()
        {
            return TerminusTypeSemiProtease;
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