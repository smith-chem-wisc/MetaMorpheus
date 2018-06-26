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

        // this parameterless constructor needs to exist to read the toml. 
        // if you can figure out a way to get rid of it, feel free...
        public DigestionParams() : this(GlobalVariables.ProteaseDictionary["trypsin"])
        {

        }

        public DigestionParams(Protease protease = null, int MaxMissedCleavages = 2, int MinPeptideLength = 7, int MaxPeptideLength=int.MaxValue, int MaxModificationIsoforms = 1024, 
            InitiatorMethionineBehavior InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable, int MaxModsForPeptides = 2, bool SemiProteaseDigestion= false, TerminusType TerminusTypeSemiProtease = TerminusType.N)
        {
            if (protease == null)
            {
                this.Protease = GlobalVariables.ProteaseDictionary["trypsin"];
            }
            else
            {
                this.Protease = protease;
            }
                
            this.MaxMissedCleavages = MaxMissedCleavages;
            this.MinPeptideLength = MinPeptideLength;
            this.MaxPeptideLength = MaxPeptideLength;
            this.MaxModificationIsoforms = MaxModificationIsoforms;
            this.InitiatorMethionineBehavior = InitiatorMethionineBehavior;
            this.MaxModsForPeptide = MaxModsForPeptides;
            this.SemiProteaseDigestion = SemiProteaseDigestion;
            this.TerminusTypeSemiProtease = TerminusTypeSemiProtease;
           
        }

        #endregion Public Constructors
        
        #region Public Methods

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

        #endregion Public Methods
    }
}