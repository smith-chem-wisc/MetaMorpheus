using EngineLayer;

namespace EngineLayer
{
    public class DigestionParams
    {
        public readonly int MaxMissedCleavages;
        public readonly InitiatorMethionineBehavior InitiatorMethionineBehavior;
        public readonly int MinPeptideLength;
        public readonly int MaxPeptideLength;
        public readonly int MaxModificationIsoforms;
        public readonly int MaxModsForPeptide;
        public readonly Protease Protease;
        public readonly bool SemiProteaseDigestion; //for nonspecific searching of proteases
        public readonly TerminusType TerminusTypeSemiProtease;

        #region Public Constructors

        public DigestionParams(string protease = "trypsin", int MaxMissedCleavages = 2, int MinPeptideLength = 7, int MaxPeptideLength = int.MaxValue, int MaxModificationIsoforms = 1024,
            InitiatorMethionineBehavior InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable, int MaxModsForPeptides = 2, bool SemiProteaseDigestion = false, TerminusType TerminusTypeSemiProtease = TerminusType.N)
        {
            this.Protease = GlobalVariables.ProteaseDictionary[protease];
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