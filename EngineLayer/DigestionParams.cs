namespace EngineLayer
{
    public class DigestionParams
    {
        #region Public Constructors

        public DigestionParams()
        {
            MaxMissedCleavages = 2;
            MinPeptideLength = 5;
            MaxPeptideLength = 50;
            MaxModificationIsoforms = 4096;
            Protease = GlobalEngineLevelSettings.ProteaseDictionary["trypsin"];
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            MaxModsForPeptide = 3;
            SemiProteaseDigestion = false;
        }

        #endregion Public Constructors

        #region Public Properties

        public int MaxMissedCleavages { get; set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }
        public int? MinPeptideLength { get; set; }
        public int? MaxPeptideLength { get; set; }
        public int MaxModificationIsoforms { get; set; }
        public int MaxModsForPeptide { get; set; }
        public Protease Protease { get; set; }
        public bool SemiProteaseDigestion { get; set; } //for nonspecific searching of proteases
        public TerminusType TerminusTypeSemiProtease { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override bool Equals(object obj)
        {
            var a = obj as DigestionParams;
            return a != null
                && MaxMissedCleavages.Equals(a.MaxMissedCleavages)
                && MinPeptideLength.Equals(a.MinPeptideLength)
                && MaxPeptideLength.Equals(a.MaxPeptideLength)
                && InitiatorMethionineBehavior.Equals(a.InitiatorMethionineBehavior)
                && MaxModificationIsoforms.Equals(a.MaxModificationIsoforms)
                && MaxModsForPeptide.Equals(a.MaxModsForPeptide)
                && Protease.Equals(a.Protease);
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