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

        // this parameterless constructor needs to exist to read the toml.
        // if you can figure out a way to get rid of it, feel free...
        public DigestionParams() : this("trypsin")
        {
        }

        public DigestionParams(string protease = "trypsin", int maxMissedCleavages = 2, int minPeptideLength = 7, int maxPeptideLength = int.MaxValue, int maxModificationIsoforms = 1024,
            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable, int maxModsForPeptides = 2, bool semiProteaseDigestion = false, TerminusType terminusTypeSemiProtease = TerminusType.N)
        {
            Protease = GlobalVariables.ProteaseDictionary[protease];
            MaxMissedCleavages = maxMissedCleavages;
            MinPeptideLength = minPeptideLength;
            MaxPeptideLength = maxPeptideLength;
            MaxModificationIsoforms = maxModificationIsoforms;
            InitiatorMethionineBehavior = initiatorMethionineBehavior;
            MaxModsForPeptide = maxModsForPeptides;
            SemiProteaseDigestion = semiProteaseDigestion;
            TerminusTypeSemiProtease = terminusTypeSemiProtease;
        }

        public override bool Equals(object obj)
        {
            DigestionParams a = obj as DigestionParams;

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
    }
}