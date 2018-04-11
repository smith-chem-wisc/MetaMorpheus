namespace EngineLayer
{
    public interface IDigestionParams
    {
        #region Public Properties

        int MaxMissedCleavages { get; }
        int MinPeptideLength { get; }
        int MaxPeptideLength { get; }
        InitiatorMethionineBehavior InitiatorMethionineBehavior { get; }
        int MaxModificationIsoforms { get; }
        int MaxModsForPeptide { get; }
        Protease Protease { get; }
        bool SemiProteaseDigestion { get; }
        TerminusType TerminusTypeSemiProtease { get; }

        #endregion Public Properties
    }
}