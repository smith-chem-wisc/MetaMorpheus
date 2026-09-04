namespace GuiFunctions.MetaDraw;

public enum BioPolymerCoverageType
{
    Unique,
    UniqueMissedCleavage,
    TandemRepeat,               // Shared within the same biopolymer
    TandemRepeatMissedCleavage,
    Shared,                     // Shared between different biopolymers
    SharedMissedCleavage,
}