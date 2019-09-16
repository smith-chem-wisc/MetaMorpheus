
namespace EngineLayer
{
    public class PercolatorFeatures
    {
        public float Intensity { get; set; }
        public float ScanPrecursorCharge { get; set; }
        public float DeltaScore { get; set; }
        public float Notch { get; set; }
        public float PsmCount { get; set; }
        public float ModsCount { get; set; }
        public float MissedCleavagesCount { get; set; }
        public float Ambiguity { get; set; }
        public float LongestFragmentIonSeries { get; set; }
        public float HydrophobicityZScore { get; set; }
        public bool IsVariantPeptide { get; set; }
    }
}
