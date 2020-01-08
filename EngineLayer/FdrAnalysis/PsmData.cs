using Microsoft.ML.Data;

namespace EngineLayer.FdrAnalysis
{
    public class PsmData
    {
        [LoadColumn(0)]
        public float Intensity { get; set; }

        [LoadColumn(1)]
        public float ScanPrecursorCharge { get; set; }

        [LoadColumn(2)]
        public float DeltaScore { get; set; }

        [LoadColumn(3)]
        public float Notch { get; set; }

        [LoadColumn(4)]
        public float PsmCount { get; set; }

        [LoadColumn(5)]
        public float ModsCount { get; set; }

        [LoadColumn(6)]
        public float MissedCleavagesCount { get; set; }

        [LoadColumn(7)]
        public float Ambiguity { get; set; }

        [LoadColumn(8)]
        public float LongestFragmentIonSeries { get; set; }

        [LoadColumn(9)]
        public float HydrophobicityZScore { get; set; }

        [LoadColumn(10)]
        public bool Label { get; set; }
    }
}