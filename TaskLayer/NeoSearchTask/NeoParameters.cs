namespace TaskLayer
{
    public class NeoParameters
    {
        public bool Calibrate { get; set; } = true;
        public double? PrecursorTolerancePPM { get; set; }
        public double? ProductTolerancePPM { get; set; }
        public bool GPTMD { get; set; } = true;
        public bool TargetSearch { get; set; } = true;
        public string TargetFilePath { get; set; }
        public bool DecoySearch { get; set; } = true;
        public string DecoyFilePath { get; set; }

        public bool SearchNTerminus { get; set; } = true;
        public string NFilePath { get; set; }
        public bool SearchCTerminus { get; set; } = true;
        public string CFilePath { get; set; }

        public int MaxMissedConsecutiveCleavages { get; set; } = 2;
        public int MaxMissedCleavages { get; set; } = 5;
        public int MaxCandidatesPerSpectrum { get; set; } = 2000;

        public int MinDistanceAllowed { get; set; } = 1;
        public int MaxDistanceAllowed { get; set; } = 25;
        public bool NormalCis { get; set; } = true;
        public bool ReverseCis { get; set; } = true;
    }
}