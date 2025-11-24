using System;

namespace MetaMorpheusGUI
{
    // Per-channel annotation (reporter tag) for a plex (Technical replicate removed).
    public class PlexAnnotation
    {
        public string Tag { get; set; } = "";
        public string SampleName { get; set; } = "";
        public string Condition { get; set; } = "";
        public int BiologicalReplicate { get; set; }
    }

    // Mapping of an input file to its plex, fraction index, and technical replicate (moved here).
    public class PlexFileEntry
    {
        public PlexFileEntry(string filePath, int fraction, string plex, int technicalReplicate)
        {
            FilePath = filePath;
            Fraction = fraction;
            Plex = plex;
            TechnicalReplicate = technicalReplicate;
        }

        public string FilePath { get; }
        public int Fraction { get; }
        public string Plex { get; }
        public int TechnicalReplicate { get; }
    }
}