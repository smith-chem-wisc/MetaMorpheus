using System;
using System.Collections.Generic;

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

    // Mapping of an input file to its plex, fraction index, and technical replicate, plus the plex annotations.
    public class PlexFileEntry
    {
        public PlexFileEntry(string filePath, int fraction, string plex, int technicalReplicate, IReadOnlyList<PlexAnnotation>? annotations = null)
        {
            FilePath = filePath;
            Fraction = fraction;
            Plex = plex;
            TechnicalReplicate = technicalReplicate;
            Annotations = annotations;
        }

        public string FilePath { get; }
        public int Fraction { get; }
        public string Plex { get; }
        public int TechnicalReplicate { get; }
        // Optional: the reporter-channel annotations associated with this plex
        public IReadOnlyList<PlexAnnotation>? Annotations { get; }
    }
}