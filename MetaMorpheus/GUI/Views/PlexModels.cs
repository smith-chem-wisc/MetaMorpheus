using System;

namespace MetaMorpheusGUI
{
    // Per-channel annotation (reporter tag) for a plex.
    public class PlexAnnotation
    {
        public string Tag { get; set; } = "";
        public string SampleName { get; set; } = "";
        public string Condition { get; set; } = "";
        public int BiologicalReplicate { get; set; } = 1;

        /// <summary>
        /// What this channel is for, bound to the Sample Type column of the annotation grid.
        /// Held as the design file's own spelling so the grid, the file and
        /// EngineLayer.TmtSampleType cannot drift apart; TmtExperimentalDesign.TryParseSampleType
        /// is the single place that interprets it.
        /// </summary>
        public string SampleType { get; set; } = "study sample";
    }
}
