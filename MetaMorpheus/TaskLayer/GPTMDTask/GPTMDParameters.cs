using System.Collections.Generic;
using System.Linq;
using EngineLayer.Gptmd;
using Nett;

namespace EngineLayer
{
    public class GptmdParameters
    {
        /// <summary>
        /// The default modifications for a G-PTM-D search can be altered here by adding to b.Modification.Type.Equals("  NAME OF MODS TO ADD TO DEFAULT ")
        /// </summary>
        public GptmdParameters()
        {
            GptmdFilters = new();
            GptmdFilterTypes = new();
            ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                b.ModificationType.Equals("Common Artifact")
                || b.ModificationType.Equals("Common Biological")
                || b.ModificationType.Equals("Metal")
                //|| b.ModificationType.Equals("Less Common")
            ).Select(b => (b.ModificationType, b.IdWithMotif)).ToList();
        }

        public List<(string, string)> ListOfModsGptmd { get; set; }

        /// <summary>
        /// Runtime list of G-PTM-D acceptance filters (e.g. <see cref="ImprovedScoreFilter"/>). This is a
        /// list of interface instances and cannot round-trip through the task toml; the toml-serializable
        /// representation is <see cref="GptmdFilterTypes"/>. The GUI sets this directly for in-memory runs.
        /// Use <see cref="GetActiveFilters"/> to obtain the effective set (this list plus any named in the toml).
        /// </summary>
        [TomlIgnore]
        public List<IGptmdFilter> GptmdFilters { get; set; }

        /// <summary>
        /// Toml-serializable filter selection by type name (e.g. "ImprovedScoreFilter"). Lets headless/CMD
        /// runs enable G-PTM-D acceptance filters that previously could only be set in the GUI. A common
        /// choice is "ImprovedScoreFilter", which only adds a modification when it improves the match score
        /// (greatly reducing the number of modifications added).
        /// </summary>
        public List<string> GptmdFilterTypes { get; set; }

        public bool WriteDecoys { get; set; } = true;

        /// <summary>
        /// The effective filter set: the in-memory <see cref="GptmdFilters"/> plus any filters named in
        /// <see cref="GptmdFilterTypes"/> (deduplicated by type). Unknown names are ignored.
        /// </summary>
        public List<IGptmdFilter> GetActiveFilters()
        {
            var active = new List<IGptmdFilter>(GptmdFilters ?? new List<IGptmdFilter>());
            foreach (var name in GptmdFilterTypes ?? new List<string>())
            {
                var filter = CreateFilter(name);
                if (filter != null && active.All(f => f.GetType() != filter.GetType()))
                    active.Add(filter);
            }
            return active;
        }

        /// <summary>Creates a filter instance from its type name; returns null for an unrecognized name.</summary>
        public static IGptmdFilter CreateFilter(string typeName) => typeName switch
        {
            nameof(ImprovedScoreFilter) => new ImprovedScoreFilter(),
            nameof(DualDirectionalIonCoverageFilter) => new DualDirectionalIonCoverageFilter(),
            nameof(UniDirectionalIonCoverageFilter) => new UniDirectionalIonCoverageFilter(),
            nameof(FlankingIonCoverageFilter) => new FlankingIonCoverageFilter(),
            _ => null
        };
    }
}
