using System.Collections.Generic;

namespace EngineLayer
{
    public class AnalyteType
    {
        public static Dictionary<string, AnalyteType> AnalyteTypes = new Dictionary<string, AnalyteType>
        {
            { "Peptide", new AnalyteType("PSM", "Peptide", "ProteinGroup") },
            { "Proteoform", new AnalyteType("PSM", "Proteoform", "ProteinGroup") },
            { "Oligo", new AnalyteType("OSM", "Oligo", "TranscriptGroup") },
        };

        public AnalyteType(string spectralMatchLabel, string uniqueFormLabel, string groupedFormLabel)
        {
            SpectralMatchLabel = spectralMatchLabel;
            UniqueFormLabel = uniqueFormLabel;
            GroupedFormLabel = groupedFormLabel;
        }

        /// <summary>
        /// Label for spectral matches (e.g. PSM)
        /// </summary>
        public string SpectralMatchLabel { get; init; }

        /// <summary>
        /// Label for unique forms (e.g. Peptide)
        /// </summary>
        public string UniqueFormLabel { get; init; }

        /// <summary>
        /// Label for grouped forms (e.g. ProteinGroup)
        /// </summary>
        public string GroupedFormLabel { get; init; }

        public override string ToString()
        {
            return UniqueFormLabel;
        }

        public static AnalyteType GetAnalyteType(string analyteTypeLabel)
        {
            if (AnalyteTypes.TryGetValue(analyteTypeLabel, out AnalyteType type))
                return type;
            else
                throw new MetaMorpheusException("Unrecognized analyte type: " + analyteTypeLabel);
        }
    }
}
