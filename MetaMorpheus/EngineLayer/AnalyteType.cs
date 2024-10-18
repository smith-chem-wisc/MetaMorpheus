using System.Collections.Generic;
using System.Dynamic;
using Microsoft.ML.Trainers;

namespace EngineLayer
{
    /// <summary>
    /// Represents an analyte type and is used to determine the output format of the analyte type.
    /// </summary>
    public class AnalyteType
    {
        public static Dictionary<string, AnalyteType> AnalyteTypes = new Dictionary<string, AnalyteType>
            {
                { "Peptide", new AnalyteType("PSM", "Peptide", "Protein", "psmtsv") },
                { "Proteoform", new AnalyteType("PSM", "Proteoform", "Protein", "psmtsv") },
                { "Oligo", new AnalyteType("OSM", "Oligo", "Transcript", "osmtsv") },
            };


        public AnalyteType(string spectralMatchLabel, string uniqueFormLabel, string bioPolymerLabel, string spectralMatchExtension)
        {
            SpectralMatchLabel = spectralMatchLabel;
            UniqueFormLabel = uniqueFormLabel;
            BioPolymerLabel = bioPolymerLabel;
            SpectralMatchExtension = spectralMatchExtension;
        }

        /// <summary>
        /// Gets or sets the label for spectral matches (e.g. PSM).
        /// </summary>
        public string SpectralMatchLabel { get; init; }

        /// <summary>
        /// Extension for spectral matches (e.g. .psmtsv).
        /// </summary>
        public string SpectralMatchExtension { get; init; }

        /// <summary>
        /// Gets or sets the label for unique forms (e.g. Peptide).
        /// </summary>
        public string UniqueFormLabel { get; init; }

        /// <summary>
        /// Gets or sets the label for grouped forms (e.g. Protein).
        /// </summary>
        public string BioPolymerLabel { get; init; }

        /// <summary>
        /// Returns the string representation of the AnalyteType.
        /// This is used to determine the output format of the AnalyteType.
        /// </summary>
        /// <returns>The unique form label.</returns>
        public override string ToString()
        {
            return UniqueFormLabel;
        }

        /// <summary>
        /// Gets the AnalyteType based on the specified analyte type label.
        /// </summary>
        /// <param name="analyteTypeLabel">The analyte type label.</param>
        /// <returns>The AnalyteType.</returns>
        public static AnalyteType GetAnalyteType(string analyteTypeLabel)
        {
            if (AnalyteTypes.TryGetValue(analyteTypeLabel, out AnalyteType type))
                return type;
            else
                throw new MetaMorpheusException("Unrecognized analyte type: " + analyteTypeLabel);
        }
    }
}
