using System.Collections.Generic;

namespace EngineLayer
{
    public enum AnalyteType
    {
        Peptide,
        Proteoform,
        Oligo
    }

    /// <summary>
    /// Accessor methods for specific information about certain analyte types
    /// </summary>
    public static class AnalyteTypeExtensions
    {
        private static readonly Dictionary<AnalyteType, AnalyteTypeData> AnalyteTypes = new()
            {
                { AnalyteType.Peptide, new AnalyteTypeData("PSM", "Peptide", "Protein", "psmtsv", "Protease", "Precursors", "MS2 Scans") },
                { AnalyteType.Proteoform, new AnalyteTypeData("PSM", "Proteoform", "Protein", "psmtsv", "Protease", "Precursor", "MS2 Scans") },
                { AnalyteType.Oligo, new AnalyteTypeData("OSM", "Oligo", "Transcript", "osmtsv", "RNase", "Precursors", "MS2 Scans") },
            };

        public static string GetSpectralMatchLabel(this AnalyteType analyteType) => AnalyteTypes[analyteType].SpectralMatchLabel;
        public static string GetSpectralMatchExtension(this AnalyteType analyteType) => AnalyteTypes[analyteType].SpectralMatchExtension;
        public static string GetUniqueFormLabel(this AnalyteType analyteType) => AnalyteTypes[analyteType].UniqueFormLabel;
        public static string GetBioPolymerLabel(this AnalyteType analyteType) => AnalyteTypes[analyteType].BioPolymerLabel;
        public static string GetDigestionAgentLabel(this AnalyteType analyteType) => AnalyteTypes[analyteType].DigestionAgentLabel;
        public static string GetPrecursorLabel(this AnalyteType analyteType) => AnalyteTypes[analyteType].PrecursorLabel;
        public static string GetMs2ScanCountLabel(this AnalyteType analyteType) => AnalyteTypes[analyteType].Ms2ScanCountLabel;
    }

    /// <summary>
    /// Represents an analyte type and is used to determine the output format of the analyte type.
    /// </summary>
    internal class AnalyteTypeData(
        string spectralMatchLabel, 
        string uniqueFormLabel, 
        string bioPolymerLabel, 
        string spectralMatchExtension, 
        string digestionAgentLabel, 
        string precursorLabel, 
        string ms2ScanCountLabel)
    {
        /// <summary>
        /// Gets or sets the label for spectral matches (e.g. PSM).
        /// </summary>
        internal string SpectralMatchLabel { get; init; } = spectralMatchLabel;

        /// <summary>
        /// Extension for spectral matches (e.g. psmtsv).
        /// </summary>
        internal string SpectralMatchExtension { get; init; } = spectralMatchExtension;

        /// <summary>
        /// Gets the label for unique forms (e.g. Peptide).
        /// </summary>
        internal string UniqueFormLabel { get; init; } = uniqueFormLabel;

        /// <summary>
        /// Gets the label for grouped forms (e.g. Protein).
        /// </summary>
        internal string BioPolymerLabel { get; init; } = bioPolymerLabel;

        /// <summary>
        /// Gets the label for the digestion agent class (e.g. Protease)
        /// </summary>
        internal string DigestionAgentLabel { get; init; } = digestionAgentLabel;
        /// <summary>
        /// Gets the label for precursor (e.g. Precursor)
        /// 
        internal string PrecursorLabel { get; init; } = precursorLabel;
        /// <summary>
        /// Gets the label for MS2 scan count (e.g. MS2 Scan Count)
        /// 
        internal string Ms2ScanCountLabel { get; init; } = ms2ScanCountLabel;
    }
}

