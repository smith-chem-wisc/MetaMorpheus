using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using Omics.Modifications;
using MassSpectrometry;
using MzLibUtil;

namespace EngineLayer
{
    public enum IsobaricMassTagType
    {
        TMT6,
        TMT10,
        TMT11,
        TMT18,
        iTRAQ4,
        iTRAQ8,
        diLeu4,
        diLeu12
    }

    /// <summary>
    /// This class contains information about the Isobaric Mass Tag (e.g., TMT), including the theoretical m/z values of the reporter ions
    /// as well as methods designed to retrieve the intensities of those reporter ions from a given MzSpectrum. 
    /// It does not store any intensity information itself. Intensity information is associated with each Ms2ScanWithSpecificMass object or SpectralMatch,
    /// in the IsobaricMassTagReporterIonIntensities property.
    /// </summary>
    public class IsobaricMassTag
    {
        /// <summary>
        /// Theoretical m/z values corresponding to each reporter ion for this isobaric mass tag type.
        /// </summary>
        public double[] ReporterIonMzs { get; private set; }
        public IsobaricMassTagType TagType { get; private set; }
        /// <summary>
        /// Hardcoded tolerance value in Daltons for isobaric mass tag reporter ions.
        /// Taken from https://doi.org/10.1021/acs.jproteome.1c00168
        /// </summary>
        public const double AbsoluteToleranceValue = 0.003; // 3 millidaltons. This is
        public AbsoluteTolerance Tolerance { get; private set; }

        /// <summary>
        /// Returns an IsobaricMassTag that corresponds to the given modification ID.
        /// ReporterIonMzs are extracted from GlobalVariables.AllModsKnown.
        /// </summary>
        /// <param name="modificationId">The type of isobaric mass tag, can be take from SearchParameters.MultiplexModId </param>
        public static IsobaricMassTag GetIsobaricMassTag(string modificationId) =>
            GetIsobaricMassTag(GetTagTypeFromModificationId(modificationId));

        /// <summary>
        /// Returns an IsobaricMassTag that corresponds to the given TagType.
        /// ReporterIonMzs are extracted from GlobalVariables.AllModsKnown.
        /// </summary>
        /// <param name="tagType">The type of isobaric mass tag (e.g., TMT10, TMT11, iTRAQ4, etc.)</param>
        public static IsobaricMassTag GetIsobaricMassTag(IsobaricMassTagType? tagType)
        {
            if (tagType == null)
            {
                return null;
            }

            // Map the enum to the modification name pattern
            string modificationPattern = GetModificationPattern((IsobaricMassTagType)tagType);

            // Find the corresponding modification in GlobalVariables.AllModsKnown
            Modification isobaricModification = GlobalVariables.AllModsKnown
                .FirstOrDefault(m => m.ModificationType == "Multiplex Label"
                                     && m.IdWithMotif.Contains(modificationPattern));

            // Verify that the modification exists and has diagnostic ions
            if (isobaricModification == null
                || !isobaricModification.DiagnosticIons.TryGetValue(DissociationType.HCD, out var reporterIonMasses) 
                || reporterIonMasses.IsNullOrEmpty()) 
            {
                return null;
            }

            // Extract the reporter ion m/z values from the first dissociation type's diagnostic ions
            // and convert them to m/z values (charge state 1)
            double[] reporterIonMzs = reporterIonMasses
                .Select(mass => mass.ToMz(1))
                .OrderBy(mz => mz)
                .ToArray();

            return new IsobaricMassTag
            {
                TagType = (IsobaricMassTagType)tagType,
                ReporterIonMzs = reporterIonMzs,
                Tolerance = new AbsoluteTolerance(AbsoluteToleranceValue)
            };
        }

        /// <summary>
        /// Takes in a spectrum, searches for the reporter ions within the specified ppm tolerance,
        /// and returns their intensities. If a reporter ion is not found, its intensity will be zero.
        /// </summary>
        public double[]? GetReporterIonIntensities(MzSpectrum spectrum)
        {
            if (spectrum.Size < 1) return null;
            double[] reporterIonIntensities = new double[ReporterIonMzs.Length];
            if (spectrum.XArray[0] > Tolerance.GetMaximumValue(ReporterIonMzs[^1]))
            {
                return reporterIonIntensities; // If the spectrum starts after the last reporter ion, return all zeros.
            }

            // Two-pointer search. The for loop iterates through the theoretical ions mzs (theoreticalIonIndex),
            // and the while loop advances through the spectrum peaks (spectrumIndex) until we find a match or surpass the theoretical mz.
            // If there are multiple peaks that match the theoretical mz within tolerance, we take the one with the highest intensity.
            int spectrumIndex = 0;
            for (int theoreticalIonIndex = 0; theoreticalIonIndex < reporterIonIntensities.Length; theoreticalIonIndex++)
            {
                double minMz = Tolerance.GetMinimumValue(ReporterIonMzs[theoreticalIonIndex]);
                double maxMz = Tolerance.GetMaximumValue(ReporterIonMzs[theoreticalIonIndex]);
                double maxIntensity = 0;
                while (spectrumIndex < spectrum.Size 
                    && spectrum.XArray[spectrumIndex] < minMz)
                {
                    spectrumIndex++;
                }
                if (spectrumIndex >= spectrum.Size)
                {
                    break;
                }
                while (spectrumIndex < spectrum.Size 
                    && minMz <= spectrum.XArray[spectrumIndex] 
                    && spectrum.XArray[spectrumIndex] <= maxMz)
                {
                    if (spectrum.YArray[spectrumIndex] > maxIntensity)
                    {
                        maxIntensity = spectrum.YArray[spectrumIndex];
                    }
                    spectrumIndex++;
                }
                reporterIonIntensities[theoreticalIonIndex] = maxIntensity;
            }

            return reporterIonIntensities;
        }

        /// <summary>
        /// Maps the IsobaricMassTagType enum to the corresponding modification name pattern
        /// </summary>
        private static string GetModificationPattern(IsobaricMassTagType tagType)
        {
            return tagType switch
            {
                IsobaricMassTagType.TMT6 => "TMT6", // Note, the full ID for TMT6 is "TMT6-plex". 
                IsobaricMassTagType.TMT10 => "TMT10",
                IsobaricMassTagType.TMT11 => "TMT11",
                IsobaricMassTagType.TMT18 => "TMT18",
                IsobaricMassTagType.iTRAQ4 => "iTRAQ-4plex",
                IsobaricMassTagType.iTRAQ8 => "iTRAQ-8plex",
                IsobaricMassTagType.diLeu4 => "DiLeu-4plex",
                IsobaricMassTagType.diLeu12 => "DiLeu-12plex",
                _ => throw new ArgumentException($"Unknown isobaric mass tag type: {tagType}")
            };
        }

        /// <summary>
        /// Returns the labels corresponding to each reporter ion for the given modificationId.
        /// </summary>
        /// <param name="modificationId"></param>
        public static List<string> GetReporterIonLabels(string modificationId)
        {
            IsobaricMassTagType? tagType = GetTagTypeFromModificationId(modificationId);
            if (tagType == null)
            {
                return null;
            }
            return GetReporterIonLabels(tagType.Value);
        }

        /// <summary>
        /// Returns the labels corresponding to each reporter ion for the given isobaric mass tag type.
        /// </summary>
        /// <param name="tagType"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static List<string> GetReporterIonLabels(IsobaricMassTagType tagType)
        {
            return tagType switch
            {
                IsobaricMassTagType.TMT6 => new List<string> { "126", "127", "128", "129", "130", "131" },
                IsobaricMassTagType.TMT10 => new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N" },
                IsobaricMassTagType.TMT11 => new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" },
                IsobaricMassTagType.TMT18 => new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" },
                IsobaricMassTagType.iTRAQ4 => new List<string> { "114", "115", "116", "117" },
                IsobaricMassTagType.iTRAQ8 => new List<string> { "113", "114", "115", "116", "117", "118", "119", "120" },
                IsobaricMassTagType.diLeu4 => new List<string> { "115", "116", "117", "118" },
                IsobaricMassTagType.diLeu12 => new List<string> { "115a", "115b", "116a", "116b", "116c", "117a", "117b", "117c", "118a", "118b", "118c", "118d" },
                _ => null
            };
        }

        /// <summary>
        /// Converts a modification ID (IdWithMotif or OriginalId) to the corresponding IsobaricMassTagType enum.
        /// </summary>
        /// <param name="modificationId">The modification ID string (e.g., "TMT10 on K", "iTRAQ-4plex on K", "DiLeu-12plex on X")</param>
        /// <returns>The corresponding IsobaricMassTagType enum value, or null if no match is found</returns>
        public static IsobaricMassTagType? GetTagTypeFromModificationId(string modificationId)
        {
            if (string.IsNullOrWhiteSpace(modificationId))
            {
                return null;
            }

            // Convert to uppercase for case-insensitive comparison
            string upperModId = modificationId.ToUpper();

            // Check for each isobaric tag type using switch expression with pattern matching
            // Order matters: check more specific patterns first (e.g., TMT18 before TMT1)
            return upperModId switch
            {
                string s when s.Contains("TMT18") => IsobaricMassTagType.TMT18,
                string s when s.Contains("TMT11") => IsobaricMassTagType.TMT11,
                string s when s.Contains("TMT10") => IsobaricMassTagType.TMT10,
                string s when s.Contains("TMT6") => IsobaricMassTagType.TMT6,
                string s when s.Contains("ITRAQ-8") || s.Contains("ITRAQ8") => IsobaricMassTagType.iTRAQ8,
                string s when s.Contains("ITRAQ-4") || s.Contains("ITRAQ4") => IsobaricMassTagType.iTRAQ4,
                string s when s.Contains("DILEU-12") || s.Contains("DILEU12") => IsobaricMassTagType.diLeu12,
                string s when s.Contains("DILEU-4") || s.Contains("DILEU4") => IsobaricMassTagType.diLeu4,
                _ => null
            };
        }
    }
}
