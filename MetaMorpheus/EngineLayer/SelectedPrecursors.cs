using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Globalization;
using Chemistry;

namespace EngineLayer
{
    public class SelectedPrecursors
    {
        private static readonly string SelectedPrecursorsHeader = "m/z\tCharge\tRT Start (min)\tRT End (min)";

        /// <summary>
        /// Represents information about a selected precursor for targeted analysis
        /// </summary>
        public class PrecursorInfo
        {
            public double Mz { get; }
            public int Charge { get; }
            public double RtStartInMinutes { get; }
            public double RtEndInMinutes { get; }

            public PrecursorInfo(double mz, int charge, double rtStartInMinutes, double rtEndInMinutes)
            {
                Mz = mz;
                Charge = charge;
                RtStartInMinutes = rtStartInMinutes;
                RtEndInMinutes = rtEndInMinutes;
            }

            public override string ToString()
            {
                return Mz.ToString(CultureInfo.InvariantCulture) + "\t" + Charge + "\t" + RtStartInMinutes.ToString(CultureInfo.InvariantCulture) + "\t" + RtEndInMinutes.ToString(CultureInfo.InvariantCulture);
            }

            public Precursor ToPrecursor()
            {
                return new Precursor(Mz, Charge, Mz.ToMass(Charge), 1.0, 1, null);
            }
        }

        /// <summary>
        /// Reads selected precursors from a tab-delimited file.
        /// </summary>
        /// <param name="precursorsFilePath">Path to the file containing selected precursors data</param>
        /// <param name="errors">List to store any errors encountered during reading</param>
        /// <returns>List of PrecursorInfo objects</returns>
        public static List<PrecursorInfo> ReadSelectedPrecursors(string precursorsFilePath, out List<string> errors)
        {
            var precursors = new List<PrecursorInfo>();
            errors = new List<string>();

            if (!File.Exists(precursorsFilePath))
            {
                errors.Add("Selected precursors file not found!");
                return precursors;
            }

            var lines = File.ReadAllLines(precursorsFilePath);

            for (int i = 1; i < lines.Length; i++) // Skip header line
            {
                var split = lines[i].Split(new char[] { '\t' });

                if (split.Length < 4)
                {
                    errors.Add($"Error: The selected precursors file was not formatted correctly. Expected 4 cells, but found {split.Length} on line {i + 1}");
                    return precursors;
                }

                string strMz = split[0];
                string strCharge = split[1];
                string strRtStart = split[2];
                string strRtEnd = split[3];

                if (!double.TryParse(strMz, NumberStyles.Any, CultureInfo.InvariantCulture, out double mz))
                {
                    errors.Add($"Error: The m/z value on line {i + 1} is not a valid number");
                    return precursors;
                }

                if (!int.TryParse(strCharge, out int charge))
                {
                    errors.Add($"Error: The charge on line {i + 1} is not a valid integer");
                    return precursors;
                }

                if (!double.TryParse(strRtStart, NumberStyles.Any, CultureInfo.InvariantCulture, out double rtStart))
                {
                    errors.Add($"Error: The RT Start value on line {i + 1} is not a valid number");
                    return precursors;
                }

                if (!double.TryParse(strRtEnd, NumberStyles.Any, CultureInfo.InvariantCulture, out double rtEnd))
                {
                    errors.Add($"Error: The RT End value on line {i + 1} is not a valid number");
                    return precursors;
                }

                // Validate the values
                if (mz <= 0)
                {
                    errors.Add($"Error: The m/z value on line {i + 1} must be greater than zero");
                    return precursors;
                }

                if (charge == 0)
                {
                    errors.Add($"Error: The charge on line {i + 1} can not be zero");
                    return precursors;
                }

                if (rtStart < 0)
                {
                    errors.Add($"Error: The RT Start value on line {i + 1} cannot be negative");
                    return precursors;
                }

                if (rtEnd < rtStart)
                {
                    errors.Add($"Error: The RT End value must be greater than or equal to RT Start on line {i + 1}");
                    return precursors;
                }

                var precursorInfo = new PrecursorInfo(mz, charge, rtStart, rtEnd);
                precursors.Add(precursorInfo);
            }

            // Check for duplicate precursors
            var duplicates = precursors
                .GroupBy(p => new { p.Mz, p.Charge })
                .Where(g => g.Count() > 1)
                .Select(g => g.Key);

            if (duplicates.Any())
            {
                var duplicatesList = string.Join(", ", duplicates.Select(d => $"m/z: {d.Mz}, charge: {d.Charge}"));
                errors.Add($"Error: Duplicate precursors found: {duplicatesList}");
            }

            return precursors;
        }

        /// <summary>
        /// Writes a list of precursor information to a tab-delimited file.
        /// </summary>
        /// <param name="precursorInfos">List of precursor information to write</param>
        /// <param name="outputDirectory">Directory where the file will be saved</param>
        /// <param name="fileName">Name of the file (without path)</param>
        /// <returns>Path to the written file</returns>
        public static string WriteSelectedPrecursorsToFile(List<PrecursorInfo> precursorInfos, string outputDirectory, string fileName = "selected_precursors.tsv")
        {
            var filePath = Path.Combine(outputDirectory, fileName);

            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine(SelectedPrecursorsHeader);

                foreach (var precursor in precursorInfos)
                {
                    output.WriteLine(
                        precursor.Mz.ToString(CultureInfo.InvariantCulture) +
                        "\t" + precursor.Charge +
                        "\t" + precursor.RtStartInMinutes.ToString(CultureInfo.InvariantCulture) +
                        "\t" + precursor.RtEndInMinutes.ToString(CultureInfo.InvariantCulture));
                }
            }

            return filePath;
        }

        /// <summary>
        /// Validates a list of precursor information for consistency and returns any errors.
        /// </summary>
        /// <param name="precursorInfos">List of precursor information to validate</param>
        /// <returns>Error message string, or null if no errors found</returns>
        public static string GetErrorsInSelectedPrecursors(List<PrecursorInfo> precursorInfos)
        {
            if (precursorInfos == null || !precursorInfos.Any())
            {
                return "No precursors defined!";
            }

            // Check for overlapping retention time windows for the same precursor (m/z and charge)
            var precursorGroups = precursorInfos.GroupBy(p => new { p.Mz, p.Charge });
            
            foreach (var group in precursorGroups)
            {
                if (group.Count() > 1)
                {
                    // Sort by RT start time
                    var sortedPrecursors = group.OrderBy(p => p.RtStartInMinutes).ToList();
                    
                    // Check for overlaps
                    for (int i = 0; i < sortedPrecursors.Count - 1; i++)
                    {
                        if (sortedPrecursors[i].RtEndInMinutes > sortedPrecursors[i + 1].RtStartInMinutes)
                        {
                            return $"Overlapping retention time windows found for m/z {group.Key.Mz}, charge {group.Key.Charge}";
                        }
                    }
                }
            }

            // Check for invalid retention time ranges
            var invalidRtRanges = precursorInfos.Where(p => p.RtStartInMinutes >= p.RtEndInMinutes);
            if (invalidRtRanges.Any())
            {
                var firstInvalid = invalidRtRanges.First();
                return $"Invalid retention time range: RT Start ({firstInvalid.RtStartInMinutes}) must be less than RT End ({firstInvalid.RtEndInMinutes}) for m/z {firstInvalid.Mz}, charge {firstInvalid.Charge}";
            }

            return null;
        }
    }
}