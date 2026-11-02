using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;

namespace EngineLayer
{
    public class ExperimentalDesign
    {
        private static string ExperimentalDesignHeader = "FileName\tCondition\tBiorep\tFraction\tTechrep";

        public static List<SpectraFileInfo> ReadExperimentalDesign(string experimentalDesignPath, List<string> fullFilePathsWithExtension, out List<string> errors)
        {
            var expDesign = new List<SpectraFileInfo>();
            errors = new List<string>();

            if (!File.Exists(experimentalDesignPath))
            {
                errors.Add("Experimental design file not found!");
                return expDesign;
            }

            var lines = File.ReadAllLines(experimentalDesignPath);

            for (int i = 1; i < lines.Length; i++)
            {
                var split = lines[i].Split(new char[] { '\t' });

                if (split.Length < 5)
                {
                    errors.Add("Error: The experimental design was not formatted correctly. Expected 5 cells, but found " + split.Length + " on line " + (i + 1));
                    return expDesign;
                }

                string fileNameWithExtension = split[0];
                string condition = split[1];
                string strBiorep = split[2];
                string strFraction = split[3];
                string strTechrep = split[4];

                if (!int.TryParse(strBiorep, out int biorep))
                {
                    errors.Add("Error: The experimental design was not formatted correctly. The biorep on line " + (i + 1) + " is not an integer");
                    return expDesign;
                }
                if (!int.TryParse(strFraction, out int fraction))
                {
                    errors.Add("Error: The experimental design was not formatted correctly. The fraction on line " + (i + 1) + " is not an integer");
                    return expDesign;
                }
                if (!int.TryParse(strTechrep, out int techrep))
                {
                    errors.Add("Error: The experimental design was not formatted correctly. The techrep on line " + (i + 1) + " is not an integer");
                    return expDesign;
                }

                var foundFilePath = fullFilePathsWithExtension.FirstOrDefault(p => Path.GetFileName(p) == fileNameWithExtension);
                if (foundFilePath == null)
                {
                    // the experimental design could include files that aren't in the spectra file list but that's ok.
                    // it's fine to have extra files defined in the experimental design as long as the remainder is valid
                    continue;
                }

                var fileInfo = new SpectraFileInfo(foundFilePath, condition, biorep - 1, techrep - 1, fraction - 1);
                expDesign.Add(fileInfo);
            }

            // check to see if there are any files missing from the experimental design
            var filesDefinedInExpDesign = expDesign.Select(p => p.FullFilePathWithExtension).ToList();
            var notDefined = fullFilePathsWithExtension.Where(p => !filesDefinedInExpDesign.Contains(p));
            if (notDefined.Any())
            {
                errors.Add("Error: The experimental design did not contain the file(s): " + string.Join(", ", notDefined));
                return expDesign;
            }

            // check to see if the design is valid
            var designError = GetErrorsInExperimentalDesign(expDesign);
            if (designError != null)
            {
                errors.Add(designError);
                return expDesign;
            }

            // all files passed in are defined in the experimental design and the exp design is valid
            return expDesign;
        }

        public static string WriteExperimentalDesignToFile(List<SpectraFileInfo> spectraFileInfos)
        {
            var dir = Directory.GetParent(spectraFileInfos.First().FullFilePathWithExtension).FullName;
            var filePath = Path.Combine(dir, GlobalVariables.ExperimentalDesignFileName);

            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine(ExperimentalDesignHeader);

                foreach (var spectraFile in spectraFileInfos)
                {
                    output.WriteLine(
                        Path.GetFileName(spectraFile.FullFilePathWithExtension) +
                        "\t" + spectraFile.Condition +
                        "\t" + (spectraFile.BiologicalReplicate + 1) +
                        "\t" + (spectraFile.Fraction + 1) +
                        "\t" + (spectraFile.TechnicalReplicate + 1));
                }
            }

            return filePath;
        }

        /// <summary>
        /// Checks for errors in the experimental design. Will return null if there are no errors.
        /// The one error is a duplicate: two files at the same condition, biorep, fraction and techrep.
        /// A gap in biorep, fraction or techrep numbers is not an error. The numbers are kept as given
        /// (a biorep number can name a subject across conditions, and a gap can be a lost sample) and
        /// reported by GetWarningsInExperimentalDesign.
        /// </summary>
        public static string GetErrorsInExperimentalDesign(List<SpectraFileInfo> spectraFileInfos)
        {
            var duplicate = spectraFileInfos
                .GroupBy(p => (p.Condition, p.BiologicalReplicate, p.Fraction, p.TechnicalReplicate))
                .Where(g => g.Count() > 1)
                .Select(g => g.Key)
                .OrderBy(k => k.Condition, StringComparer.Ordinal).ThenBy(k => k.BiologicalReplicate).ThenBy(k => k.Fraction).ThenBy(k => k.TechnicalReplicate)
                .FirstOrDefault();

            if (duplicate != default)
            {
                return "Duplicates are not allowed:\n" +
                    "Condition \"" + duplicate.Condition + "\" biorep " + (duplicate.BiologicalReplicate + 1) +
                    " fraction " + (duplicate.Fraction + 1) + " techrep " + (duplicate.TechnicalReplicate + 1);
            }

            return null;
        }

        /// <summary>
        /// One warning for each place the numbering has a gap, naming the missing numbers: a condition's
        /// bioreps, a biorep's fractions, and a fraction's techreps, each expected to run 1..N. Such a design
        /// is quantified as numbered; the warning is there because a missing number may be a sample or
        /// file that was lost or not searched. Empty when there is nothing to say.
        /// </summary>
        public static List<string> GetWarningsInExperimentalDesign(List<SpectraFileInfo> spectraFileInfos)
        {
            var warnings = new List<string>();

            foreach (var condition in spectraFileInfos.GroupBy(p => p.Condition).OrderBy(c => c.Key, StringComparer.Ordinal))
            {
                string conditionLabel = "Condition \"" + condition.Key + "\"";
                AddGapWarning(warnings, conditionLabel, "biological replicate", condition.Select(p => p.BiologicalReplicate));

                foreach (var biorep in condition.GroupBy(p => p.BiologicalReplicate).OrderBy(b => b.Key))
                {
                    string biorepLabel = conditionLabel + " biorep " + (biorep.Key + 1);
                    AddGapWarning(warnings, biorepLabel, "fraction", biorep.Select(p => p.Fraction));

                    foreach (var fraction in biorep.GroupBy(p => p.Fraction).OrderBy(f => f.Key))
                    {
                        AddGapWarning(warnings, biorepLabel + " fraction " + (fraction.Key + 1), "technical replicate", fraction.Select(p => p.TechnicalReplicate));
                    }
                }
            }

            return warnings;
        }

        /// <summary>
        /// Adds a warning when the given zero-based numbers, read one-based, are not 1..N. Worded as mzLib's
        /// SdrfLabelFreeDesign notes the same gap ("Condition 'A' biorep 1: fractions 1, 3, kept as the SDRF
        /// numbers them."), so the SDRF report and the run say one thing.
        /// </summary>
        private static void AddGapWarning(List<string> warnings, string owner, string level, IEnumerable<int> zeroBasedNumbers)
        {
            var present = zeroBasedNumbers.Select(n => n + 1).Distinct().OrderBy(n => n).ToList();
            var missing = Enumerable.Range(1, Math.Max(0, present.Max())).Except(present).ToList();

            if (missing.Any())
            {
                warnings.Add(owner + ": " + level + "s " + string.Join(", ", present) + ", quantified as numbered; " +
                    level + (missing.Count == 1 ? " " : "s ") + string.Join(", ", missing) + (missing.Count == 1 ? " is" : " are") + " not in the design. " +
                    "A missing number may be a sample or file that was lost or not searched.");
            }
        }
    }
}
