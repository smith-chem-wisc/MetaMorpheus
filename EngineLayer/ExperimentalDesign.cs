using FlashLFQ;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

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
        /// </summary>
        public static string GetErrorsInExperimentalDesign(List<SpectraFileInfo> spectraFileInfos)
        {
            // check for correct iteration of integer values and duplicates
            var conditions = spectraFileInfos.GroupBy(p => p.Condition);

            foreach (var condition in conditions)
            {
                var temp = condition.OrderBy(p => p.BiologicalReplicate).ThenBy(p => p.Fraction).ThenBy(p => p.TechnicalReplicate);
                int numB = temp.Max(p => p.BiologicalReplicate + 1);

                // check bioreps are in order
                for (int b = 0; b < numB; b++)
                {
                    var biorepFiles = temp.Where(p => p.BiologicalReplicate == b);

                    if (!biorepFiles.Any())
                    {
                        return "Condition \"" + condition.Key + "\" biorep " + (b + 1) + " is missing!";
                    }

                    // check fractions are in order
                    int numF = biorepFiles.Max(p => p.Fraction + 1);

                    for (int f = 0; f < numF; f++)
                    {
                        var fractionFiles = biorepFiles.Where(p => p.Fraction == f);

                        if (!fractionFiles.Any())
                        {
                            return "Condition \"" + condition.Key + "\" biorep " + (b + 1) + " fraction " + (f + 1) + " is missing!";
                        }

                        // check techreps are in order
                        int numT = fractionFiles.Max(p => p.TechnicalReplicate + 1);

                        for (int t = 0; t < numT; t++)
                        {
                            var techrepFiles = fractionFiles.Where(p => p.TechnicalReplicate == t);

                            if (!techrepFiles.Any())
                            {
                                return "Condition \"" + condition.Key + "\" biorep " + (b + 1) + " fraction " + (f + 1) + " techrep " + (t + 1) + " is missing!";
                            }

                            if (techrepFiles.Count() > 1)
                            {
                                return "Duplicates are not allowed:\n" +
                                    "Condition \"" + condition.Key + "\" biorep " + (b + 1) + " fraction " + (f + 1) + " techrep " + (t + 1);
                            }
                        }
                    }
                }
            }

            return null;
        }
    }
}
