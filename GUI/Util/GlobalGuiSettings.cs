using System;
using System.Collections.Generic;
using System.Globalization;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Provides filters and error handling for GUI forms
    /// </summary>
    public class GlobalGuiSettings
    {
        /// <summary>
        /// Checks the validity of each setting passed from the GUI task windows
        /// </summary>
        public static bool CheckTaskSettingsValidity(string precursorMassTolerance,
            string productMassTolerance,
            string maxMissedCleavages,
            string maxModificationIsoforms,
            string minPeptideLength,
            string maxPeptideLength,
            string maxThreads,
            string minScore,
            string peakFindingTolerance,
            string histogramBinWidth,
            string deconMaxAssumedCharge,
            string numPeaks,
            string minRatio,
            string numberOfDatabaseSearches,
            string maxModsPerPeptide,
            string maxFragmentMass,
            string qValueFilter
            )
        {
            maxMissedCleavages = MaxValueConversion(maxMissedCleavages);
            maxPeptideLength = MaxValueConversion(maxPeptideLength);

            List<bool> results = new List<bool>();
            results.Add((CheckPrecursorMassTolerance(precursorMassTolerance)));
            results.Add((CheckProductMassTolerance(productMassTolerance)));
            results.Add((CheckMaxMissedCleavages(maxMissedCleavages)));
            results.Add((CheckMaxModificationIsoForms(maxModificationIsoforms)));
            results.Add((CheckPeptideLength(minPeptideLength, maxPeptideLength)));
            results.Add((CheckMaxThreads(maxThreads)));
            results.Add((CheckMinScoreAllowed(minScore)));
            results.Add((CheckPeakFindingTolerance(peakFindingTolerance)));
            results.Add((CheckHistogramBinWidth(histogramBinWidth)));
            results.Add((CheckDeconvolutionMaxAssumedChargeState(deconMaxAssumedCharge)));
            results.Add((CheckTopNPeaks(numPeaks)));
            results.Add((CheckMinRatio(minRatio)));
            results.Add((CheckNumberOfDatabasePartitions(numberOfDatabaseSearches)));
            results.Add((CheckMaxModsPerPeptide(maxModsPerPeptide)));
            results.Add((CheckMaxFragementMass(maxFragmentMass)));
            results.Add((CheckQValueFilter(qValueFilter)));

            if (results.Contains(false))
            {
                return false;
            }
            return true;
        }

        /// <summary>
        /// Checks to see if the given text contains non-numerical characters (letters, etc.)
        /// </summary>
        public static bool CheckIsNumber(string text)
        {
            bool result = true;
            foreach (var character in text)
            {
                if (!Char.IsDigit(character) && !(character == '.') && !(character == '-'))
                {
                    result = false;
                }
            }
            return result;
        }

        #region Check Task Validity

        public static string MaxValueConversion(string text)
        {
            if (string.IsNullOrEmpty(text))
            {
                text = int.MaxValue.ToString();
            }
            return text;
        }

        public static bool CheckDeconvolutionMaxAssumedChargeState(string text)
        {
            if (!int.TryParse(text, out int deconMaxAssumedCharge) || deconMaxAssumedCharge < 1)
            {
                MessageBox.Show("The maximum assumed charge state for deconvolution is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive, non-zero number.");
                return false;
            }
            return true;
        }

        public static bool CheckTopNPeaks(string text)
        {
            if (text.Length == 0)
            {
                text = int.MaxValue.ToString();
            }

            if (!int.TryParse(text, out int numPeaks) || numPeaks < 1)
            {
                MessageBox.Show("The Top N Peaks to be retained must be greater than zero. \n You entered " + '"' + text + '"' + "\n Please enter a positive number.");
                return false;
            }
            return true;
        }

        public static bool CheckMinRatio(string text)
        {
            if (!double.TryParse(text, NumberStyles.Any, CultureInfo.InvariantCulture, out double minRatio) || minRatio < 0 || minRatio > 1)
            {
                MessageBox.Show("The minimum intensity ratio must be between zero and one. \n You entered " + '"' + text + '"');
                return false;
            }
            return true;
        }

        public static bool CheckPrecursorMassTolerance(string text)
        {

            if (!double.TryParse(text, NumberStyles.Any, CultureInfo.InvariantCulture, out double precursorMassTolerance) || precursorMassTolerance <= 0)
            {
                MessageBox.Show("The precursor mass tolerance is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive number.");
                return false;
            }
            return true;
        }

        public static bool CheckProductMassTolerance(string text)
        {
            if (!double.TryParse(text, NumberStyles.Any, CultureInfo.InvariantCulture, out double productMassTolerance) || productMassTolerance <= 0)
            {
                MessageBox.Show("The product mass tolerance is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive number.");
                return false;
            }
            return true;
        }

        public static bool CheckNumberOfDatabasePartitions(string text)
        {
            if (!int.TryParse(text, out int numberOfDatabaseSearches) || numberOfDatabaseSearches <= 0)
            {
                MessageBox.Show("The number of database partitions is invalid. At least one database is required for searching.");
                return false;
            }
            return true;
        }

        public static bool CheckMaxMissedCleavages(string text)
        {
            if (!int.TryParse(text, out int maxMissedCleavages) || maxMissedCleavages < 0)
            {
                MessageBox.Show("The number of missed cleavages is invalid. Please enter an integer zero or greater.");
                return false;
            }

            return true;
        }

        public static bool CheckMaxModificationIsoForms(string text)
        {
            if (!int.TryParse(text, out int maxModificationIsoforms) || maxModificationIsoforms < 1)
            {
                MessageBox.Show("The maximum number of modification isoforms is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive, non-zero number.");
                return false;
            }
            return true;
        }

        public static bool CheckPeptideLength(string min, string max)
        {
            if (!int.TryParse(min, out int minPeptideLength) || minPeptideLength < 1)
            {
                MessageBox.Show("The minimum peptide length must be a positive integer");
                return false;
            }

            if (!int.TryParse(max, out int maxPeptideLength) || maxPeptideLength < 1)
            {
                MessageBox.Show("The maximum peptide length must be a positive integer");
                return false;
            }

            if (Convert.ToInt32(min) > Convert.ToInt32(max))
            {
                MessageBox.Show("The maximum peptide length must be greater than or equal to the minimum peptide length.");
                return false;
            }
            return true;
        }

        public static bool CheckMaxModsPerPeptide(string text)
        {
            if (!int.TryParse(text, out int maxModsPerPeptide) || maxModsPerPeptide < 0)
            {
                MessageBox.Show("The max mods per peptide allowed is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a number greater than or equal to zero.");
                return false;
            }
            return true;
        }

        public static bool CheckMaxFragementMass(string text)
        {
            if (!int.TryParse(text, out int maxFragmentMass) || maxFragmentMass < 0)
            {
                MessageBox.Show("The max fragment mass is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive, non-zero number.");
                return false;
            }
            return true;
        }

        public static bool CheckMaxThreads(string text)
        {
            if (!int.TryParse(text, out int maxThreads) || maxThreads > Environment.ProcessorCount || maxThreads < 1)
            {
                MessageBox.Show("Your current device has " + Environment.ProcessorCount + " processors. \n Please select a positive value less than or equal to this number.");
                return false;
            }
            return true;
        }

        public static bool CheckMinScoreAllowed(string text)
        {
            if (!double.TryParse(text, NumberStyles.Any, CultureInfo.InvariantCulture, out double minScore) || minScore < 1)
            {
                MessageBox.Show("The minimum score allowed is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive, non-zero number.");
                return false;
            }
            return true;
        }

        public static bool CheckPeakFindingTolerance(string text)
        {
            if (!double.TryParse(text, NumberStyles.Any, CultureInfo.InvariantCulture, out double peakFindingTolerance) || peakFindingTolerance <= 0)
            {
                MessageBox.Show("The peak finding tolerance is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive number.");
                return false;
            }
            return true;
        }

        public static bool CheckHistogramBinWidth(string text)
        {
            if (!float.TryParse(text, out float binWidth) || binWidth < 0 || binWidth > 1)
            {
                MessageBox.Show("The histogram bin width must be between zero and one Daltons. \n You entered " + '"' + text + '"');
                return false;
            }
            return true;
        }

        public static bool CheckQValueFilter(string text)
        {
            if (!double.TryParse(text, NumberStyles.Any, CultureInfo.InvariantCulture, out double qValue) || qValue < 0 || qValue > 1)
            {
                MessageBox.Show("The q-value cutoff must be a number between 0 and 1");
                return false;
            }
            return true;
        }

        public static bool VariableModCheck(List<(string, string)> listOfModsVariable)
        {
            if (listOfModsVariable.Count > 1)
            {
                var dialogResult = MessageBox.Show("More than 1 modification has been selected as variable. Using the GPTMD task to discover modifications is recommended instead. \n\nContinue anyway?", "Multiple Variable Mods Detected", MessageBoxButton.OKCancel);
                if (dialogResult == MessageBoxResult.Cancel)
                {
                    return false;
                }
            }
            return true;
        }

        #endregion
    }
}
