using System;
using System.Collections.Generic;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Provides Filters and error handling for GUI forms
    /// </summary>
    public class GlobalGuiSettings
    {
        /// <summary>
        /// Conducts all checks for fields to all forms
        /// </summary>
        /// <param name="precursorMassTolerance"></param>
        /// <param name="productMassTolerance"></param>
        /// <param name="maxMissedCleavages"></param>
        /// <param name="maxModificationIsoforms"></param>
        /// <param name="minPeptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <param name="maxThreads"></param>
        /// <param name="minScore"></param>
        /// <returns></returns>
        public static bool CheckGeneralFilters(string precursorMassTolerance,
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
            string maxFragmentMass
            )
        {
            List<string> results = new List<string>();
            results.Add((CheckPrecursorMassTolerance(precursorMassTolerance)).ToString());
            results.Add((CheckProductMassTolerance(productMassTolerance)).ToString());
            results.Add((CheckMaxMissedCleavages(maxMissedCleavages)).ToString());
            results.Add((CheckMaxModificationIsoForms(maxModificationIsoforms)).ToString());
            results.Add((CheckMinPeptideLength(minPeptideLength)).ToString());
            results.Add((CheckMaxPeptideLength(maxPeptideLength)).ToString());
            results.Add((CheckMaxThreads(maxThreads)).ToString());
            results.Add((CheckMinScoreAllowed(minScore)).ToString());
            results.Add((CheckPeakFindingTolerance(peakFindingTolerance)).ToString());
            results.Add((CheckHistogramBinWidth(histogramBinWidth)).ToString());
            results.Add((CheckDeconvolutionMaxAssumedChargeState(deconMaxAssumedCharge)).ToString());
            results.Add((CheckTopNPeaks(numPeaks)).ToString());
            results.Add((CheckMinRatio(minRatio)).ToString());
            results.Add((CheckNumberOfDatabasePartitions(numberOfDatabaseSearches)).ToString());
            results.Add((CheckMaxModsPerPeptide(maxModsPerPeptide)).ToString());
            results.Add((CheckMaxFragementMass(maxFragmentMass)).ToString());

            if (results.Contains("-1"))
            {
                return false;
            }
            if (CheckMinPeptideLength(minPeptideLength) > CheckMaxPeptideLength(maxPeptideLength))
            {
                MessageBox.Show("The maximum peptide length must be greater than or equal to the minimum peptide length.");
                return false;
            }
            return true;
        }

        /// <summary>
        /// Filters out illgal characters
        /// </summary>
        /// <param name="text"></param>
        /// <returns></returns>
        public static bool CheckIsNumber(string text)
        {
           bool result = true;
           foreach (var character in text)
            {
                if (!Char.IsDigit(character) && !(character=='.'))
                {
                    result = false;
                }
            }
            return result;
        }

        #region Check Task Validity

        public static int CheckDeconvolutionMaxAssumedChargeState(string text)
        {
            if (!double.TryParse(text, out double deconMaxAssumedCharge) || deconMaxAssumedCharge < 1)
            {
                MessageBox.Show("The maximum assumed charge state for deconvolution is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive, non-zero number.");
                return -1;
            }
            return (int)deconMaxAssumedCharge;
        }

        public static int CheckTopNPeaks(string text)
        {            
            if (text.Length == 0)
            {
                text = int.MaxValue.ToString();
            }

            if (!double.TryParse(text, out double numPeaks) || numPeaks < 1)
            {
                MessageBox.Show("The Top N Peaks to be retained must be greater than zero. \n You entered " + '"' + text + '"' + "\n Please enter a positive number.");
                return 1;
            }
            return (int)numPeaks;
        }

        public static double CheckMinRatio(string text)
        {
            if (!double.TryParse(text, out double minRatio) || minRatio < 0 || minRatio > 1)
            {
                MessageBox.Show("The minimum ratio was not set to a number between zero and one. \n You entered " + '"' + text + '"');
                return -1;
            }
            return minRatio;
        }

        public static double CheckPrecursorMassTolerance(string text)
        {
          
            if (!double.TryParse(text, out double precursorMassTolerance) || precursorMassTolerance <= 0)
            {
                MessageBox.Show("The precursor mass tolerance is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive number.");
                return -1;
            }
            return precursorMassTolerance;
        }

        public static double CheckProductMassTolerance(string text)
        {
            if (!double.TryParse(text, out double productMassTolerance) || productMassTolerance <= 0)
            {
                MessageBox.Show("The product mass tolerance is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive number.");
                return -1;
            }
            return productMassTolerance;
        }

        public static int CheckNumberOfDatabasePartitions(string text)
        {
            if (!double.TryParse(text, out double numberOfDatabaseSearches) || numberOfDatabaseSearches <= 0)
            {
                MessageBox.Show("The number of database partitions is invalid. At least one database is required for searching.");
                return -1;
            }
            return (int)numberOfDatabaseSearches;
        }

        public static int CheckMaxMissedCleavages(string text)
        {
            if (string.IsNullOrEmpty(text))
            {
                return int.MaxValue;
            }

            if (!double.TryParse(text, out double maxMissedCleavages) || maxMissedCleavages < 0)
            {
                MessageBox.Show("The number of missed cleavages is invalid. Please enter an integer zero or greater.");
                return -1;
            }

            return (int)maxMissedCleavages;
        }

        public static int CheckMaxModificationIsoForms(string text)
        {
            if (!double.TryParse(text, out double maxModificationIsoforms) || maxModificationIsoforms < 1)
            {
                MessageBox.Show("The maximum number of modification isoforms is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive, non-zero number.");
                return -1;
            }
            return (int)maxModificationIsoforms;
        }

        public static int CheckMinPeptideLength(string text)
        {
            if (!double.TryParse(text, out double minPeptideLength) || minPeptideLength < 1)
            {
                MessageBox.Show("The minimum peptide length must be a positive integer");
                return -1;
            }
            return (int)minPeptideLength;
        }

        public static int CheckMaxPeptideLength(string text)
        {
            if (string.IsNullOrEmpty(text))
            {
                return int.MaxValue;
            }
            if (!double.TryParse(text, out double maxPeptideLength) || maxPeptideLength < 1)
            {
                MessageBox.Show("The minimum peptide length must be a positive integer");
                return -1;
            }
            return (int)maxPeptideLength;

            // TO DO:
            // may have to stay in each class for task
            //
            //if (maxPeptideLength < minPeptideLength)
            //{
            //    MessageBox.Show("The maximum peptide length must be greater than or equal to the minimum peptide length.");
            //    return;
            //}
        }

        public static int CheckMaxModsPerPeptide(string text)
        {
            if (!double.TryParse(text, out double maxModsPerPeptide) || maxModsPerPeptide < 1)
            {
                MessageBox.Show("The mods per peptide allowed is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive, non-zero number.");
                return -1;
            }
            return (int)maxModsPerPeptide;
        }

        public static double CheckMaxFragementMass(string text)
        {
            if (!double.TryParse(text, out double maxFragmentMass) || maxFragmentMass < 0)
            {
                MessageBox.Show("The fragment mass is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive, non-zero number.");
                return -1;
            }
            return maxFragmentMass;
        }

        public static int CheckMaxThreads(string text)
        {
            if (!double.TryParse(text, out double maxThreads) || maxThreads > Environment.ProcessorCount || maxThreads < 1)
            {
                MessageBox.Show("Your current device has " + Environment.ProcessorCount + " processors. \n Please select a positive value less than or equal to this number.");
                return -1;
            }
            return (int)maxThreads;
        }

        public static double CheckMinScoreAllowed(string text)
        {
            if (!double.TryParse(text, out double minScore) || minScore < 1)
            {
                MessageBox.Show("The minimum score allowed is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive, non-zero number.");
                return -1;
            }
            return minScore;
        }

        public static double CheckPeakFindingTolerance(string text)
        {
            if (!double.TryParse(text, out double peakFindingTolerance) || peakFindingTolerance <= 0)
            {
                MessageBox.Show("The peak finding tolerance is invalid. \n You entered " + '"' + text + '"' + "\n Please enter a positive number.");
                return -1;
            }
            return peakFindingTolerance;
        }

        public static float CheckHistogramBinWidth(string text)
        {
            if (!float.TryParse(text, out float binWidth) || binWidth < 0 || binWidth > 1)
            {
                MessageBox.Show("The precursor mass tolerance was not set to a number between zero and one. \n You entered " + '"' + text + '"' );
                return -1;
            }
            return binWidth;
        }

        //public static int Check()
        //{

        //}
        #endregion
    }
}
