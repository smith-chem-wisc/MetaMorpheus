using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public abstract class PsmParent
    {

        #region Public Fields

        public double[] quantIntensity;
        public double apexMz;
        public double quantRT;
        public double mostAbundantMass;

        #endregion Public Fields

        #region Protected Constructors

        protected PsmParent(int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan)
        {
            this.Notch = notch;
            this.Score = score;
            this.ScanIndex = scanIndex;
            this.FileName = scan.FileNameWithoutExtension;
            this.ScanNumber = scan.TheScan.OneBasedScanNumber;
            this.PrecursorScanNumber = scan.TheScan.OneBasedPrecursorScanNumber;
            this.ScanRetentionTime = scan.TheScan.RetentionTime;
            this.ScanExperimentalPeaks = scan.TheScan.MassSpectrum.Size;
            this.TotalIonCurrent = scan.TheScan.TotalIonCurrent;
            this.ScanPrecursorCharge = scan.PrecursorCharge;
            this.ScanPrecursorMonoisotopicPeak = scan.PrecursorMonoisotopicPeak;
            this.ScanPrecursorMass = scan.PrecursorMass;
            quantIntensity = new double[1];
        }

        #endregion Protected Constructors

        #region Public Properties

        public int Notch { get; }
        public double Score { get; }
        public int ScanNumber { get; }
        public int PrecursorScanNumber { get; }
        public double ScanRetentionTime { get; }
        public int ScanExperimentalPeaks { get; }
        public double TotalIonCurrent { get; }
        public int ScanPrecursorCharge { get; }
        public IMzPeak ScanPrecursorMonoisotopicPeak { get; }
        public double ScanPrecursorMass { get; }
        public string FileName { get; }

        public int ScanIndex { get; }

        public int NumAmbiguous { get; internal set; }

        public ProteinLevelInfo Pli { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public static double MatchIons(IMsDataScan<IMzSpectrum<IMzPeak>> thisScan, Tolerance productMassTolerance, double[] sorted_theoretical_product_masses_for_this_peptide, double[] matchedIonMassesListPositiveIsMatch)
        {
            var TotalProductsHere = sorted_theoretical_product_masses_for_this_peptide.Length;
            if (TotalProductsHere == 0)
                return 0;
            int MatchingProductsHere = 0;
            double MatchingIntensityHere = 0;

            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.XArray;
            double[] experimental_intensities = thisScan.MassSpectrum.YArray;
            int num_experimental_peaks = experimental_mzs.Length;

            int currentTheoreticalIndex = -1;
            double currentTheoreticalMass;
            do
            {
                currentTheoreticalIndex++;
                currentTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[currentTheoreticalIndex];
            } while (double.IsNaN(currentTheoreticalMass) && currentTheoreticalIndex < sorted_theoretical_product_masses_for_this_peptide.Length - 1);

            if (double.IsNaN(currentTheoreticalMass))
                return 0;

            double currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;

            int testTheoreticalIndex;
            double testTheoreticalMZ;
            double testTheoreticalMass;
            // Loop over all experimenal indices
            for (int experimentalIndex = 0; experimentalIndex < num_experimental_peaks; experimentalIndex++)
            {
                double currentExperimentalMZ = experimental_mzs[experimentalIndex];
                // If found match
                if (productMassTolerance.Within(currentExperimentalMZ, currentTheoreticalMz))
                {
                    MatchingProductsHere++;
                    MatchingIntensityHere += experimental_intensities[experimentalIndex];
                    matchedIonMassesListPositiveIsMatch[currentTheoreticalIndex] = currentTheoreticalMass;
                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                        break;
                    currentTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;
                }
                // Else if for sure did not reach the next theoretical yet, move to next experimental
                else if (currentExperimentalMZ < currentTheoreticalMz)
                    continue;
                // Else if for sure passed a theoretical
                else
                {
                    // Mark the theoretical as missed
                    matchedIonMassesListPositiveIsMatch[currentTheoreticalIndex] = -currentTheoreticalMass;

                    // Move on to next index and never come back!
                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                        break;
                    currentTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;

                    // Start with the current ones
                    testTheoreticalIndex = currentTheoreticalIndex;
                    testTheoreticalMZ = currentTheoreticalMz;
                    testTheoreticalMass = currentTheoreticalMass;
                    // Mark the skipped theoreticals as not found. The last one is not for sure, might be flipped!
                    while (currentExperimentalMZ > testTheoreticalMZ)
                    {
                        matchedIonMassesListPositiveIsMatch[testTheoreticalIndex] = -currentTheoreticalMass;
                        // Store old info for possible reuse
                        currentTheoreticalMass = testTheoreticalMass;
                        currentTheoreticalMz = testTheoreticalMZ;
                        currentTheoreticalIndex = testTheoreticalIndex;

                        // Update test stuff!
                        testTheoreticalIndex++;
                        if (testTheoreticalIndex == TotalProductsHere)
                            break;
                        testTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[testTheoreticalIndex];
                        testTheoreticalMZ = testTheoreticalMass + Constants.protonMass;
                    }

                    experimentalIndex--;
                }
            }
            return MatchingProductsHere + MatchingIntensityHere / thisScan.TotalIonCurrent;
        }

        public abstract CompactPeptide GetCompactPeptide(Dictionary<ModificationWithMass, ushort> modsDictionary);

        public void ComputeProteinLevelInfo(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> matching, Tolerance fragmentTolerance, Ms2ScanWithSpecificMass theScan, List<ProductType> lp, Dictionary<Proteomics.ModificationWithMass, ushort> modsDictionary)
        {
            Pli = new ProteinLevelInfo(matching[GetCompactPeptide(modsDictionary)], fragmentTolerance, theScan, lp);
        }

        #endregion Public Methods

        #region Internal Methods

        internal static string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("fileName" + '\t');
            sb.Append("scanNumber" + '\t');
            sb.Append("scanRetentionTime" + '\t');
            sb.Append("scanExperimentalPeaks" + '\t');
            sb.Append("totalIonCurrent" + '\t');
            sb.Append("precursorScanNumber" + '\t');
            sb.Append("scanPrecursorCharge" + '\t');
            sb.Append("scanPrecursorMZ" + '\t');
            sb.Append("scanPrecursorMass" + '\t');
            sb.Append("score" + '\t');
            sb.Append("notch" + '\t');
            sb.Append("quantificationIntensity" + '\t');
            sb.Append("quantificationRT" + '\t');
            sb.Append("matched ions" + '\t');
            sb.Append("matched ion counts" + '\t');
            sb.Append("localized scores" + '\t');
            sb.Append("improvement" + '\t');
            sb.Append("terminal localization");
            sb.Append("Protein Accession" + '\t');
            sb.Append("Protein FullName" + '\t');
            sb.Append("Peptide Description" + '\t');
            sb.Append("Start and End ResidueInProtein" + '\t');
            sb.Append("PreviousAminoAcid" + '\t');
            sb.Append("NextAminoAcid" + '\t');
            sb.Append("BaseSequence" + '\t');
            sb.Append("FullSequence" + '\t');
            sb.Append("numVariableMods" + '\t');
            sb.Append("MissedCleavages" + '\t');
            sb.Append("PeptideMonoisotopicMass" + '\t');
            sb.Append("MassDiff (Da)" + '\t');
            sb.Append("MassDiff (ppm)" + '\t');
            sb.Append("Decoy/Contaminant/Target");
            return sb.ToString();
        }

        #endregion Internal Methods

    }
}