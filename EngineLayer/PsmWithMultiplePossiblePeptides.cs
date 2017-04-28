using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class PsmWithMultiplePossiblePeptides
    {

        #region Public Fields

        public PsmParent newPsm;

        #endregion Public Fields

        #region Public Constructors

        public PsmWithMultiplePossiblePeptides(PsmParent newPsm, HashSet<PeptideWithSetModifications> peptidesWithSetModifications, Tolerance fragmentTolerance, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, List<ProductType> lp)
        {
            this.newPsm = newPsm;
            IsDecoy = peptidesWithSetModifications.Any(b => b.Protein.IsDecoy);
            IsContaminant = peptidesWithSetModifications.Any(b => b.Protein.IsContaminant);
            this.PeptidesWithSetModifications = peptidesWithSetModifications;

            var representative = peptidesWithSetModifications.First();

            IMsDataScan<IMzSpectrum<IMzPeak>> theScan;
            if (myMsDataFile != null && newPsm.matchedIonsListPositiveIsMatch == null)
            {
                theScan = myMsDataFile.GetOneBasedScan(newPsm.scanNumber);

                var MatchedIonDictPositiveIsMatch = new Dictionary<ProductType, double[]>();
                foreach (var huh in lp)
                {
                    var df = representative.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { huh });
                    Array.Sort(df);
                    double[] matchedIonMassesListPositiveIsMatch = new double[df.Length];
                    MatchIons(theScan, fragmentTolerance, df, matchedIonMassesListPositiveIsMatch);
                    MatchedIonDictPositiveIsMatch.Add(huh, matchedIonMassesListPositiveIsMatch);
                }

                newPsm.matchedIonsListPositiveIsMatch = MatchedIonDictPositiveIsMatch;
            }

            if (myMsDataFile != null && newPsm.LocalizedScores == null)
            {
                theScan = myMsDataFile.GetOneBasedScan(newPsm.scanNumber);
                var localizedScores = new List<double>();
                for (int indexToLocalize = 0; indexToLocalize < representative.Length; indexToLocalize++)
                {
                    PeptideWithSetModifications localizedPeptide = representative.Localize(indexToLocalize, ScanPrecursorMass - representative.MonoisotopicMass);

                    var gg = localizedPeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                    Array.Sort(gg);
                    double[] matchedIonMassesListPositiveIsMatch = new double[gg.Length];
                    var score = MatchIons(theScan, fragmentTolerance, gg, matchedIonMassesListPositiveIsMatch);
                    localizedScores.Add(score);
                }
                newPsm.LocalizedScores = localizedScores;
            }

            PeptideMonoisotopicMass = representative.MonoisotopicMass;

            // Look for better match in MS1 spectrum!!!
            if (myMsDataFile != null && !newPsm.precursorScanBestMass.HasValue)
            {
                var precursorScan = myMsDataFile.GetOneBasedScan(newPsm.precursorScanNumber);
                newPsm.precursorScanBestMass = precursorScan.MassSpectrum.GetClosestPeakXvalue(this.PeptideMonoisotopicMass.ToMz(this.newPsm.scanPrecursorCharge)).ToMass(this.newPsm.scanPrecursorCharge);
            }

            FullSequence = representative.Sequence;
            BaseSequence = representative.BaseSequence;
            MissedCleavages = representative.MissedCleavages;
            NumVariableMods = representative.NumMods - representative.numFixedMods;
        }

        #endregion Public Constructors

        #region Public Properties

        public HashSet<PeptideWithSetModifications> PeptidesWithSetModifications { get; private set; }
        public bool IsDecoy { get; private set; }
        public bool IsContaminant { get; private set; }

        public double Score
        {
            get
            {
                return newPsm.score;
            }
        }

        public double ScanPrecursorMass
        {
            get
            {
                return newPsm.scanPrecursorMass;
            }
        }

        public List<double> LocalizedScores
        {
            get
            {
                return newPsm.LocalizedScores;
            }
        }

        public double PeptideMonoisotopicMass { get; private set; }

        public string FullSequence { get; private set; }

        public string BaseSequence { get; private set; }

        public int MissedCleavages { get; private set; }

        public int NumVariableMods { get; private set; }

        public string SequenceWithChemicalFormulas
        {
            get
            {
                return PeptidesWithSetModifications.First().SequenceWithChemicalFormulas;
            }
        }

        #endregion Public Properties

        #region Internal Properties

        internal static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append(PsmParent.GetTabSeparatedHeader() + '\t');
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
                sb.Append("BestMassInPrecursor" + '\t');
                sb.Append("MassDiffToBestMass (Da)" + '\t');
                sb.Append("MassDiffToBestMass (ppm)" + '\t');
                sb.Append("Decoy/Contaminant/Target");
                return sb.ToString();
            }
        }

        #endregion Internal Properties

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

        public override string ToString()
        {
            var sb = new StringBuilder();

            sb.Append(newPsm.ToString() + '\t');

            sb.Append(string.Join(" or ", PeptidesWithSetModifications.Select(b => b.Protein.Accession)) + "\t");
            sb.Append(string.Join(" or ", PeptidesWithSetModifications.Select(b => b.Protein.FullName)) + "\t");
            sb.Append(string.Join(" or ", PeptidesWithSetModifications.Select(b => b.PeptideDescription)) + "\t");
            sb.Append(string.Join(" or ", PeptidesWithSetModifications.Select(b => "[" + b.OneBasedStartResidueInProtein + " to " + b.OneBasedEndResidueInProtein + "]")) + "\t");
            sb.Append(string.Join(" or ", PeptidesWithSetModifications.Select(b => b.PreviousAminoAcid)) + "\t");
            sb.Append(string.Join(" or ", PeptidesWithSetModifications.Select(b => b.NextAminoAcid)) + "\t");
            sb.Append(BaseSequence.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(FullSequence.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(NumVariableMods.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(MissedCleavages.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(PeptideMonoisotopicMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append((ScanPrecursorMass - PeptideMonoisotopicMass).ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(((ScanPrecursorMass - PeptideMonoisotopicMass) / PeptideMonoisotopicMass * 1e6).ToString("F5", CultureInfo.InvariantCulture) + '\t');

            sb.Append(newPsm.precursorScanBestMass.Value.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append((newPsm.precursorScanBestMass - PeptideMonoisotopicMass).Value.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append(((newPsm.precursorScanBestMass - PeptideMonoisotopicMass) / PeptideMonoisotopicMass * 1e6).Value.ToString("F5", CultureInfo.InvariantCulture) + '\t');

            if (IsDecoy)
                sb.Append("D");
            else if (IsContaminant)
                sb.Append("C");
            else
                sb.Append("T");

            return sb.ToString();
        }

        #endregion Public Methods

    }
}