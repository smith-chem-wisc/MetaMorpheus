using System.Collections.Generic;
using Chemistry;
using System;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System.Linq;

namespace EngineLayer.CrosslinkSearch
{
    public class PsmCross : Psm
    {
        #region Private Fields

        private CompactPeptide compactPeptide;

        #endregion Private Fields

        #region Public Constructors

        public PsmCross(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(theBestPeptide, notch, score, scanIndex, scan)
        {
            compactPeptide = theBestPeptide;
        }

        #endregion Public Constructors

        public CompactPeptide CompactPeptide { get { return compactPeptide; } set { compactPeptide = value; } }

        //public ProductMassesMightHave pmmh { get; set; }
        public double peptideMass { get; set; }
        public double XLBestScore { get; set; }
        public MatchedIonInfo matchedIonInfo { get; set; }
        public double NScore { get; set; }
        public double XLTotalScore { get; set; }
        

        //Compute ProductMassesMightHave: the theoritical masses of psmCross
        #region ProductMassesMightHave

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;     

        public ProductMassesMightHave ProductMassesMightHaveDuplicatesAndNaNs(List<ProductType> productTypes)
        {
            int massLen = 0;
            bool containsAdot = productTypes.Contains(ProductType.Adot);
            bool containsB = productTypes.Contains(ProductType.B);
            bool containsC = productTypes.Contains(ProductType.C);
            bool containsX = productTypes.Contains(ProductType.X);
            bool containsY = productTypes.Contains(ProductType.Y);
            bool containsZdot = productTypes.Contains(ProductType.Zdot);

            if (containsAdot)
                throw new NotImplementedException();
            if (containsB)
                massLen += CompactPeptide.NTerminalMasses.Length - 1;
            if (containsC)
                massLen += CompactPeptide.NTerminalMasses.Length;
            if (containsX)
                throw new NotImplementedException();
            if (containsY)
                massLen += CompactPeptide.CTerminalMasses.Length;
            if (containsZdot)
                massLen += CompactPeptide.CTerminalMasses.Length;

            ProductMassesMightHave productMassMightHave = new ProductMassesMightHave(massLen);
            int i = 0;
            for (int j = 0; j < compactPeptide.NTerminalMasses.Length; j++)
            {
                var hm = compactPeptide.NTerminalMasses[j];
                if (containsB)
                {
                    if (j > 0)
                    {
                        productMassMightHave.ProductMz[i] = hm;
                        productMassMightHave.ProductName[i] = "b" + i.ToString();
                        i++;
                    }
                }
                if (containsC)
                {
                    productMassMightHave.ProductMz[i] = hm + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "c" + i.ToString();
                    i++;
                }
            }
            for (int j = 0; j < compactPeptide.CTerminalMasses.Length; j++)
            {
                var hm = compactPeptide.CTerminalMasses[j];
                if (containsY)
                {
                    productMassMightHave.ProductMz[i] = hm + waterMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "y" + (compactPeptide.CTerminalMasses.Length + 2 - i).ToString();
                    i++;
                }
                if (containsZdot)
                {
                    productMassMightHave.ProductMz[i] = hm + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "z" + (compactPeptide.CTerminalMasses.Length + 2 - i).ToString();
                    i++;
                }
            }
            return productMassMightHave;
        }

        //Compute matched ions
        public static double XLMatchIons(IMsDataScan<IMzSpectrum<IMzPeak>> thisScan, Tolerance productMassTolerance, double[] sorted_theoretical_product_masses_for_this_peptide, string[] sorted_theoretical_product_name_for_this_peptide, MatchedIonInfo matchedIonMassesListPositiveIsMatch)
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
                    matchedIonMassesListPositiveIsMatch.MatchedIonMz[currentTheoreticalIndex] = currentTheoreticalMass;
                    matchedIonMassesListPositiveIsMatch.MatchedIonIntensity[currentTheoreticalIndex] = experimental_intensities[experimentalIndex];
                    matchedIonMassesListPositiveIsMatch.MatchedIonName[currentTheoreticalIndex] = sorted_theoretical_product_name_for_this_peptide[currentTheoreticalIndex];
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
                    matchedIonMassesListPositiveIsMatch.MatchedIonMz[currentTheoreticalIndex] = -currentTheoreticalMass;

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
                        matchedIonMassesListPositiveIsMatch.MatchedIonMz[testTheoreticalIndex] = -currentTheoreticalMass;
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

        //Compute if crosslink amino acid exist and return its position based on compactPeptide
        public static List<int> xlPosCal(CompactPeptide compactPeptide, CrosslinkerTypeClass crosslinker)
        {
            Tolerance tolerance = new PpmTolerance(0.001);
            List<int> xlpos = new List<int>();
            if (tolerance.Within( compactPeptide.NTerminalMasses[0] , Residue.GetResidue(crosslinker.CrosslinkerModSite).MonoisotopicMass))
            {
                xlpos.Add(0);
            }
            for (int i = 1; i < compactPeptide.NTerminalMasses.Length; i++)
            {

                if (tolerance.Within(compactPeptide.NTerminalMasses[i] - compactPeptide.NTerminalMasses[i-1] , Residue.GetResidue(crosslinker.CrosslinkerModSite).MonoisotopicMass))
                {
                    xlpos.Add(i);
                }
            }
            if (tolerance.Within(compactPeptide.CTerminalMasses[0] , Residue.GetResidue(crosslinker.CrosslinkerModSite).MonoisotopicMass))
            {
                xlpos.Add(compactPeptide.NTerminalMasses.Length);
            }
            return xlpos;
        }
        #endregion

    }
}
