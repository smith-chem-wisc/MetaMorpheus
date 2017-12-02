using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.CrosslinkSearch
{
    public class PsmCross : Psm
    {
        #region Private Fields

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        //Compute ProductMassesMightHave: the theoritical masses of psmCross
        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;

        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        public CompactPeptide compactPeptide;

        #endregion Private Fields

        #region Public Constructors

        public PsmCross(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(theBestPeptide, notch, score, scanIndex, scan)
        {
            compactPeptide = theBestPeptide;
        }

        #endregion Public Constructors

        #region Public Properties

        //public CompactPeptide CompactPeptide { get { return compactPeptide; } set { compactPeptide = value; } }

        //public ProductMassesMightHave pmmh { get; set; }
        public double peptideMass { get; set; }

        public double XLBestScore { get; set; }
        public MatchedIonInfo matchedIonInfo { get; set; }
        public double XLTotalScore { get; set; }
        public double XLQvalueTotalScore { get; set; }
        public int XlPos { get; set; }
        public int[] XlRank { get; set; }
        public string ParentIonExist { get; set; }
        public int Charge2IonExist { get; set; }
        public PsmCross BetaPsmCross { get; set; }
        public PsmCrossType CrossType { get; set; }
        public double DScore { get; set; }

        #endregion Public Properties

        #region Public Methods

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
            Tolerance tolerance = new PpmTolerance(1);
            List<int> xlpos = new List<int>();
            if (tolerance.Within(compactPeptide.NTerminalMasses[0], Residue.GetResidue(crosslinker.CrosslinkerModSite).MonoisotopicMass))
            {
                xlpos.Add(0);
            }
            for (int i = 1; i < compactPeptide.NTerminalMasses.Length; i++)
            {
                if (tolerance.Within(compactPeptide.NTerminalMasses[i] - compactPeptide.NTerminalMasses[i - 1], Residue.GetResidue(crosslinker.CrosslinkerModSite).MonoisotopicMass))
                {
                    xlpos.Add(i);
                }
            }
            if (tolerance.Within(compactPeptide.CTerminalMasses[0], Residue.GetResidue(crosslinker.CrosslinkerModSite).MonoisotopicMass))
            {
                xlpos.Add(compactPeptide.NTerminalMasses.Length);
            }
            return xlpos;
        }

        public static void XLCalculateTotalProductMassesMightHave(Ms2ScanWithSpecificMass theScan, PsmCross psmCross, CrosslinkerTypeClass crosslinker, List<ProductType> lp, Tolerance fragmentTolerance)
        {
            bool CalculateHighCharge = false;
            var modMass = theScan.PrecursorMass - psmCross.compactPeptide.MonoisotopicMassIncludingFixedMods - crosslinker.TotalMass;
            int length = psmCross.compactPeptide.NTerminalMasses.Length;
            var pmmh = psmCross.ProductMassesMightHaveDuplicatesAndNaNs(lp);
            ProductMassesMightHave pmmhTop = new ProductMassesMightHave();

            List<ProductMassesMightHave> pmmhList = new List<ProductMassesMightHave>();

            var linkPos = PsmCross.xlPosCal(psmCross.compactPeptide, crosslinker);

            foreach (var ipos in linkPos)
            {
                ProductMassesMightHave pmmhCurr = new ProductMassesMightHave();
                pmmhCurr.XlPos = ipos;
                List<double> x = new List<double>();
                List<string> y = new List<string>();
                if (crosslinker.Cleavable)
                {
                    x.Add(theScan.PrecursorMass - modMass - crosslinker.CleaveMassLong);
                    y.Add("PepS");
                    x.Add(theScan.PrecursorMass - modMass - crosslinker.CleaveMassShort);
                    y.Add("PepL");
                }
                for (int i = 0; i < pmmh.ProductMz.Length; i++)
                {
                    var cr = pmmh.ProductName[i][0];
                    //get the position of amino acid
                    var nm = Int32.Parse(System.Text.RegularExpressions.Regex.Match(pmmh.ProductName[i], @"\d+").Value);
                    if ((cr == 'b' || cr == 'c') && nm < ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i]);
                        y.Add(pmmh.ProductName[i]);
                    }
                    if ((cr == 'y' || cr == 'z') && nm < length - ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i]);
                        y.Add(pmmh.ProductName[i]);
                    }
                    if (cr == 'b' && nm >= ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i] + modMass + crosslinker.TotalMass);
                        y.Add("t1b" + nm.ToString());
                        if (CalculateHighCharge)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 2);
                            y.Add("t2b" + nm.ToString());
                            //x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 3);
                            //y.Add("t3b" + nm.ToString());
                        }
                        if (crosslinker.Cleavable)
                        {
                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                            y.Add("sb" + nm.ToString());

                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                            y.Add("lb" + nm.ToString());
                        }
                    }

                    if (cr == 'c' && nm >= ipos)
                    {
                        x.Add(pmmh.ProductMz[i] + modMass + crosslinker.TotalMass);
                        y.Add("t1c" + nm.ToString());
                        if (CalculateHighCharge)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 2);
                            y.Add("t2c" + nm.ToString());
                            //x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 3);
                            //y.Add("t3c" + nm.ToString());
                        }

                        if (crosslinker.Cleavable)
                        {
                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                            y.Add("sc" + nm.ToString());

                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                            y.Add("lc" + nm.ToString());
                        }
                    }

                    if (cr == 'y' && (nm >= length - ipos + 1))
                    {
                        x.Add(pmmh.ProductMz[i] + modMass + crosslinker.TotalMass);
                        y.Add("t1y" + nm.ToString());

                        if (CalculateHighCharge)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 2);
                            y.Add("t2y" + nm.ToString());
                            //x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 3);
                            //y.Add("t3y" + nm.ToString());
                        }
                        if (crosslinker.Cleavable)
                        {
                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                            y.Add("sy" + nm.ToString());

                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                            y.Add("ly" + nm.ToString());
                        }
                    }

                    if (cr == 'z' && (nm >= length - ipos + 1))
                    {
                        x.Add(pmmh.ProductMz[i] + modMass + crosslinker.TotalMass);
                        y.Add("t1z" + nm.ToString());
                        if (CalculateHighCharge)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 2);
                            y.Add("t2z" + nm.ToString());
                            //x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 3);
                            //y.Add("t3z" + nm.ToString());
                        }
                        if (crosslinker.Cleavable)
                        {
                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                            y.Add("sz" + nm.ToString());

                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                            y.Add("lz" + nm.ToString());
                        }
                    }
                }
                pmmhCurr.ProductMz = x.ToArray();
                pmmhCurr.ProductName = y.ToArray();
                Array.Sort(pmmhCurr.ProductMz, pmmhCurr.ProductName);
                pmmhList.Add(pmmhCurr);
            }

            List<double> scoreList = new List<double>();
            List<MatchedIonInfo> miil = new List<MatchedIonInfo>();
            foreach (var pmm in pmmhList)
            {
                var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(pmm.ProductMz.Length);
                double pmmScore = PsmCross.XLMatchIons(theScan.TheScan, fragmentTolerance, pmm.ProductMz, pmm.ProductName, matchedIonMassesListPositiveIsMatch);
                miil.Add(matchedIonMassesListPositiveIsMatch);
                scoreList.Add(pmmScore);
            }

            psmCross.XLBestScore = scoreList.Max();
            psmCross.matchedIonInfo = miil[scoreList.IndexOf(scoreList.Max())];
            psmCross.XlPos = pmmhList[scoreList.IndexOf(scoreList.Max())].XlPos + 1;
            if (crosslinker.Cleavable)
            {
                if (psmCross.matchedIonInfo.MatchedIonName.Contains("PepS"))
                {
                    psmCross.ParentIonExist += "PepS";
                }
                if (psmCross.matchedIonInfo.MatchedIonName.Contains("PepL"))
                {
                    psmCross.ParentIonExist += "PepL";
                }
            }
            if (CalculateHighCharge)
            {
                int Charge2IonExist = 0;
                for (int i = 0; i < psmCross.matchedIonInfo.MatchedIonName.Length; i++)
                {
                    if (psmCross.matchedIonInfo.MatchedIonName[i] != null && (psmCross.matchedIonInfo.MatchedIonName[i].Contains("t2") || psmCross.matchedIonInfo.MatchedIonName[i].Contains("t3")))
                    {
                        Charge2IonExist++;
                    }
                }
                psmCross.Charge2IonExist = Charge2IonExist;
            }
        }

        public static void XLCalculateTotalProductMassesMightHaveDeadend(Ms2ScanWithSpecificMass theScan, PsmCross psmCross, CrosslinkerTypeClass crosslinker, List<ProductType> lp, Tolerance fragmentTolerance, Double modMass)
        {
            bool CalculateHighCharge = false;
            int length = psmCross.compactPeptide.NTerminalMasses.Length;
            var pmmh = psmCross.ProductMassesMightHaveDuplicatesAndNaNs(lp);
            ProductMassesMightHave pmmhTop = new ProductMassesMightHave();

            List<ProductMassesMightHave> pmmhList = new List<ProductMassesMightHave>();

            var linkPos = PsmCross.xlPosCal(psmCross.compactPeptide, crosslinker);

            foreach (var ipos in linkPos)
            {
                ProductMassesMightHave pmmhCurr = new ProductMassesMightHave();
                pmmhCurr.XlPos = ipos;
                List<double> x = new List<double>();
                List<string> y = new List<string>();
                if (crosslinker.Cleavable)
                {
                    x.Add(psmCross.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.CleaveMassLong);
                    y.Add("PepL");
                    x.Add(psmCross.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.CleaveMassShort);
                    y.Add("PepS");
                }
                for (int i = 0; i < pmmh.ProductMz.Length; i++)
                {
                    var cr = pmmh.ProductName[i][0];
                    //get the position of amino acid
                    var nm = Int32.Parse(System.Text.RegularExpressions.Regex.Match(pmmh.ProductName[i], @"\d+").Value);
                    if ((cr == 'b' || cr == 'c') && nm < ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i]);
                        y.Add(pmmh.ProductName[i]);
                    }
                    if ((cr == 'y' || cr == 'z') && nm < length - ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i]);
                        y.Add(pmmh.ProductName[i]);
                    }
                    if (cr == 'b' && nm >= ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i] + modMass);
                        y.Add("t1b" + nm.ToString());
                        if (CalculateHighCharge)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass) / 2);
                            y.Add("t2b" + nm.ToString());
                            //x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 3);
                            //y.Add("t3b" + nm.ToString());
                        }
                        if (crosslinker.Cleavable)
                        {
                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                            y.Add("sb" + nm.ToString());

                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                            y.Add("lb" + nm.ToString());
                        }
                    }

                    if (cr == 'c' && nm >= ipos)
                    {
                        x.Add(pmmh.ProductMz[i] + modMass);
                        y.Add("t1c" + nm.ToString());
                        if (CalculateHighCharge)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass) / 2);
                            y.Add("t2c" + nm.ToString());
                            //x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 3);
                            //y.Add("t3c" + nm.ToString());
                        }

                        if (crosslinker.Cleavable)
                        {
                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                            y.Add("sc" + nm.ToString());

                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                            y.Add("lc" + nm.ToString());
                        }
                    }

                    if (cr == 'y' && (nm >= length - ipos + 1))
                    {
                        x.Add(pmmh.ProductMz[i] + modMass);
                        y.Add("t1y" + nm.ToString());

                        if (CalculateHighCharge)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass) / 2);
                            y.Add("t2y" + nm.ToString());
                            //x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 3);
                            //y.Add("t3y" + nm.ToString());
                        }
                        if (crosslinker.Cleavable)
                        {
                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                            y.Add("sy" + nm.ToString());

                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                            y.Add("ly" + nm.ToString());
                        }
                    }

                    if (cr == 'z' && (nm >= length - ipos + 1))
                    {
                        x.Add(pmmh.ProductMz[i] + modMass);
                        y.Add("t1z" + nm.ToString());
                        if (CalculateHighCharge)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass) / 2);
                            y.Add("t2z" + nm.ToString());
                            //x.Add((pmmh.ProductMz[i] + modMass + crosslinker.TotalMass) / 3);
                            //y.Add("t3z" + nm.ToString());
                        }
                        if (crosslinker.Cleavable)
                        {
                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                            y.Add("sz" + nm.ToString());

                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                            y.Add("lz" + nm.ToString());
                        }
                    }
                }
                pmmhCurr.ProductMz = x.ToArray();
                pmmhCurr.ProductName = y.ToArray();
                Array.Sort(pmmhCurr.ProductMz, pmmhCurr.ProductName);
                pmmhList.Add(pmmhCurr);
            }

            List<double> scoreList = new List<double>();
            List<MatchedIonInfo> miil = new List<MatchedIonInfo>();
            foreach (var pmm in pmmhList)
            {
                var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(pmm.ProductMz.Length);
                double pmmScore = PsmCross.XLMatchIons(theScan.TheScan, fragmentTolerance, pmm.ProductMz, pmm.ProductName, matchedIonMassesListPositiveIsMatch);
                miil.Add(matchedIonMassesListPositiveIsMatch);
                scoreList.Add(pmmScore);
            }

            psmCross.XLBestScore = scoreList.Max();
            psmCross.matchedIonInfo = miil[scoreList.IndexOf(scoreList.Max())];
            psmCross.XlPos = pmmhList[scoreList.IndexOf(scoreList.Max())].XlPos + 1;
            if (crosslinker.Cleavable)
            {
                if (psmCross.matchedIonInfo.MatchedIonName.Contains("PepS"))
                {
                    psmCross.ParentIonExist += "PepS";
                }
                if (psmCross.matchedIonInfo.MatchedIonName.Contains("PepL"))
                {
                    psmCross.ParentIonExist += "PepL";
                }
            }
            if (CalculateHighCharge)
            {
                int Charge2IonExist = 0;
                for (int i = 0; i < psmCross.matchedIonInfo.MatchedIonName.Length; i++)
                {
                    if (psmCross.matchedIonInfo.MatchedIonName[i] != null && (psmCross.matchedIonInfo.MatchedIonName[i].Contains("t2") || psmCross.matchedIonInfo.MatchedIonName[i].Contains("t3")))
                    {
                        Charge2IonExist++;
                    }
                }
                psmCross.Charge2IonExist = Charge2IonExist;
            }
        }

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
                massLen += compactPeptide.NTerminalMasses.Length - 1;
            if (containsC)
                massLen += compactPeptide.NTerminalMasses.Length;
            if (containsX)
                throw new NotImplementedException();
            if (containsY)
                massLen += compactPeptide.CTerminalMasses.Length;
            if (containsZdot)
                massLen += compactPeptide.CTerminalMasses.Length;

            ProductMassesMightHave productMassMightHave = new ProductMassesMightHave(massLen);
            int i = 0;
            int ib = 0;
            int ic = 0;
            for (int j = 0; j < compactPeptide.NTerminalMasses.Length; j++)
            {
                var hm = compactPeptide.NTerminalMasses[j];
                if (containsB)
                {
                    if (j > 0)
                    {
                        productMassMightHave.ProductMz[i] = hm;
                        productMassMightHave.ProductName[i] = "b" + (ib + 2).ToString();
                        i++;
                        ib++;
                    }
                }
                if (containsC)
                {
                    productMassMightHave.ProductMz[i] = hm + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "c" + (ic + 1).ToString();
                    i++;
                    ic++;
                }
            }
            int iy = compactPeptide.CTerminalMasses.Length - 1;
            int iz = compactPeptide.CTerminalMasses.Length - 1;
            for (int j = 0; j < compactPeptide.CTerminalMasses.Length; j++)
            {
                var hm = compactPeptide.CTerminalMasses[j];
                if (containsY)
                {
                    productMassMightHave.ProductMz[i] = hm + waterMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "y" + (compactPeptide.CTerminalMasses.Length - iy).ToString();
                    i++;
                    iy--;
                }
                if (containsZdot)
                {
                    productMassMightHave.ProductMz[i] = hm + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "z" + (compactPeptide.CTerminalMasses.Length - iz).ToString();
                    i++;
                    iz--;
                }
            }
            return productMassMightHave;
        }

        #endregion Public Methods
    }
}