using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.CrosslinkSearch
{
    public class PsmCross : PeptideSpectralMatch
    {
        public CompactPeptide CompactPeptide;

        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double NitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double OxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double HydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;

        public PsmCross(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, DigestionParams digestionParams)
            : base(theBestPeptide, notch, score, scanIndex, scan, digestionParams)
        {
            CompactPeptide = theBestPeptide;
        }

        public double XLBestScore { get; set; }
        public MatchedIonInfo MatchedIonInfo { get; set; }
        public double XLTotalScore { get; set; }
        public double XLQvalueTotalScore { get; set; }
        public int XlPos { get; set; }
        public int XlPos2 { get; set; }
        public int XlProteinPos { get; set; }
        public int[] XlRank { get; set; }
        public string ParentIonExist { get; set; }
        public int ParentIonExistNum { get; set; }
        public List<int> ParentIonMaxIntensityRanks { get; set; }
        public int Charge2IonExist { get; set; }
        public PsmCross BetaPsmCross { get; set; }
        public PsmCrossType CrossType { get; set; }
        public double DScore { get; set; }

        //Calculate score based on Product Masses.
        public static double XlMatchIons(MsDataScan thisScan, Tolerance productMassTolerance, double[] sorted_theoretical_product_masses_for_this_peptide, string[] sorted_theoretical_product_name_for_this_peptide, MatchedIonInfo matchedIonMassesListPositiveIsMatch)
        {
            var TotalProductsHere = sorted_theoretical_product_masses_for_this_peptide.Length;
            if (TotalProductsHere == 0)
                return 0;
            int MatchingProductsHere = 0;
            double MatchingIntensityHere = 0;

            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.XArray;
            double[] experimental_intensities = thisScan.MassSpectrum.YArray;
            int[] experimental_intensities_rank = GenerateIntensityRanks(experimental_intensities);
            int num_experimental_peaks = experimental_mzs.Length;

            int currentTheoreticalIndex = -1;
            double currentTheoreticalMass;
            do
            {
                currentTheoreticalIndex++;
                currentTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[currentTheoreticalIndex];
            } while (double.IsNaN(currentTheoreticalMass) && currentTheoreticalIndex < sorted_theoretical_product_masses_for_this_peptide.Length - 1);

            if (double.IsNaN(currentTheoreticalMass))
            {
                return 0;
            }

            double currentTheoreticalMz = currentTheoreticalMass + Constants.ProtonMass;

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
                    matchedIonMassesListPositiveIsMatch.MatchedIonIntensityRank[currentTheoreticalIndex] = experimental_intensities_rank[experimentalIndex];
                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                    {
                        break;
                    }
                    currentTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.ProtonMass;
                }
                // Else if for sure did not reach the next theoretical yet, move to next experimental
                else if (currentExperimentalMZ < currentTheoreticalMz)
                {
                    continue;
                }
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
                    currentTheoreticalMz = currentTheoreticalMass + Constants.ProtonMass;

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
                        {
                            break;
                        }
                        testTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[testTheoreticalIndex];
                        testTheoreticalMZ = testTheoreticalMass + Constants.ProtonMass;
                    }

                    experimentalIndex--;
                }
            }
            return MatchingProductsHere + MatchingIntensityHere / thisScan.TotalIonCurrent;
        }

        //Calculate if crosslink amino acid exist and return its position based on compactPeptide
        public static List<int> XlPosCal(CompactPeptide compactPeptide, string crosslinkerModSites)
        {
            Tolerance tolerance = new PpmTolerance(1);
            List<int> xlpos = new List<int>();
            foreach (char item in crosslinkerModSites)
            {
                if (tolerance.Within(compactPeptide.NTerminalMasses[0], Residue.GetResidue(item).MonoisotopicMass))
                {
                    xlpos.Add(0);
                }
                for (int i = 1; i < compactPeptide.NTerminalMasses.Length; i++)
                {
                    if (tolerance.Within(compactPeptide.NTerminalMasses[i] - compactPeptide.NTerminalMasses[i - 1], Residue.GetResidue(item).MonoisotopicMass))
                    {
                        xlpos.Add(i);
                    }
                }
                if (tolerance.Within(compactPeptide.CTerminalMasses[0], Residue.GetResidue(item).MonoisotopicMass))
                {
                    xlpos.Add(compactPeptide.NTerminalMasses.Length);
                }
            }
            xlpos.Sort();
            return xlpos;
        }

        //Calculate All possible Products Masses based on ModMass and linkPos
        public static List<ProductMassesMightHave> XlCalculateTotalProductMasses(PsmCross psmCross, double modMass,
            CrosslinkerTypeClass crosslinker, List<ProductType> lp, bool Charge_2_3, bool Charge_2_3_PrimeFragment, List<int> linkPos)
        {
            int length = psmCross.CompactPeptide.NTerminalMasses.Length;
            var pmmh = psmCross.ProductMassesMightHaveDuplicatesAndNaNs(lp);
            ProductMassesMightHave pmmhTop = new ProductMassesMightHave();

            List<ProductMassesMightHave> pmmhList = new List<ProductMassesMightHave>();

            foreach (var ipos in linkPos)
            {
                var pmmhCurr = new ProductMassesMightHave();
                pmmhCurr.XlPos = ipos;
                List<double> x = new List<double>();
                List<string> y = new List<string>();
                if (crosslinker.Cleavable)
                {
                    x.Add((double)psmCross.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.CleaveMassShort);
                    y.Add("PepS");
                    x.Add(((double)psmCross.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.CleaveMassShort) / 2);
                    y.Add("PepS2");
                    x.Add((double)psmCross.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.CleaveMassLong);
                    y.Add("PepL");
                    x.Add(((double)psmCross.CompactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.CleaveMassLong) / 2);
                    y.Add("PepL2");
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
                        if (Charge_2_3_PrimeFragment && cr == 'b')
                        {
                            x.Add(pmmh.ProductMz[i] / 2);
                            y.Add("t2b" + nm.ToString());
                            x.Add(pmmh.ProductMz[i] / 3);
                            y.Add("t3b" + nm.ToString());
                        }
                        if (Charge_2_3_PrimeFragment && cr == 'c')
                        {
                            x.Add(pmmh.ProductMz[i] / 2);
                            y.Add("t2c" + nm.ToString());
                            x.Add(pmmh.ProductMz[i] / 3);
                            y.Add("t3c" + nm.ToString());
                        }
                    }
                    if ((cr == 'y' || cr == 'z') && nm < length - ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i]);
                        y.Add(pmmh.ProductName[i]);
                        if (Charge_2_3_PrimeFragment && cr == 'y')
                        {
                            x.Add(pmmh.ProductMz[i] / 2);
                            y.Add("t2y" + nm.ToString());
                            x.Add(pmmh.ProductMz[i] / 3);
                            y.Add("t3y" + nm.ToString());
                        }
                        if (Charge_2_3_PrimeFragment && cr == 'z')
                        {
                            x.Add(pmmh.ProductMz[i] / 2);
                            y.Add("t2z" + nm.ToString());
                            x.Add(pmmh.ProductMz[i] / 3);
                            y.Add("t3z" + nm.ToString());
                        }
                    }
                    if (cr == 'b' && nm >= ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i] + modMass);
                        y.Add("t1b" + nm.ToString());
                        if (Charge_2_3)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass) / 2);
                            y.Add("t2b" + nm.ToString());
                            x.Add((pmmh.ProductMz[i] + modMass) / 3);
                            y.Add("t3b" + nm.ToString());
                        }
                        if (crosslinker.Cleavable)
                        {
                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassShort);
                            y.Add("sb" + nm.ToString());

                            x.Add(pmmh.ProductMz[i] + crosslinker.CleaveMassLong);
                            y.Add("lb" + nm.ToString());
                        }
                    }

                    if (cr == 'c' && nm >= ipos + 1)
                    {
                        x.Add(pmmh.ProductMz[i] + modMass);
                        y.Add("t1c" + nm.ToString());
                        if (Charge_2_3)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass) / 2);
                            y.Add("t2c" + nm.ToString());
                            x.Add((pmmh.ProductMz[i] + modMass) / 3);
                            y.Add("t3c" + nm.ToString());
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

                        if (Charge_2_3)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass) / 2);
                            y.Add("t2y" + nm.ToString());
                            x.Add((pmmh.ProductMz[i] + modMass) / 3);
                            y.Add("t3y" + nm.ToString());
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
                        if (Charge_2_3)
                        {
                            x.Add((pmmh.ProductMz[i] + modMass) / 2);
                            y.Add("t2z" + nm.ToString());
                            x.Add((pmmh.ProductMz[i] + modMass) / 3);
                            y.Add("t3z" + nm.ToString());
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

            return pmmhList;
        }

        //Calculate score based on All possible Products Masses for inter- or intra- crosslinks and deadend.
        public static void XlLocalization(Ms2ScanWithSpecificMass theScan, PsmCross psmCross, double modMass,
            CrosslinkerTypeClass crosslinker, List<ProductType> lp, Tolerance fragmentTolerance, bool Charge_2_3, bool Charge_2_3_PrimeFragment, List<int> linkPos)
        {
            var pmmhList = PsmCross.XlCalculateTotalProductMasses(psmCross, modMass, crosslinker, lp, Charge_2_3, Charge_2_3_PrimeFragment, linkPos);

            List<double> scoreList = new List<double>();
            List<MatchedIonInfo> miil = new List<MatchedIonInfo>();
            foreach (var pmm in pmmhList)
            {
                var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(pmm.ProductMz.Length);
                double pmmScore = PsmCross.XlMatchIons(theScan.TheScan, fragmentTolerance, pmm.ProductMz, pmm.ProductName, matchedIonMassesListPositiveIsMatch);
                miil.Add(matchedIonMassesListPositiveIsMatch);
                scoreList.Add(pmmScore);
            }

            psmCross.XLBestScore = scoreList.Max();
            psmCross.MatchedIonInfo = miil[scoreList.IndexOf(scoreList.Max())];
            psmCross.XlPos = pmmhList[scoreList.IndexOf(scoreList.Max())].XlPos + 1;
            if (crosslinker.Cleavable)
            {
                psmCross.ParentIonMaxIntensityRanks = new List<int>();
                if (psmCross.MatchedIonInfo.MatchedIonName.Any(p => p != null && p.Contains("PepS")))
                {
                    psmCross.ParentIonExist += "PepS";
                    psmCross.ParentIonExistNum += 1;
                }
                if (psmCross.MatchedIonInfo.MatchedIonName.Any(p => p != null && p.Contains("PepL")))
                {
                    psmCross.ParentIonExist += "PepL";
                    psmCross.ParentIonExistNum += 1;
                }
                //if (psmCross.MatchedIonInfo.MatchedIonName.Any(p => p != null && p.Equals("PepS2")))
                //{
                //    psmCross.ParentIonExist += "PepS2";
                //    psmCross.ParentIonExistNum += 1;
                //}
                //if (psmCross.MatchedIonInfo.MatchedIonName.Any(p => p != null && p.Equals("PepL2")))
                //{
                //    psmCross.ParentIonExist += "PepL2";
                //    psmCross.ParentIonExistNum += 1;
                //}
                for (int i = 0; i < psmCross.MatchedIonInfo.MatchedIonName.Length; i++)
                {
                    if (psmCross.MatchedIonInfo.MatchedIonName[i] != null)
                    {
                        if (psmCross.MatchedIonInfo.MatchedIonName[i].Contains("Pep"))
                        {
                            psmCross.ParentIonMaxIntensityRanks.Add(psmCross.MatchedIonInfo.MatchedIonIntensityRank[i]);
                        }
                    }
                }
            }
            if (Charge_2_3 || Charge_2_3_PrimeFragment)
            {
                int Charge2IonExist = 0;
                for (int i = 0; i < psmCross.MatchedIonInfo.MatchedIonName.Length; i++)
                {
                    if (psmCross.MatchedIonInfo.MatchedIonName[i] != null && (psmCross.MatchedIonInfo.MatchedIonName[i].Contains("t2") || psmCross.MatchedIonInfo.MatchedIonName[i].Contains("t3")))
                    {
                        Charge2IonExist++;
                    }
                }
                psmCross.Charge2IonExist = Charge2IonExist;
            }
        }

        //Calculate All possible Products Masses based on ModMass and linkPos
        public static List<ProductMassesMightHave> XlCalculateTotalProductMassesForLoopCrosslink(PsmCross psmCross, double modMass, CrosslinkerTypeClass crosslinker, List<ProductType> lp, List<int> linkPos)
        {
            int length = psmCross.CompactPeptide.NTerminalMasses.Length;
            var pmmh = psmCross.ProductMassesMightHaveDuplicatesAndNaNs(lp);
            ProductMassesMightHave pmmhTop = new ProductMassesMightHave();

            List<ProductMassesMightHave> pmmhList = new List<ProductMassesMightHave>();

            if (linkPos.Count >= 2)
            {
                for (int ipos = 0; ipos < linkPos.Count - 1; ipos++)
                {
                    for (int jpos = ipos + 1; jpos < linkPos.Count; jpos++)
                    {
                        var pmmhCurr = new ProductMassesMightHave();
                        pmmhCurr.XlPos = linkPos[ipos];
                        pmmhCurr.XlPos2 = linkPos[jpos];
                        List<double> x = new List<double>();
                        List<string> y = new List<string>();

                        for (int i = 0; i < pmmh.ProductMz.Length; i++)
                        {
                            var cr = pmmh.ProductName[i][0];
                            //get the position of amino acid
                            var nm = Int32.Parse(System.Text.RegularExpressions.Regex.Match(pmmh.ProductName[i], @"\d+").Value);
                            if ((cr == 'b' || cr == 'c') && nm < linkPos[ipos] + 1)
                            {
                                x.Add(pmmh.ProductMz[i]);
                                y.Add(pmmh.ProductName[i]);
                            }
                            if ((cr == 'y' || cr == 'z') && nm < length - linkPos[jpos] + 1)
                            {
                                x.Add(pmmh.ProductMz[i]);
                                y.Add(pmmh.ProductName[i]);
                            }
                            if (cr == 'b' && nm >= linkPos[jpos] + 1)
                            {
                                x.Add(pmmh.ProductMz[i] + modMass);
                                y.Add("t1b" + nm.ToString());
                            }

                            if (cr == 'c' && nm >= linkPos[jpos] + 1)
                            {
                                x.Add(pmmh.ProductMz[i] + modMass);
                                y.Add("t1c" + nm.ToString());
                            }

                            if (cr == 'y' && (nm >= length - linkPos[ipos] + 1))
                            {
                                x.Add(pmmh.ProductMz[i] + modMass);
                                y.Add("t1y" + nm.ToString());
                            }

                            if (cr == 'z' && (nm >= length - linkPos[ipos] + 1))
                            {
                                x.Add(pmmh.ProductMz[i] + modMass);
                                y.Add("t1z" + nm.ToString());
                            }
                        }
                        pmmhCurr.ProductMz = x.ToArray();
                        pmmhCurr.ProductName = y.ToArray();
                        Array.Sort(pmmhCurr.ProductMz, pmmhCurr.ProductName);
                        pmmhList.Add(pmmhCurr);
                    }
                }
            }
            return pmmhList;
        }

        //Calculate score based on All possible Products Masses for inter- or intra- crosslinks and deadend.
        public static void XlLocalizationForLoopCrosslink(Ms2ScanWithSpecificMass theScan, PsmCross psmCross, double modMass, CrosslinkerTypeClass crosslinker, List<ProductType> lp, Tolerance fragmentTolerance, List<int> linkPos)
        {
            var pmmhList = PsmCross.XlCalculateTotalProductMassesForLoopCrosslink(psmCross, modMass, crosslinker, lp, linkPos);

            List<double> scoreList = new List<double>();
            List<MatchedIonInfo> miil = new List<MatchedIonInfo>();
            foreach (var pmm in pmmhList)
            {
                var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(pmm.ProductMz.Length);
                double pmmScore = PsmCross.XlMatchIons(theScan.TheScan, fragmentTolerance, pmm.ProductMz, pmm.ProductName, matchedIonMassesListPositiveIsMatch);
                miil.Add(matchedIonMassesListPositiveIsMatch);
                scoreList.Add(pmmScore);
            }

            psmCross.XLBestScore = scoreList.Max();
            psmCross.MatchedIonInfo = miil[scoreList.IndexOf(scoreList.Max())];
            psmCross.XlPos = pmmhList[scoreList.IndexOf(scoreList.Max())].XlPos + 1;
            psmCross.XlPos2 = pmmhList[scoreList.IndexOf(scoreList.Max())].XlPos2 + 1;
        }

        public static int[] GenerateIntensityRanks(double[] experimental_intensities)
        {
            var y = experimental_intensities.ToArray();
            var x = Enumerable.Range(1, y.Length).OrderBy(p => p).ToArray();
            Array.Sort(y, x);
            var experimental_intensities_rank = Enumerable.Range(1, y.Length).OrderByDescending(p => p).ToArray();
            Array.Sort(x, experimental_intensities_rank);
            return experimental_intensities_rank;
        }

        public ProductMassesMightHave ProductMassesMightHaveDuplicatesAndNaNs(List<ProductType> productTypes)
        {
            int massLen = 0;
            bool containsAdot = productTypes.Contains(ProductType.Adot);
            bool containsB = productTypes.Contains(ProductType.B);
            bool containsBnoB1 = productTypes.Contains(ProductType.BnoB1ions);
            bool containsC = productTypes.Contains(ProductType.C);
            bool containsX = productTypes.Contains(ProductType.X);
            bool containsY = productTypes.Contains(ProductType.Y);
            bool containsZdot = productTypes.Contains(ProductType.Zdot);

            if (containsAdot)
                throw new NotImplementedException();
            if (containsBnoB1)
                massLen += CompactPeptide.NTerminalMasses.Length - 1;
            if (containsB)
                massLen += CompactPeptide.NTerminalMasses.Length;
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
            int ib = 0;
            int ic = 0;
            for (int j = 0; j < CompactPeptide.NTerminalMasses.Length; j++)
            {
                var hm = CompactPeptide.NTerminalMasses[j];
                if (containsBnoB1)
                {
                    if (j > 0)
                    {
                        productMassMightHave.ProductMz[i] = hm;
                        productMassMightHave.ProductName[i] = "b" + (ib + 2).ToString();
                        i++;
                        ib++;
                    }
                }
                if (containsB)
                {
                    productMassMightHave.ProductMz[i] = hm;
                    productMassMightHave.ProductName[i] = "b" + (ib + 1).ToString();
                    i++;
                    ib++;
                }
                if (containsC)
                {
                    productMassMightHave.ProductMz[i] = hm + NitrogenAtomMonoisotopicMass + 3 * HydrogenAtomMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "c" + (ic + 1).ToString();
                    i++;
                    ic++;
                }
            }
            int iy = CompactPeptide.CTerminalMasses.Length - 1;
            int iz = CompactPeptide.CTerminalMasses.Length - 1;
            for (int j = 0; j < CompactPeptide.CTerminalMasses.Length; j++)
            {
                var hm = CompactPeptide.CTerminalMasses[j];
                if (containsY)
                {
                    productMassMightHave.ProductMz[i] = hm + WaterMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "y" + (CompactPeptide.CTerminalMasses.Length - iy).ToString();
                    i++;
                    iy--;
                }
                if (containsZdot)
                {
                    productMassMightHave.ProductMz[i] = hm + OxygenAtomMonoisotopicMass - NitrogenAtomMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "z" + (CompactPeptide.CTerminalMasses.Length - iz).ToString();
                    i++;
                    iz--;
                }
            }
            return productMassMightHave;
        }
    }
}