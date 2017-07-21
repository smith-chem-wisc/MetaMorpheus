using Proteomics;
using System.Collections.Generic;
using Chemistry;
using System;

namespace EngineLayer.CrosslinkSearch
{
    public class PsmCross : Psm
    {
        #region Private Fields

        private CompactPeptide compactPeptide;
        private PeptideWithSetModifications ps;

        #endregion Private Fields

        #region Public Constructors

        public PsmCross(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(theBestPeptide, notch, score, scanIndex, scan)
        {
            compactPeptide = theBestPeptide;
        }

        public PsmCross(PeptideWithSetModifications ps, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(ps.CompactPeptide, notch, score, scanIndex, scan)
        {
            this.ps = ps;
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
        public Dictionary<int, ModificationWithMass> allModsOneIsNterminus;      

        public ProductMassesMightHave ProductMassesMightHaveDuplicatesAndNaNs(List<ProductType> productTypes)
        {
            //PeptideFragmentMasses p = new PeptideFragmentMasses();
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

        #endregion

    }
}
