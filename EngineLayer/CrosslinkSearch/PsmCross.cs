using System.Linq;
using Proteomics;
using System.Collections.Generic;
using Chemistry;

namespace EngineLayer.CrosslinkSearch
{
    public class PsmCross : PsmParent
    {
        #region Private Fields

        private CompactPeptide compactPeptide;
        private PeptideWithSetModifications ps;
        //private ProductMassesMightHave pmmh;

        #endregion Private Fields

        #region Public Constructors

        public PsmCross(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(notch, score, scanIndex, scan)
        {
            compactPeptide = theBestPeptide;
        }

        public PsmCross(PeptideWithSetModifications ps, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(notch, score, scanIndex, scan)
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

            //if (containsAdot)
            //    throw new NotImplementedException();
            if (containsB)
                massLen += CompactPeptide.NTerminalMasses.Length - 1;
            if (containsC)
                massLen += CompactPeptide.NTerminalMasses.Length;
            //if (containsX)
            //    throw new NotImplementedException();
            if (containsY)
                massLen += CompactPeptide.CTerminalMasses.Length;
            if (containsZdot)
                massLen += CompactPeptide.CTerminalMasses.Length;


            ProductMassesMightHave indexesMassesToReturn = new ProductMassesMightHave(massLen);

            //int i = 0;
            //foreach (var hm in CompactPeptide.NTerminalMasses)
            //{
            //    if (hm.index > 1 && containsB)
            //    {
            //        indexesMassesToReturn.ProductMz[i] = hm.mass;
            //        indexesMassesToReturn.ProductName[i] = "b" + hm.index.ToString();
            //        i++;
            //    }
            //    if (containsC)
            //    {
            //        indexesMassesToReturn.ProductMz[i] = hm.mass + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
            //        indexesMassesToReturn.ProductName[i] = "c" + hm.index.ToString();
            //        i++;
            //    }
            //}
            //foreach (var hm in p.cTerminalMasses)
            //{
            //    if (containsY)
            //    {
            //        indexesMassesToReturn.ProductMz[i] = hm.mass + waterMonoisotopicMass;
            //        indexesMassesToReturn.ProductName[i] = "y" + (compactPeptide.CTerminalMasses.Length+1 - hm.index + 1).ToString();
            //        i++;
            //    }
            //    if (containsZdot)
            //    {
            //        indexesMassesToReturn.ProductMz[i] = hm.mass + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass;
            //        indexesMassesToReturn.ProductName[i] = "z" + (compactPeptide.CTerminalMasses.Length + 1 - hm.index + 1).ToString();
            //        i++;
            //    }
            //}
            return indexesMassesToReturn;
        }
        

        public class ProductMassesMightHave
        {
            public double[] ProductMz { get; set; }
            public string[] ProductName { get; set; }

            public ProductMassesMightHave(int length)
            {
                ProductMz = new double[length];
                ProductName = new string[length];
            }

            public ProductMassesMightHave()
            {
            }
        }

        #endregion

    }
}
