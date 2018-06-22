using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using Chemistry;


namespace MetaDrawGUI
{
    //public class PsmDraw : PeptideSpectralMatch
    public class PsmDraw
    {
        private CompactPeptide compactPeptide;

        //public PsmDraw(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, DigestionParams digestionParams) : base(theBestPeptide, notch, score, scanIndex, scan, digestionParams)
        //{
        //    compactPeptide = theBestPeptide;
        //}
        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;

        public PsmDraw() { }

        public MatchedIonInfo MatchedIonInfo { get; set; }
        public  int ScanNumber { get; set; }
        public  string BaseSequence { get; set; }
        public  string FullSequence { get; set; }
        public  bool IsDecoy { get; set; }
        public  double? PeptideMonisotopicMass { get; set; }
        public  int ScanPrecursorCharge { get; set; }
        public CompactPeptide CompactPeptide { get { return compactPeptide; } set { this.compactPeptide = value; } }
        public CompactPeptideDraw CompactPeptideDraw { get; set;}

        //Calculate All possible Products Masses based on ModMass and linkPos
        public static ProductMassesMightHave XlCalculateTotalProductMassesForSingle(PsmDraw psmCross, List<ProductType> lp, bool Charge_2_3_PrimeFragment)
        {
            int length = psmCross.CompactPeptideDraw.NTerminalMasses.Length;
            var pmmh = psmCross.ProductMassesMightHaveDuplicatesAndNaNs(lp);

            var pmmhCurr = new ProductMassesMightHave();

            List<double> x = new List<double>();
            List<string> y = new List<string>();

            for (int i = 0; i < pmmh.ProductMz.Length; i++)
            {
                var cr = pmmh.ProductName[i][0];
                //get the position of amino acid
                var nm = Int32.Parse(System.Text.RegularExpressions.Regex.Match(pmmh.ProductName[i], @"\d+").Value);
                if ((cr == 'b' || cr == 'c') && nm <= length + 1)
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
                if ((cr == 'y' || cr == 'z') && nm <= length + 1)
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
            }
            pmmhCurr.ProductMz = x.ToArray();
            pmmhCurr.ProductName = y.ToArray();

            Array.Sort(pmmhCurr.ProductMz, pmmhCurr.ProductName);


            return pmmhCurr;
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
                massLen += CompactPeptideDraw.NTerminalMasses.Length - 1;
            if (containsB)
                massLen += CompactPeptideDraw.NTerminalMasses.Length;
            if (containsC)
                massLen += CompactPeptideDraw.NTerminalMasses.Length;
            if (containsX)
                throw new NotImplementedException();
            if (containsY)
                massLen += CompactPeptideDraw.CTerminalMasses.Length;
            if (containsZdot)
                massLen += CompactPeptideDraw.CTerminalMasses.Length;

            ProductMassesMightHave productMassMightHave = new ProductMassesMightHave(massLen);
            int i = 0;
            int ib = 0;
            int ic = 0;
            for (int j = 0; j < CompactPeptideDraw.NTerminalMasses.Length; j++)
            {
                var hm = CompactPeptideDraw.NTerminalMasses[j];
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
                    productMassMightHave.ProductMz[i] = hm + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "c" + (ic + 1).ToString();
                    i++;
                    ic++;
                }
            }
            int iy = CompactPeptideDraw.CTerminalMasses.Length - 1;
            int iz = CompactPeptideDraw.CTerminalMasses.Length - 1;
            for (int j = 0; j < CompactPeptideDraw.CTerminalMasses.Length; j++)
            {
                var hm = CompactPeptideDraw.CTerminalMasses[j];
                if (containsY)
                {
                    productMassMightHave.ProductMz[i] = hm + waterMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "y" + (CompactPeptideDraw.CTerminalMasses.Length - iy).ToString();
                    i++;
                    iy--;
                }
                if (containsZdot)
                {
                    productMassMightHave.ProductMz[i] = hm + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass;
                    productMassMightHave.ProductName[i] = "z" + (CompactPeptideDraw.CTerminalMasses.Length - iz).ToString();
                    i++;
                    iz--;
                }
            }
            return productMassMightHave;
        }
    }
}
