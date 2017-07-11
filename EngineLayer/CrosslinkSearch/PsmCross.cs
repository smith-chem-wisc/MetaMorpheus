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

        public PsmCross(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(notch, score, scanIndex, scan, theBestPeptide.MonoisotopicMassIncludingFixedMods)
        {
            compactPeptide = theBestPeptide;
        }

        public PsmCross(PeptideWithSetModifications ps, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(notch, score, scanIndex, scan, ps.MonoisotopicMass)
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


        #region Public Methods

        public override CompactPeptide GetCompactPeptide(Dictionary<ModificationWithMass, ushort> modsDictionary)
        {
            return compactPeptide;

        }
        public CompactPeptide GetCompactPeptidePs(Dictionary<ModificationWithMass, ushort> modsDictionary)
        {
            if (compactPeptide == null)
                compactPeptide = new CompactPeptide(ps, modsDictionary);
            return compactPeptide;
        }

        #endregion Public Methods

        //Compute ProductMassesMightHave: the theoritical masses of psmCross
        #region ProductMassesMightHave

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        public Dictionary<int, ModificationWithMass> allModsOneIsNterminus;
        
        private PeptideFragmentMasses p;

        private void computeAllModsOneIsNterminus(Dictionary<ModificationWithMass, ushort> modsDictionary)
        {
            allModsOneIsNterminus = new Dictionary<int, ModificationWithMass>();
            ModificationWithMass ok;
            int i;
            if (compactPeptide.varMod1Type != 0)
            {
                i = (int)compactPeptide.varMod1Loc;
                ok = modsDictionary.FirstOrDefault(x => x.Value == compactPeptide.varMod1Type).Key;
                allModsOneIsNterminus.Add(i, ok);
            }
            if (compactPeptide.varMod2Type != 0)
            {
                i = (int)compactPeptide.varMod1Loc;
                ok = modsDictionary.FirstOrDefault(x => x.Value == compactPeptide.varMod1Type).Key;
                allModsOneIsNterminus.Add(i, ok);
            }
            if (compactPeptide.varMod3Type != 0)
            {
                i = (int)compactPeptide.varMod1Loc;
                ok = modsDictionary.FirstOrDefault(x => x.Value == compactPeptide.varMod1Type).Key;
                allModsOneIsNterminus.Add(i, ok);
            }
        }

        public ProductMassesMightHave ProductMassesMightHaveDuplicatesAndNaNs(List<ProductType> productTypes, Dictionary<ModificationWithMass, ushort> modsDictionary)
        {
            computeAllModsOneIsNterminus(modsDictionary);
            //PeptideFragmentMasses p = new PeptideFragmentMasses();
            int massLen = 0;
            bool containsAdot = productTypes.Contains(ProductType.Adot);
            bool containsB = productTypes.Contains(ProductType.B);
            bool containsC = productTypes.Contains(ProductType.C);
            bool containsX = productTypes.Contains(ProductType.X);
            bool containsY = productTypes.Contains(ProductType.Y);
            bool containsZdot = productTypes.Contains(ProductType.Zdot);

            if (p == null)
                ComputeFragmentMasses();

            //if (containsAdot)
            //    throw new NotImplementedException();
            if (containsB)
                massLen += p.nTerminalMasses.Count(b => b.index > 1);
            if (containsC)
                massLen += p.nTerminalMasses.Count;
            //if (containsX)
            //    throw new NotImplementedException();
            if (containsY)
                massLen += p.cTerminalMasses.Count;
            if (containsZdot)
                massLen += p.cTerminalMasses.Count;


            ProductMassesMightHave indexesMassesToReturn = new ProductMassesMightHave(massLen);

            int i = 0;
            foreach (var hm in p.nTerminalMasses)
            {
                if (hm.index > 1 && containsB)
                {
                    indexesMassesToReturn.ProductMz[i] = hm.mass;
                    indexesMassesToReturn.ProductName[i] = "b" + hm.index.ToString();
                    i++;
                }
                if (containsC)
                {
                    indexesMassesToReturn.ProductMz[i] = hm.mass + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
                    indexesMassesToReturn.ProductName[i] = "c" + hm.index.ToString();
                    i++;
                }
            }
            foreach (var hm in p.cTerminalMasses)
            {
                if (containsY)
                {
                    indexesMassesToReturn.ProductMz[i] = hm.mass + waterMonoisotopicMass;
                    indexesMassesToReturn.ProductName[i] = "y" + (compactPeptide.BaseSequence.Length - hm.index + 1).ToString();
                    i++;
                }
                if (containsZdot)
                {
                    indexesMassesToReturn.ProductMz[i] = hm.mass + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass;
                    indexesMassesToReturn.ProductName[i] = "z" + (compactPeptide.BaseSequence.Length - hm.index + 1).ToString();
                    i++;
                }
            }
            return indexesMassesToReturn;
        }

        private void ComputeFragmentMasses()
        {
            p = new PeptideFragmentMasses();

            ModificationWithMass pep_n_term_variable_mod;
            double theMass = 0;
            if (allModsOneIsNterminus.TryGetValue(1, out pep_n_term_variable_mod))
                foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                    theMass = pep_n_term_variable_mod.monoisotopicMass - nl;
            else
                theMass = 0;
            p.nTerminalMasses = ComputeFollowingFragmentMasses(theMass, 1, 1).ToList();

            ModificationWithMass pep_c_term_variable_mod;
            theMass = 0;
            if (allModsOneIsNterminus.TryGetValue(compactPeptide.BaseSequence.Length + 2, out pep_c_term_variable_mod))
                foreach (double nl in pep_c_term_variable_mod.neutralLosses)
                    theMass = pep_c_term_variable_mod.monoisotopicMass - nl;
            else
                theMass = 0;
            p.cTerminalMasses = ComputeFollowingFragmentMasses(theMass, compactPeptide.BaseSequence.Length, -1).ToList();
        }

        private IEnumerable<MetaMorpheusFragment> ComputeFollowingFragmentMasses(double prevMass, int oneBasedIndexToLookAt, int direction)
        {
            ModificationWithMass residue_variable_mod = null;
            do
            {
                prevMass += Residue.ResidueMonoisotopicMass[compactPeptide.BaseSequence[oneBasedIndexToLookAt - 1]];

                allModsOneIsNterminus.TryGetValue(oneBasedIndexToLookAt + 1, out residue_variable_mod);
                if (residue_variable_mod == null)
                {
                    var theFrag = new MetaMorpheusFragment()
                    {
                        mass = prevMass,
                        index = oneBasedIndexToLookAt
                    };
                    yield return theFrag;
                }
                else if (residue_variable_mod.neutralLosses.Count() == 1)
                {
                    prevMass += residue_variable_mod.monoisotopicMass - residue_variable_mod.neutralLosses.First();
                    var theFrag = new MetaMorpheusFragment()
                    {
                        mass = prevMass,
                        index = oneBasedIndexToLookAt
                    };
                    yield return theFrag;
                }
                else
                {
                    foreach (double nl in residue_variable_mod.neutralLosses)
                    {
                        var theMass = prevMass + residue_variable_mod.monoisotopicMass - nl;
                        var theFrag = new MetaMorpheusFragment()
                        {
                            mass = theMass,
                            index = oneBasedIndexToLookAt
                        };
                        yield return theFrag;
                        if ((direction == 1 && oneBasedIndexToLookAt + direction < compactPeptide.BaseSequence.Length) ||
                            (direction == -1 && oneBasedIndexToLookAt + direction > 1))
                            foreach (var nextMass in ComputeFollowingFragmentMasses(theMass, oneBasedIndexToLookAt + direction, direction))
                                yield return nextMass;
                    }
                    break;
                }
                oneBasedIndexToLookAt += direction;
            } while ((oneBasedIndexToLookAt > 1 && direction == -1) || (oneBasedIndexToLookAt < compactPeptide.BaseSequence.Length && direction == 1));
        }

        private class PeptideFragmentMasses
        {

            #region Internal Fields

            internal List<MetaMorpheusFragment> cTerminalMasses;
            internal List<MetaMorpheusFragment> nTerminalMasses;

            #endregion Internal Fields

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
