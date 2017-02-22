using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class PeptideWithSetModifications : Peptide
    {

        #region Public Fields

        public Dictionary<int, ModificationWithMass> allModsOneIsNterminus;

        #endregion Public Fields

        #region Private Fields

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        private static readonly double protonMass = Constants.protonMass;
        private static readonly double electronMass = Constants.electronMass;

        private readonly PeptideWithPossibleModifications modPep;
        private double? monoisotopicMass;

        private string sequence;

        private PeptideFragmentMasses p;

        private bool? hasChemicalFormulas;

        private string sequenceWithChemicalFormulas;

        #endregion Private Fields

        #region Internal Constructors

        internal PeptideWithSetModifications(PeptideWithPossibleModifications modPep, Dictionary<int, ModificationWithMass> allModsOneIsNterminus)
                                                                                            : base(modPep.Protein, modPep.OneBasedStartResidueInProtein, modPep.OneBasedEndResidueInProtein)
        {
            this.modPep = modPep;
            this.allModsOneIsNterminus = allModsOneIsNterminus;
        }

        #endregion Internal Constructors

        #region Public Properties

        public double MonoisotopicMass
        {
            get
            {
                if (!monoisotopicMass.HasValue)
                {
                    monoisotopicMass = waterMonoisotopicMass;
                    foreach (var mod in allModsOneIsNterminus.Values)
                        monoisotopicMass += mod.monoisotopicMass;
                    monoisotopicMass += BaseSequence.Select(b => Residue.ResidueMonoisotopicMass[b]).Sum();
                }
                return monoisotopicMass.Value;
            }
        }

        public virtual string Sequence
        {
            get
            {
                if (sequence == null)
                {
                    var sbsequence = new StringBuilder();

                    // variable modification on peptide N-terminus
                    ModificationWithMass pep_n_term_variable_mod;
                    if (allModsOneIsNterminus.TryGetValue(1, out pep_n_term_variable_mod))
                        sbsequence.Append('[' + pep_n_term_variable_mod.database + ":" + pep_n_term_variable_mod.id + ']');

                    for (int r = 0; r < Length; r++)
                    {
                        sbsequence.Append(this[r]);
                        // variable modification on this residue
                        ModificationWithMass residue_variable_mod;
                        if (allModsOneIsNterminus.TryGetValue(r + 2, out residue_variable_mod))
                            sbsequence.Append('[' + residue_variable_mod.database + ":" + residue_variable_mod.id + ']');
                    }

                    // variable modification on peptide C-terminus
                    ModificationWithMass pep_c_term_variable_mod;
                    if (allModsOneIsNterminus.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                        sbsequence.Append('[' + pep_c_term_variable_mod.database + ":" + pep_c_term_variable_mod.id + ']');

                    sequence = sbsequence.ToString();
                }
                return sequence;
            }
        }

        public int NumMods
        {
            get
            {
                return allModsOneIsNterminus.Count;
            }
        }

        public int MissedCleavages
        {
            get { return modPep.MissedCleavages; }
        }

        public override string PeptideDescription
        {
            get { return modPep.PeptideDescription; }
        }

        public string SequenceWithChemicalFormulas
        {
            get
            {
                if (!hasChemicalFormulas.HasValue)
                {
                    hasChemicalFormulas = true;
                    var sbsequence = new StringBuilder();

                    // variable modification on peptide N-terminus
                    ModificationWithMass pep_n_term_variable_mod;
                    if (allModsOneIsNterminus.TryGetValue(1, out pep_n_term_variable_mod))
                    {
                        var jj = pep_n_term_variable_mod as ModificationWithMassAndCf;
                        if (jj != null)
                            sbsequence.Append('[' + jj.chemicalFormula.Formula + ']');
                        else
                            return null;
                    }

                    for (int r = 0; r < Length; r++)
                    {
                        sbsequence.Append(this[r]);
                        // variable modification on this residue
                        ModificationWithMass residue_variable_mod;
                        if (allModsOneIsNterminus.TryGetValue(r + 2, out residue_variable_mod))
                        {
                            var jj = residue_variable_mod as ModificationWithMassAndCf;
                            if (jj != null)
                                sbsequence.Append('[' + jj.chemicalFormula.Formula + ']');
                            else
                                return null;
                        }
                    }

                    // variable modification on peptide C-terminus
                    ModificationWithMass pep_c_term_variable_mod;
                    if (allModsOneIsNterminus.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                    {
                        var jj = pep_c_term_variable_mod as ModificationWithMassAndCf;
                        if (jj != null)
                            sbsequence.Append('[' + jj.chemicalFormula.Formula + ']');
                        else
                            return null;
                    }

                    sequenceWithChemicalFormulas = sbsequence.ToString();
                }
                return sequenceWithChemicalFormulas;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public PeptideWithSetModifications Localize(int j, double v)
        {
            var vvv = new Dictionary<int, ModificationWithMass>(allModsOneIsNterminus);
            ModificationWithMass existingMod;
            double massInMs2OfExistingMod = 0;
            if (vvv.TryGetValue(j + 2, out existingMod))
            {
                massInMs2OfExistingMod = existingMod.monoisotopicMass - existingMod.neutralLoss;
                vvv.Remove(j + 2);
            }
            vvv.Add(j + 2, new ModificationWithMass(null, null, null, ModificationSites.Any, v + massInMs2OfExistingMod, null, 0, new List<double> { v + massInMs2OfExistingMod }, null, null));
            var hm = new PeptideWithSetModifications(modPep, vvv);
            return hm;
        }

        public double[] FastSortedProductMasses(List<ProductType> productTypes)
        {
            if (p == null)
                ComputeFragmentMasses();

            double[] productsB = null;
            double[] productsY = null;
            double[] productsC = null;
            double[] productsZ = null;

            if (productTypes.Contains(ProductType.B))
                productsB = new double[Length - 2];
            if (productTypes.Contains(ProductType.Y))
                productsY = new double[Length - 1];
            if (productTypes.Contains(ProductType.C))
                productsC = new double[Length - 1];
            if (productTypes.Contains(ProductType.Zdot))
                productsZ = new double[Length - 1];

            int iB = 0;
            int iY = 0;
            int iC = 0;
            int iZ = 0;

            for (int r = 1; r < Length; r++)
            {
                foreach (ProductType product_type in productTypes)
                {
                    if (!(product_type == ProductType.C && r < Length && this[r] == 'P') &&
                       !(product_type == ProductType.Zdot && Length - r < Length && this[Length - r] == 'P') &&
                       !(product_type == ProductType.B && r == 1))
                    {
                        switch (product_type)
                        {
                            // p.cumulativeNTerminalMass[r] refers to a generic fragment; ion type masses are determined 
                            // by adding/subtracting a constant value, depending on the fragment type
                            case ProductType.Adot:
                                throw new NotImplementedException();

                            case ProductType.B:
                                productsB[iB] = p.cumulativeNTerminalMass[r];
                                iB++;
                                break;

                            case ProductType.C:
                                productsC[iC] = p.cumulativeNTerminalMass[r] + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
                                iC++;
                                break;

                            case ProductType.X:
                                throw new NotImplementedException();

                            case ProductType.Y:
                                productsY[iY] = p.cumulativeCTerminalMass[r] + waterMonoisotopicMass;
                                iY++;
                                break;

                            case ProductType.Zdot:
                                productsZ[iZ] = p.cumulativeCTerminalMass[r] + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass - hydrogenAtomMonoisotopicMass + protonMass + electronMass;
                                iZ++;
                                break;
                        }
                    }
                }
            }
            iB = 0;
            iY = 0;
            iC = 0;
            iZ = 0;

            IEnumerable<double> enumProducts = new double[0];

            if (productsB != null)
                enumProducts = enumProducts.Concat(productsB);
            if (productsY != null)
                enumProducts = enumProducts.Concat(productsY);
            if (productsC != null)
                enumProducts = enumProducts.Concat(productsC);
            if (productsZ != null)
                enumProducts = enumProducts.Concat(productsZ);

            double[] products = enumProducts.Where(f => !double.IsNaN(f)).ToArray();
            Array.Sort(products);
            
            return products;
        }

        public override bool Equals(object obj)
        {
            var q = obj as PeptideWithSetModifications;
            return q != null && q.Sequence.Equals(Sequence) && q.OneBasedStartResidueInProtein == OneBasedStartResidueInProtein && q.Protein.Equals(Protein);
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode();
        }

        #endregion Public Methods

        #region Private Methods

        private void ComputeFragmentMasses()
        {
            p = new PeptideFragmentMasses();

            p.cumulativeNTerminalMass = new double[Length];

            // N-terminus
            ModificationWithMass pep_n_term_variable_mod;
            if (allModsOneIsNterminus.TryGetValue(1, out pep_n_term_variable_mod))
            {
                p.cumulativeNTerminalMass[0] = pep_n_term_variable_mod.monoisotopicMass - pep_n_term_variable_mod.neutralLoss;
            }

            // Loop for cumulative n mass
            for (int r = 1; r < Length; r++)
            {
                p.cumulativeNTerminalMass[r] = p.cumulativeNTerminalMass[r - 1] + Residue.ResidueMonoisotopicMass[this[r - 1]];
                ModificationWithMass residue_variable_mod;
                if (allModsOneIsNterminus.TryGetValue(r + 1, out residue_variable_mod))
                {
                    p.cumulativeNTerminalMass[r] += residue_variable_mod.monoisotopicMass - residue_variable_mod.neutralLoss;
                }
            }

            p.cumulativeCTerminalMass = new double[Length];

            // C-terminus
            ModificationWithMass pep_c_term_variable_mod;
            if (allModsOneIsNterminus.TryGetValue(Length + 2, out pep_c_term_variable_mod))
            {
                p.cumulativeCTerminalMass[0] = pep_c_term_variable_mod.monoisotopicMass - pep_c_term_variable_mod.neutralLoss;
            }

            // Loop for cumulative c terminal mass
            for (int r = 1; r < Length; r++)
            {
                p.cumulativeCTerminalMass[r] = p.cumulativeCTerminalMass[r - 1] + Residue.ResidueMonoisotopicMass[this[Length - r]];
                ModificationWithMass residue_variable_mod;
                if (allModsOneIsNterminus.TryGetValue(Length - r + 2, out residue_variable_mod))
                {
                    p.cumulativeCTerminalMass[r] += residue_variable_mod.monoisotopicMass - residue_variable_mod.neutralLoss;
                }
            }
        }

        #endregion Private Methods

        #region Private Classes

        private class PeptideFragmentMasses
        {

            #region Internal Fields

            internal double[] cumulativeCTerminalMass;
            internal double[] cumulativeNTerminalMass;

            #endregion Internal Fields

        }

        #endregion Private Classes

    }
}