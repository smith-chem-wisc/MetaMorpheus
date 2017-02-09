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

        private readonly PeptideWithPossibleModifications modPep;
        private double? monoisotopicMass;

        private string sequence;

        private PeptideFragmentMasses p;

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
                        sbsequence.Append('[' + pep_n_term_variable_mod.id + ']');

                    for (int r = 0; r < Length; r++)
                    {
                        sbsequence.Append(this[r]);
                        // variable modification on this residue
                        ModificationWithMass residue_variable_mod;
                        if (allModsOneIsNterminus.TryGetValue(r + 2, out residue_variable_mod))
                            sbsequence.Append('[' + residue_variable_mod.id + ']');
                    }

                    // variable modification on peptide C-terminus
                    ModificationWithMass pep_c_term_variable_mod;
                    if (allModsOneIsNterminus.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                        sbsequence.Append('[' + pep_c_term_variable_mod.id + ']');

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

        bool? hasChemicalFormulas;
        private string sequenceWithChemicalFormulas;
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
            vvv.Add(j + 2, new ModificationWithMass(null, null, null, ModificationSites.Any, v + massInMs2OfExistingMod, 0, null));
            var hm = new PeptideWithSetModifications(modPep, vvv);
            return hm;
        }

        public double[] FastSortedProductMasses(List<ProductType> productTypes)
        {
            if (p == null)
                p = ComputeFragmentMasses();

            double[] products1 = null;
            double[] products2 = null;
            if (productTypes.Contains(ProductType.B))
                products1 = new double[Length - 2];
            if (productTypes.Contains(ProductType.Y))
                products2 = new double[Length - 1];

            int i1 = 0;
            int i2 = 0;
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
                            case ProductType.Adot:
                                throw new NotImplementedException();
                            case ProductType.B:
                                products1[i1] = p.cumulativeNTerminalMass[r];
                                i1++;
                                break;

                            case ProductType.C:
                                throw new NotImplementedException();

                            case ProductType.X:
                                throw new NotImplementedException();
                            case ProductType.Y:
                                products2[i2] = p.cumulativeCTerminalMass[r] + waterMonoisotopicMass;
                                i2++;
                                break;

                            case ProductType.Zdot:
                                throw new NotImplementedException();
                        }
                    }
                }
            }
            i1 = 0;
            i2 = 0;
            int len = (productTypes.Contains(ProductType.B) ? Length - 2 : 0) +
                      (productTypes.Contains(ProductType.Y) ? Length - 1 : 0);
            double[] products = new double[len];
            for (int i = 0; i < len; i++)
            {
                if (products1 != null && (products2 == null || (i1 != products1.Length && (i2 == products2.Length || products1[i1] <= products2[i2]))))
                {
                    products[i] = products1[i1];
                    i1++;
                }
                else
                {
                    products[i] = products2[i2];
                    i2++;
                }
            }
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

        private PeptideFragmentMasses ComputeFragmentMasses()
        {
            var p = new PeptideFragmentMasses();

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

            return p;
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