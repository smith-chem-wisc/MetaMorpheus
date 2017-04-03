using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    internal struct MetaMorpheusFragment
    {

        #region Public Fields

        public int index;
        public double mass;

        #endregion Public Fields

    }

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
                        sbsequence.Append('[' + pep_n_term_variable_mod.modificationType + ":" + pep_n_term_variable_mod.id + ']');

                    for (int r = 0; r < Length; r++)
                    {
                        sbsequence.Append(this[r]);
                        // variable modification on this residue
                        ModificationWithMass residue_variable_mod;
                        if (allModsOneIsNterminus.TryGetValue(r + 2, out residue_variable_mod))
                            sbsequence.Append('[' + residue_variable_mod.modificationType + ":" + residue_variable_mod.id + ']');
                    }

                    // variable modification on peptide C-terminus
                    ModificationWithMass pep_c_term_variable_mod;
                    if (allModsOneIsNterminus.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                        sbsequence.Append('[' + pep_c_term_variable_mod.modificationType + ":" + pep_c_term_variable_mod.id + ']');

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

        public PeptideWithSetModifications Localize(int j, double massToLocalize)
        {
            var vvv = new Dictionary<int, ModificationWithMass>(allModsOneIsNterminus);
            ModificationWithMass modToReplace;
            double massOfExistingMod = 0;
            if (vvv.TryGetValue(j + 2, out modToReplace))
            {
                massOfExistingMod = modToReplace.monoisotopicMass;
                vvv.Remove(j + 2);
            }
            vvv.Add(j + 2, new ModificationWithMass(null, null, null, ModificationSites.Any, massToLocalize + massOfExistingMod, null, new List<double> { 0 }, new List<double> { massToLocalize + massOfExistingMod }, null, null));
            var hm = new PeptideWithSetModifications(modPep, vvv);
            return hm;
        }

        public double[] SortedProductMassesMightNotBeUnique(List<ProductType> productTypes)
        {
            if (p == null)
                ComputeFragmentMasses();

            IEnumerable<double> allMasses = new List<double>();

            // TODO: THIS
            //if (!(product_type == ProductType.C && r < Length && this[r] == 'P') &&
            //!(product_type == ProductType.Zdot && Length - r < Length && this[Length - r] == 'P')
            if (productTypes.Contains(ProductType.Adot))
                throw new NotImplementedException();
            if (productTypes.Contains(ProductType.B))
                allMasses = allMasses.Concat(p.nTerminalMasses.Where(b => b.index > 1).Select(b => b.mass));
            if (productTypes.Contains(ProductType.C))
                allMasses = allMasses.Concat(p.nTerminalMasses.Select(b => b.mass + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass));
            if (productTypes.Contains(ProductType.X))
                throw new NotImplementedException();
            if (productTypes.Contains(ProductType.Y))
                allMasses = allMasses.Concat(p.cTerminalMasses.Select(b => b.mass + waterMonoisotopicMass));
            if (productTypes.Contains(ProductType.Zdot))
                allMasses = allMasses.Concat(p.cTerminalMasses.Select(b => b.mass + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass));

            return allMasses.Where(f => !double.IsNaN(f)).OrderBy(b => b).ToArray();
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

        /// <summary>
        /// Computes fragment masses from both termini. Neutral losses are accounted for! Does NOT include the last residue!
        /// </summary>
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
            if (allModsOneIsNterminus.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                foreach (double nl in pep_c_term_variable_mod.neutralLosses)
                    theMass = pep_c_term_variable_mod.monoisotopicMass - nl;
            else
                theMass = 0;
            p.cTerminalMasses = ComputeFollowingFragmentMasses(theMass, Length, -1).ToList();
        }

        private IEnumerable<MetaMorpheusFragment> ComputeFollowingFragmentMasses(double prevMass, int oneBasedIndexToLookAt, int direction)
        {
            ModificationWithMass residue_variable_mod = null;
            do
            {
                prevMass += Residue.ResidueMonoisotopicMass[this[oneBasedIndexToLookAt - 1]];

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
                oneBasedIndexToLookAt += direction;
            } while (((oneBasedIndexToLookAt > 1 && direction == -1) || (oneBasedIndexToLookAt < Length && direction == 1)) && (residue_variable_mod == null || residue_variable_mod.neutralLosses.Count() == 1));
            if (residue_variable_mod != null && residue_variable_mod.neutralLosses.Count() > 1)
                foreach (double nl in residue_variable_mod.neutralLosses)
                {
                    var theMass = prevMass + residue_variable_mod.monoisotopicMass - nl;
                    var theFrag = new MetaMorpheusFragment()
                    {
                        mass = theMass,
                        index = oneBasedIndexToLookAt
                    };
                    yield return theFrag;
                    if ((direction == 1 && oneBasedIndexToLookAt + direction < Length) ||
                        (direction == -1 && oneBasedIndexToLookAt + direction > 1))
                        foreach (var nextMass in ComputeFollowingFragmentMasses(theMass, oneBasedIndexToLookAt + direction, direction))
                            yield return nextMass;
                }
        }

        #endregion Private Methods

        #region Private Classes

        private class PeptideFragmentMasses
        {

            #region Internal Fields

            internal List<MetaMorpheusFragment> cTerminalMasses;
            internal List<MetaMorpheusFragment> nTerminalMasses;

            #endregion Internal Fields

        }

        #endregion Private Classes

    }
}