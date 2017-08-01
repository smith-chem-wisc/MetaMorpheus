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

        public readonly int numFixedMods;
        public readonly Dictionary<int, ModificationWithMass> allModsOneIsNterminus;
        public readonly PeptideWithPossibleModifications modPep;

        #endregion Public Fields

        #region Private Fields

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private string sequence;
        private bool? hasChemicalFormulas;
        private string sequenceWithChemicalFormulas;
        private object lockObj = new object();

        private CompactPeptide compactPeptide;
        private double? monoisotopicMass;

        #endregion Private Fields

        #region Internal Constructors

        internal PeptideWithSetModifications(PeptideWithPossibleModifications modPep, Dictionary<int, ModificationWithMass> allModsOneIsNterminus, int numFixedMods, bool addCompIons)
                                                                                                            : base(modPep.Protein, modPep.OneBasedStartResidueInProtein, modPep.OneBasedEndResidueInProtein, addCompIons)
        {
            this.modPep = modPep;
            this.allModsOneIsNterminus = allModsOneIsNterminus;
            this.numFixedMods = numFixedMods;
        }

        internal PeptideWithSetModifications(PeptideWithSetModifications modsFromThisOne, PeptideWithSetModifications everythingElseFromThisOne, bool addCompIons) : base(everythingElseFromThisOne.Protein, everythingElseFromThisOne.OneBasedStartResidueInProtein, everythingElseFromThisOne.OneBasedEndResidueInProtein, addCompIons)
        {
            this.modPep = everythingElseFromThisOne.modPep;
            this.allModsOneIsNterminus = modsFromThisOne.allModsOneIsNterminus;
            this.numFixedMods = modsFromThisOne.numFixedMods;
        }

        #endregion Internal Constructors

        #region Public Properties

        public CompactPeptide CompactPeptide
        {
            get
            {
                if (compactPeptide == null)
                {
                    compactPeptide = new CompactPeptide(this, addCompIons);
                }
                return compactPeptide;
            }
        }

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
                        if (jj != null && Math.Abs(jj.chemicalFormula.MonoisotopicMass - jj.monoisotopicMass) < 1e-5)
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
                            if (jj != null && Math.Abs(jj.chemicalFormula.MonoisotopicMass - jj.monoisotopicMass) < 1e-5)
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
                        if (jj != null && Math.Abs(jj.chemicalFormula.MonoisotopicMass - jj.monoisotopicMass) < 1e-5)
                            sbsequence.Append('[' + jj.chemicalFormula.Formula + ']');
                        else
                            return null;
                    }

                    sequenceWithChemicalFormulas = sbsequence.ToString();
                }
                return sequenceWithChemicalFormulas;
            }
        }

        public int NumVariableMods { get { return this.NumMods - this.numFixedMods; } }

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
            vvv.Add(j + 2, new ModificationWithMass(null, null, null, ModificationSites.Any, massToLocalize + massOfExistingMod, null, new List<double> { 0 }, new List<double> { massToLocalize + massOfExistingMod }, null));
            var hm = new PeptideWithSetModifications(modPep, vvv, numFixedMods, addCompIons);
            return hm;
        }

        public override bool Equals(object obj)
        {
            var q = obj as PeptideWithSetModifications;
            //foreach (var ok in q.allModsOneIsNterminus)
            //    if (!allModsOneIsNterminus.Contains(ok))
            //        return false;
            return q != null
                && q.Sequence.Equals(Sequence)
                && q.OneBasedStartResidueInProtein == OneBasedStartResidueInProtein
                && q.Protein.Equals(Protein);
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode();
        }

        #endregion Public Methods

    }
}