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
        public readonly int? missedCleavages;

        #endregion Public Fields

        #region Private Fields

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private readonly Dictionary<TerminusType, CompactPeptide> compactPeptides = new Dictionary<TerminusType, CompactPeptide>();
        private string sequence;
        private bool? hasChemicalFormulas;
        private string sequenceWithChemicalFormulas;
        private object lockObj = new object();
        private double? monoisotopicMass;

        #endregion Private Fields

        #region Public Constructors

        public PeptideWithSetModifications(PeptideWithPossibleModifications modPep, Dictionary<int, ModificationWithMass> allModsOneIsNterminus, int numFixedMods)
            : base(modPep.Protein, modPep.OneBasedStartResidueInProtein, modPep.OneBasedEndResidueInProtein, modPep.PeptideDescription)
        {
            this.missedCleavages = modPep.MissedCleavages;
            this.allModsOneIsNterminus = allModsOneIsNterminus;
            this.numFixedMods = numFixedMods;
        }

        public PeptideWithSetModifications(PeptideWithSetModifications modsFromThisOne, PeptideWithSetModifications everythingElseFromThisOne) : base(everythingElseFromThisOne.Protein, everythingElseFromThisOne.OneBasedStartResidueInProtein, everythingElseFromThisOne.OneBasedEndResidueInProtein, everythingElseFromThisOne.PeptideDescription)
        {
            this.missedCleavages = everythingElseFromThisOne.missedCleavages;
            this.allModsOneIsNterminus = modsFromThisOne.allModsOneIsNterminus;
            this.numFixedMods = modsFromThisOne.numFixedMods;
        }

        public PeptideWithSetModifications(PeptideWithSetModifications modsFromThisOne, int proteinOneBasedStart, int proteinOneBasedEnd) : base(modsFromThisOne.Protein, proteinOneBasedStart, proteinOneBasedEnd,modsFromThisOne.PeptideDescription)
        {
            this.allModsOneIsNterminus = modsFromThisOne.allModsOneIsNterminus.Where(b => b.Key > (1 + proteinOneBasedStart - modsFromThisOne.OneBasedStartResidueInProtein) && b.Key <= (2 + proteinOneBasedEnd - modsFromThisOne.OneBasedStartResidueInProtein)).ToDictionary(b => (b.Key+modsFromThisOne.OneBasedStartResidueInProtein-proteinOneBasedStart), b => b.Value);
        }

        public PeptideWithSetModifications(Dictionary<int, ModificationWithMass> allModsOneIsNterminus, int numFixedMods, Protein protein, int proteinOneBasedStart, int proteinOneBasedEnd, int? missedCleavages) : base(protein, proteinOneBasedStart, proteinOneBasedEnd)
        {
            this.missedCleavages = missedCleavages;
            this.numFixedMods = numFixedMods;
            this.allModsOneIsNterminus = allModsOneIsNterminus;
        }

        #endregion Public Constructors

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
                    if (allModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                        sbsequence.Append('[' + pep_n_term_variable_mod.modificationType + ":" + pep_n_term_variable_mod.id + ']');

                    for (int r = 0; r < Length; r++)
                    {
                        sbsequence.Append(this[r]);
                        // variable modification on this residue
                        if (allModsOneIsNterminus.TryGetValue(r + 2, out ModificationWithMass residue_variable_mod))
                            sbsequence.Append('[' + residue_variable_mod.modificationType + ":" + residue_variable_mod.id + ']');
                    }

                    // variable modification on peptide C-terminus
                    if (allModsOneIsNterminus.TryGetValue(Length + 2, out ModificationWithMass pep_c_term_variable_mod))
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

        public string SequenceWithChemicalFormulas
        {
            get
            {
                if (!hasChemicalFormulas.HasValue)
                {
                    hasChemicalFormulas = true;
                    var sbsequence = new StringBuilder();

                    // variable modification on peptide N-terminus
                    if (allModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                    {
                        if (pep_n_term_variable_mod is ModificationWithMassAndCf jj && Math.Abs(jj.chemicalFormula.MonoisotopicMass - jj.monoisotopicMass) < 1e-5)
                            sbsequence.Append('[' + jj.chemicalFormula.Formula + ']');
                        else
                            return null;
                    }

                    for (int r = 0; r < Length; r++)
                    {
                        sbsequence.Append(this[r]);
                        // variable modification on this residue
                        if (allModsOneIsNterminus.TryGetValue(r + 2, out ModificationWithMass residue_variable_mod))
                        {
                            if (residue_variable_mod is ModificationWithMassAndCf jj && Math.Abs(jj.chemicalFormula.MonoisotopicMass - jj.monoisotopicMass) < 1e-5)
                                sbsequence.Append('[' + jj.chemicalFormula.Formula + ']');
                            else
                                return null;
                        }
                    }

                    // variable modification on peptide C-terminus
                    if (allModsOneIsNterminus.TryGetValue(Length + 2, out ModificationWithMass pep_c_term_variable_mod))
                    {
                        if (pep_c_term_variable_mod is ModificationWithMassAndCf jj && Math.Abs(jj.chemicalFormula.MonoisotopicMass - jj.monoisotopicMass) < 1e-5)
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

        public CompactPeptide CompactPeptide(TerminusType terminusType)
        {
            if (compactPeptides.TryGetValue(terminusType, out CompactPeptide compactPeptide))
            {
                return compactPeptide;
            }
            else
            {
                CompactPeptide cp = new CompactPeptide(this, terminusType);
                compactPeptides.Add(terminusType, cp);
                return cp;
            }
        }

        public PeptideWithSetModifications Localize(int j, double massToLocalize)
        {
            var vvv = new Dictionary<int, ModificationWithMass>(allModsOneIsNterminus);
            double massOfExistingMod = 0;
            if (vvv.TryGetValue(j + 2, out ModificationWithMass modToReplace))
            {
                massOfExistingMod = modToReplace.monoisotopicMass;
                vvv.Remove(j + 2);
            }

            vvv.Add(j + 2, new ModificationWithMass(null, null, null, TerminusLocalization.Any, massToLocalize + massOfExistingMod, null, new List<double> { 0 }, new List<double> { massToLocalize + massOfExistingMod }, null));
            var hm = new PeptideWithSetModifications(vvv, numFixedMods, this.Protein, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein, this.missedCleavages);

            return hm;
        }

        public override bool Equals(object obj)
        {
            var q = obj as PeptideWithSetModifications;
            return q != null
                && q.Sequence.Equals(Sequence)
                && q.OneBasedStartResidueInProtein == OneBasedStartResidueInProtein
                && q.Protein.Accession.Equals(Protein.Accession);
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode();
        }

        #endregion Public Methods
    }
}