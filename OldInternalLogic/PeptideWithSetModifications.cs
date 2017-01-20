using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Text;

namespace OldInternalLogic
{
    public class PeptideWithSetModifications : Peptide
    {

        #region Public Fields

        public Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss;

        #endregion Public Fields

        #region Private Fields

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        private readonly PeptideWithPossibleModifications modPep;
        private double monoisotopicMass = double.NaN;

        private string extendedSequence;

        private string sequence;

        #endregion Private Fields

        #region Public Constructors

        public PeptideWithSetModifications(PeptideWithPossibleModifications modPep, Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss)
                                                            : base(modPep.Protein, modPep.OneBasedStartResidueInProtein, modPep.OneBasedEndResidueInProtein)
        {
            this.modPep = modPep;
            this.twoBasedVariableAndLocalizeableModificationss = twoBasedVariableAndLocalizeableModificationss;
        }

        #endregion Public Constructors

        #region Public Properties

        public double MonoisotopicMass
        {
            get
            {
                if (double.IsNaN(monoisotopicMass))
                    ComputeFragmentMasses();
                return monoisotopicMass;
            }
        }

        public virtual string ExtendedSequence
        {
            get
            {
                if (extendedSequence == null)
                {
                    var exSeq = new StringBuilder();
                    exSeq.Append(PreviousAminoAcid);
                    exSeq.Append(".");
                    exSeq.Append(Sequence);
                    exSeq.Append(".");
                    exSeq.Append(NextAminoAcid);
                    extendedSequence = exSeq.ToString();
                }
                return extendedSequence;
            }
        }

        public virtual string Sequence
        {
            get
            {
                if (sequence == null)
                {
                    var sbsequence = new StringBuilder();
                    List<MorpheusModification> value;
                    // fixed modifications on protein N-terminus
                    if (modPep.twoBasedFixedModificationss.TryGetValue(0, out value))
                        foreach (var fixed_modification in value)
                            sbsequence.Append('[' + fixed_modification.Description + ']');
                    // variable modification on protein N-terminus
                    MorpheusModification prot_n_term_variable_mod;
                    if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(0, out prot_n_term_variable_mod))
                        sbsequence.Append('(' + prot_n_term_variable_mod.Description + ')');

                    // fixed modifications on peptide N-terminus
                    if (modPep.twoBasedFixedModificationss.TryGetValue(1, out value))
                        foreach (var fixed_modification in value)
                            sbsequence.Append('[' + fixed_modification.Description + ']');

                    // variable modification on peptide N-terminus
                    MorpheusModification pep_n_term_variable_mod;
                    if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(1, out pep_n_term_variable_mod))
                        sbsequence.Append('(' + pep_n_term_variable_mod.Description + ')');

                    for (int r = 0; r < Length; r++)
                    {
                        sbsequence.Append(this[r]);
                        // fixed modifications on this residue
                        if (modPep.twoBasedFixedModificationss.TryGetValue(r + 2, out value))
                            foreach (var fixed_modification in value)
                                sbsequence.Append('[' + fixed_modification.Description + ']');
                        // variable modification on this residue
                        MorpheusModification residue_variable_mod;
                        if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(r + 2, out residue_variable_mod))
                            sbsequence.Append('(' + residue_variable_mod.Description + ')');
                    }

                    // fixed modifications on peptide C-terminus
                    if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 2, out value))
                        foreach (var fixed_modification in value)
                            sbsequence.Append('[' + fixed_modification.Description + ']');

                    // variable modification on peptide C-terminus
                    MorpheusModification pep_c_term_variable_mod;
                    if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                        sbsequence.Append('(' + pep_c_term_variable_mod.Description + ')');

                    // fixed modifications on protein C-terminus
                    if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 3, out value))
                        foreach (var fixed_modification in value)
                            sbsequence.Append('[' + fixed_modification.Description + ']');
                    // variable modification on protein C-terminus
                    MorpheusModification prot_c_term_variable_mod;
                    if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 3, out prot_c_term_variable_mod))
                        sbsequence.Append('(' + prot_c_term_variable_mod.Description + ')');
                    sequence = sbsequence.ToString();
                }
                return sequence;
            }
        }

        public int NumVariableMods
        {
            get
            {
                return twoBasedVariableAndLocalizeableModificationss.Count;
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

        public IEnumerable<KeyValuePair<int, MorpheusModification>> Modifications
        {
            get
            {
                foreach (var ye in twoBasedVariableAndLocalizeableModificationss)
                    yield return ye;
                foreach (var aass in modPep.twoBasedFixedModificationss)
                {
                    int key = aass.Key;
                    foreach (var fff in aass.Value)
                        yield return new KeyValuePair<int, MorpheusModification>(key, fff);
                }
            }
        }

        public string SequenceWithChemicalFormulas
        {
            get
            {
                var sbsequence = new StringBuilder();
                List<MorpheusModification> value;
                // fixed modifications on protein N-terminus
                if (modPep.twoBasedFixedModificationss.TryGetValue(0, out value))
                    foreach (var fixed_modification in value)
                        sbsequence.Append('[' + fixed_modification.ChemicalFormula.Formula + ']');
                // variable modification on protein N-terminus
                MorpheusModification prot_n_term_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(0, out prot_n_term_variable_mod))
                    sbsequence.Append('[' + prot_n_term_variable_mod.ChemicalFormula.Formula + ']');

                // fixed modifications on peptide N-terminus
                if (modPep.twoBasedFixedModificationss.TryGetValue(1, out value))
                    foreach (var fixed_modification in value)
                        sbsequence.Append('[' + fixed_modification.ChemicalFormula.Formula + ']');

                // variable modification on peptide N-terminus
                MorpheusModification pep_n_term_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(1, out pep_n_term_variable_mod))
                    sbsequence.Append('[' + pep_n_term_variable_mod.ChemicalFormula.Formula + ']');

                for (int r = 0; r < Length; r++)
                {
                    sbsequence.Append(this[r]);
                    // fixed modifications on this residue
                    if (modPep.twoBasedFixedModificationss.TryGetValue(r + 2, out value))
                        foreach (var fixed_modification in value)
                            sbsequence.Append('[' + fixed_modification.ChemicalFormula.Formula + ']');
                    // variable modification on this residue
                    MorpheusModification residue_variable_mod;
                    if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(r + 2, out residue_variable_mod))
                        sbsequence.Append('[' + residue_variable_mod.ChemicalFormula.Formula + ']');
                }

                // fixed modifications on peptide C-terminus
                if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 2, out value))
                    foreach (var fixed_modification in value)
                        sbsequence.Append('[' + fixed_modification.ChemicalFormula.Formula + ']');

                // variable modification on peptide C-terminus
                MorpheusModification pep_c_term_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                    sbsequence.Append('[' + pep_c_term_variable_mod.ChemicalFormula.Formula + ']');

                // fixed modifications on protein C-terminus
                if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 3, out value))
                    foreach (var fixed_modification in value)
                        sbsequence.Append('[' + fixed_modification.ChemicalFormula.Formula + ']');
                // variable modification on protein C-terminus
                MorpheusModification prot_c_term_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 3, out prot_c_term_variable_mod))
                    sbsequence.Append('[' + prot_c_term_variable_mod.ChemicalFormula.Formula + ']');
                return sbsequence.ToString();
            }
        }

        #endregion Public Properties

        #region Public Methods

        public PeptideWithSetModifications Localize(int j, double v)
        {
            var vvv = new Dictionary<int, MorpheusModification>(twoBasedVariableAndLocalizeableModificationss);
            MorpheusModification existingMod;
            double massInMs2OfExistingMod = 0;
            if (vvv.TryGetValue(j + 2, out existingMod))
            {
                massInMs2OfExistingMod = existingMod.Labile ? 0 : existingMod.MonoisotopicMassShift;
                vvv.Remove(j + 2);
            }
            vvv.Add(j + 2, new MorpheusModification(null, ModificationType.AminoAcidResidue, '\0', v + massInMs2OfExistingMod, null, null, '\0', double.NaN, false, null));
            var hm = new PeptideWithSetModifications(modPep, vvv);
            return hm;
        }

        public double[] FastSortedProductMasses(List<ProductType> productTypes)
        {
            PeptideFragmentMasses p = ComputeFragmentMasses();
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
            return ExtendedSequence.GetHashCode();
        }

        #endregion Public Methods

        #region Private Methods

        private PeptideFragmentMasses ComputeFragmentMasses()
        {
            var p = new PeptideFragmentMasses();

            monoisotopicMass = waterMonoisotopicMass;

            double mass_shift;
            p.cumulativeNTerminalMass = new double[Length];
            mass_shift = 0.0f;
            List<MorpheusModification> modificationList;
            // fixed modifications on protein N-terminus
            if (OneBasedEndResidueInProtein == 1 || OneBasedEndResidueInProtein == 2)
                if (modPep.twoBasedFixedModificationss.TryGetValue(0, out modificationList))
                    foreach (var fixed_modification in modificationList)
                    {
                        if (!fixed_modification.Labile)
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                        monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                    }
            // variable modification on the protein N-terminus
            MorpheusModification prot_n_term_var_mod;
            if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(0, out prot_n_term_var_mod))
            {
                if (!prot_n_term_var_mod.Labile)
                    mass_shift += prot_n_term_var_mod.MonoisotopicMassShift;
                monoisotopicMass += prot_n_term_var_mod.MonoisotopicMassShift;
            }
            // fixed modifications on peptide N-terminus
            if (modPep.twoBasedFixedModificationss.TryGetValue(1, out modificationList))
                foreach (var fixed_modification in modificationList)
                {
                    if (!fixed_modification.Labile)
                        mass_shift += fixed_modification.MonoisotopicMassShift;
                    monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                }
            // variable modification on peptide N-terminus
            MorpheusModification pep_n_term_variable_mod;
            if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(1, out pep_n_term_variable_mod))
            {
                if (!pep_n_term_variable_mod.Labile)
                    mass_shift += pep_n_term_variable_mod.MonoisotopicMassShift;
                monoisotopicMass += pep_n_term_variable_mod.MonoisotopicMassShift;
            }
            p.cumulativeNTerminalMass[0] = mass_shift;

            for (int r = 1; r < Length; r++)
            {
                mass_shift = 0.0f;
                // fixed modifications on this residue
                if (modPep.twoBasedFixedModificationss.TryGetValue(r + 1, out modificationList))
                    foreach (var fixed_modification in modificationList)
                    {
                        if (!fixed_modification.Labile)
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                        monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                    }
                // variable modification on this residue
                MorpheusModification residue_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(r + 1, out residue_variable_mod))
                {
                    if (!residue_variable_mod.Labile)
                        mass_shift += residue_variable_mod.MonoisotopicMassShift;
                    monoisotopicMass += residue_variable_mod.MonoisotopicMassShift;
                }
                p.cumulativeNTerminalMass[r] = p.cumulativeNTerminalMass[r - 1] + Residue.ResidueMonoisotopicMass[this[r - 1]] + mass_shift;
            }

            p.cumulativeCTerminalMass = new double[Length];

            mass_shift = 0.0f;
            // fixed modifications on peptide C-terminus
            if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 2, out modificationList))
                foreach (var fixed_modification in modificationList)
                {
                    if (!fixed_modification.Labile)
                        mass_shift += fixed_modification.MonoisotopicMassShift;
                    monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                }
            // variable modification on peptide C-terminus
            MorpheusModification pep_c_term_variable_mod;
            if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 2, out pep_c_term_variable_mod))
            {
                if (!pep_c_term_variable_mod.Labile)
                    mass_shift += pep_c_term_variable_mod.MonoisotopicMassShift;
                monoisotopicMass += pep_c_term_variable_mod.MonoisotopicMassShift;
            }
            // fixed modifications on protein C-terminus
            if (OneBasedEndResidueInProtein == Protein.Length)
                if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 3, out modificationList))
                    foreach (var fixed_modification in modificationList)
                    {
                        if (!fixed_modification.Labile)
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                        monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                    }
            // variable modification on protein C-terminus
            MorpheusModification prot_c_term_variable_mod;
            if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 3, out prot_c_term_variable_mod))
            {
                if (!prot_c_term_variable_mod.Labile)
                    mass_shift += prot_c_term_variable_mod.MonoisotopicMassShift;
                monoisotopicMass += prot_c_term_variable_mod.MonoisotopicMassShift;
            }

            p.cumulativeCTerminalMass[0] = mass_shift;
            monoisotopicMass += Residue.ResidueMonoisotopicMass[BaseSequence[0]];

            for (int r = 1; r < Length; r++)
            {
                mass_shift = 0.0f;
                monoisotopicMass += Residue.ResidueMonoisotopicMass[BaseSequence[r]];
                // fixed modifications on this residue
                if (modPep.twoBasedFixedModificationss.TryGetValue(Length - r + 2, out modificationList))
                    foreach (var fixed_modification in modificationList)
                        if (!fixed_modification.Labile)
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                // variable modification on this residue
                MorpheusModification residue_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length - r + 2, out residue_variable_mod))
                    if (!residue_variable_mod.Labile)
                        mass_shift += residue_variable_mod.MonoisotopicMassShift;
                p.cumulativeCTerminalMass[r] = p.cumulativeCTerminalMass[r - 1] + Residue.ResidueMonoisotopicMass[this[Length - r]] + mass_shift;
            }

            return p;
        }

        #endregion Private Methods

    }
}