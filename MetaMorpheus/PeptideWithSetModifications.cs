using System;
using System.Collections.Generic;
using System.Text;

namespace MetaMorpheus
{
    public class PeptideWithSetModifications : Peptide
    {
        public PeptideWithSetModifications(PeptideWithPossibleModifications modPep, Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss)
            : base(modPep.protein, modPep.OneBasedStartResidueInProtein, modPep.OneBasedEndResidueInProtein)
        {
            this.modPep = modPep;
            this.twoBasedVariableAndLocalizeableModificationss = twoBasedVariableAndLocalizeableModificationss;
        }

        private double monoisotopicMass = double.NaN;

        public double MonoisotopicMass
        {
            get
            {
                if (double.IsNaN(monoisotopicMass))
                    computeFragmentMasses();
                return monoisotopicMass;
            }
            set
            {
                monoisotopicMass = value;
            }
        }

        private string extendedSequence;

        public virtual string ExtendedSequence
        {
            get
            {
                if (extendedSequence == null)
                {
                    StringBuilder exSeq = new StringBuilder();
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

        private string sequence;

        public virtual string Sequence
        {
            get
            {
                if (sequence == null)
                {
                    StringBuilder sbsequence = new StringBuilder();
                    List<MorpheusModification> value = null;
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

        public PeptideWithSetModifications Localize(int j, double v)
        {
            Dictionary<int, MorpheusModification> vvv = new Dictionary<int, MorpheusModification>(twoBasedVariableAndLocalizeableModificationss);
            MorpheusModification existingMod;
            double massInMs2OfExistingMod = 0;
            if (vvv.TryGetValue(j + 2, out existingMod))
            {
                massInMs2OfExistingMod = existingMod.labile ? 0 : existingMod.MonoisotopicMassShift;
                vvv.Remove(j + 2);
            }
            vvv.Add(j + 2, new MorpheusModification(null, ModificationType.AminoAcidResidue, '\0', v + massInMs2OfExistingMod, null, null, '\0', double.NaN, false, null));
            var hm = new PeptideWithSetModifications(modPep, vvv);
            return hm;
        }

        public PeptideFragmentMasses computeFragmentMasses()
        {
            PeptideFragmentMasses p = new PeptideFragmentMasses();

            monoisotopicMass = Constants.WATER_MONOISOTOPIC_MASS;

            double mass_shift;
            p.cumulativeNTerminalMass = new double[Length];
            mass_shift = 0.0f;
            List<MorpheusModification> modificationList;
            // fixed modifications on protein N-terminus
            if (OneBasedEndResidueInProtein == 1 || OneBasedEndResidueInProtein == 2)
                if (modPep.twoBasedFixedModificationss.TryGetValue(0, out modificationList))
                    foreach (var fixed_modification in modificationList)
                    {
                        if (!fixed_modification.labile)
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                        monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                    }
            // variable modification on the protein N-terminus
            MorpheusModification prot_n_term_var_mod;
            if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(0, out prot_n_term_var_mod))
            {
                if (!prot_n_term_var_mod.labile)
                    mass_shift += prot_n_term_var_mod.MonoisotopicMassShift;
                monoisotopicMass += prot_n_term_var_mod.MonoisotopicMassShift;
            }
            // fixed modifications on peptide N-terminus
            if (modPep.twoBasedFixedModificationss.TryGetValue(1, out modificationList))
                foreach (var fixed_modification in modificationList)
                {
                    if (!fixed_modification.labile)
                        mass_shift += fixed_modification.MonoisotopicMassShift;
                    monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                }
            // variable modification on peptide N-terminus
            MorpheusModification pep_n_term_variable_mod;
            if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(1, out pep_n_term_variable_mod))
            {
                if (!pep_n_term_variable_mod.labile)
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
                        if (!fixed_modification.labile)
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                        monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                    }
                // variable modification on this residue
                MorpheusModification residue_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(r + 1, out residue_variable_mod))
                {
                    if (!residue_variable_mod.labile)
                        mass_shift += residue_variable_mod.MonoisotopicMassShift;
                    monoisotopicMass += residue_variable_mod.MonoisotopicMassShift;
                }
                p.cumulativeNTerminalMass[r] = p.cumulativeNTerminalMass[r - 1] + AminoAcidMasses.GetMonoisotopicMass(this[r - 1]) + mass_shift;
            }

            p.cumulativeCTerminalMass = new double[Length];

            mass_shift = 0.0f;
            // fixed modifications on peptide C-terminus
            if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 2, out modificationList))
                foreach (var fixed_modification in modificationList)
                {
                    if (!fixed_modification.labile)
                        mass_shift += fixed_modification.MonoisotopicMassShift;
                    monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                }
            // variable modification on peptide C-terminus
            MorpheusModification pep_c_term_variable_mod;
            if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 2, out pep_c_term_variable_mod))
            {
                if (!pep_c_term_variable_mod.labile)
                    mass_shift += pep_c_term_variable_mod.MonoisotopicMassShift;
                monoisotopicMass += pep_c_term_variable_mod.MonoisotopicMassShift;
            }
            // fixed modifications on protein C-terminus
            if (OneBasedEndResidueInProtein == protein.Length)
                if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 3, out modificationList))
                    foreach (var fixed_modification in modificationList)
                    {
                        if (!fixed_modification.labile)
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                        monoisotopicMass += fixed_modification.MonoisotopicMassShift;
                    }
            // variable modification on protein C-terminus
            MorpheusModification prot_c_term_variable_mod;
            if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 3, out prot_c_term_variable_mod))
            {
                if (!prot_c_term_variable_mod.labile)
                    mass_shift += prot_c_term_variable_mod.MonoisotopicMassShift;
                monoisotopicMass += prot_c_term_variable_mod.MonoisotopicMassShift;
            }

            p.cumulativeCTerminalMass[0] = mass_shift;
            monoisotopicMass += AminoAcidMasses.GetMonoisotopicMass(BaseSequence[0]);

            for (int r = 1; r < Length; r++)
            {
                mass_shift = 0.0f;
                monoisotopicMass += AminoAcidMasses.GetMonoisotopicMass(BaseSequence[r]);
                // fixed modifications on this residue
                if (modPep.twoBasedFixedModificationss.TryGetValue(Length - r + 2, out modificationList))
                    foreach (var fixed_modification in modificationList)
                        if (!fixed_modification.labile)
                            mass_shift += fixed_modification.MonoisotopicMassShift;
                // variable modification on this residue
                MorpheusModification residue_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length - r + 2, out residue_variable_mod))
                    if (!residue_variable_mod.labile)
                        mass_shift += residue_variable_mod.MonoisotopicMassShift;
                p.cumulativeCTerminalMass[r] = p.cumulativeCTerminalMass[r - 1] + AminoAcidMasses.GetMonoisotopicMass(this[Length - r]) + mass_shift;
            }

            return p;
        }

        public int numVariableMods
        {
            get
            {
                return twoBasedVariableAndLocalizeableModificationss.Count;
            }
        }

        public Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss;
        private PeptideWithPossibleModifications modPep;

        public int MissedCleavages
        {
            get { return modPep.missedCleavages; }
        }

        public override string PeptideDescription
        {
            get { return modPep.PeptideDescription; }
        }

        public double[] FastSortedProductMasses(List<ProductType> productTypes)
        {
            PeptideFragmentMasses p = computeFragmentMasses();
            var PRODUCT_CAPS = ProductCaps.Instance;
            double[] products1 = null;
            double[] products2 = null;
            if (productTypes.Contains(ProductType.b))
                products1 = new double[Length - 2];
            if (productTypes.Contains(ProductType.y))
                products2 = new double[Length - 1];

            int i1 = 0;
            int i2 = 0;
            for (int r = 1; r < Length; r++)
            {
                foreach (ProductType product_type in productTypes)
                {
                    if (!(product_type == ProductType.c && r < Length && this[r] == 'P') &&
                       !(product_type == ProductType.zdot && Length - r < Length && this[Length - r] == 'P') &&
                       !(product_type == ProductType.b && r == 1))
                    {
                        switch (product_type)
                        {
                            case ProductType.adot:
                                throw new NotImplementedException();
                            case ProductType.b:
                                products1[i1] = p.cumulativeNTerminalMass[r] + PRODUCT_CAPS[product_type, MassType.Monoisotopic];
                                i1++;
                                break;

                            case ProductType.c:
                                throw new NotImplementedException();

                            case ProductType.x:
                                throw new NotImplementedException();
                            case ProductType.y:
                                products2[i2] = p.cumulativeCTerminalMass[r] + PRODUCT_CAPS[product_type, MassType.Monoisotopic];
                                i2++;
                                break;

                            case ProductType.zdot:
                                throw new NotImplementedException();
                        }
                    }
                }
            }
            i1 = 0;
            i2 = 0;
            int len = (productTypes.Contains(ProductType.b) ? Length - 2 : 0) +
                      (productTypes.Contains(ProductType.y) ? Length - 1 : 0);
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
                StringBuilder sbsequence = new StringBuilder();
                List<MorpheusModification> value = null;
                // fixed modifications on protein N-terminus
                if (modPep.twoBasedFixedModificationss.TryGetValue(0, out value))
                    foreach (var fixed_modification in value)
                        sbsequence.Append('[' + fixed_modification.cf.Formula + ']');
                // variable modification on protein N-terminus
                MorpheusModification prot_n_term_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(0, out prot_n_term_variable_mod))
                    sbsequence.Append('[' + prot_n_term_variable_mod.cf.Formula + ']');

                // fixed modifications on peptide N-terminus
                if (modPep.twoBasedFixedModificationss.TryGetValue(1, out value))
                    foreach (var fixed_modification in value)
                        sbsequence.Append('[' + fixed_modification.cf.Formula + ']');

                // variable modification on peptide N-terminus
                MorpheusModification pep_n_term_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(1, out pep_n_term_variable_mod))
                    sbsequence.Append('[' + pep_n_term_variable_mod.cf.Formula + ']');

                for (int r = 0; r < Length; r++)
                {
                    sbsequence.Append(this[r]);
                    // fixed modifications on this residue
                    if (modPep.twoBasedFixedModificationss.TryGetValue(r + 2, out value))
                        foreach (var fixed_modification in value)
                            sbsequence.Append('[' + fixed_modification.cf.Formula + ']');
                    // variable modification on this residue
                    MorpheusModification residue_variable_mod;
                    if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(r + 2, out residue_variable_mod))
                        sbsequence.Append('[' + residue_variable_mod.cf.Formula + ']');
                }

                // fixed modifications on peptide C-terminus
                if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 2, out value))
                    foreach (var fixed_modification in value)
                        sbsequence.Append('[' + fixed_modification.cf.Formula + ']');

                // variable modification on peptide C-terminus
                MorpheusModification pep_c_term_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                    sbsequence.Append('[' + pep_c_term_variable_mod.cf.Formula + ']');

                // fixed modifications on protein C-terminus
                if (modPep.twoBasedFixedModificationss.TryGetValue(Length + 3, out value))
                    foreach (var fixed_modification in value)
                        sbsequence.Append('[' + fixed_modification.cf.Formula + ']');
                // variable modification on protein C-terminus
                MorpheusModification prot_c_term_variable_mod;
                if (twoBasedVariableAndLocalizeableModificationss.TryGetValue(Length + 3, out prot_c_term_variable_mod))
                    sbsequence.Append('[' + prot_c_term_variable_mod.cf.Formula + ']');
                return sbsequence.ToString();
            }
        }

        public override bool Equals(object obj)
        {
            PeptideWithSetModifications q = obj as PeptideWithSetModifications;
            return q != null && q.Sequence.Equals(Sequence) && q.OneBasedStartResidueInProtein == OneBasedStartResidueInProtein && q.protein.Equals(protein);
        }

        public override int GetHashCode()
        {
            return ExtendedSequence.GetHashCode();
        }
    }
}