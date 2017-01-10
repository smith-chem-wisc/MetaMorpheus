//using System.Collections.Generic;
//using System.Linq;
//using System.Text;

//namespace MetaMorpheus
//{
//    public class PeptideWithSetModificationsHeavy : PeptideWithSetModifications
//    {
//private PeptideFragmentMasses peptideFragmentMasses;
//        public PeptideWithSetModificationsHeavy(PeptideWithPossibleModifications modPep, Dictionary<int, List<Modification>> twoBasedfixedModifications, Dictionary<int, Modification> twoBasedVariableModificationss)
//            : base(modPep.protein, modPep.OneBasedStartResidueInProtein, modPep.OneBasedEndResidueInProtein, modPep.PeptideDescription + "set", modPep.MissedCleavages)
//        {
//            this.twoBasedfixedModifications = twoBasedfixedModifications;
//            this.twoBasedVariableModificationss = twoBasedVariableModificationss;

//            double mass_shift;
//            cumulativeNTerminalMass = new double[Length];
//            mass_shift = 0.0;
//             fixed modifications on protein N-terminus
//            if (twoBasedfixedModifications != null)
//            {
//                List<Modification> prot_n_term_fixed_mods;
//                if (twoBasedfixedModifications.TryGetValue(0, out prot_n_term_fixed_mods))
//                {
//                    foreach (Modification fixed_modification in prot_n_term_fixed_mods)
//                    {
//                        if (!fixed_modification.labile)
//                            mass_shift += fixed_modification.MonoisotopicMassShift;
//                    }
//                }
//            }
//             variable modification on the protein N-terminus
//            if (twoBasedVariableModificationss != null)
//            {
//                Modification prot_n_term_var_mod;
//                if (twoBasedVariableModificationss.TryGetValue(0, out prot_n_term_var_mod))
//                {
//                    if (!prot_n_term_var_mod.labile)
//                        mass_shift += prot_n_term_var_mod.MonoisotopicMassShift;
//                }
//            }
//             fixed modifications on peptide N-terminus
//            if (twoBasedfixedModifications != null)
//            {
//                List<Modification> pep_n_term_fixed_mods;
//                if (twoBasedfixedModifications.TryGetValue(1, out pep_n_term_fixed_mods))
//                {
//                    foreach (Modification fixed_modification in pep_n_term_fixed_mods)
//                    {
//                        if (!fixed_modification.labile)
//                            mass_shift += fixed_modification.MonoisotopicMassShift;
//                    }
//                }
//            }
//             variable modification on peptide N-terminus
//            if (twoBasedVariableModificationss != null)
//            {
//                Modification pep_n_term_variable_mod;
//                if (twoBasedVariableModificationss.TryGetValue(1, out pep_n_term_variable_mod))
//                {
//                    if (!pep_n_term_variable_mod.labile)
//                        mass_shift += pep_n_term_variable_mod.MonoisotopicMassShift;
//                }
//            }
//            cumulativeNTerminalMass[0] = mass_shift;

//            for (int r = 1; r < Length; r++)
//            {
//                mass_shift = 0.0;
//                 fixed modifications on this residue
//                if (twoBasedfixedModifications != null)
//                {
//                    List<Modification> residue_fixed_mods;
//                    if (twoBasedfixedModifications.TryGetValue(r + 1, out residue_fixed_mods))
//                    {
//                        foreach (Modification fixed_modification in residue_fixed_mods)
//                        {
//                            if (!fixed_modification.labile)
//                                mass_shift += fixed_modification.MonoisotopicMassShift;
//                        }
//                    }
//                }
//                 variable modification on this residue
//                if (twoBasedVariableModificationss != null)
//                {
//                    Modification residue_variable_mod;
//                    if (twoBasedVariableModificationss.TryGetValue(r + 1, out residue_variable_mod))
//                    {
//                        if (!residue_variable_mod.labile)
//                            mass_shift += residue_variable_mod.MonoisotopicMassShift;
//                    }
//                }
//                cumulativeNTerminalMass[r] = cumulativeNTerminalMass[r - 1] + AminoAcidMasses.GetMonoisotopicMass(this[r - 1]) + mass_shift;
//            }

//            cumulativeCTerminalMass = new double[Length];

//            mass_shift = 0.0;
//             fixed modifications on protein C-terminus
//            if (twoBasedfixedModifications != null)
//            {
//                List<Modification> prot_c_term_fixed_mods;
//                if (twoBasedfixedModifications.TryGetValue(Length + 3, out prot_c_term_fixed_mods))
//                {
//                    foreach (Modification fixed_modification in prot_c_term_fixed_mods)
//                    {
//                        if (!fixed_modification.labile)
//                            mass_shift += fixed_modification.MonoisotopicMassShift;
//                    }
//                }
//            }
//             variable modification on protein C-terminus
//            if (twoBasedVariableModificationss != null)
//            {
//                Modification prot_c_term_variable_mod;
//                if (twoBasedVariableModificationss.TryGetValue(Length + 3, out prot_c_term_variable_mod))
//                {
//                    if (!prot_c_term_variable_mod.labile)
//                        mass_shift += prot_c_term_variable_mod.MonoisotopicMassShift;
//                }
//            }
//             fixed modifications on peptide C-terminus
//            if (twoBasedfixedModifications != null)
//            {
//                List<Modification> pep_c_term_fixed_mods;
//                if (twoBasedfixedModifications.TryGetValue(Length + 2, out pep_c_term_fixed_mods))
//                {
//                    foreach (Modification fixed_modification in pep_c_term_fixed_mods)
//                    {
//                        if (!fixed_modification.labile)
//                            mass_shift += fixed_modification.MonoisotopicMassShift;
//                    }
//                }
//            }
//             variable modification on peptide C-terminus
//            if (twoBasedVariableModificationss != null)
//            {
//                Modification pep_c_term_variable_mod;
//                if (twoBasedVariableModificationss.TryGetValue(Length + 2, out pep_c_term_variable_mod))
//                {
//                    if (!pep_c_term_variable_mod.labile)
//                        mass_shift += pep_c_term_variable_mod.MonoisotopicMassShift;
//                }
//            }
//            cumulativeCTerminalMass[0] = mass_shift;

//            for (int r = 1; r < Length; r++)
//            {
//                mass_shift = 0.0;
//                 fixed modifications on this residue
//                if (twoBasedfixedModifications != null)
//                {
//                    List<Modification> residue_fixed_mods;
//                    if (twoBasedfixedModifications.TryGetValue(Length - r + 2, out residue_fixed_mods))
//                    {
//                        foreach (Modification fixed_modification in residue_fixed_mods)
//                        {
//                            if (!fixed_modification.labile)
//                                mass_shift += fixed_modification.MonoisotopicMassShift;
//                        }
//                    }
//                }
//                 variable modification on this residue
//                if (twoBasedVariableModificationss != null)
//                {
//                    Modification residue_variable_mod;
//                    if (twoBasedVariableModificationss.TryGetValue(Length - r + 2, out residue_variable_mod))
//                    {
//                        if (!residue_variable_mod.labile)
//                            mass_shift += residue_variable_mod.MonoisotopicMassShift;
//                    }
//                }
//                cumulativeCTerminalMass[r] = cumulativeCTerminalMass[r - 1] + AminoAcidMasses.GetMonoisotopicMass(this[Length - r]) + mass_shift;
//            }

//             ADD EVEN THE LABILE ONES!!!
//            MonoisotopicMass = Constants.WATER_MONOISOTOPIC_MASS;

//            foreach (char amino_acid in BaseSequence)
//                MonoisotopicMass += AminoAcidMasses.GetMonoisotopicMass(amino_acid);
//            if (twoBasedfixedModifications != null)
//                foreach (List<Modification> fixed_modifications in twoBasedfixedModifications.Values)
//                    foreach (Modification fixed_modification in fixed_modifications)
//                        MonoisotopicMass += fixed_modification.MonoisotopicMassShift;
//            if (twoBasedVariableModificationss != null)
//                foreach (Modification variable_modification in twoBasedVariableModificationss.Values)
//                    MonoisotopicMass += variable_modification.MonoisotopicMassShift;
//        }

//        private string extendedSequence;

//        public string ExtendedSequence
//        {
//            get
//            {
//                if (extendedSequence == null)
//                {
//                    StringBuilder sequence = new StringBuilder();
//                    sequence.Append(PreviousAminoAcid);
//                    sequence.Append(".");
//                    sequence.Append(Sequence);
//                    sequence.Append(".");
//                    sequence.Append(NextAminoAcid);
//                    extendedSequence = sequence.ToString();
//                }
//                return extendedSequence;
//            }
//        }

//        private string sequence;

//        public string Sequence
//        {
//            get
//            {
//                if (sequence == null)
//                {
//                    StringBuilder sbsequence = new StringBuilder();
//                     fixed modifications on protein N-terminus
//                    if (twoBasedfixedModifications != null)
//                    {
//                        List<Modification> prot_n_term_fixed_mods;
//                        if (twoBasedfixedModifications.TryGetValue(0, out prot_n_term_fixed_mods))
//                        {
//                            foreach (Modification fixed_modification in prot_n_term_fixed_mods)
//                            {
//                                sbsequence.Append('[' + fixed_modification.Description + ']');
//                            }
//                        }
//                    }
//                     variable modification on protein N-terminus
//                    if (twoBasedVariableModificationss != null)
//                    {
//                        Modification prot_n_term_variable_mod;
//                        if (twoBasedVariableModificationss.TryGetValue(0, out prot_n_term_variable_mod))
//                        {
//                            sbsequence.Append('(' + prot_n_term_variable_mod.Description + ')');
//                        }
//                    }

//                     fixed modifications on peptide N-terminus
//                    if (twoBasedfixedModifications != null)
//                    {
//                        List<Modification> pep_n_term_fixed_mods;
//                        if (twoBasedfixedModifications.TryGetValue(1, out pep_n_term_fixed_mods))
//                        {
//                            foreach (Modification fixed_modification in pep_n_term_fixed_mods)
//                            {
//                                sbsequence.Append('[' + fixed_modification.Description + ']');
//                            }
//                        }
//                    }
//                     variable modification on peptide N-terminus
//                    if (twoBasedVariableModificationss != null)
//                    {
//                        Modification pep_n_term_variable_mod;
//                        if (twoBasedVariableModificationss.TryGetValue(1, out pep_n_term_variable_mod))
//                        {
//                            sbsequence.Append('(' + pep_n_term_variable_mod.Description + ')');
//                        }
//                    }

//                    for (int r = 0; r < Length; r++)
//                    {
//                        sbsequence.Append(this[r]);
//                         fixed modifications on this residue
//                        if (twoBasedfixedModifications != null)
//                        {
//                            List<Modification> residue_fixed_mods;
//                            if (twoBasedfixedModifications.TryGetValue(r + 2, out residue_fixed_mods))
//                            {
//                                foreach (Modification fixed_modification in residue_fixed_mods)
//                                {
//                                    sbsequence.Append('[' + fixed_modification.Description + ']');
//                                }
//                            }
//                        }
//                         variable modification on this residue
//                        if (twoBasedVariableModificationss != null)
//                        {
//                            Modification residue_variable_mod;
//                            if (twoBasedVariableModificationss.TryGetValue(r + 2, out residue_variable_mod))
//                            {
//                                sbsequence.Append('(' + residue_variable_mod.Description + ')');
//                            }
//                        }
//                    }

//                     fixed modifications on peptide C-terminus
//                    if (twoBasedfixedModifications != null)
//                    {
//                        List<Modification> pep_c_term_fixed_mods;
//                        if (twoBasedfixedModifications.TryGetValue(Length + 2, out pep_c_term_fixed_mods))
//                        {
//                            foreach (Modification fixed_modification in pep_c_term_fixed_mods)
//                            {
//                                sbsequence.Append('[' + fixed_modification.Description + ']');
//                            }
//                        }
//                    }
//                     variable modification on peptide C-terminus
//                    if (twoBasedVariableModificationss != null)
//                    {
//                        Modification pep_c_term_variable_mod;
//                        if (twoBasedVariableModificationss.TryGetValue(Length + 2, out pep_c_term_variable_mod))
//                        {
//                            sbsequence.Append('(' + pep_c_term_variable_mod.Description + ')');
//                        }
//                    }

//                     fixed modifications on protein C-terminus
//                    if (twoBasedfixedModifications != null)
//                    {
//                        List<Modification> prot_c_term_fixed_mods;
//                        if (twoBasedfixedModifications.TryGetValue(Length + 3, out prot_c_term_fixed_mods))
//                        {
//                            foreach (Modification fixed_modification in prot_c_term_fixed_mods)
//                            {
//                                sbsequence.Append('[' + fixed_modification.Description + ']');
//                            }
//                        }
//                    }
//                     variable modification on protein C-terminus
//                    if (twoBasedVariableModificationss != null)
//                    {
//                        Modification prot_c_term_variable_mod;
//                        if (twoBasedVariableModificationss.TryGetValue(Length + 3, out prot_c_term_variable_mod))
//                        {
//                            sbsequence.Append('(' + prot_c_term_variable_mod.Description + ')');
//                        }
//                    }
//                    sequence = sbsequence.ToString();
//                }
//                return sequence;
//            }
//        }

//        public double MonoisotopicMass { get; private set; }

//        public IEnumerable<KeyValuePair<int, Modification>> Modifications
//        {
//            get
//            {
//                foreach (var ye in twoBasedVariableModificationss)
//                    yield return ye;
//                if (twoBasedfixedModifications != null)
//                    foreach (var heh in twoBasedfixedModifications)
//                        foreach (var hm in heh.Value)
//                            yield return new KeyValuePair<int, Modification>(heh.Key, hm);
//            }
//        }

//        public int numVariableMods
//        {
//            get
//            {
//                return twoBasedVariableModificationss.Count;
//            }
//        }

//        private double[] cumulativeNTerminalMass;
//        private double[] cumulativeCTerminalMass;

//        private static readonly ProductCaps PRODUCT_CAPS = ProductCaps.Instance;
//        private Dictionary<int, List<Modification>> twoBasedfixedModifications;
//        private Dictionary<int, Modification> twoBasedVariableModificationss;

//        private Dictionary<List<ProductType>, double[]> productMassesSingleArray = new Dictionary<List<ProductType>, double[]>();

//        public IEnumerable<double> FastProductMasses(List<ProductType> productTypes)
//        {
//            int len = (productTypes.Contains(ProductType.b) ? Length - 1 : 0) +
//                      (productTypes.Contains(ProductType.y) ? Length - 1 : 0);

//            foreach (ProductType product_type in productTypes)
//            {
//                for (int r = 1; r < Length; r++)
//                {
//                    if (!(product_type == ProductType.c && r < Length && this[r] == 'P') &&
//                       !(product_type == ProductType.zdot && Length - r < Length && this[Length - r] == 'P'))
//                    {
//                        switch (product_type)
//                        {
//                            case ProductType.adot:
//                            case ProductType.b:
//                            case ProductType.c:
//                                yield return cumulativeNTerminalMass[r] + PRODUCT_CAPS[product_type, MassType.Monoisotopic];
//                                break;

//                            case ProductType.x:
//                            case ProductType.y:
//                            case ProductType.zdot:
//                                yield return cumulativeCTerminalMass[r] + PRODUCT_CAPS[product_type, MassType.Monoisotopic];
//                                break;
//                        }
//                    }
//                }
//            }
//        }

//        public double[] GetProductMassesSingleArray(List<ProductType> productTypes)
//        {
//            if (!productMassesSingleArray.ContainsKey(productTypes))
//            {
//                int len = (productTypes.Contains(ProductType.b) ? Length - 1 : 0) +
//                          (productTypes.Contains(ProductType.y) ? Length - 1 : 0);

//                List<double> product_masses = new List<double>();

//                foreach (ProductType product_type in productTypes)
//                {
//                    double[] product_masses_for_product_type = new double[Length - 1];
//                    for (int r = 1; r < Length; r++)
//                    {
//                        if (!(product_type == ProductType.c && r < Length && this[r] == 'P') &&
//                           !(product_type == ProductType.zdot && Length - r < Length && this[Length - r] == 'P'))
//                        {
//                            switch (product_type)
//                            {
//                                case ProductType.adot:
//                                case ProductType.b:
//                                case ProductType.c:
//                                    product_masses.Add(cumulativeNTerminalMass[r] + PRODUCT_CAPS[product_type, MassType.Monoisotopic]);
//                                    break;

//                                case ProductType.x:
//                                case ProductType.y:
//                                case ProductType.zdot:
//                                    product_masses.Add(cumulativeCTerminalMass[r] + PRODUCT_CAPS[product_type, MassType.Monoisotopic]);
//                                    break;
//                            }
//                        }
//                    }
//                }
//                IEqualityComparer<double> comparer = new MyDoubleEqualityComparer(0.01);
//                productMassesSingleArray[productTypes] = product_masses.OrderBy(b => b).Distinct(comparer).ToArray();
//            }
//            return productMassesSingleArray[productTypes];
//        }

//        public PeptideWithSetModifications(PeptideWithSetModifications unmodified, int oneBasedModIndex, double modMass)
//            : base(unmodified.protein, unmodified.OneBasedStartResidueInProtein, unmodified.OneBasedEndResidueInProtein, unmodified.PeptideDescription + "ind" + oneBasedModIndex + "mass" + modMass, unmodified.MissedCleavages)
//        {
//             ONLY NEED TO SET CORRECT cumulativeNTerminalMass and cumulativeCTerminalMass!!!
//            cumulativeCTerminalMass = new double[Length];
//            cumulativeNTerminalMass = new double[Length];

//             Always at least one here!!!
//            for (int i = 0; i < oneBasedModIndex; i++)
//            {
//                cumulativeNTerminalMass[i] = unmodified.cumulativeNTerminalMass[i];
//            }
//             Might not have any here if oneBasedModIndex==Length
//            for (int i = oneBasedModIndex; i < Length; i++)
//            {
//                cumulativeNTerminalMass[i] = unmodified.cumulativeNTerminalMass[i] + modMass;
//            }
//             Always at least one here!!!
//            for (int i = 0; i < Length - oneBasedModIndex + 1; i++)
//            {
//                cumulativeCTerminalMass[i] = unmodified.cumulativeCTerminalMass[i];
//            }
//             Might not have any here if oneBasedModIndex==1
//            for (int i = Length - oneBasedModIndex + 1; i < Length; i++)
//            {
//                cumulativeCTerminalMass[i] = unmodified.cumulativeCTerminalMass[i] + modMass;
//            }
//        }

//        public PeptideWithSetModificationsHeavy(string v, string b) : base(null, 0, 0, null, 0)
//        {
//            this.extendedSequence = v;
//            this.baseSequence = b;
//        }

//public PeptideWithSetModifications(PeptideWithSetModifications unmodified, int oneBasedModIndex, double modMass)
//            : base(unmodified.protein, unmodified.OneBasedStartResidueInProtein, unmodified.OneBasedEndResidueInProtein, unmodified.PeptideDescription + "ind" + oneBasedModIndex + "mass" + modMass, unmodified.MissedCleavages)
//        {
//    // ONLY NEED TO SET CORRECT cumulativeNTerminalMass and cumulativeCTerminalMass!!!
//    cumulativeCTerminalMass = new double[Length];
//    cumulativeNTerminalMass = new double[Length];

//    // Always at least one here!!!
//    for (int i = 0; i < oneBasedModIndex; i++)
//    {
//        cumulativeNTerminalMass[i] = unmodified.cumulativeNTerminalMass[i];
//    }
//    // Might not have any here if oneBasedModIndex==Length
//    for (int i = oneBasedModIndex; i < Length; i++)
//    {
//        cumulativeNTerminalMass[i] = unmodified.cumulativeNTerminalMass[i] + modMass;
//    }
//    // Always at least one here!!!
//    for (int i = 0; i < Length - oneBasedModIndex + 1; i++)
//    {
//        cumulativeCTerminalMass[i] = unmodified.cumulativeCTerminalMass[i];
//    }
//    // Might not have any here if oneBasedModIndex==1
//    for (int i = Length - oneBasedModIndex + 1; i < Length; i++)
//    {
//        cumulativeCTerminalMass[i] = unmodified.cumulativeCTerminalMass[i] + modMass;
//    }
//}

//public double[] GetProductMassesSingleArray(List<ProductType> productTypes)
//{
//    if (!productMassesSingleArray.ContainsKey(productTypes))
//    {
//        int len = (productTypes.Contains(ProductType.b) ? Length - 1 : 0) +
//                  (productTypes.Contains(ProductType.y) ? Length - 1 : 0);

//        List<double> product_masses = new List<double>();

//        foreach (ProductType product_type in productTypes)
//        {
//            double[] product_masses_for_product_type = new double[Length - 1];
//            for (int r = 1; r < Length; r++)
//            {
//                if (!(product_type == ProductType.c && r < Length && this[r] == 'P') &&
//                   !(product_type == ProductType.zdot && Length - r < Length && this[Length - r] == 'P'))
//                {
//                    switch (product_type)
//                    {
//                        case ProductType.adot:
//                        case ProductType.b:
//                        case ProductType.c:
//                            product_masses.Add(cumulativeNTerminalMass[r] + PRODUCT_CAPS[product_type, MassType.Monoisotopic]);
//                            break;

//                        case ProductType.x:
//                        case ProductType.y:
//                        case ProductType.zdot:
//                            product_masses.Add(cumulativeCTerminalMass[r] + PRODUCT_CAPS[product_type, MassType.Monoisotopic]);
//                            break;
//                    }
//                }
//            }
//        }
//        IEqualityComparer<double> comparer = new MyDoubleEqualityComparer(0.01);
//        productMassesSingleArray[productTypes] = product_masses.OrderBy(b => b).Distinct(comparer).ToArray();
//    }
//    return productMassesSingleArray[productTypes];
////}

//public PeptideWithSetModifications(PeptideWithSetModifications unmodified, int oneBasedModIndex, double modMass)
//                : base(unmodified.protein, unmodified.OneBasedStartResidueInProtein, unmodified.OneBasedEndResidueInProtein, unmodified.PeptideDescription + "ind" + oneBasedModIndex + "mass" + modMass, unmodified.MissedCleavages)
//        {
//    // ONLY NEED TO SET CORRECT cumulativeNTerminalMass and cumulativeCTerminalMass!!!
//    cumulativeCTerminalMass = new double[Length];
//    cumulativeNTerminalMass = new double[Length];

//    // Always at least one here!!!
//    for (int i = 0; i < oneBasedModIndex; i++)
//    {
//        cumulativeNTerminalMass[i] = unmodified.cumulativeNTerminalMass[i];
//    }
//    // Might not have any here if oneBasedModIndex==Length
//    for (int i = oneBasedModIndex; i < Length; i++)
//    {
//        cumulativeNTerminalMass[i] = unmodified.cumulativeNTerminalMass[i] + modMass;
//    }
//    // Always at least one here!!!
//    for (int i = 0; i < Length - oneBasedModIndex + 1; i++)
//    {
//        cumulativeCTerminalMass[i] = unmodified.cumulativeCTerminalMass[i];
//    }
//    // Might not have any here if oneBasedModIndex==1
//    for (int i = Length - oneBasedModIndex + 1; i < Length; i++)
//    {
//        cumulativeCTerminalMass[i] = unmodified.cumulativeCTerminalMass[i] + modMass;
//    }
//}
//    }
//}