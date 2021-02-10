using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;

namespace EngineLayer.GlycoSearch
{
    public static class GlycoPeptides
    {
        //Motif should be writen with required form
        public static List<int> GetPossibleModSites(PeptideWithSetModifications peptide, string[] motifs)
        {
            List<int> possibleModSites = new List<int>();

            List<Modification> modifications = new List<Modification>();

            foreach (var mtf in motifs)
            {
                if (ModificationMotif.TryGetMotif(mtf, out ModificationMotif aMotif))
                {
                    Modification modWithMotif = new Modification(_target: aMotif, _locationRestriction: "Anywhere.");
                    modifications.Add(modWithMotif);
                }
            }

            foreach (var modWithMotif in modifications)
            {
                for (int r = 0; r < peptide.Length; r++)
                {
                    if (peptide.AllModsOneIsNterminus.Keys.Contains(r + 2))
                    {
                        continue;
                    }

                    //FullSequence is used here to avoid duplicated modification on same sites?
                    if (ModificationLocalization.ModFits(modWithMotif, peptide.BaseSequence, r + 1, peptide.Length, r + 1))
                    {
                        possibleModSites.Add(r + 2);
                    }
                }
            }

            return possibleModSites;
        }

        public static bool MotifExist(string baseSeq, string[] motifs)
        {
            List<Modification> modifications = new List<Modification>();

            foreach (var mtf in motifs)
            {
                if (ModificationMotif.TryGetMotif(mtf, out ModificationMotif aMotif))
                {
                    Modification modWithMotif = new Modification(_target: aMotif, _locationRestriction: "Anywhere.");
                    modifications.Add(modWithMotif);
                }
            }

            foreach (var modWithMotif in modifications)
            {
                for (int r = 0; r < baseSeq.Length; r++)
                {
                    //Modification is not considered.                  
                    if (ModificationLocalization.ModFits(modWithMotif, baseSeq, r + 1, baseSeq.Length, r + 1))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        public static double[] ScanOxoniumIonFilter(Ms2ScanWithSpecificMass theScan, MassDiffAcceptor massDiffAcceptor, DissociationType dissociationType)
        {
            double[] oxoniumIonsintensities = new double[Glycan.AllOxoniumIons.Length];

            if (dissociationType != DissociationType.HCD && dissociationType != DissociationType.CID && dissociationType != DissociationType.EThcD)
            {
                return oxoniumIonsintensities;
            }

            for (int i = 0; i < Glycan.AllOxoniumIons.Length; i++)
            {
                var oxoMass = ((double)Glycan.AllOxoniumIons[i] / 1E5).ToMass(1);
                var envelope = theScan.GetClosestExperimentalIsotopicEnvelope(oxoMass);
                if (massDiffAcceptor.Accepts(envelope.MonoisotopicMass, oxoMass) >= 0)
                {
                    oxoniumIonsintensities[i] = envelope.TotalIntensity;
                }

            }

            //Normalize by 204. What will happen if 204 is 0.
            if (oxoniumIonsintensities[9] != 0)
            {
                var x204 = oxoniumIonsintensities[9];
                for (int i = 0; i < Glycan.AllOxoniumIons.Length; i++)
                {

                    oxoniumIonsintensities[i] = oxoniumIonsintensities[i] / x204;
                }
            }

            return oxoniumIonsintensities;
        }

        //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons. 
        //Rules are coded in the function.    
        public static bool OxoniumIonsAnalysis(double[] oxoniumIonsintensities, GlycanBox glycanBox)
        {
            //If a glycopeptide spectrum does not have 292.1027 or 274.0921, then remove all glycans that have NeuAc from the search.
            if (oxoniumIonsintensities[10] <= 0 && oxoniumIonsintensities[12] <= 0)
            {
                if (glycanBox.Kind[2] != 0 )
                {
                    return false;
                }
            }

            //If a glycopeptide spectrum does not have 290.0870 or 308.0976, then remove all glycans that have NeuGc from the search.
            if (oxoniumIonsintensities[11] <= 0 && oxoniumIonsintensities[13] <= 0)
            {
                if (glycanBox.Kind[3] != 0)
                {
                    return false;
                }
            }

            //pGlyco3 paper also mentioned the ion 243.0264304315 which is important for phospho-Hex (phoH) exist.

            //If a spectrum has 366.1395, remove glycans that do not have HexNAc(1)Hex(1) or more. Here use the total glycan of glycanBox to calculate. 
            if (oxoniumIonsintensities[14] > 0)
            {
                if (glycanBox.Kind[0] < 1 && glycanBox.Kind[1] < 1)
                {
                    return false;
                }
            }

            //Other rules:
            //A spectrum needs to have 204.0867 to be considered as a glycopeptide.              
            //Ratio of 138.055 to 144.0655 can seperate O/N glycan.

            return true;
        }

        //Get specific glycan Y ions such 'NN' or 'NH'.
        public static Product GetIndicatorYIon(double peptideMonomassWithNoGlycan, string glycanString)
        {
            Product product = new Product(ProductType.M, FragmentationTerminus.Both, peptideMonomassWithNoGlycan + (double)Glycan.GetMass(glycanString) / 1E5, 0, 0, 0);
            return product;
        }

        public static bool MatchIndicatorYIon(Ms2ScanWithSpecificMass scan, Product theoreticalProduct, CommonParameters commonParameters)
        {
            List<Product> products = new List<Product>();
            products.Add(theoreticalProduct);
            var x = MetaMorpheusEngine.MatchFragmentIons(scan, products, commonParameters).Count();

            foreach (var childScan in scan.ChildScans)
            {
                x += MetaMorpheusEngine.MatchFragmentIons(childScan, products, commonParameters).Count();
            }

            return x > 0;
        }

        //Find Glycan index or glycanBox index.
        public static int BinarySearchGetIndex(double[] massArray, double targetMass)
        {
            var iD = Array.BinarySearch(massArray, targetMass);
            if (iD < 0) { iD = ~iD; }
            else
            {
                while (iD - 1 >= 0 && massArray[iD - 1] >= targetMass - 0.00000001)
                {
                    iD--;
                }
            }
            return iD;
        }

        public static bool DissociationTypeContainETD(DissociationType dissociationType)
        {
            if (dissociationType == DissociationType.ETD || dissociationType == DissociationType.EThcD)
            {
                return true;
            }

            return false;
        }

        public static bool DissociationTypeContainHCD(DissociationType dissociationType)
        {
            if (dissociationType == DissociationType.HCD || dissociationType == DissociationType.CID || dissociationType == DissociationType.EThcD)
            {
                return true;
            }

            return false;
        }

        #region Glycopeptide and fragmentation functions

        //TO THINK: if one peptide modified with >1 nglycans, how Yions will be generated. Here we just simplified to treat it as two seperate generation.
        public static List<Product> GetGlycanYIons(PeptideWithSetModifications peptide, Glycan[] glycans)
        {
            double possiblePeptideMass = peptide.MonoisotopicMass;
            List<Product> YIons = new List<Product>();
            HashSet<int> masses = new HashSet<int>();

            foreach (var glycan in glycans)
            {
                foreach (var ion in glycan.Ions)
                {
                    if (masses.Contains(ion.IonMass))
                    {
                        continue;
                    }
                    Product product = new Product(ProductType.M, FragmentationTerminus.Both, possiblePeptideMass + (double)ion.IonMass / 1E5, 0, 0, (double)ion.IonMass / 1E5);
                    YIons.Add(product);
                    masses.Add(ion.IonMass);
                }
            }
            return YIons;
        }

        public static List<Product> GetOxoniumIons(Glycan[] glycans)
        {
            List<Product> OxoIons = new List<Product>();

            HashSet<int> masses = new HashSet<int>();

            foreach (var g in glycans)
            {
                foreach (var d in g.DiagnosticIons)
                {
                    if (masses.Contains(d))
                    {
                        continue;
                    }
                    masses.Add(d);

                    var dmass = (double)d / 1E5;
                    Product D = new Product(ProductType.D, FragmentationTerminus.None, dmass, (int)Math.Round(dmass.ToMz(1), 0), 0, 0);
                    OxoIons.Add(D);
                }
            }
            return OxoIons;
        }

        public static PeptideWithSetModifications GlyGetTheoreticalPeptide(int[] modPos, PeptideWithSetModifications peptide, Modification[] modifications)
        {
            if (modPos.Length == 0)
            {
                return null;
            }

            Dictionary<int, Modification> testMods = new Dictionary<int, Modification>();
            for (int i = 0; i < modPos.Length; i++)
            {
                testMods.Add(modPos[i], modifications[i]);
            }

            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                testMods.Add(mod.Key, mod.Value);
            }

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

            return testPeptide;
        }

        public static PeptideWithSetModifications GlyGetTheoreticalPeptide(int[] modPos, PeptideWithSetModifications peptide, GlycanBox glycanBox, Modification[] globalMods)
        {
            Modification[] modifications = new Modification[glycanBox.ModCount];
            for (int i = 0; i < glycanBox.ModCount; i++)
            {
                modifications[i] = globalMods[glycanBox.ModIds.ElementAt(i)];
            }

            return GlyGetTheoreticalPeptide(modPos, peptide, modifications);
        }

        public static PeptideWithSetModifications GlyGetTheoreticalPeptide(Route theModPositions, PeptideWithSetModifications peptide, Modification[] globalMods)
        {
            Modification[] modifications = new Modification[theModPositions.Mods.Count];
            for (int i = 0; i < theModPositions.Mods.Count; i++)
            {
                modifications[i] = globalMods[theModPositions.Mods[i].GlycanID];
            }

            return GlyGetTheoreticalPeptide(theModPositions.Mods.Select(p=>p.ModSite).ToArray(), peptide, modifications);
        }

        public static List<Product> GlyGetTheoreticalFragments(GlycoType glycanType, DissociationType dissociationType, PeptideWithSetModifications peptide, PeptideWithSetModifications modPeptide, List<int> NPos, Glycan[] glycans)
        {
            List<Product> theoreticalProducts = new List<Product>();
            HashSet<double> masses = new HashSet<double>();
            List<Product> products = new List<Product>();

            if (dissociationType == DissociationType.HCD || dissociationType == DissociationType.CID)
            {
                peptide.Fragment(dissociationType, FragmentationTerminus.Both, products);

                if (glycanType == GlycoType.NGlycoPep || glycanType == GlycoType.MixedGlycoPep)
                {
                    List<Product> shiftProducts = new List<Product>();

                    foreach (var pd in products)
                    {
                        if (pd.ProductType == ProductType.b)
                        {
                            int count = NPos.Count(p => pd.AminoAcidPosition >= p - 1 );
                            if (count > 0)
                            {
                                Product b = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass + count * Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, count * Glycan.HexNAcMass);
                                shiftProducts.Add(b);
                            }
                        }
                        else if (pd.ProductType == ProductType.y)
                        {
                            int count = NPos.Count(p => pd.AminoAcidPosition < p );
                            if (count > 0)
                            {
                                Product y = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass + count * Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, count * Glycan.HexNAcMass);
                                shiftProducts.Add(y);
                            }
                        }
                    }

                    products.AddRange(shiftProducts);
                }

                products.AddRange(GetOxoniumIons(glycans));

                products.AddRange(GetGlycanYIons(peptide, glycans));

            }
            else if (dissociationType == DissociationType.ETD)
            {
                modPeptide.Fragment(dissociationType, FragmentationTerminus.Both, products);
            }
            else if (dissociationType == DissociationType.EThcD)
            {
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);

                if (glycanType == GlycoType.NGlycoPep || glycanType == GlycoType.MixedGlycoPep)
                {
                    List<Product> shiftProducts = new List<Product>();

                    foreach (var pd in products)
                    {
                        if (pd.ProductType == ProductType.b)
                        {
                            int count = NPos.Count(p => p > pd.AminoAcidPosition);
                            if (count > 0)
                            {
                                Product b = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass + count * Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, count * Glycan.HexNAcMass);
                                shiftProducts.Add(b);
                            }
                        }
                        else if (pd.ProductType == ProductType.y)
                        {
                            int count = NPos.Count(p => p < pd.AminoAcidPosition);
                            if (count > 0)
                            {
                                Product y = new Product(pd.ProductType, pd.Terminus, pd.NeutralMass + count * Glycan.HexNAcMass, pd.FragmentNumber, pd.AminoAcidPosition, count * Glycan.HexNAcMass);
                                shiftProducts.Add(y);
                            }
                        }
                    }

                    products.AddRange(shiftProducts);
                }

                List<Product> etdProducts = new List<Product>();

                modPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, etdProducts);

                products.AddRange(etdProducts.Where(p => p.ProductType != ProductType.y));

                products.AddRange(GetOxoniumIons(glycans));

                products.AddRange(GetGlycanYIons(peptide, glycans));

            }          

            foreach (var fragment in products)
            {
                if (!masses.Contains(fragment.NeutralMass))
                {
                    masses.Add(fragment.NeutralMass);
                    theoreticalProducts.Add(fragment);
                }
            }

            return theoreticalProducts;
        }

        //Find FragmentMass for the fragments that contain localization Information.
        public static List<double> GetLocalFragmentGlycan(List<Product> products, int modInd, int childBoxInd, LocalizationGraph localizationGraph)
        {
            List<double> newFragments = new List<double>();

            var local_c_fragments = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] - 1 && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1] - 1).ToList();
            foreach (var c in local_c_fragments)
            {
                var newMass = c.NeutralMass + localizationGraph.ChildModBoxes[childBoxInd].Mass;
                newFragments.Add(newMass);
            }

            var local_z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1]).ToList();

            foreach (var z in local_z_fragments)
            {
                var newMass = z.NeutralMass + (localizationGraph.ModBox.Mass - localizationGraph.ChildModBoxes[childBoxInd].Mass);
                newFragments.Add(newMass);
            }

            var local_b_fragments = products.Where(p => p.ProductType == ProductType.b && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] - 1 && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1] - 1).ToList();
            foreach (var b in local_b_fragments)
            {
                newFragments.Add(b.NeutralMass);
                var newMass = b.NeutralMass + ((GlycanBox)localizationGraph.ChildModBoxes[childBoxInd]).NGlycanCount * Glycan.HexNAcMass;
                newFragments.Add(newMass);
            }

            var local_y_fragments = products.Where(p => p.ProductType == ProductType.y && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1]).ToList();
            foreach (var y in local_y_fragments)
            {
                newFragments.Add(y.NeutralMass);
                var newMass = y.NeutralMass + (((GlycanBox)localizationGraph.ModBox).NGlycanCount - ((GlycanBox)localizationGraph.ChildModBoxes[childBoxInd]).NGlycanCount) * Glycan.HexNAcMass;
                newFragments.Add(newMass);
            }

            return newFragments;
        }

        //Find FragmentMass for the fragments that doesn't contain localization Information. For example, "A|TAABBS|B", c1 and c7, z1 and z7, z8 ion don't contain localization information.
        public static List<double> GetUnlocalFragmentGlycan(List<Product> products, int[] modPoses, ModBox modBox)
        {
            var mass = modBox.Mass;

            List<double> newFragments = new List<double>();

            var c_fragments = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition < modPoses.First() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(c_fragments);

            var c_fragments_shift = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition >= modPoses.Last() - 1).Select(p => p.NeutralMass);

            foreach (var c in c_fragments_shift)
            {
                var newMass = c + mass;
                newFragments.Add(newMass);
            }

            var z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition > modPoses.Last() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(z_fragments);

            var z_fragments_shift = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition <= modPoses.First() - 1).Select(p => p.NeutralMass);

            foreach (var z in z_fragments_shift)
            {
                var newMass = z + mass;
                newFragments.Add(newMass);
            }

            var b_fragments = products.Where(p => p.ProductType == ProductType.b && p.AminoAcidPosition < modPoses.First() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(b_fragments);

            var b_fragments_shift = products.Where(p => p.ProductType == ProductType.b && p.AminoAcidPosition >= modPoses.Last() - 1).Select(p => p.NeutralMass);
            double nlocalMass = ((GlycanBox)modBox).NGlycanCount * Glycan.HexNAcMass;
            foreach (var b in b_fragments_shift)
            {
                newFragments.Add(b);
                var newMass = b + nlocalMass;
                newFragments.Add(newMass);
            }

            var y_fragments = products.Where(p => p.ProductType == ProductType.y && p.AminoAcidPosition > modPoses.Last() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(y_fragments);

            var y_fragments_shift = products.Where(p => p.ProductType == ProductType.y && p.AminoAcidPosition <= modPoses.First() - 1).Select(p => p.NeutralMass);

            foreach (var y in y_fragments_shift)
            {
                newFragments.Add(y);
                var newMass = y + nlocalMass;
                newFragments.Add(newMass);
            }

            return newFragments;
        }

        #endregion

        #region N-Glyco related functions

        public static Dictionary<int, double> ScanGetTrimannosylCore(List<MatchedFragmentIon> matchedFragmentIons, Glycan glycan)
        {
            Dictionary<int, double> cores = new Dictionary<int, double>();

            foreach (var fragment in matchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M))
            {
                if (Glycan.TrimannosylCores.ContainsKey((int)fragment.NeutralTheoreticalProduct.NeutralLoss))
                {
                    var pair = Glycan.TrimannosylCores.Where(p => p.Key == (int)fragment.NeutralTheoreticalProduct.NeutralLoss).FirstOrDefault();
                    if (!cores.ContainsKey(pair.Key))
                    {
                        cores.Add(pair.Key, pair.Value);
                    }
                }

            }
            return cores;
        }

        public static bool ScanTrimannosylCoreFilter(List<MatchedFragmentIon> matchedFragmentIons, Glycan glycan)
        {
            Dictionary<int, double> cores = ScanGetTrimannosylCore(matchedFragmentIons, glycan);
            if (cores.Count > 2)
            {
                return true;
            }
            else if (cores.Keys.Contains(83) && cores.Keys.Contains(203))
            {
                return true;
            }
            return false;
        }     

        #endregion

        #region O-Glyco related functions, plan to be deprecated.

        //TO THINK: filter reasonable fragments here. The final solution is to change mzLib.Proteomics.PeptideWithSetModifications.Fragment
        public static List<Product> OGlyGetTheoreticalFragments(DissociationType dissociationType, PeptideWithSetModifications peptide, PeptideWithSetModifications modPeptide)
        {
            List<Product> theoreticalProducts = new List<Product>();        
            HashSet<double> masses = new HashSet<double>();

            List<Product> products = new List<Product>();
            if (dissociationType == DissociationType.HCD || dissociationType == DissociationType.CID)
            {
                List<Product> diag = new List<Product>();
                modPeptide.Fragment(dissociationType, FragmentationTerminus.Both, diag);
                peptide.Fragment(dissociationType, FragmentationTerminus.Both, products);
                products = products.Concat(diag.Where(p => p.ProductType != ProductType.b && p.ProductType != ProductType.y)).ToList();
            }
            else if(dissociationType == DissociationType.ETD)
            {
                modPeptide.Fragment(dissociationType, FragmentationTerminus.Both, products);
            }
            else if(dissociationType == DissociationType.EThcD)
            {
                List<Product> diag = new List<Product>();
                modPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, diag);
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
                products = products.Concat(diag.Where(p => p.ProductType != ProductType.b && p.ProductType != ProductType.y)).ToList();

                List<Product> etdProducts = new List<Product>();
                modPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, etdProducts);
                products = products.Concat(etdProducts.Where(p => p.ProductType != ProductType.y)).ToList();
            }

            foreach (var fragment in products)
            {
                if (!masses.Contains(fragment.NeutralMass))
                {
                    masses.Add(fragment.NeutralMass);
                    theoreticalProducts.Add(fragment);
                }           
            }

            return theoreticalProducts;
        }

        //Find FragmentMass for the fragments that contain localization Information.
        public static List<double> GetLocalFragment(List<Product> products, int modInd, int childBoxInd, LocalizationGraph localizationGraph)
        {
            List<double> newFragments = new List<double>();
            var local_c_fragments = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] - 1 && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1] - 1).ToList();

            foreach (var c in local_c_fragments)
            {
                var newMass = c.NeutralMass + localizationGraph.ChildModBoxes[childBoxInd].Mass;
                newFragments.Add(newMass);
            }

            var local_z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition >= localizationGraph.ModPos[modInd] && p.AminoAcidPosition < localizationGraph.ModPos[modInd + 1]).ToList();

            foreach (var z in local_z_fragments)
            {
                var newMass = z.NeutralMass + (localizationGraph.ModBox.Mass - localizationGraph.ChildModBoxes[childBoxInd].Mass);
                newFragments.Add(newMass);
            }

            return newFragments;
        }

        //Find FragmentMass for the fragments that doesn't contain localization Information. For example, "A|TAABBS|B", c1 and c7, z1 and z7, z8 ion don't contain localization information.
        public static List<double> GetUnlocalFragment(List<Product> products, int[] modPoses, ModBox OGlycanBox)
        {
            var mass = OGlycanBox.Mass;

            List<double> newFragments = new List<double>();
            var c_fragments = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition < modPoses.First() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(c_fragments);

            var c_fragments_shift = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition >= modPoses.Last() - 1).Select(p => p.NeutralMass);

            foreach (var c in c_fragments_shift)
            {
                var newMass = c + mass;
                newFragments.Add(newMass);
            }

            var z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition > modPoses.Last() - 1).Select(p => p.NeutralMass);
            newFragments.AddRange(z_fragments);

            var z_fragments_shift = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition <= modPoses.First() - 1).Select(p => p.NeutralMass);

            foreach (var z in z_fragments_shift)
            {
                var newMass = z + mass;
                newFragments.Add(newMass);
            }

            return newFragments;
        }

        #endregion

        #region Functions are not used now, could be useful in the future. 

        //NGlycopeptide usually contain Y ions with different charge states, especially in sceHCD data. 
        //The purpose of this function is to try match all Y ion with different charges. The usage of this function requires further investigation. 
        //Not sure about OGlycopeptide. 
        public static List<MatchedFragmentIon> GlyMatchOriginFragmentIons(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();

            // if the spectrum has no peaks
            if (scan.ExperimentalFragments != null && !scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            // search for ions in the spectrum

            for (int id = 0; id < theoreticalProducts.Count; id++)
            {
                var product = theoreticalProducts[id];
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                if (product.ProductType == ProductType.M)
                {
                    for (int i = 1; i <= scan.PrecursorCharge; i++)
                    {

                        var closestExperimentalMz = scan.GetClosestExperimentalFragmentMz(product.NeutralMass.ToMz(i), out double? intensity);

                        if (closestExperimentalMz.HasValue && commonParameters.ProductMassTolerance.Within(closestExperimentalMz.Value, product.NeutralMass.ToMz(i)))
                        {
                            matchedFragmentIons.Add(new MatchedFragmentIon(ref product, closestExperimentalMz.Value, intensity.Value, i));
                        }
                    }
                }

                else
                {
                    // get the closest peak in the spectrum to the theoretical peak
                    var closestExperimentalMass = scan.GetClosestExperimentalIsotopicEnvelope(product.NeutralMass);

                    // is the mass error acceptable?
                    if (closestExperimentalMass != null && commonParameters.ProductMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, product.NeutralMass) && closestExperimentalMass.Charge <= scan.PrecursorCharge)
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(ref product, closestExperimentalMass.MonoisotopicMass.ToMz(closestExperimentalMass.Charge),
                            closestExperimentalMass.Peaks.First().intensity, closestExperimentalMass.Charge));
                    }
                }
            }

            return matchedFragmentIons;
        }


        //The function here is to calculate permutation localization which could be used to compare with Graph-Localization.
        public static List<int[]> GetPermutations(List<int> allModPos, int[] glycanBoxId)
        {
            var length = glycanBoxId.Length;
            var indexes = Enumerable.Range(0, length).ToArray();
            int[] orderGlycan = new int[length];

            List<int[]> permutateModPositions = new List<int[]>();

            var combinations = Glycan.GetKCombs(allModPos, length);

            foreach (var com in combinations)
            {
                var permutation = Glycan.GetPermutations(com, length);

                HashSet<string> keys = new HashSet<string>();

                foreach (var per in permutation)
                {
                    Array.Sort(indexes);

                    var orderedPer = per.ToArray();
                    Array.Sort(orderedPer, indexes);

                    for (int i = 0; i < length; i++)
                    {
                        orderGlycan[i] = glycanBoxId[indexes[i]];
                    }
                    var key = string.Join(",", orderGlycan.Select(p => p.ToString()));
                    if (!keys.Contains(key))
                    {
                        keys.Add(key);
                        permutateModPositions.Add(per.ToArray());
                    }
                }
            }

            return permutateModPositions;
        }

        #endregion

    }
}
