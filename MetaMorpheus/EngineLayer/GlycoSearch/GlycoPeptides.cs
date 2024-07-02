using MassSpectrometry;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using Omics.Modifications;

namespace EngineLayer.GlycoSearch
{
    public static class GlycoPeptides 
    {
        /// <summary>
        /// Generate a list of isotopic intesitry of the oxonium ions
        /// </summary>
        /// <param name="theScan"> The MS2 Scan</param>
        /// <param name="massDiffAcceptor"></param>
        /// <returns> int[], The intensity list </returns>
        public static double[] ScanOxoniumIonFilter(Ms2ScanWithSpecificMass theScan, MassDiffAcceptor massDiffAcceptor)
        {
            double[] oxoniumIonsintensities = new double[Glycan.AllOxoniumIons.Length];

            for (int i = 0; i < Glycan.AllOxoniumIons.Length; i++)
            {
                var oxoMass = ((double)Glycan.AllOxoniumIons[i] / 1E5).ToMass(1);
                var envelope = theScan.GetClosestExperimentalIsotopicEnvelope(oxoMass);
                if (massDiffAcceptor.Accepts(envelope.MonoisotopicMass, oxoMass) >= 0)
                {
                    oxoniumIonsintensities[i] = envelope.TotalIntensity;
                }

            }

            return oxoniumIonsintensities;
        }

        #region N-Glyco related functions

        public static Dictionary<int, double> ScanGetTrimannosylCore(List<MatchedFragmentIon> matchedFragmentIons, Glycan glycan)
        {
            Dictionary<int, double> cores = new Dictionary<int, double>();

            foreach (var fragment in matchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M))
            {
                if (Glycan.TrimannosylCores.ContainsKey((int)((double)glycan.Mass / 1E5 - fragment.NeutralTheoreticalProduct.NeutralLoss)))
                {
                    var pair = Glycan.TrimannosylCores.Where(p => p.Key == (int)((double)glycan.Mass / 1E5 - fragment.NeutralTheoreticalProduct.NeutralLoss)).FirstOrDefault();
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

        public static List<Product> GetGlycanYIons(double precursorMass, Glycan glycan)
        {
            List<Product> YIons = new List<Product>();
            YIons.Add(new Product(ProductType.M, FragmentationTerminus.Both, precursorMass - (double)glycan.Mass / 1E5, 0, 0, (double)glycan.Mass/1E5)); //Y0 ion. Glycan totally loss.
            foreach (var ion in glycan.Ions)
            {
                Product product = new Product(ProductType.M, FragmentationTerminus.Both, precursorMass - (double)ion.LossIonMass / 1E5, 0, 0, (double)ion.LossIonMass/1E5);
                YIons.Add(product);
            }
            return YIons;
        }

        public static List<Product> GetGlycanYIons(PeptideWithSetModifications peptide, Glycan glycan)
        {
            double possiblePeptideMass = peptide.MonoisotopicMass;
            List<Product> YIons = new List<Product>();
            YIons.Add(new Product(ProductType.M, FragmentationTerminus.Both, possiblePeptideMass + (double)glycan.Mass/1E5, 0, 0, (double)glycan.Mass/1E5));
            foreach (var ion in glycan.Ions)
            {
                Product product = new Product(ProductType.M, FragmentationTerminus.Both, possiblePeptideMass + (double)glycan.Mass/1E5, 0, 0, (double)ion.LossIonMass/1E5);
                YIons.Add(product);
            }
            return YIons;
        }

        public static Tuple<int, double, double>[] MatchBestGlycan(Ms2ScanWithSpecificMass theScan, Glycan[] glycans, CommonParameters commonParameters)
        {
            Tuple<int, double, double>[] tuples = new Tuple<int, double, double>[glycans.Length]; //Tuple<id, Yion matched score, glycan mass> 
            //TO DO: Parallel this?
            for (int i = 0; i < glycans.Length; i++)
            {
                if (theScan.PrecursorMass - (double)glycans[i].Mass/1E5 < 350) //Filter large glycans
                {
                    continue;
                }
                List<Product> YIons = GetGlycanYIons(theScan.PrecursorMass, glycans[i]);
                List<MatchedFragmentIon> matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(theScan, YIons, commonParameters);
                if (ScanTrimannosylCoreFilter(matchedFragmentIons, glycans[i]))
                {
                    var score = MetaMorpheusEngine.CalculatePeptideScore(theScan.TheScan, matchedFragmentIons);
                    tuples[i] = new Tuple<int, double, double>(i, score, (double)glycans[i].Mass/1E5);
                }
            }

            return tuples;
        }

        public static int BinarySearchGetIndex(double[] massArray, double targetMass)
        {
            var iD = Array.BinarySearch(massArray, targetMass);
            if (iD < 0) { iD = ~iD; }
            else
            {
                while (iD - 1 >= 0 && massArray[iD-1] >= targetMass - 0.00000001)
                {
                    iD--;
                }
            }
            return iD;
        }

        public static double CalculateGlycoPeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons, double maximumMassThatFragmentIonScoreIsDoubled)
        {
            double score = 0;

            foreach (var fragment in matchedFragmentIons)
            {
                if (fragment.NeutralTheoreticalProduct.ProductType != ProductType.M && fragment.NeutralTheoreticalProduct.ProductType != ProductType.D)
                {
                    double fragmentScore = 1 + (fragment.Intensity / thisScan.TotalIonCurrent);
                    score += fragmentScore;

                    if (fragment.NeutralTheoreticalProduct.NeutralMass <= maximumMassThatFragmentIonScoreIsDoubled)
                    {
                        score += fragmentScore;
                    }
                }
            }

            return score;
        }

        public static PeptideWithSetModifications GenerateGlycopeptide(int position, PeptideWithSetModifications peptide, Glycan glycan)
        {
            Modification modification = Glycan.NGlycanToModification(glycan);


            Dictionary<int, Modification> testMods = new Dictionary<int, Modification> { { position, modification } };

            if (!peptide.AllModsOneIsNterminus.Keys.Contains(position))
            {
                foreach (var mod in peptide.AllModsOneIsNterminus)
                {
                    testMods.Add(mod.Key, mod.Value);
                }
            }

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidue,
                peptide.OneBasedEndResidue, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);
            
            return testPeptide;

        }

        #endregion

        #region O-Glyco related functions

        public static bool DissociationTypeContainETD(DissociationType dissociationType, List<ProductType> customIons)
        {
            if (dissociationType == DissociationType.ETD || dissociationType == DissociationType.EThcD)
            {
                return true;
            }

            if (dissociationType == DissociationType.Custom ) //Use the fragment type to determine the dissociation type.
            {
                if (customIons.Contains(ProductType.zDot) || customIons.Contains(ProductType.c))
                {
                    return true;
                } 
            }

            return false;
        }

        //TO THINK: filter reasonable fragments here. The final solution is to change mzLib.Proteomics.PeptideWithSetModifications.Fragment

        /// <summary>
        /// Get the theoretical fragments of the peptide with the glycan modification. With different dissociation type, the fragment ions are different.
        /// </summary>
        /// <param name="dissociationType"></param>
        /// <param name="customIons"></param>
        /// <param name="peptide"></param>
        /// <param name="modPeptide"></param>
        /// <returns> product[], Fragments list</returns>
        public static List<Product> OGlyGetTheoreticalFragments(DissociationType dissociationType, List<ProductType> customIons, PeptideWithSetModifications peptide, PeptideWithSetModifications modPeptide)
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
            else if (dissociationType == DissociationType.Custom)
            {
                List<Product> diag = new List<Product>();
                modPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, diag);

                List<Product> hcdProducts = new List<Product>();
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, hcdProducts);

                List<Product> etdProducts = new List<Product>();
                modPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, etdProducts);

                products = products.Concat(diag.Where(p => p.ProductType != ProductType.b && p.ProductType != ProductType.y)).ToList();
                foreach (var pt in customIons)
                {
                    products = products.Concat(hcdProducts.Where(p => p.ProductType == pt)).ToList();
                    products = products.Concat(etdProducts.Where(p => p.ProductType != ProductType.y && p.ProductType == pt)).ToList();
                }

            }

            foreach (var fragment in products) //this part just for the unique fragment ions. (filter the fragment with the same neturalMass)
            {
                if (!masses.Contains(fragment.NeutralMass))
                {
                    masses.Add(fragment.NeutralMass);
                    theoreticalProducts.Add(fragment);
                }           
            }

            return theoreticalProducts;
        }

        
        /// <summary>
        /// Generate the theroertical glycan modified peptide. With the glycanBox, modPos, and the peptide.
        /// </summary>
        /// <param name="theModPositions"></param>
        /// <param name="peptide"></param>
        /// <param name="glycanBox"></param>
        /// <returns> A modfiied peptide.</returns>
        public static PeptideWithSetModifications OGlyGetTheoreticalPeptide(int[] theModPositions, PeptideWithSetModifications peptide, GlycanBox glycanBox)
        {
            Modification[] modifications = new Modification[glycanBox.NumberOfMods];
            for (int i = 0; i < glycanBox.NumberOfMods; i++)
            {
                modifications[i] = GlycanBox.GlobalOGlycanModifications[glycanBox.ModIds.ElementAt(i)]; // transfer the glycanBox information to a new list.
            }

            Dictionary<int, Modification> testMods = new Dictionary<int, Modification>();
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                testMods.Add(mod.Key, mod.Value);   // transfer the AllMod information to a new list.
            }

            for (int i = 0; i < theModPositions.Count(); i++)
            {
                testMods.Add(theModPositions.ElementAt(i), modifications[i]);  //combine the glycanBox information to the AllMod list
            }

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidue,
                peptide.OneBasedEndResidue, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

            return testPeptide;
        }

        /// <summary>
        /// Generate the theroertical glycan modified peptide. With the route the peptide. Because the route contains the glycanBox and modPos information.
        /// </summary>
        /// <param name="theModPositions"></param>
        /// <param name="peptide"></param>
        /// <returns> A modfiied peptide </returns>
        public static PeptideWithSetModifications OGlyGetTheoreticalPeptide(Route theModPositions, PeptideWithSetModifications peptide)
        {
            Modification[] modifications = new Modification[theModPositions.Mods.Count];
            for (int i = 0; i < theModPositions.Mods.Count; i++)
            {
                modifications[i] = GlycanBox.GlobalOGlycanModifications[theModPositions.Mods[i].Item2];
            }

            Dictionary<int, Modification> testMods = new Dictionary<int, Modification>();
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                testMods.Add(mod.Key, mod.Value);
            }

            for (int i = 0; i < theModPositions.Mods.Count; i++)
            {
                testMods.Add(theModPositions.Mods[i].Item1, modifications[i]);
            }

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidue,
                peptide.OneBasedEndResidue, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

            return testPeptide;
        }

        //Should be revised for easier understanding.
        /// <summary>
        /// Generate all possible glycosite for the glycan set. Supposed we will put the glycan on the glycosite in sequence.
        /// </summary>
        /// <param name="allModPos"> Ex. [3,5,2,7]</param>
        /// <param name="glycanBoxId"> Ex. [2,2,3] means id2 + id2 + id3 </param>
        /// <returns> A glycosite set collection. Ex. ([2,5,7],[3,5,7]...), each one list means the glcosites for glycanBox. 
        /// [2,5,7] means we will put the glycan on position 2, 5, 7. </returns>
        /// </returns>
        public static List<int[]> GetPermutations(List<int> allModPos, int[] glycanBoxId)
        {
            var length = glycanBoxId.Length;
            var indexes = Enumerable.Range(0, length).ToArray();  // just the index for the glycanBoxId to keep the order.
            int[] orderGlycan = new int[length];

            List<int[]> permutateModPositions = new List<int[]>(); //The list to store all possible permutation localization.

            var combinations = Glycan.GetKCombs(allModPos, length); //Get all possible combinations of the mod sites. ex. four site[1,2,3,4], length:3 -> combination [1,2,3], [1,2,4], [1,3,4], [2,3,4]
        
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
                    if (!keys.Contains(key)) //Remove the duplicate permutation localization.
                    {
                        keys.Add(key);
                        permutateModPositions.Add(per.ToArray());
                    }
                }
            }

            return permutateModPositions;
        }


        /// <summary>
        /// Generate the new fragment list, we add the glycan mass to the c ions and z ions from the peptide fragment list
        /// </summary>
        /// <param name="products"></param>
        /// <param name="keyValuePair"></param>
        /// <param name="OGlycanBoxes"></param>
        /// <param name="FragmentBinsPerDalton"></param>
        /// <returns></returns>
        public static int[] GetFragmentHash(List<Product> products, Tuple<int, int[]> keyValuePair, GlycanBox[] OGlycanBoxes, int FragmentBinsPerDalton)
        {
            double[] newFragments = products.OrderBy(p=>p.ProductType).ThenBy(p=>p.FragmentNumber).Select(p => p.NeutralMass).ToArray(); // store the fragment mass in the order of c1, c2, c3, y1, y2, y3, z1, z2, z3
            var len = products.Count / 3;
            if (keyValuePair.Item2!=null)
            {
                for (int i = 0; i < keyValuePair.Item2.Length; i++) // we want to add the glycan mass to the c ions and z ions that contain the glycan.
                {                                                   // y ions didn't change in EThcD for O-glyco, so we just need to deal with c ions and z ions.
                    var j = keyValuePair.Item2[i];
                    while (j <= len + 1) // for c ions
                    {
                        newFragments[j - 2] += (double)GlycanBox.GlobalOGlycans[OGlycanBoxes[keyValuePair.Item1].ModIds[i]].Mass/1E5;
                        j++;
                    }
                    j = keyValuePair.Item2[i]; // reset the j to the position of the glycan
                    while (j >= 3)             // for z ions
                    {
                        newFragments[len * 3 - j + 2] += (double)GlycanBox.GlobalOGlycans[OGlycanBoxes[keyValuePair.Item1].ModIds[i]].Mass/1E5;
                        j--;
                    }
                }
            }


            int[] fragmentHash = new int[products.Count]; // store the fragment mass in the order of c1, c2, c3, y1, y2, y3, z1, z2, z3 and with the umit of FragmentBinsPerDalton
            for (int i = 0; i < products.Count; i++)
            {
                fragmentHash[i] = (int)Math.Round(newFragments[i] * FragmentBinsPerDalton);
            }
            return fragmentHash;
        }

 
        /// <summary>
        /// Generate the fragment list with the specific childBox located on specific modPos. At here, the ModInd is the index for modPos. Not used in the current version.
        /// </summary>
        /// <param name="products"></param>
        /// <param name="modPoses"> ModPos list </param>
        /// <param name="modInd"> Specific ModPos, index in ModPos</param>
        /// <param name="OGlycanBox"> Whole glycanBox</param>
        /// <param name="localOGlycanBox">Partial glycanBox, at here is the childBox</param>
        /// <returns></returns>
        public static List<double> GetLocalFragment(List<Product> products, int[] modPoses, int modInd, ModBox OGlycanBox, ModBox localOGlycanBox)
        {
            List<double> newFragments = new List<double>();
            var local_c_fragments = products.Where(p => p.ProductType == ProductType.c && p.AminoAcidPosition >= modPoses[modInd] - 1 && p.AminoAcidPosition < modPoses[modInd + 1] - 1).ToList();

            foreach (var c in local_c_fragments)
            {
                var newMass = c.NeutralMass + localOGlycanBox.Mass;
                newFragments.Add(newMass);
            }

            var local_z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition >= modPoses[modInd] && p.AminoAcidPosition < modPoses[modInd + 1]).ToList();

            foreach (var z in local_z_fragments)
            {
                var newMass = z.NeutralMass + (OGlycanBox.Mass - localOGlycanBox.Mass);
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

            var z_fragments_shift = products.Where(p => p.ProductType == ProductType.zDot && p.AminoAcidPosition < modPoses.First() - 1).Select(p => p.NeutralMass);

            foreach (var z in z_fragments_shift)
            {
                var newMass = z + mass;
                newFragments.Add(newMass);
            }

            return newFragments;
        }


        /// <summary>
        /// The oxoniumIonIntensities is related with Glycan.AllOxoniumIons. Filter the invalid data.
        /// </summary>
        /// <param name="oxoniumIonsintensities"></param>
        /// <param name="glycanBox"></param>
        /// <returns >True : is Oglycan and pass the filter, False: isn't Oglycan and not pass the filter</returns>
        public static bool OxoniumIonsAnalysis(double[] oxoniumIonsintensities, GlycanBox glycanBox)
        {
            //If a glycopeptide spectrum does not have 292.1027 or 274.0921, then remove all glycans that have sialic acids from the search.
            if (oxoniumIonsintensities[10] <= 0 && oxoniumIonsintensities[12] <= 0)
            {
                if (glycanBox.Kind[2] != 0 || glycanBox.Kind[3] != 0)
                {
                    return false;
                }
            }

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
            // use some other oxonium ions to determine the glycan type.

            return true;
        }

        #endregion
    }
}
