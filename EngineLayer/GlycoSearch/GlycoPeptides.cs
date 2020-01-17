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
    public class GlycoPeptides
    {
        public static double[] ScanOxoniumIonFilter(Ms2ScanWithSpecificMass theScan, MassDiffAcceptor massDiffAcceptor, DissociationType dissociationType)
        {
            double[] oxoniumIonsintensities = new double[Glycan.AllOxoniumIons.Length];

            if (dissociationType != DissociationType.HCD && dissociationType != DissociationType.CID && dissociationType != DissociationType.EThcD)
            {
                return oxoniumIonsintensities;
            }

            for (int i = 0; i < Glycan.AllOxoniumIons.Length; i++)
            {
                double? matchedIntensity;
                var matchedMass = theScan.GetClosestExperimentalFragmentMz((double)Glycan.AllOxoniumIons[i] / 1E5, out matchedIntensity);
                if (matchedMass.HasValue && massDiffAcceptor.Accepts(matchedMass.Value, (double)Glycan.AllOxoniumIons[i] / 1E5) >= 0)
                {
                    oxoniumIonsintensities[i] = matchedIntensity.Value;
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
            double possiblePeptideMass = precursorMass - (double)glycan.Mass/1E5;
            List<Product> YIons = new List<Product>();
            YIons.Add(new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, precursorMass, 0, 0), (double)glycan.Mass/1E5)); //Y0 ion. Glycan totally loss.
            foreach (var ion in glycan.Ions)
            {
                Product product = new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, precursorMass, 0, 0), (double)ion.LossIonMass/1E5);
                YIons.Add(product);
            }
            return YIons;
        }

        public static List<Product> GetGlycanYIons(PeptideWithSetModifications peptide, Glycan glycan)
        {
            double possiblePeptideMass = peptide.MonoisotopicMass;
            List<Product> YIons = new List<Product>();
            YIons.Add(new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, possiblePeptideMass + (double)glycan.Mass/1E5, 0, 0), (double)glycan.Mass/1E5));
            foreach (var ion in glycan.Ions)
            {
                Product product = new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, possiblePeptideMass + (double)glycan.Mass/1E5, 0, 0), (double)ion.LossIonMass/1E5);
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
                while (iD - 1 >= 0 && massArray[iD-1] == targetMass)
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

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);
            
            return testPeptide;

        }

        #endregion

        #region O-Glyco related functions

        public static bool DissociationTypeContainETD(DissociationType dissociationType)
        {
            if (dissociationType == DissociationType.ETD || dissociationType == DissociationType.EThcD)
            {
                return true;
            }

            return false;
        }

        //TO THINK: filter reasonable fragments here. The final solution is to change mzLib.Proteomics.PeptideWithSetModifications.Fragment
        public static List<Product> OGlyGetTheoreticalFragments(DissociationType dissociationType, PeptideWithSetModifications peptide, PeptideWithSetModifications modPeptide)
        {
            List<Product> theoreticalProducts = new List<Product>();

            if (dissociationType == DissociationType.HCD || dissociationType == DissociationType.CID)
            {   
                var diag = modPeptide.Fragment(dissociationType, FragmentationTerminus.Both).Where(p=>p.ProductType != ProductType.b && p.ProductType != ProductType.y);
                theoreticalProducts = peptide.Fragment(dissociationType, FragmentationTerminus.Both).Concat(diag).ToList();
            }
            else if(dissociationType == DissociationType.ETD)
            {
                theoreticalProducts = modPeptide.Fragment(dissociationType, FragmentationTerminus.Both).ToList();
            }
            else if(dissociationType == DissociationType.EThcD)
            {
                var diag = modPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both).Where(p => p.ProductType != ProductType.b && p.ProductType != ProductType.y);
                theoreticalProducts = peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both).Concat(diag).ToList();
                HashSet<double> masses = new HashSet<double>();

                foreach (var fragment in theoreticalProducts)
                {
                    masses.Add(fragment.NeutralMass);
                }

                foreach (var fragment in modPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both).Where(p=>p.ProductType!=ProductType.y))
                {
                    if (!masses.Contains(fragment.NeutralMass))
                    {
                        theoreticalProducts.Add(fragment);
                        masses.Add(fragment.NeutralMass);
                    }
                }

            }

            return theoreticalProducts;
        }

        public static PeptideWithSetModifications OGlyGetTheoreticalPeptide(int[] theModPositions, PeptideWithSetModifications peptide, GlycanBox glycanBox)
        {
            Modification[] modifications = new Modification[glycanBox.NumberOfMods];
            for (int i = 0; i < glycanBox.NumberOfMods; i++)
            {
                modifications[i] = GlycanBox.GlobalOGlycanModifications[glycanBox.ModIds.ElementAt(i)];
            }

            Dictionary<int, Modification> testMods = new Dictionary<int, Modification>();
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                testMods.Add(mod.Key, mod.Value);
            }

            for (int i = 0; i < theModPositions.Count(); i++)
            {
                testMods.Add(theModPositions.ElementAt(i), modifications[i]);
            }

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

            return testPeptide;
        }

        public static PeptideWithSetModifications OGlyGetTheoreticalPeptide(Tuple<int, int>[] theModPositions, PeptideWithSetModifications peptide)
        {
            Modification[] modifications = new Modification[theModPositions.Length];
            for (int i = 0; i < theModPositions.Length; i++)
            {
                modifications[i] = GlycanBox.GlobalOGlycanModifications[theModPositions[i].Item2];
            }

            Dictionary<int, Modification> testMods = new Dictionary<int, Modification>();
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                testMods.Add(mod.Key, mod.Value);
            }

            for (int i = 0; i < theModPositions.Length; i++)
            {
                testMods.Add(theModPositions[i].Item1, modifications[i]);
            }

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

            return testPeptide;
        }

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

        //The purpose of the funtion is to generate hash fragment ions without generate the PeptideWithMod. keyValuePair key:GlycanBoxId, Value:mod sites
        public static int[] GetFragmentHash(List<Product> products, Tuple<int, int[]> keyValuePair, GlycanBox[] OGlycanBoxes, int FragmentBinsPerDalton)
        {
            double[] newFragments = products.Select(p => p.NeutralMass).ToArray();
            var len = products.Count / 3;
            if (keyValuePair.Item2!=null)
            {
                for (int i = 0; i < keyValuePair.Item2.Length; i++)
                {
                    var j = keyValuePair.Item2[i];
                    while (j <= len + 1)
                    {
                        newFragments[j - 2] += GlycanBox.GlobalOGlycans[OGlycanBoxes[keyValuePair.Item1].ModIds[i]].Mass;
                        j++;
                    }
                    j = keyValuePair.Item2[i];
                    while (j >= 3)
                    {
                        newFragments[len * 2 - j + 2] += GlycanBox.GlobalOGlycans[OGlycanBoxes[keyValuePair.Item1].ModIds[i]].Mass;
                        newFragments[len * 3 - j + 2] += GlycanBox.GlobalOGlycans[OGlycanBoxes[keyValuePair.Item1].ModIds[i]].Mass;
                        j--;
                    }
                }
            }


            int[] fragmentHash = new int[products.Count];
            for (int i = 0; i < products.Count; i++)
            {
                fragmentHash[i] = (int)Math.Round(newFragments[i] * FragmentBinsPerDalton);
            }
            return fragmentHash;
        }

        //Find FragmentHsh for current box at modInd. 
        //TO DO: How about y-ions from ETD?
        public static int[] GetLocalFragmentHash(List<Product> products, int peptideLength, int[] modPoses, int modInd, GlycanBox OGlycanBox, GlycanBox localOGlycanBox, int FragmentBinsPerDalton)
        {
            List<double> newFragments = new List<double>();
            var local_c_fragments = products.Where(p => p.ProductType == ProductType.c && p.TerminusFragment.AminoAcidPosition >= modPoses[modInd]-1 && p.TerminusFragment.AminoAcidPosition < modPoses[modInd+1]-1).ToList();

            foreach (var c in local_c_fragments)
            {
                var newMass = c.NeutralMass + localOGlycanBox.Mass;
                newFragments.Add(newMass);
            }

            var local_z_fragments = products.Where(p => p.ProductType == ProductType.zDot && p.TerminusFragment.AminoAcidPosition >= modPoses[modInd]  && p.TerminusFragment.AminoAcidPosition < modPoses[modInd + 1] ).ToList();

            foreach (var z in local_z_fragments)
            {
                var newMass = z.NeutralMass + (OGlycanBox.Mass - localOGlycanBox.Mass);
                newFragments.Add(newMass);
            }


            int[] fragmentHash = new int[newFragments.Count];
            for (int i = 0; i < newFragments.Count; i++)
            {
                fragmentHash[i] = (int)Math.Round(newFragments[i] * FragmentBinsPerDalton);
            }
            return fragmentHash;
        }

        //The oxoniumIonIntensities is related with Glycan.AllOxoniumIons. 
        //Rules are coded in the function.    
        public static bool OxoniumIonsAnalysis(double[] oxoniumIonsintensities, GlycanBox glycanBox)
        {
            //If a glycopeptide spectrum does not have 292.1027 or 274.0921, then remove all glycans that have sialic acids from the search.
            if (oxoniumIonsintensities[10] == 0 && oxoniumIonsintensities[12] == 0)
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

            return true;
        }

        #endregion


    }
}
