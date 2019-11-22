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
        public static double[] ScanGetOxoniumIons(Ms2ScanWithSpecificMass theScan, MassDiffAcceptor massDiffAcceptor)
        {
            double[] oxoniumIonsIntensities = new double[Glycan.allOxoniumIons.Count()];

            for (int i = 0; i < Glycan.allOxoniumIons.Count(); i++)
            {
                var ioxo = (double)Glycan.allOxoniumIons[i]/1E5;
                //Match oxoniumIons from original scan
                int matchedPeakIndex = theScan.TheScan.MassSpectrum.GetClosestPeakIndex(ioxo).Value;

                //Another way to match oxoniumIons
                //int diagnosticIonLabel = (int)Math.Round(ioxo.ToMz(1), 0);
                //HashSet<Product> diagnosticIons = new HashSet<Product>();
                //diagnosticIons.Add(new Product(ProductType.D, new NeutralTerminusFragment(FragmentationTerminus.Both, ioxo, diagnosticIonLabel, 0), 0));
                //var matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(theScan, diagnosticIons.ToList(), commonParameters);

                if (massDiffAcceptor.Accepts(theScan.TheScan.MassSpectrum.XArray[matchedPeakIndex], ioxo) >= 0)
                {
                    oxoniumIonsIntensities[i] = theScan.TheScan.MassSpectrum.YArray[matchedPeakIndex];
                }
                else
                {
                    oxoniumIonsIntensities[i] = 0;
                }
            }

         var largest = oxoniumIonsIntensities.Max();
         return oxoniumIonsIntensities.Select(p => p / largest).ToArray();
        }

        public static void OxoniumIonsAnalysis(double[] vs)
        {

        }

        public static int ScanOxoniumIonFilter(Ms2ScanWithSpecificMass theScan, MassDiffAcceptor massDiffAcceptor, DissociationType dissociationType)
        {
            if (dissociationType != DissociationType.HCD && dissociationType != DissociationType.CID && dissociationType != DissociationType.EThcD)
            {
                return 1;
            }

            int totalNum = 0;

            foreach (var ioxo in Glycan.oxoniumIons)
            {
                int matchedPeakIndex = theScan.TheScan.MassSpectrum.GetClosestPeakIndex((double)ioxo/1E5).Value;
                if (massDiffAcceptor.Accepts(theScan.TheScan.MassSpectrum.XArray[matchedPeakIndex], (double)ioxo/1E5) >= 0)
                {
                    totalNum++;
                }
            }

            return totalNum;
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

                foreach (var fragment in modPeptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both))
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
            Modification[] modifications = new Modification[glycanBox.NumberOfGlycans];
            for (int i = 0; i < glycanBox.NumberOfGlycans; i++)
            {
                modifications[i] = GlycanBox.GlobalOGlycanModifications[glycanBox.GlycanIds.ElementAt(i)];
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
