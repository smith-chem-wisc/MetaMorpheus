using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.CrosslinkSearch
{
    public class GlycoPeptides
    {
        public static bool ScanOxoniumIonFilter(Ms2ScanWithSpecificMass theScan, DissociationType dissociationType)
        {
            if (dissociationType == DissociationType.ETD)
            {
                return true;
            }
            var massDiffAcceptor = new SinglePpmAroundZeroSearchMode(10);

            int intensity138 = 101;
            int intensity204 = 101;
            int intensity366 = 101;
            int totalNum = 0;

            double[] experimental_intensities = theScan.TheScan.MassSpectrum.YArray;
            int[] experimental_intensities_rank = CrosslinkSpectralMatch.GenerateIntensityRanks(experimental_intensities);

            foreach (var ioxo in Glycan.oxoniumIons)
            {
                int matchedPeakIndex = theScan.TheScan.MassSpectrum.GetClosestPeakIndex(ioxo.Value).Value;
                if (massDiffAcceptor.Accepts(theScan.TheScan.MassSpectrum.XArray[matchedPeakIndex], ioxo.Value) >= 0)
                {
                    totalNum++;
                    if (totalNum > 1)
                    {
                        return true;
                    }
                    if (ioxo.Key == 138)
                        intensity138 = experimental_intensities_rank[matchedPeakIndex];
                    if (ioxo.Key == 204)
                        intensity204 = experimental_intensities_rank[matchedPeakIndex];
                    if (ioxo.Key == 366)
                        intensity366 = experimental_intensities_rank[matchedPeakIndex];
                }
            }

            if (intensity138 < 101 || intensity204 < 101 || intensity366 < 101 || totalNum > 1)
            {
                return true;
            }

            return false;
        }

        public static Dictionary<int, double> ScanGetTrimannosylCore(List<MatchedFragmentIon> matchedFragmentIons, Glycan glycan)
        {
            Dictionary<int, double> cores = new Dictionary<int, double>();

            foreach (var fragment in matchedFragmentIons.Where(p=>p.NeutralTheoreticalProduct.ProductType == ProductType.M))
            {
                if (Glycan.TrimannosylCores.ContainsKey((int)(glycan.Mass - fragment.NeutralTheoreticalProduct.NeutralLoss)))
                {
                    var pair = Glycan.TrimannosylCores.Where(p=>p.Key == (int)(glycan.Mass - fragment.NeutralTheoreticalProduct.NeutralLoss)).FirstOrDefault();
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

        public static List<Product> GetGlycanYIons(Ms2ScanWithSpecificMass theScan, Glycan glycan)
        {
            double possiblePeptideMass = theScan.PrecursorMass - glycan.Mass;
            List<Product> YIons = new List<Product>();
            YIons.Add(new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, theScan.PrecursorMass, 0, 0), glycan.Mass)); //Y0 ion. Glycan totally loss.
            foreach (var ion in glycan.Ions)
            {
                Product product = new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, theScan.PrecursorMass, 0, 0), ion.LossIonMass);
                YIons.Add(product);
            }
            return YIons;
        }

        public static List<Product> GetGlycanYIons(PeptideWithSetModifications peptide, Glycan glycan)
        {
            double possiblePeptideMass = peptide.MonoisotopicMass;
            List<Product> YIons = new List<Product>();
            YIons.Add(new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, possiblePeptideMass + glycan.Mass, 0, 0), glycan.Mass));
            foreach (var ion in glycan.Ions)
            {
                Product product = new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, possiblePeptideMass + glycan.Mass, 0, 0), ion.LossIonMass);
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
                if (theScan.PrecursorMass - glycans[i].Mass < 350) //Filter large glycans
                {
                    continue;
                }
                List<Product> YIons = GetGlycanYIons(theScan, glycans[i]);
                List<MatchedFragmentIon> matchedFragmentIons = MetaMorpheusEngine.MatchFragmentIons(theScan, YIons, commonParameters);
                if (ScanTrimannosylCoreFilter(matchedFragmentIons, glycans[i]))
                {
                    var score = MetaMorpheusEngine.CalculatePeptideScore(theScan.TheScan, matchedFragmentIons, 0);
                    tuples[i] = new Tuple<int, double, double>(i, score, glycans[i].Mass);
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
                while (massArray[iD-1] == targetMass)
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

        public static IEnumerable<Tuple<int, List<Product>>> NGlyGetTheoreticalFragments(Ms2ScanWithSpecificMass theScan, DissociationType dissociationType, 
            List<int> possibleModPositions, PeptideWithSetModifications peptide, Glycan glycan)
        {
            Modification modification = GlycanToModification(glycan);

            foreach (var position in possibleModPositions)
            {               
                Dictionary<int, Modification> testMods = new Dictionary<int, Modification> { { position + 1, modification } };

                foreach (var mod in peptide.AllModsOneIsNterminus)
                {
                    testMods.Append(mod);
                }

                var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                    peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

                List<Product> theoreticalProducts = testPeptide.Fragment(dissociationType, FragmentationTerminus.Both).Where(p=>p.ProductType!= ProductType.M).ToList();
                theoreticalProducts.AddRange(GetGlycanYIons(theScan, glycan));

                yield return new Tuple<int, List<Product>>(position, theoreticalProducts);
            }
        }

        public static Modification GlycanToModification(Glycan glycan)
        {          
            Dictionary<DissociationType, List<double>> neutralLosses = new Dictionary<DissociationType, List<double>>();
            List<double> lossMasses = glycan.Ions.Where(p=>p.IonMass < 570).Select(p => p.LossIonMass).OrderBy(p => p).ToList(); //570 is a cutoff for glycan ion size 2N1H, which will generate fragment ions. 
            neutralLosses.Add(DissociationType.HCD, lossMasses);
            neutralLosses.Add(DissociationType.CID, lossMasses);
            
            Dictionary<DissociationType, List<double>> diagnosticIons = new Dictionary<DissociationType, List<double>>();
            diagnosticIons.Add(DissociationType.HCD, glycan.GetDiagnosticIons().Values.ToList());
            diagnosticIons.Add(DissociationType.CID, glycan.GetDiagnosticIons().Values.ToList());

            Modification modification = new Modification(_originalId:"Glycan", _monoisotopicMass: glycan.Mass, _neutralLosses: neutralLosses, _diagnosticIons : diagnosticIons);
            return modification;
        }

        public static IEnumerable<Tuple<int[] , Tuple<int[], List<Product>>>> OGlyGetTheoreticalFragments(DissociationType dissociationType, 
            List<int> possibleModPositions, PeptideWithSetModifications peptide, GlycanBox glycanBox)
        {
            Modification[] modifications = new Modification[glycanBox.glycans.Count];

            for (int i = 0; i < glycanBox.glycans.Count; i++)
            {
                modifications[i] = GlycanToModification(glycanBox.glycans[i]);
            }

            foreach (var modcombine in Glycan.GetPermutations(Enumerable.Range(0, glycanBox.glycans.Count), glycanBox.glycans.Count))
            {
                foreach (var combine in Glycan.GetKCombs(possibleModPositions, glycanBox.glycans.Count))
                {
                    Dictionary<int, Modification> testMods = new Dictionary<int, Modification>();

                    for (int i = 0; i < glycanBox.glycans.Count; i++)
                    {
                        testMods.Add(combine.ElementAt(i), modifications[modcombine.ElementAt(i)]);
                    }

                    var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                    peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

                    List<Product> theoreticalProducts = testPeptide.Fragment(dissociationType, FragmentationTerminus.Both).ToList();

                    yield return new Tuple<int[], Tuple<int[], List<Product>>>(combine.ToArray(), new Tuple<int[], List<Product>>(modcombine.ToArray(), theoreticalProducts));
                }
            }


        }

    }
}
