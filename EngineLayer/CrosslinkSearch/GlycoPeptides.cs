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

        private static Dictionary<int, double> ScanGetTrimannosylCore(List<MatchedFragmentIon> matchedFragmentIons, Glycan glycan)
        {
            Dictionary<int, double> cores = new Dictionary<int, double>();
            MassDiffAcceptor massDiffAcceptor = new SinglePpmAroundZeroSearchMode(5);
            foreach (var fragment in matchedFragmentIons.Where(p=>p.NeutralTheoreticalProduct.ProductType == ProductType.M))
            {
                if (Glycan.TrimannosylCores.ContainsKey((int)(glycan.Mass - fragment.NeutralTheoreticalProduct.NeutralLoss)))
                {
                    var pair = Glycan.TrimannosylCores.Where(p=>p.Key == (int)(glycan.Mass - fragment.NeutralTheoreticalProduct.NeutralLoss)).FirstOrDefault();
                    cores.Add(pair.Key, pair.Value);
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

        public static IEnumerable<Tuple<int, List<Product>>> NGlyGetTheoreticalFragments(DissociationType dissociationType, 
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

                List<Product> theoreticalProducts = testPeptide.Fragment(dissociationType, FragmentationTerminus.Both).ToList();

                yield return new Tuple<int, List<Product>>(position, theoreticalProducts);
            }
        }

        public static Modification GlycanToModification(Glycan glycan)
        {          
            Dictionary<DissociationType, List<double>> neutralLosses = new Dictionary<DissociationType, List<double>>();
            List<double> lossMasses = glycan.Ions.Where(p => p.IonMass <= 1000).Select(p => glycan.Mass - p.IonMass).OrderBy(p => p).ToList();
            lossMasses.Add(glycan.Mass);
            neutralLosses.Add(DissociationType.HCD, lossMasses);
            neutralLosses.Add(DissociationType.CID, lossMasses);
            //TO DO: add diagnosticIons as a property for Glycan
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
