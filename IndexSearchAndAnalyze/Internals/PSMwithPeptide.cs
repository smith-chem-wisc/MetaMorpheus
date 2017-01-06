using Chemistry;
using MassSpectrometry;
using MetaMorpheus;
using Spectra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace IndexSearchAndAnalyze
{
    public class PSMwithPeptide
    {
        public NewPsm newPsm;
        public PeptideWithSetModifications peptideWithSetModifications;

        public bool isDecoy
        {
            get
            {
                return peptideWithSetModifications.protein.isDecoy;
            }
        }

        public double ScoreFromSearch
        {
            get
            {
                return newPsm.ScoreFromSearch;
            }
        }

        public double scanPrecursorMass
        {
            get
            {
                return newPsm.scanPrecursorMass;
            }
        }

        public List<double> LocalizedScores
        {
            get
            {
                return newPsm.LocalizedScores;
            }
        }

        public double PeptideMonoisotopicMass
        {
            get
            {
                return peptideWithSetModifications.MonoisotopicMass;
            }
        }

        public string FullSequence
        {
            get
            {
                return peptideWithSetModifications.Sequence;
            }
        }

        public string BaseSequence
        {
            get
            {
                return peptideWithSetModifications.BaseSequence;
            }
        }

        public PSMwithPeptide(NewPsm newPsm, PeptideWithSetModifications peptideWithSetModifications, Tolerance fragmentTolerance, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile)
        {
            this.newPsm = newPsm;
            this.peptideWithSetModifications = peptideWithSetModifications;

            var allProductTypes = new List<ProductType>() { ProductType.b, ProductType.y };
            IMsDataScan<IMzSpectrum<MzPeak>> theScan;
            if (newPsm.matchedIonsList == null)
            {
                theScan = myMsDataFile.GetOneBasedScan(newPsm.scanNumber);
                double selectedMZ;
                int selectedCharge;
                theScan.TryGetSelectedIonGuessMonoisotopicMZ(out selectedMZ);
                theScan.TryGetSelectedIonGuessChargeStateGuess(out selectedCharge);
                var scanPrecursorMass = selectedMZ.ToMass(selectedCharge);

                Dictionary<ProductType, double[]> MatchedIonDict = new Dictionary<ProductType, double[]>();
                double score = 0;
                foreach (var huh in allProductTypes)
                {
                    var df = peptideWithSetModifications.FastSortedProductMasses(new List<ProductType>() { huh });
                    double[] matchedIonList = new double[df.Length];
                    score += MatchIons(theScan, fragmentTolerance, df, matchedIonList);
                    MatchedIonDict.Add(huh, matchedIonList);
                }

                newPsm.matchedIonsList = MatchedIonDict;
                newPsm.ScoreFromMatch = score;
            }

            if (newPsm.LocalizedScores == null)
            {
                theScan = myMsDataFile.GetOneBasedScan(newPsm.scanNumber);
                List<double> localizedScores = new List<double>();
                for (int indexToLocalize = 0; indexToLocalize < peptideWithSetModifications.Length; indexToLocalize++)
                {
                    PeptideWithSetModifications localizedPeptide = peptideWithSetModifications.Localize(indexToLocalize, scanPrecursorMass - peptideWithSetModifications.MonoisotopicMass);

                    var gg = localizedPeptide.FastSortedProductMasses(allProductTypes);
                    double[] matchedIonList = new double[gg.Length];
                    var score = MatchIons(theScan, fragmentTolerance, gg, matchedIonList);
                    localizedScores.Add(score);
                }
                newPsm.LocalizedScores = localizedScores;
            }
        }

        internal static double MatchIons(IMsDataScan<IMzSpectrum<MzPeak>> thisScan, Tolerance product_mass_tolerance_value, double[] sorted_theoretical_product_masses_for_this_peptide, double[] matchedIonsList)
        {
            var TotalProductsHere = sorted_theoretical_product_masses_for_this_peptide.Length;
            if (TotalProductsHere == 0)
                return 0;
            int MatchingProductsHere = 0;
            double MatchingIntensityHere = 0;

            int theoreticalLeft = TotalProductsHere;
            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.xArray;
            double[] experimental_intensities = thisScan.MassSpectrum.yArray;
            int num_experimental_peaks = experimental_mzs.Length;

            int theoreticalIndex = 0;
            double nextTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[0];
            double nextTheoreticalMZ = nextTheoreticalMass + 1.007276466879;

            double currentExperimentalMZ;
            for (int i = 0; i < num_experimental_peaks; i++)
            {
                currentExperimentalMZ = experimental_mzs[i];
                if (product_mass_tolerance_value.Within(currentExperimentalMZ, nextTheoreticalMZ))
                {
                    MatchingProductsHere++;
                    MatchingIntensityHere += experimental_intensities[i];
                    matchedIonsList[theoreticalIndex] = nextTheoreticalMass;
                }
                else if (currentExperimentalMZ < nextTheoreticalMZ)
                    continue;
                else
                    matchedIonsList[theoreticalIndex] = -nextTheoreticalMass;
                i--;
                // Passed a theoretical! Move counter forward
                theoreticalIndex++;
                if (theoreticalIndex == TotalProductsHere)
                    break;
                nextTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[theoreticalIndex];
                nextTheoreticalMZ = nextTheoreticalMass + 1.007276466879;
            }
            return MatchingProductsHere + MatchingIntensityHere / thisScan.TotalIonCurrent;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(newPsm.ToString() + '\t');

            sb.Append(peptideWithSetModifications.PreviousAminoAcid.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(peptideWithSetModifications.Sequence.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(peptideWithSetModifications.NextAminoAcid.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(peptideWithSetModifications.numVariableMods.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(peptideWithSetModifications.OneBasedStartResidueInProtein.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(peptideWithSetModifications.OneBasedEndResidueInProtein.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(peptideWithSetModifications.PeptideDescription + '\t');
            sb.Append(peptideWithSetModifications.MissedCleavages.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(peptideWithSetModifications.MonoisotopicMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');

            sb.Append(peptideWithSetModifications.protein.FullDescription.ToString(CultureInfo.InvariantCulture) + '\t');

            sb.Append((scanPrecursorMass - peptideWithSetModifications.MonoisotopicMass).ToString("F5", CultureInfo.InvariantCulture) + '\t');

            return sb.ToString();
        }
    }
}