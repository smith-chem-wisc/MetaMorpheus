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
    internal class PSMwithPeptide
    {
        private NewPsm newPsm;
        private PeptideWithSetModifications peptideWithSetModifications;

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

        public PSMwithPeptide(NewPsm newPsm, PeptideWithSetModifications peptideWithSetModifications, double fragmentTolerance, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile)
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

                Dictionary<ProductType, Dictionary<double, double>> MatchedIonDict = new Dictionary<ProductType, Dictionary<double, double>>();
                foreach (var huh in allProductTypes)
                    MatchedIonDict.Add(huh, MatchIons(theScan, fragmentTolerance, peptideWithSetModifications.FastUnsortedProductMasses(new List<ProductType>() { huh })));

                newPsm.matchedIonsList = MatchedIonDict;
            }

            if (newPsm.LocalizedScores == null)
            {
                theScan = myMsDataFile.GetOneBasedScan(newPsm.scanNumber);
                List<double> localizedScores = new List<double>();
                for (int indexToLocalize = 0; indexToLocalize < peptideWithSetModifications.Length; indexToLocalize++)
                {
                    PeptideWithSetModifications localizedPeptide = peptideWithSetModifications.Localize(indexToLocalize, scanPrecursorMass - peptideWithSetModifications.MonoisotopicMass);
                    var hm = MatchIons(theScan, fragmentTolerance, localizedPeptide.FastUnsortedProductMasses(allProductTypes));
                    localizedScores.Add(hm.Count + hm.Select(b => b.Value / theScan.TotalIonCurrent).Sum());
                }
                newPsm.LocalizedScores = localizedScores;
            }
        }

        private static Dictionary<double, double> MatchIons(IMsDataScan<IMzSpectrum<MzPeak>> thisScan, double product_mass_tolerance_value, double[] theoretical_product_masses_for_this_peptide)
        {
            Dictionary<double, double> matchedIonsList = new Dictionary<double, double>();
            if (theoretical_product_masses_for_this_peptide.Length == 0)
                return matchedIonsList;
            int MatchingProductsHere = 0;
            double MatchingIntensityHere = 0;

            List<Tuple<double, double>> listOfMzSandMasses = theoretical_product_masses_for_this_peptide.Select(b => new Tuple<double, double>(b.ToMassToChargeRatio(1), b)).ToList();

            var theoretical_product_mzs_for_this_peptide = listOfMzSandMasses.OrderBy(b => b.Item1).Select(b => b.Item1).ToArray();
            var masses_all = listOfMzSandMasses.OrderBy(b => b.Item1).Select(b => b.Item2).ToArray();
            var TotalProductsHere = theoretical_product_mzs_for_this_peptide.Count();
            int theoreticalLeft = TotalProductsHere;
            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.xArray;
            double[] experimental_intensities = thisScan.MassSpectrum.yArray;
            int num_experimental_peaks = experimental_mzs.Length;

            int theoreticalIndex = 0;
            double nextTheoreticalMZ;
            double mass_difference;
            double currentExperimentalMZ;
            double nextTheoreticalMass;
            for (int i = 0; i < num_experimental_peaks; i++)
            {
                currentExperimentalMZ = experimental_mzs[i];
                nextTheoreticalMZ = theoretical_product_mzs_for_this_peptide[theoreticalIndex];
                nextTheoreticalMass = masses_all[theoreticalIndex];
                mass_difference = currentExperimentalMZ - nextTheoreticalMZ;
                if (Math.Abs(mass_difference) <= product_mass_tolerance_value)
                {
                    MatchingProductsHere++;
                    MatchingIntensityHere += experimental_intensities[i];
                    matchedIonsList.Add(nextTheoreticalMass, experimental_intensities[i]);
                }
                else if (currentExperimentalMZ < nextTheoreticalMZ)
                    continue;
                i--;
                // Passed a theoretical! Move counter forward
                theoreticalIndex++;
                if (theoreticalIndex == TotalProductsHere)
                    break;
            }
            return matchedIonsList;
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