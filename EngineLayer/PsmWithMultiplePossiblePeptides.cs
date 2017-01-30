using Chemistry;
using MassSpectrometry;

using Spectra;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class PsmWithMultiplePossiblePeptides
    {

        #region Public Fields

        public PsmParent newPsm;

        #endregion Public Fields

        #region Public Constructors

        public PsmWithMultiplePossiblePeptides(PsmParent newPsm, HashSet<PeptideWithSetModifications> peptidesWithSetModifications, Tolerance fragmentTolerance, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, List<ProductType> lp)
        {
            this.newPsm = newPsm;
            IsDecoy = peptidesWithSetModifications.Any(b => b.Protein.IsDecoy);
            IsContaminant = peptidesWithSetModifications.Any(b => b.Protein.IsContaminant);
            this.peptidesWithSetModifications = peptidesWithSetModifications;

            var representative = peptidesWithSetModifications.First();

            IMsDataScan<IMzSpectrum<MzPeak>> theScan;
            if (myMsDataFile != null && newPsm.matchedIonsList == null)
            {
                theScan = myMsDataFile.GetOneBasedScan(newPsm.scanNumber);

                var MatchedIonDict = new Dictionary<ProductType, double[]>();
                foreach (var huh in lp)
                {
                    var df = representative.FastSortedProductMasses(new List<ProductType> { huh });
                    double[] matchedIonList = new double[df.Length];
                    MatchIons(theScan, fragmentTolerance, df, matchedIonList);
                    MatchedIonDict.Add(huh, matchedIonList);
                }

                newPsm.matchedIonsList = MatchedIonDict;
            }

            if (myMsDataFile != null && newPsm.LocalizedScores == null)
            {
                theScan = myMsDataFile.GetOneBasedScan(newPsm.scanNumber);
                var localizedScores = new List<double>();
                for (int indexToLocalize = 0; indexToLocalize < representative.Length; indexToLocalize++)
                {
                    PeptideWithSetModifications localizedPeptide = representative.Localize(indexToLocalize, ScanPrecursorMass - representative.MonoisotopicMass);

                    var gg = localizedPeptide.FastSortedProductMasses(lp);
                    double[] matchedIonList = new double[gg.Length];
                    var score = MatchIons(theScan, fragmentTolerance, gg, matchedIonList);
                    localizedScores.Add(score);
                }
                newPsm.LocalizedScores = localizedScores;
            }

            PeptideMonoisotopicMass = representative.MonoisotopicMass;

            FullSequence = representative.Sequence;
            BaseSequence = representative.BaseSequence;
            MissedCleavages = representative.MissedCleavages;
            NumVariableMods = representative.NumVariableMods;
        }

        #endregion Public Constructors

        #region Public Properties

        public HashSet<PeptideWithSetModifications> peptidesWithSetModifications { get; private set; }
        public bool IsDecoy { get; private set; }
        public bool IsContaminant { get; private set; }

        public double Score
        {
            get
            {
                return newPsm.score;
            }
        }

        public double ScanPrecursorMass
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

        public double PeptideMonoisotopicMass { get; private set; }

        public string FullSequence { get; private set; }

        public string BaseSequence { get; private set; }

        public int MissedCleavages { get; private set; }

        public int NumVariableMods { get; private set; }

        public string SequenceWithChemicalFormulas
        {
            get
            {
                return peptidesWithSetModifications.First().SequenceWithChemicalFormulas;
            }
        }

        #endregion Public Properties

        #region Internal Properties

        internal static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append(PsmParent.GetTabSeparatedHeader() + '\t');
                sb.Append("Protein Accession" + '\t');
                sb.Append("Protein FullName" + '\t');
                sb.Append("Peptide Description" + '\t');
                sb.Append("Start and End ResidueInProtein" + '\t');
                sb.Append("BaseSequence" + '\t');
                sb.Append("FullSequence" + '\t');
                sb.Append("numVariableMods" + '\t');
                sb.Append("MissedCleavages" + '\t');
                sb.Append("PeptideMonoisotopicMass" + '\t');
                sb.Append("MassDiff (Da)" + '\t');
                sb.Append("Decoy/Contaminant");
                return sb.ToString();
            }
        }

        #endregion Internal Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();

            sb.Append(newPsm.ToString() + '\t');

            sb.Append(string.Join(" or ", peptidesWithSetModifications.Select(b => b.Protein.Accession)) + "\t");
            sb.Append(string.Join(" or ", peptidesWithSetModifications.Select(b => b.Protein.FullName)) + "\t");
            sb.Append(string.Join(" or ", peptidesWithSetModifications.Select(b => b.PeptideDescription)) + "\t");
            sb.Append(string.Join(" or ", peptidesWithSetModifications.Select(b => "[" + b.OneBasedStartResidueInProtein + " to " + b.OneBasedEndResidueInProtein + "]")) + "\t");
            sb.Append(BaseSequence.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(FullSequence.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(NumVariableMods.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(MissedCleavages.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(PeptideMonoisotopicMass.ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append((ScanPrecursorMass - PeptideMonoisotopicMass).ToString("F5", CultureInfo.InvariantCulture) + '\t');
            sb.Append("" + (IsDecoy ? "D" : " ") + (IsContaminant ? "C" : " "));

            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal static double MatchIons(IMsDataScan<IMzSpectrum<MzPeak>> thisScan, Tolerance product_mass_tolerance_value, double[] sorted_theoretical_product_masses_for_this_peptide, double[] matchedIonMassesList)
        {
            var TotalProductsHere = sorted_theoretical_product_masses_for_this_peptide.Length;
            if (TotalProductsHere == 0)
                return 0;
            int MatchingProductsHere = 0;
            double MatchingIntensityHere = 0;

            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.XArray;
            double[] experimental_intensities = thisScan.MassSpectrum.YArray;
            int num_experimental_peaks = experimental_mzs.Length;

            int theoreticalIndex = 0;
            double nextTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[0];
            double nextTheoreticalMZ = nextTheoreticalMass + Constants.ProtonMass;

            double currentExperimentalMZ;
            for (int i = 0; i < num_experimental_peaks; i++)
            {
                currentExperimentalMZ = experimental_mzs[i];
                if (product_mass_tolerance_value.Within(currentExperimentalMZ, nextTheoreticalMZ))
                {
                    MatchingProductsHere++;
                    MatchingIntensityHere += experimental_intensities[i];
                    matchedIonMassesList[theoreticalIndex] = nextTheoreticalMass;
                }
                else if (currentExperimentalMZ < nextTheoreticalMZ)
                    continue;
                else
                    matchedIonMassesList[theoreticalIndex] = -nextTheoreticalMass;
                i--;
                // Passed a theoretical! Move counter forward
                theoreticalIndex++;
                if (theoreticalIndex == TotalProductsHere)
                    break;
                nextTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[theoreticalIndex];
                nextTheoreticalMZ = nextTheoreticalMass + Constants.ProtonMass;
            }
            return MatchingProductsHere + MatchingIntensityHere / thisScan.TotalIonCurrent;
        }

        #endregion Internal Methods

    }
}