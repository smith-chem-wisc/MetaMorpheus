using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace MetaMorpheus
{
    public class PeptideSpectrumMatch
    {
        public PeptideWithSetModifications Peptide { get; private set; }

        public double PrecursorMassErrorDa { get; private set; }

        public double PrecursorMassErrorPpm { get; private set; }

        // THE SIX OUTPUTS IN PSMS.TSV
        public int MatchingProducts { get; private set; }

        public int TotalProducts { get; private set; }
        public double MatchingProductsFraction { get; private set; }
        public double MatchingIntensity { get; private set; }
        public double MatchingIntensityFraction { get; private set; }
        public double MetaMorpheusScore { get; private set; }

        //private MassTolerance productMassTolerance;

        public Dictionary<ProductType, List<double>> MatchedIonsList = new Dictionary<ProductType, List<double>>();

        public bool isDecoy { get; private set; }

        public string spectrumFilename { get; private set; }
        public string spectrumID { get; private set; }
        public int spectrumPrecursorCharge { get; private set; }
        public double spectrumPrecursorMZ { get; private set; }
        public int SpectrumNumber { get; private set; }
        public int SpectrumIndexHere { get; internal set; }

        private string SpectrumTitle;
        private double RetentionTimeMinutes;
        private double PrecursorIntensity;
        private double PrecursorMass;
        private int totalExperimentalPeaks;
        private double TotalIntensity;

        //public PeptideSpectrumMatch(TandemMassSpectrum Spectrum, PeptideWithSetModifications Peptide, MassTolerance productMassTolerance, int needToMatch, List<ProductType> productTypes, int indexHere)
        //{
        //    this.Peptide = Peptide;

        //    spectrumFilename = Spectrum.Filename;
        //    spectrumID = Spectrum.SpectrumId;
        //    spectrumPrecursorCharge = Spectrum.PrecursorCharge;
        //    spectrumPrecursorMZ = Spectrum.PrecursorMZ;
        //    SpectrumNumber = Spectrum.SpectrumNumber;

        //    this.SpectrumIndexHere = indexHere;

        //    PrecursorMassErrorDa = Spectrum.PrecursorMass - Peptide.MonoisotopicMass;
        //    PrecursorMassErrorPpm = PrecursorMassErrorDa / Peptide.MonoisotopicMass * 1e6;
        //    this.productMassTolerance = productMassTolerance;
        //    this.productTypes = productTypes;

        //    SpectrumTitle = Spectrum.SpectrumTitle;
        //    RetentionTimeMinutes = Spectrum.RetentionTimeMinutes;
        //    PrecursorIntensity = Spectrum.PrecursorIntensity;
        //    PrecursorMass = Spectrum.PrecursorMass;
        //    totalExperimentalPeaks = Spectrum.mzs.Count();
        //    TotalIntensity = Spectrum.TotalIntensity;

        //    //var hm = Peptide.GetProductMassesSingleArray(productTypes);

        //    //var ok = ScoreMatch(hm, needToMatch, ProductType.none, Spectrum);

        //    //if (ok != null)
        //    //{
        //    //    //TotalProducts = hm.Count();
        //    //    MatchingProducts = ok.Item1;
        //    //    MatchingProductsFraction = (double)MatchingProducts / TotalProducts;
        //    //    MatchingIntensity = ok.Item2;
        //    //    MatchingIntensityFraction = MatchingIntensity / Spectrum.TotalIntensity;
        //    //    MetaMorpheusScore = MatchingProducts + MatchingIntensityFraction;
        //    //}
        //}

        public PeptideSpectrumMatch(bool isDecoy, double PrecursorMassErrorDa, List<double> LocalizedScores, double MetaMorpheusScore, PeptideWithSetModifications peptide)
        {
            this.isDecoy = isDecoy;
            this.PrecursorMassErrorDa = PrecursorMassErrorDa;
            this.MetaMorpheusScore = MetaMorpheusScore;
            this.LocalizedScores = LocalizedScores;
            this.Peptide = peptide;
        }

        //private Tuple<int, double> ScoreMatch(double[] theoretical_product_mzs_for_this_peptide, int needToMatch, ProductType productTypeToDetailList, TandemMassSpectrum Spectrum)
        //{
        //    int MatchingProductsHere = 0;
        //    double MatchingIntensityHere = 0;

        //    var TotalProductsHere = theoretical_product_mzs_for_this_peptide.Count();
        //    if (TotalProductsHere < needToMatch)
        //        return null;
        //    int theoreticalLeft = TotalProductsHere;
        //    // speed optimizations
        //    List<double> experimental_masses = Spectrum.mzs;
        //    double[] experimental_intensities = Spectrum.Intensities;
        //    int num_experimental_peaks = experimental_masses.Count;
        //    double product_mass_tolerance_value = productMassTolerance.Value;
        //    MassToleranceUnits product_mass_tolerance_units = productMassTolerance.Units;
        //    if (productTypeToDetailList != ProductType.none)
        //        MatchedIonsList[productTypeToDetailList] = new List<double>();

        //    int theoreticalIndex = 0;
        //    double nextTheoretical;
        //    double mass_difference;
        //    double currentExperimentalMass;
        //    for (int i = 0; i < num_experimental_peaks; i++)
        //    {
        //        currentExperimentalMass = experimental_masses[i];
        //        nextTheoretical = theoretical_product_mzs_for_this_peptide[theoreticalIndex];
        //        mass_difference = currentExperimentalMass - nextTheoretical;
        //        if (product_mass_tolerance_units == MassToleranceUnits.ppm)
        //            mass_difference = mass_difference / nextTheoretical * 1e6;
        //        if (Math.Abs(mass_difference) <= product_mass_tolerance_value)
        //        {
        //            MatchingProductsHere++;
        //            MatchingIntensityHere += experimental_intensities[i];
        //            if (productTypeToDetailList != ProductType.none)
        //                MatchedIonsList[productTypeToDetailList].Add(-nextTheoretical);
        //            //Console.WriteLine(theoreticalLeft + " " + nextTheoretical + " + ");
        //        }
        //        else if (currentExperimentalMass < nextTheoretical)
        //            continue;
        //        else if (productTypeToDetailList != ProductType.none)
        //            MatchedIonsList[productTypeToDetailList].Add(nextTheoretical);
        //        i--;
        //        theoreticalLeft--;
        //        if (needToMatch > theoreticalLeft + MatchingProductsHere)
        //            return null;
        //        // Passed a theoretical! Move counter forward
        //        theoreticalIndex++;
        //        if (theoreticalIndex == TotalProductsHere)
        //            break;
        //    }
        //    return new Tuple<int, double>(MatchingProductsHere, MatchingIntensityHere);
        //}

        // left is new, right is current best
        public static int DescendingMetaMorpheusScoreComparison(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            // Higher score is better, so have a negative
            int comparison = -(left.MetaMorpheusScore.CompareTo(right.MetaMorpheusScore));
            if (comparison != 0)
                return comparison;
            // Low precursor mass is better, so no negative
            if (Math.Abs(left.PrecursorMassErrorDa) <= 0.5)
                comparison = Math.Abs(left.PrecursorMassErrorDa).CompareTo(Math.Abs(right.PrecursorMassErrorDa));
            if (comparison != 0)
                return comparison;
            // Low number of variable ptms is better, so no negative
            comparison = left.Peptide.numVariableMods.CompareTo(right.Peptide.numVariableMods);
            if (comparison != 0)
                return comparison;
            return right.isDecoy.CompareTo(left.isDecoy);
        }

        public static readonly string Header = "Filename\tSpectrum Number\tSpectrum ID\tSpectrum Title\tRetention Time (minutes)\tPrecursor m/z\tPrecursor Intensity\tPrecursor Charge\tPrecursor Mass (Da)\tExperimental Peaks\tTotal Intensity"
            + "\tPeptide Sequence\tBase Peptide Sequence\tProtein Description\tPeptide Description\tStart Residue Number\tStop Residue Number\tMissed Cleavages"
            + "\tTheoretical Mass (Da)\tPrecursor Mass Error (Da)\tPrecursor Mass Error (ppm)"
            + "\tMatching Products\tTotal Products\tRatio of Matching Products\tMatching Intensity\tFraction of Intensity Matching\tMetaMorpheus Score\tIon Matches\tIon Counts\tExtended Scores\tImprovement\tImprovementResidue\tImprovementTerminus";

        private List<ProductType> productTypes;
        public List<double> LocalizedScores;

        //private List<double> theoreticaklMzMatches;

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(spectrumFilename.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(SpectrumNumber.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(spectrumID.ToString(CultureInfo.InvariantCulture) + '\t');
            if (SpectrumTitle == null)
                sb.Append('\t');
            else
                sb.Append(SpectrumTitle.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(RetentionTimeMinutes.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(spectrumPrecursorMZ.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(PrecursorIntensity.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(spectrumPrecursorCharge.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(PrecursorMass.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(totalExperimentalPeaks.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(TotalIntensity.ToString(CultureInfo.InvariantCulture) + '\t');

            sb.Append(Peptide.ExtendedSequence + '\t');
            sb.Append(Peptide.BaseSequence + '\t');
            sb.Append(Peptide.protein.FullDescription + '\t');
            sb.Append(Peptide.PeptideDescription + '\t');
            sb.Append(Peptide.OneBasedStartResidueInProtein.ToString() + '\t');
            sb.Append(Peptide.OneBasedEndResidueInProtein.ToString() + '\t');
            //sb.Append(Peptide.MissedCleavages.ToString() + '\t');
            sb.Append(Peptide.MonoisotopicMass.ToString(CultureInfo.InvariantCulture) + '\t');

            sb.Append(PrecursorMassErrorDa.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(PrecursorMassErrorPpm.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(MatchingProducts.ToString() + '\t');
            sb.Append(TotalProducts.ToString() + '\t');
            sb.Append(MatchingProductsFraction.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(MatchingIntensity.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(MatchingIntensityFraction.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append(MetaMorpheusScore.ToString(CultureInfo.InvariantCulture) + '\t');
            sb.Append("[");
            foreach (var kvp in MatchedIonsList)
                sb.Append("[" + string.Join(",", kvp.Value) + "];");
            sb.Append("]" + '\t');
            sb.Append(string.Join(";", MatchedIonsList.Select(b => b.Value.Where(c => c < 0).Count())) + '\t');
            sb.Append("[" + string.Join(",", LocalizedScores) + "]" + '\t');
            sb.Append((LocalizedScores.Max() - MetaMorpheusScore).ToString() + '\t');
            sb.Append(Peptide[LocalizedScores.IndexOf(LocalizedScores.Max())].ToString() + '\t');
            if (LocalizedScores.IndexOf(LocalizedScores.Max()) == 0)
                sb.Append("N");
            else if (LocalizedScores.Max() == LocalizedScores.Last())
                sb.Append("C");
            else
                sb.Append("");

            return sb.ToString();
        }

        //internal void ComputeIonMatchesAndCounts(TandemMassSpectrum Spectrum)
        //{
        //    // foreach (var pt in productTypes)
        //    //    ScoreMatch(Peptide.GetProductMassesSingleArray(new List<ProductType>() { pt }), 1, pt, Spectrum);
        //}

        //internal void GetLocalizedScores(TandemMassSpectrum Spectrum)
        //{
        //    LocalizedScores = new List<double>();
        //    for (int i = 1; i <= Peptide.Length; i++)
        //    {
        //        //PeptideWithSetModifications ok = new PeptideWithSetModifications(Peptide, i, PrecursorMassErrorDa);
        //        //PeptideSpectrumMatch psm = new PeptideSpectrumMatch(Spectrum, ok, productMassTolerance, 1, productTypes, SpectrumIndexHere);
        //        //LocalizedScores.Add(psm.MetaMorpheusScore);
        //    }
        //}
    }
}