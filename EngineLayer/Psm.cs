using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class Psm
    {

        #region Private Fields

        private const double tolInDaForPreferringHavingMods = 0.03;
        private Dictionary<CompactPeptide, Tuple<int, HashSet<PeptideWithSetModifications>>> compactPeptides = new Dictionary<CompactPeptide, Tuple<int, HashSet<PeptideWithSetModifications>>>();

        #endregion Private Fields

        #region Public Constructors

        public Psm(CompactPeptide peptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan)
        {
            this.Score = score;
            this.ScanIndex = scanIndex;
            this.FullFilePath = scan.FullFilePath;
            this.ScanNumber = scan.TheScan.OneBasedScanNumber;
            this.PrecursorScanNumber = scan.TheScan.OneBasedPrecursorScanNumber;
            this.ScanRetentionTime = scan.TheScan.RetentionTime;
            this.ScanExperimentalPeaks = scan.TheScan.MassSpectrum.Size;
            this.TotalIonCurrent = scan.TheScan.TotalIonCurrent;
            this.ScanPrecursorCharge = scan.PrecursorCharge;
            this.ScanPrecursorMonoisotopicPeak = scan.PrecursorMonoisotopicPeak;
            this.ScanPrecursorMass = scan.PrecursorMass;
            Add(peptide, notch);
        }

        #endregion Public Constructors

        #region Public Properties

        public int ScanNumber { get; }
        public int PrecursorScanNumber { get; }
        public double ScanRetentionTime { get; }
        public int ScanExperimentalPeaks { get; }
        public double TotalIonCurrent { get; }
        public int ScanPrecursorCharge { get; }
        public IMzPeak ScanPrecursorMonoisotopicPeak { get; }
        public double ScanPrecursorMass { get; }
        public string FullFilePath { get; }
        public int ScanIndex { get; }

        public IEnumerable<KeyValuePair<CompactPeptide, Tuple<int, HashSet<PeptideWithSetModifications>>>> CompactPeptides { get { return compactPeptides.AsEnumerable(); } }

        public int NumAmbiguous { get { return compactPeptides.Count; } }

        public FdrInfo FdrInfo { get; private set; }
        public double Score { get; private set; }

        public LocalizationResults LocalizationResults { get; internal set; }
        public double QuantIntensity { get; set; }

        public ProteinLinkedInfo MostProbableProteinInfo { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public static double MatchIons(IMsDataScan<IMzSpectrum<IMzPeak>> thisScan, Tolerance productMassTolerance, double[] sorted_theoretical_product_masses_for_this_peptide, double[] matchedIonMassesListPositiveIsMatch, bool addComp, double precursorMass, List<ProductType> lp)
        {
            var TotalProductsHere = sorted_theoretical_product_masses_for_this_peptide.Length;
            if (TotalProductsHere == 0)
                return 0;
            int MatchingProductsHere = 0;
            double MatchingIntensityHere = 0;

            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.XArray;
            double[] experimental_intensities = thisScan.MassSpectrum.YArray;

            if (addComp)
            {
                List<MzPeak> complementaryPeaks = new List<MzPeak>();

                //keep original peeks
                for (int i = 0; i < experimental_mzs.Length; i++)
                {
                    complementaryPeaks.Add(new MzPeak(experimental_mzs[i], experimental_intensities[i]));
                }
                //If HCD
                if (lp.Contains(ProductType.B) | lp.Contains(ProductType.Y))
                {
                    for (int i = 0; i < experimental_mzs.Length; i++)
                    {
                        complementaryPeaks.Add(new MzPeak((precursorMass - experimental_mzs[i] + Constants.protonMass * 2), (experimental_intensities[i] / 100)));
                    }
                }
                //If ETD
                if (lp.Contains(ProductType.C) | lp.Contains(ProductType.Zdot))
                {
                    for (int i = 0; i < experimental_mzs.Length; i++)
                    {
                        complementaryPeaks.Add(new MzPeak((precursorMass - experimental_mzs[i] + Constants.protonMass * 3), (experimental_intensities[i] / 100)));
                    }
                }

                IEnumerable<MzPeak> sortedPeaksMZ = complementaryPeaks.OrderBy(x => x.Mz);
                experimental_mzs = new double[sortedPeaksMZ.Count()];
                experimental_intensities = new double[sortedPeaksMZ.Count()];
                int index = 0;
                foreach (MzPeak peak in sortedPeaksMZ)
                {
                    experimental_mzs[index] = peak.Mz;
                    experimental_intensities[index] = peak.Intensity;
                    index++;
                }
            }

            int num_experimental_peaks = experimental_mzs.Length;

            int currentTheoreticalIndex = -1;
            double currentTheoreticalMass;
            do
            {
                currentTheoreticalIndex++;
                currentTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[currentTheoreticalIndex];
            } while (double.IsNaN(currentTheoreticalMass) && currentTheoreticalIndex < sorted_theoretical_product_masses_for_this_peptide.Length - 1);

            if (double.IsNaN(currentTheoreticalMass))
                return 0;

            double currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;

            int testTheoreticalIndex;
            double testTheoreticalMZ;
            double testTheoreticalMass;
            // Loop over all experimenal indices
            for (int experimentalIndex = 0; experimentalIndex < num_experimental_peaks; experimentalIndex++)
            {
                double currentExperimentalMZ = experimental_mzs[experimentalIndex];
                // If found match
                if (productMassTolerance.Within(currentExperimentalMZ, currentTheoreticalMz))
                {
                    MatchingProductsHere++;
                    MatchingIntensityHere += experimental_intensities[experimentalIndex];
                    matchedIonMassesListPositiveIsMatch[currentTheoreticalIndex] = currentTheoreticalMass;
                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                        break;
                    currentTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;
                }
                // Else if for sure did not reach the next theoretical yet, move to next experimental
                else if (currentExperimentalMZ < currentTheoreticalMz)
                    continue;
                // Else if for sure passed a theoretical
                else
                {
                    // Mark the theoretical as missed
                    matchedIonMassesListPositiveIsMatch[currentTheoreticalIndex] = -currentTheoreticalMass;

                    // Move on to next index and never come back!
                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                        break;
                    currentTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;

                    // Start with the current ones
                    testTheoreticalIndex = currentTheoreticalIndex;
                    testTheoreticalMZ = currentTheoreticalMz;
                    testTheoreticalMass = currentTheoreticalMass;
                    // Mark the skipped theoreticals as not found. The last one is not for sure, might be flipped!
                    while (currentExperimentalMZ > testTheoreticalMZ)
                    {
                        matchedIonMassesListPositiveIsMatch[testTheoreticalIndex] = -currentTheoreticalMass;
                        // Store old info for possible reuse
                        currentTheoreticalMass = testTheoreticalMass;
                        currentTheoreticalMz = testTheoreticalMZ;
                        currentTheoreticalIndex = testTheoreticalIndex;

                        // Update test stuff!
                        testTheoreticalIndex++;
                        if (testTheoreticalIndex == TotalProductsHere)
                            break;
                        testTheoreticalMass = sorted_theoretical_product_masses_for_this_peptide[testTheoreticalIndex];
                        testTheoreticalMZ = testTheoreticalMass + Constants.protonMass;
                    }

                    experimentalIndex--;
                }
            }
            return MatchingProductsHere + MatchingIntensityHere / thisScan.TotalIonCurrent;
        }

        public static string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("File Name");
            sb.Append('\t' + "Scan Number");
            sb.Append('\t' + "Scan Retention Time");
            sb.Append('\t' + "Num Experimental Peaks");
            sb.Append('\t' + "Total Ion Current");
            sb.Append('\t' + "Precursor Scan Number");
            sb.Append('\t' + "Precursor Charge");
            sb.Append('\t' + "Precursor MZ");
            sb.Append('\t' + "Precursor Intensity");
            sb.Append('\t' + "Precursor Mass");
            sb.Append('\t' + "Score");
            sb.Append('\t' + "Compact Peptide Matches");

            // Single info, common for all peptides/proteins in most cases...
            sb.Append('\t' + "Notch");
            sb.Append('\t' + "Base Sequence");
            sb.Append('\t' + "Full Sequence");
            sb.Append('\t' + "Num Variable Mods");
            sb.Append('\t' + "Missed Cleavages");
            sb.Append('\t' + "Peptide Monoisotopic Mass");
            sb.Append('\t' + "Mass Diff (Da)");
            sb.Append('\t' + "Mass Diff (ppm)");
            sb.Append('\t' + "Decoy/Contaminant/Target");

            // Could have MANY options
            sb.Append('\t' + "Identical Sequence Ambiguity");
            sb.Append('\t' + "Protein Accession");
            sb.Append('\t' + "Protein Name");
            sb.Append('\t' + "Gene Name");
            sb.Append('\t' + "Peptide Description");
            sb.Append('\t' + "Start and End Residues In Protein");
            sb.Append('\t' + "Previous Amino Acid");
            sb.Append('\t' + "Next Amino Acid");

            sb.Append('\t' + LocalizationResults.GetTabSeparatedHeader());
            sb.Append('\t' + "Improvement Possible");

            sb.Append('\t' + "Cumulative Target");
            sb.Append('\t' + "Cumulative Decoy");
            sb.Append('\t' + "QValue");
            sb.Append('\t' + "Cumulative Target Notch");
            sb.Append('\t' + "Cumulative Decoy Notch");
            sb.Append('\t' + "QValue Notch");

            return sb.ToString();
        }

        public void Replace(CompactPeptide correspondingCompactPeptide, double score, int v)
        {
            compactPeptides = new Dictionary<CompactPeptide, Tuple<int, HashSet<PeptideWithSetModifications>>>
            {
                { correspondingCompactPeptide, new  Tuple<int, HashSet<PeptideWithSetModifications>>(v,null)}
            };
            Score = score;
        }

        public void ResolveProteinsAndMostProbablePeptide(Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> matching)
        {
            foreach (var ok in compactPeptides.Keys.ToList())
            {
                compactPeptides[ok] = new Tuple<int, HashSet<PeptideWithSetModifications>>(compactPeptides[ok].Item1, matching[ok]);

                var candidatePli = new ProteinLinkedInfo(matching[ok], compactPeptides[ok].Item1);
                if (MostProbableProteinInfo == null || FirstIsPreferable(candidatePli, MostProbableProteinInfo))
                    MostProbableProteinInfo = candidatePli;
            }
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            sb.Append(Path.GetFileNameWithoutExtension(FullFilePath));
            sb.Append('\t' + ScanNumber.ToString(CultureInfo.InvariantCulture));
            sb.Append('\t' + ScanRetentionTime.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + ScanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + TotalIonCurrent.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + PrecursorScanNumber.ToString(CultureInfo.InvariantCulture));
            sb.Append('\t' + ScanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + ScanPrecursorMonoisotopicPeak.Mz.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + ScanPrecursorMonoisotopicPeak.Intensity.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + ScanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + Score.ToString("F3", CultureInfo.InvariantCulture));
            sb.Append('\t' + NumAmbiguous.ToString("F5", CultureInfo.InvariantCulture));

            var firstNotch = compactPeptides.First().Value.Item1;
            if (compactPeptides.All(b => b.Value.Item1.Equals(firstNotch)))
                sb.Append("\t" + firstNotch.ToString(CultureInfo.InvariantCulture));
            else
            {
                var s = string.Join(" or ", compactPeptides.Select(b => b.Value.Item1.ToString(CultureInfo.InvariantCulture)));
                if (s.Length > 32000)
                    s = "too many";
                sb.Append("\t" + s);
            }

            // These assume that the peptides int each hashset share identical infos
            if (compactPeptides.First().Value.Item2 != null)
            {
                {
                    var firstBaseSeq = compactPeptides.First().Value.Item2.First().BaseSequence;
                    if (compactPeptides.All(b => b.Value.Item2.All(c => c.BaseSequence.Equals(firstBaseSeq))))
                        sb.Append("\t" + firstBaseSeq);
                    else
                    {
                        var s = string.Join(" or ", compactPeptides.Select(b => b.Value.Item2.First().BaseSequence));
                        if (s.Length > 32000)
                            s = "too many";
                        sb.Append("\t" + s);
                    }
                }
                {
                    var firstFullSequence = compactPeptides.First().Value.Item2.First().Sequence;
                    if (compactPeptides.All(b => b.Value.Item2.All(c => c.Sequence.Equals(firstFullSequence))))
                        sb.Append("\t" + firstFullSequence);
                    else
                    {
                        var s = string.Join(" or ", compactPeptides.Select(b => b.Value.Item2.First().Sequence));
                        if (s.Length > 32000)
                            s = "too many";
                        sb.Append("\t" + s);
                    }
                }
                {
                    var firstNumVariableMods = compactPeptides.First().Value.Item2.First().NumVariableMods;
                    if (compactPeptides.All(b => b.Value.Item2.All(c => c.NumVariableMods.Equals(firstNumVariableMods))))
                        sb.Append("\t" + firstNumVariableMods.ToString(CultureInfo.InvariantCulture));
                    else
                    {
                        var s = string.Join(" or ", compactPeptides.Select(b => b.Value.Item2.First().NumVariableMods.ToString(CultureInfo.InvariantCulture)));
                        if (s.Length > 32000)
                            s = "too many";
                        sb.Append("\t" + s);
                    }
                }
                {
                    var firstMissedCleavages = compactPeptides.First().Value.Item2.First().MissedCleavages;
                    if (compactPeptides.All(b => b.Value.Item2.All(c => c.MissedCleavages.Equals(firstMissedCleavages))))
                        sb.Append("\t" + firstMissedCleavages.ToString(CultureInfo.InvariantCulture));
                    else
                    {
                        var s = string.Join(" or ", compactPeptides.Select(b => b.Value.Item2.First().MissedCleavages.ToString(CultureInfo.InvariantCulture)));
                        if (s.Length > 32000)
                            s = "too many";
                        sb.Append("\t" + s);
                    }
                }
                {
                    var firstMass = compactPeptides.First().Value.Item2.First().MonoisotopicMass;
                    if (compactPeptides.All(b => b.Value.Item2.All(c => c.MonoisotopicMass.Equals(firstMass))))
                    {
                        sb.Append("\t" + firstMass.ToString("F5", CultureInfo.InvariantCulture));
                        sb.Append("\t" + (ScanPrecursorMass - firstMass).ToString("F5", CultureInfo.InvariantCulture));
                        sb.Append("\t" + ((ScanPrecursorMass - firstMass) / firstMass * 1e6).ToString("F5", CultureInfo.InvariantCulture));
                    }
                    else
                    {
                        var s = string.Join(" or ", compactPeptides.Select(b => b.Value.Item2.First().MonoisotopicMass.ToString(CultureInfo.InvariantCulture)));
                        if (s.Length > 32000)
                            s = "too many";
                        sb.Append("\t" + s);
                        s = string.Join(" or ", compactPeptides.Select(b => (ScanPrecursorMass - b.Value.Item2.First().MonoisotopicMass).ToString("F5", CultureInfo.InvariantCulture)));
                        if (s.Length > 32000)
                            s = "too many";
                        sb.Append("\t" + s);
                        s = string.Join(" or ", compactPeptides.Select(b => ((ScanPrecursorMass - b.Value.Item2.First().MonoisotopicMass) / b.Value.Item2.First().MonoisotopicMass * 1e6).ToString("F5", CultureInfo.InvariantCulture)));
                        if (s.Length > 32000)
                            s = "too many";
                        sb.Append("\t" + s);
                    }
                }

                // Unambiguous
                if (compactPeptides.Any(b => b.Value.Item2.Any(c => c.Protein.IsDecoy)))
                    sb.Append("\t" + "D");
                else if (compactPeptides.Any(b => b.Value.Item2.Any(c => c.Protein.IsContaminant)))
                    sb.Append("\t" + "C");
                else
                    sb.Append("\t" + "T");

                sb.Append("\t" + string.Join(" or ", compactPeptides.Select(b => b.Value.Item2.Count.ToString(CultureInfo.InvariantCulture))));
            }
            else
            {
                sb.Append('\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " ");
            }

            // These outputs require a single compact peptide. These we expect may be different for the same compact peptide
            if (compactPeptides.Count == 1 && compactPeptides.First().Value.Item2 != null)
            {
                var PeptidesWithSetModifications = compactPeptides.First().Value.Item2;

                {
                    var first = PeptidesWithSetModifications.First().Protein.Accession;
                    if (PeptidesWithSetModifications.All(b => b.Protein.Accession.Equals(first)))
                    {
                        sb.Append("\t" + first);
                    }
                    else
                    {
                        var st = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.Protein.Accession));
                        if (st.Length > 32000)
                            st = "too many";
                        sb.Append("\t" + st);
                    }
                }
                {
                    var first = PeptidesWithSetModifications.First().Protein.FullName;
                    if ((first == null && PeptidesWithSetModifications.All(b => b.Protein.FullName == null)) || (first != null && PeptidesWithSetModifications.All(b => first.Equals(b.Protein.FullName))))
                    {
                        sb.Append("\t" + first);
                    }
                    else
                    {
                        var st = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.Protein.FullName));
                        if (st.Length > 32000)
                            st = "too many";
                        sb.Append("\t" + st);
                    }
                }
                {
                    var first = string.Join(",", PeptidesWithSetModifications.First().Protein.GeneNames.Select(c => c.Item1 + ":" + c.Item2));
                    if (PeptidesWithSetModifications.All(b => first.Equals(string.Join(",", b.Protein.GeneNames.Select(c => c.Item1 + ":" + c.Item2)))))
                    {
                        sb.Append("\t" + first);
                    }
                    else
                    {
                        var st = string.Join(" or ", PeptidesWithSetModifications.Select(b => string.Join(",", b.Protein.GeneNames.Select(c => c.Item1 + ":" + c.Item2))));
                        if (st.Length > 32000)
                            st = "too many";
                        sb.Append("\t" + st);
                    }
                }
                {
                    var first = PeptidesWithSetModifications.First().PeptideDescription;
                    if (PeptidesWithSetModifications.All(b => b.PeptideDescription.Equals(first)))
                    {
                        sb.Append("\t" + first);
                    }
                    else
                    {
                        var st = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.PeptideDescription));
                        if (st.Length > 32000)
                            st = "too many";
                        sb.Append("\t" + st);
                    }
                }
                {
                    var first = PeptidesWithSetModifications.First().OneBasedStartResidueInProtein.ToString(CultureInfo.InvariantCulture) + " to " + PeptidesWithSetModifications.First().OneBasedEndResidueInProtein.ToString(CultureInfo.InvariantCulture) + "]";
                    if (PeptidesWithSetModifications.All(b => first.Equals("[" + b.OneBasedStartResidueInProtein.ToString(CultureInfo.InvariantCulture) + " to " + b.OneBasedEndResidueInProtein.ToString(CultureInfo.InvariantCulture) + "]")))
                    {
                        sb.Append("\t" + first);
                    }
                    else
                    {
                        var st = string.Join(" or ", PeptidesWithSetModifications.Select(b => ("[" + b.OneBasedStartResidueInProtein.ToString(CultureInfo.InvariantCulture) + " to " + b.OneBasedEndResidueInProtein.ToString(CultureInfo.InvariantCulture) + "]")));
                        if (st.Length > 32000)
                            st = "too many";
                        sb.Append("\t" + st);
                    }
                }
                {
                    var first = PeptidesWithSetModifications.First().PreviousAminoAcid;
                    if (PeptidesWithSetModifications.All(b => b.PreviousAminoAcid.Equals(first)))
                    {
                        sb.Append("\t" + first);
                    }
                    else
                    {
                        var st = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.PreviousAminoAcid));
                        if (st.Length > 32000)
                            st = "too many";
                        sb.Append("\t" + st);
                    }
                }
                {
                    var first = PeptidesWithSetModifications.First().NextAminoAcid;
                    if (PeptidesWithSetModifications.All(b => b.NextAminoAcid.Equals(first)))
                    {
                        sb.Append("\t" + first);
                    }
                    else
                    {
                        var st = string.Join(" or ", PeptidesWithSetModifications.Select(b => b.NextAminoAcid));
                        if (st.Length > 32000)
                            st = "too many";
                        sb.Append("\t" + st);
                    }
                }
            }
            else
            {
                sb.Append('\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " ");
            }

            if (LocalizationResults != null)
            {
                sb.Append('\t' + LocalizationResults.ToString());
                sb.Append('\t' + (LocalizationResults.LocalizedScores.Max() - Score).ToString("F3", CultureInfo.InvariantCulture));
            }
            else
            {
                sb.Append('\t' + " " + '\t' + " " + '\t' + " " + '\t' + " ");
            }

            if (FdrInfo != null)
            {
                sb.Append('\t' + FdrInfo.cumulativeTarget.ToString(CultureInfo.InvariantCulture));
                sb.Append('\t' + FdrInfo.cumulativeDecoy.ToString(CultureInfo.InvariantCulture));
                sb.Append('\t' + FdrInfo.QValue.ToString("F6", CultureInfo.InvariantCulture));
                sb.Append('\t' + FdrInfo.cumulativeTargetNotch.ToString(CultureInfo.InvariantCulture));
                sb.Append('\t' + FdrInfo.cumulativeDecoyNotch.ToString(CultureInfo.InvariantCulture));
                sb.Append('\t' + FdrInfo.QValueNotch.ToString("F6", CultureInfo.InvariantCulture));
            }
            else
                sb.Append('\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " ");

            return sb.ToString();
        }

        public void SetFdrValues(int cumulativeTarget, int cumulativeDecoy, double tempQValue, int cumulativeTargetNotch, int cumulativeDecoyNotch, double tempQValueNotch)
        {
            FdrInfo = new FdrInfo()
            {
                cumulativeTarget = cumulativeTarget,
                cumulativeDecoy = cumulativeDecoy,
                QValue = tempQValue,
                cumulativeTargetNotch = cumulativeTargetNotch,
                cumulativeDecoyNotch = cumulativeDecoyNotch,
                QValueNotch = tempQValueNotch
            };
        }

        #endregion Public Methods

        #region Internal Methods

        internal void Add(CompactPeptide compactPeptide, int v)
        {
            compactPeptides[compactPeptide] = new Tuple<int, HashSet<PeptideWithSetModifications>>(v, null);
        }

        internal void Add(Psm psmParent)
        {
            foreach (var kvp in psmParent.compactPeptides)
                Add(kvp.Key, kvp.Value.Item1);
        }

        #endregion Internal Methods

        #region Private Methods

        private bool FirstIsPreferable(ProteinLinkedInfo firstPli, ProteinLinkedInfo secondPli)
        {
            if (firstPli.IsDecoy && !secondPli.IsDecoy)
                return true;
            if (!firstPli.IsDecoy && secondPli.IsDecoy)
                return false;

            // Score is same, need to see if accepts and if prefer the new one
            var first = firstPli.PeptidesWithSetModifications.First();
            var second = secondPli.PeptidesWithSetModifications.First();

            // Prefer to be at zero rather than fewer modifications
            if ((Math.Abs(first.MonoisotopicMass - ScanPrecursorMass) < tolInDaForPreferringHavingMods)
                && (Math.Abs(second.MonoisotopicMass - ScanPrecursorMass) > tolInDaForPreferringHavingMods))
                return true;
            if ((Math.Abs(second.MonoisotopicMass - ScanPrecursorMass) < tolInDaForPreferringHavingMods)
                && (Math.Abs(first.MonoisotopicMass - ScanPrecursorMass) > tolInDaForPreferringHavingMods))
                return false;

            //int firstVarMods = first.allModsOneIsNterminus.Count(b => variableMods.Contains(b.Value));
            //int secondVarMods = second.allModsOneIsNterminus.Count(b => variableMods.Contains(b.Value));

            // Want the lowest number of localizeable mods!!! Even at the expense of many variable and fixed mods.
            //if (first.NumMods - firstVarMods - first.numFixedMods < second.NumMods - secondVarMods - second.numFixedMods)
            //    return true;
            //if (first.NumMods - firstVarMods - first.numFixedMods > second.NumMods - secondVarMods - second.numFixedMods)
            //    return false;

            // If have same number of localizeable mods, pick the lowest number of localizeable + variable mods
            if (first.NumMods - first.numFixedMods < second.NumMods - second.numFixedMods)
                return true;
            if (first.NumMods - first.numFixedMods > second.NumMods - second.numFixedMods)
                return false;

            // If have same numbers of localizeable and variable mods, prefer not to have substitutions and removals!
            int firstNumRazor = first.allModsOneIsNterminus.Count(b => b.Value.modificationType.Equals("substitution") || b.Value.modificationType.Equals("missing") || b.Value.modificationType.Equals("trickySubstitution"));
            int secondNumRazor = second.allModsOneIsNterminus.Count(b => b.Value.modificationType.Equals("substitution") || b.Value.modificationType.Equals("missing") || b.Value.modificationType.Equals("trickySubstitution"));

            if (firstNumRazor < secondNumRazor)
                return true;
            if (firstNumRazor > secondNumRazor)
                return false;

            return true;
        }

        #endregion Private Methods

    }
}