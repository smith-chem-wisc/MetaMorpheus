using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
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

        private Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>> compactPeptides = new Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>();

        #endregion Private Fields

        #region Public Constructors

        public Psm(CompactPeptideBase peptide, int notch, double score, int scanIndex, IScan scan)
        {
            this.ScanIndex = scanIndex;
            this.FullFilePath = scan.FullFilePath;
            this.ScanNumber = scan.OneBasedScanNumber;
            this.PrecursorScanNumber = scan.OneBasedPrecursorScanNumber;
            this.ScanRetentionTime = scan.RetentionTime;
            this.ScanExperimentalPeaks = scan.NumPeaks;
            this.TotalIonCurrent = scan.TotalIonCurrent;
            this.ScanPrecursorCharge = scan.PrecursorCharge;
            this.ScanPrecursorMonoisotopicPeak = scan.PrecursorMonoisotopicPeak;
            this.ScanPrecursorMass = scan.PrecursorMass;
            AddOrReplace(peptide, score, notch);
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

        public IEnumerable<KeyValuePair<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>> CompactPeptides { get { return compactPeptides.AsEnumerable(); } }

        public int NumDifferentCompactPeptides { get { return compactPeptides.Count; } }

        public FdrInfo FdrInfo { get; private set; }
        public double Score { get; private set; }

        public double QuantIntensity { get; set; }

        public ProteinLinkedInfo MostProbableProteinInfo { get; private set; }
        public bool IsDecoy { get; private set; }
        public string FullSequence { get; private set; }
        public int? Notch { get; private set; }
        public string BaseSequence { get; private set; }
        public int? PeptideLength { get; private set; }
        public int? OneBasedStartResidueInProtein { get; private set; }
        public int? OneBasedEndResidueInProtein { get; private set; }
        public double? PeptideMonisotopicMass { get; private set; }
        public int? ProteinLength { get; private set; }
        public List<double> LocalizedScores { get; internal set; }
        public MatchedIonMassesListPositiveIsMatch MatchedIonDictPositiveIsMatch { get; internal set; }
        public string ProteinAccesion { get; private set; }
        public Dictionary<string, int> ModsIdentified { get; private set; }

        #endregion Public Properties

        // id and number of occurences

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
                if (lp.Contains(ProductType.B) || lp.Contains(ProductType.Y))
                {
                    for (int i = 0; i < experimental_mzs.Length; i++)
                    {
                        complementaryPeaks.Add(new MzPeak((precursorMass - experimental_mzs[i] + Constants.protonMass * 2), (experimental_intensities[i] / 100)));
                    }
                }
                //If ETD
                if (lp.Contains(ProductType.C) || lp.Contains(ProductType.Zdot))
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
            sb.Append('\t' + "Notch");
            sb.Append('\t' + "Different Peak Matches");

            sb.Append('\t' + "Peptides Sharing Same Peaks");
            sb.Append('\t' + "Base Sequence");
            sb.Append('\t' + "Full Sequence");
            sb.Append('\t' + "Mods");
            sb.Append('\t' + "Num Variable Mods");
            sb.Append('\t' + "Missed Cleavages");
            sb.Append('\t' + "Peptide Monoisotopic Mass");
            sb.Append('\t' + "Mass Diff (Da)");
            sb.Append('\t' + "Mass Diff (ppm)");
            sb.Append('\t' + "Protein Accession");
            sb.Append('\t' + "Protein Name");
            sb.Append('\t' + "Gene Name");
            sb.Append('\t' + "Contaminant");
            sb.Append('\t' + "Decoy");
            sb.Append('\t' + "Peptide Description");
            sb.Append('\t' + "Start and End Residues In Protein");
            sb.Append('\t' + "Previous Amino Acid");
            sb.Append('\t' + "Next Amino Acid");
            sb.Append('\t' + "Decoy/Contaminant/Target");

            sb.Append('\t' + "Matched Ion Counts");
            sb.Append('\t' + "Matched Ion Masses");

            sb.Append('\t' + "Localized Scores");
            sb.Append('\t' + "Improvement Possible");

            sb.Append('\t' + "Cumulative Target");
            sb.Append('\t' + "Cumulative Decoy");
            sb.Append('\t' + "QValue");
            sb.Append('\t' + "Cumulative Target Notch");
            sb.Append('\t' + "Cumulative Decoy Notch");
            sb.Append('\t' + "QValue Notch");

            return sb.ToString();
        }

        public void AddOrReplace(CompactPeptideBase compactPeptide, double score, int notch)
        {
            if (score - Score > 1e-9)
            {
                compactPeptides = new Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>
                {
                    { compactPeptide, new  Tuple<int, HashSet<PeptideWithSetModifications>>(notch,null)}
                };
                Score = score;
            }
            else if (score - Score > -1e-9)
            {
                compactPeptides[compactPeptide] = new Tuple<int, HashSet<PeptideWithSetModifications>>(notch, null);
            }
        }

        public void CompactCompactPeptides()
        {
            List<Tuple<CompactPeptideBase, int>> cps = new List<Tuple<CompactPeptideBase, int>>();
            foreach (KeyValuePair<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>> kvp in compactPeptides)
            {
                //Change CPWM to reflect actual CP
                    Tuple<CompactPeptideBase, int> tempTuple = new Tuple<CompactPeptideBase, int>(kvp.Key, kvp.Value.Item1);
                    if (!cps.Contains(tempTuple))
                        cps.Add(tempTuple);
            }
            compactPeptides = new Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>();
            foreach (Tuple<CompactPeptideBase, int> cp in cps)
                compactPeptides[cp.Item1] = new Tuple<int, HashSet<PeptideWithSetModifications>>(cp.Item2, null);
        }

        public void MatchToProteinLinkedPeptides(Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching)
        {
            foreach (var cpKey in compactPeptides.Keys.ToList())
            {
                compactPeptides[cpKey] = new Tuple<int, HashSet<PeptideWithSetModifications>>(compactPeptides[cpKey].Item1, matching[cpKey]);
                var candidatePli = new ProteinLinkedInfo(matching[cpKey]);
                if (MostProbableProteinInfo == null || FirstIsPreferable(candidatePli, MostProbableProteinInfo))
                    MostProbableProteinInfo = candidatePli;
            }

            IsDecoy = compactPeptides.Any(b => b.Value.Item2.Any(c => c.Protein.IsDecoy));

            FullSequence = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Sequence)).Item2;

            BaseSequence = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.BaseSequence)).Item2;

            PeptideLength = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Length)).Item2;

            OneBasedStartResidueInProtein = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.OneBasedStartResidueInProtein)).Item2;

            OneBasedEndResidueInProtein = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.OneBasedEndResidueInProtein)).Item2;

            ProteinLength = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.Length)).Item2;

            PeptideMonisotopicMass = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.MonoisotopicMass)).Item2;

            ProteinAccesion = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.Accession)).Item2;

            ModsIdentified = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.allModsOneIsNterminus)).Item2;

            Notch = Resolve(compactPeptides.Select(b => b.Value.Item1)).Item2;
        }

        public bool CompactPeptidesContainsKey(CompactPeptideBase key)
        {
            return compactPeptides.ContainsKey(key);
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
            sb.Append("\t" + Resolve(compactPeptides.Select(b => b.Value.Item1)).Item1); // Notch
            sb.Append('\t' + NumDifferentCompactPeptides.ToString("F5", CultureInfo.InvariantCulture));

            if (compactPeptides.First().Value.Item2 != null)
            {
                sb.Append("\t" + string.Join(" or ", compactPeptides.Select(b => b.Value.Item2.Count.ToString(CultureInfo.InvariantCulture))));

                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.BaseSequence)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Sequence)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.allModsOneIsNterminus)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.NumVariableMods)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.missedCleavages.HasValue ? b.missedCleavages.Value.ToString(CultureInfo.InvariantCulture) : "unknown")).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.MonoisotopicMass)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => ScanPrecursorMass - b.MonoisotopicMass)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => ((ScanPrecursorMass - b.MonoisotopicMass) / b.MonoisotopicMass * 1e6))).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.Accession)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.FullName)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => string.Join(", ", b.Protein.GeneNames.Select(d => d.Item1 + ":" + d.Item2)))).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.IsContaminant ? "Y" : "N")).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.IsDecoy ? "Y" : "N")).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.PeptideDescription)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => ("[" + b.OneBasedStartResidueInProtein.ToString(CultureInfo.InvariantCulture) + " to " + b.OneBasedEndResidueInProtein.ToString(CultureInfo.InvariantCulture) + "]"))).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.PreviousAminoAcid.ToString())).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.NextAminoAcid.ToString())).Item1);

                // Unambiguous
                if (IsDecoy)
                    sb.Append("\t" + "D");
                else if (compactPeptides.Any(b => b.Value.Item2.Any(c => c.Protein.IsContaminant)))
                    sb.Append("\t" + "C");
                else
                    sb.Append("\t" + "T");
            }
            else
            {
                sb.Append('\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " ");
            }

            if (MatchedIonDictPositiveIsMatch != null)
            {
                sb.Append('\t' + string.Join(";", MatchedIonDictPositiveIsMatch.Select(b => b.Value.Count(c => c > 0))));

                sb.Append('\t' + "[");
                foreach (var kvp in MatchedIonDictPositiveIsMatch)
                    sb.Append("[" + string.Join(",", kvp.Value.Where(b => b > 0).Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
                sb.Append("]");
            }
            else
            {
                sb.Append('\t' + " " + '\t' + " ");
            }

            if (LocalizedScores != null)
            {
                sb.Append('\t' + "[" + string.Join(",", LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]");
                sb.Append('\t' + (LocalizedScores.Max() - Score).ToString("F3", CultureInfo.InvariantCulture));
            }
            else
            {
                sb.Append('\t' + " " + '\t' + " ");
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

        internal void Add(Psm psmParent)
        {
            foreach (var kvp in psmParent.compactPeptides)
                AddOrReplace(kvp.Key, psmParent.Score, kvp.Value.Item1);
        }

        #endregion Internal Methods

        #region Private Methods

        private Tuple<string, Dictionary<string, int>> Resolve(IEnumerable<Dictionary<int, ModificationWithMass>> enumerable)
        {
            Dictionary<string, int> ok = enumerable.First().Values.OrderBy(b => b.id).GroupBy(b => b.id).ToDictionary(b => b.Key, b => b.Count());
            bool notEqual = false;
            foreach (var ha in enumerable)
            {
                Dictionary<string, int> okTest = ha.Values.OrderBy(b => b.id).GroupBy(b => b.id).ToDictionary(b => b.Key, b => b.Count());
                if (!ok.SequenceEqual(okTest))
                {
                    notEqual = true;
                    break;
                }
            }
            if (notEqual)
            {
                var possibleReturn = string.Join(" or ", enumerable.Select(b => string.Join(" ", b.Values.Select(c => c.id).OrderBy(c => c))));
                if (possibleReturn.Length > 32000)
                    return new Tuple<string, Dictionary<string, int>>("too many", null);
                else
                    return new Tuple<string, Dictionary<string, int>>(possibleReturn, null);
            }
            else
            {
                return new Tuple<string, Dictionary<string, int>>(string.Join(" ", ok), ok);
            }
        }

        private Tuple<string, double?> Resolve(IEnumerable<double> enumerable)
        {
            double doubleTol = 1e-6;
            var list = enumerable.ToList();
            if (list.Max() - list.Min() < doubleTol)
            {
                return new Tuple<string, double?>(list.Average().ToString("F5", CultureInfo.InvariantCulture), list.Average());
            }
            else
            {
                var possibleReturn = string.Join(" or ", list.Select(b => b.ToString("F5", CultureInfo.InvariantCulture)));
                if (possibleReturn.Length > 32000)
                    return new Tuple<string, double?>("too many", null);
                else
                    return new Tuple<string, double?>(possibleReturn, null);
            }
        }

        private Tuple<string, int?> Resolve(IEnumerable<int> enumerable)
        {
            var list = enumerable.ToList();
            var first = list[0];
            if (list.All(b => first.Equals(b)))
            {
                return new Tuple<string, int?>(first.ToString(CultureInfo.InvariantCulture), first);
            }
            else
            {
                var possibleReturn = string.Join(" or ", list.Select(b => b.ToString(CultureInfo.InvariantCulture)));
                if (possibleReturn.Length > 32000)
                    return new Tuple<string, int?>("too many", null);
                else
                    return new Tuple<string, int?>(possibleReturn, null);
            }
        }

        private Tuple<string, string> Resolve(IEnumerable<string> enumerable)
        {
            var list = enumerable.ToList();
            var first = list.FirstOrDefault(b => b != null);
            // Only first if list is either all null or all equal to the first
            if (list.All(b => b == null) || list.All(b => first.Equals(b)))
            {
                return new Tuple<string, string>(first, first);
            }
            else
            {
                var possibleReturn = string.Join(" or ", list);
                if (possibleReturn.Length > 32000)
                    return new Tuple<string, string>("too many", null);
                else
                    return new Tuple<string, string>(possibleReturn, null);
            }
        }

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