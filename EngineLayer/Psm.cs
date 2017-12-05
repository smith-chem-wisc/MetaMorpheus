using Chemistry;
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

        private const double tolForDoubleResolution = 1e-6;
        private const double tolForScoreDifferentiation = 1e-9;

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
            this.ScanPrecursorMonoisotopicPeakMz = scan.PrecursorMonoisotopicPeakMz;
            this.ScanPrecursorMass = scan.PrecursorMass;
            AddOrReplace(peptide, score, notch, true);
            this.ExcelCompatible = true;
        }

        public Psm(CompactPeptideBase peptide, int notch, double score, int scanIndex, IScan scan, bool excelCompatible) : this(peptide, notch, score, scanIndex, scan)
        {
            this.ExcelCompatible = excelCompatible;
        }

        #endregion Public Constructors

        #region Public Properties

        public ChemicalFormula ModsChemicalFormula { get; private set; }
        public int ScanNumber { get; }
        public int? PrecursorScanNumber { get; }
        public double ScanRetentionTime { get; }
        public int ScanExperimentalPeaks { get; }
        public double TotalIonCurrent { get; }
        public int ScanPrecursorCharge { get; }
        public double ScanPrecursorMonoisotopicPeakMz { get; }
        public double ScanPrecursorMass { get; }
        public string FullFilePath { get; }
        public int ScanIndex { get; }
        public IEnumerable<KeyValuePair<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>> CompactPeptides { get { return compactPeptides.AsEnumerable(); } }
        public int NumDifferentCompactPeptides { get { return compactPeptides.Count; } }
        public FdrInfo FdrInfo { get; private set; }
        public double Score { get; private set; }

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
        public Dictionary<ProductType, double[]> MatchedIonDictOnlyMatches { get; internal set; }
        public string ProteinAccesion { get; private set; }
        public Dictionary<string, int> ModsIdentified { get; private set; }
        public Dictionary<ProductType, double[]> ProductMassErrorDa { get; internal set; }
        public Dictionary<ProductType, double[]> ProductMassErrorPpm { get; internal set; }

        #endregion Public Properties

        #region Private Properties

        private bool ExcelCompatible { get; set; }

        #endregion Private Properties

        #region Public Methods

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
            sb.Append('\t' + "Precursor Mass");
            sb.Append('\t' + "Score");
            sb.Append('\t' + "Notch");
            sb.Append('\t' + "Different Peak Matches");

            sb.Append('\t' + "Peptides Sharing Same Peaks");
            sb.Append('\t' + "Base Sequence");
            sb.Append('\t' + "Full Sequence");
            sb.Append('\t' + "Mods");
            sb.Append('\t' + "Mods Chemical Formula");
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
            sb.Append('\t' + "Matched Ion Mass Diff (Da)");
            sb.Append('\t' + "Matched Ion Mass Diff (Ppm)");

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

        public void AddOrReplace(CompactPeptideBase compactPeptide, double score, int notch, bool reportAllAmbiguity)
        {
            if (score - Score > tolForScoreDifferentiation) //if new score beat the old score, overwrite it
            {
                compactPeptides = new Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>
                {
                    { compactPeptide, new  Tuple<int, HashSet<PeptideWithSetModifications>>(notch,null)}
                };
                Score = score;
            }
            else if (score - Score > -tolForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
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
                compactPeptides[cpKey] = new Tuple<int, HashSet<PeptideWithSetModifications>>(compactPeptides[cpKey].Item1, matching[cpKey]);

            IsDecoy = compactPeptides.Any(b => b.Value.Item2.All(c => c.Protein.IsDecoy));

            FullSequence = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Sequence)).Item2;

            BaseSequence = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.BaseSequence)).Item2;

            PeptideLength = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Length)).Item2;

            OneBasedStartResidueInProtein = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.OneBasedStartResidueInProtein)).Item2;

            OneBasedEndResidueInProtein = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.OneBasedEndResidueInProtein)).Item2;

            ProteinLength = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.Length)).Item2;

            PeptideMonisotopicMass = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.MonoisotopicMass)).Item2;

            ProteinAccesion = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.Accession)).Item2;

            ModsIdentified = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.allModsOneIsNterminus)).Item2;

            ModsChemicalFormula = Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.allModsOneIsNterminus.Select(c => (c.Value as ModificationWithMassAndCf)))).Item2;

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
            sb.Append('\t' + (PrecursorScanNumber.HasValue ? PrecursorScanNumber.Value.ToString(CultureInfo.InvariantCulture) : "unknown"));
            sb.Append('\t' + ScanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + ScanPrecursorMonoisotopicPeakMz.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + ScanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture));
            sb.Append('\t' + Score.ToString("F3", CultureInfo.InvariantCulture));
            sb.Append("\t" + Resolve(compactPeptides.Select(b => b.Value.Item1)).Item1); // Notch
            sb.Append('\t' + NumDifferentCompactPeptides.ToString("F5", CultureInfo.InvariantCulture));

            if (compactPeptides.First().Value.Item2 != null)
            {
                sb.Append("\t" + TrimStringForExcel(string.Join(" or ", compactPeptides.Select(b => b.Value.Item2.Count.ToString(CultureInfo.InvariantCulture)))));

                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.BaseSequence)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Sequence)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.allModsOneIsNterminus)).Item1);
                sb.Append('\t' + Resolve(compactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.allModsOneIsNterminus.Select(c => (c.Value as ModificationWithMassAndCf)))).Item1);
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
                sb.Append('\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " " + '\t' + " ");
            }

            if (MatchedIonDictOnlyMatches != null)
            {
                //Count
                sb.Append('\t' + string.Join(";", MatchedIonDictOnlyMatches.Select(b => b.Value.Count(c => c > 0))));

                //Masses
                sb.Append('\t' + "[");
                foreach (var kvp in MatchedIonDictOnlyMatches)
                    sb.Append("[" + string.Join(",", kvp.Value.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
                sb.Append("]");

                //Mass error Da
                sb.Append('\t' + "[");
                foreach (var kvp in ProductMassErrorDa)
                    sb.Append("[" + string.Join(",", kvp.Value.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
                sb.Append("]");

                //Mass error ppm
                sb.Append('\t' + "[");
                foreach (var kvp in ProductMassErrorPpm)
                    sb.Append("[" + string.Join(",", kvp.Value.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
                sb.Append("]");
            }
            else
            {
                sb.Append('\t' + " " + '\t' + " " + '\t' + " " + '\t' + " ");
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
            FdrInfo = new FdrInfo
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

        internal void AddOrReplace(Psm psmParent, bool reportAllAmbiguity)
        {
            foreach (var kvp in psmParent.compactPeptides)
                AddOrReplace(kvp.Key, psmParent.Score, kvp.Value.Item1, reportAllAmbiguity);
        }

        #endregion Internal Methods

        #region Private Methods

        private static string TrimStringForExcel(string s)
        {
            return s.Length > 32000 ? "too many" : s;
        }

        private static (string, ChemicalFormula) Resolve(IEnumerable<IEnumerable<ModificationWithMassAndCf>> enumerable)
        {
            ChemicalFormula f = new ChemicalFormula();
            {
                var firstEnum = enumerable.First();
                foreach (var mod in firstEnum)
                {
                    if (mod == null)
                        return ("unknown", null);
                    f.Add(mod.chemicalFormula);
                }
            }
            bool equals = true;
            List<ChemicalFormula> formulas = new List<ChemicalFormula>();
            foreach (var anEnum in enumerable)
            {
                ChemicalFormula fhere = new ChemicalFormula();
                foreach (var mod in anEnum)
                {
                    if (mod == null)
                        return ("unknown", null);
                    fhere.Add(mod.chemicalFormula);
                }
                if (!f.Equals(fhere))
                    equals = false;
                formulas.Add(fhere);
            }
            if (!equals)
            {
                var possibleReturn = string.Join(" or ", formulas.Select(b => b.Formula));
                if (possibleReturn.Length > 32000)
                    return ("too many", null);
                else
                    return (possibleReturn, null);
            }
            else
            {
                return (f.Formula, f);
            }
        }

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
                return (ExcelCompatible && possibleReturn.Length > 32000) ? new Tuple<string, Dictionary<string, int>>("(too many)", null) : new Tuple<string, Dictionary<string, int>>(possibleReturn, null);
            }
            else
            {
                return new Tuple<string, Dictionary<string, int>>(string.Join(" ", enumerable.First().Values.Select(c => c.id).OrderBy(c => c)), ok);
            }
        }

        private Tuple<string, double?> Resolve(IEnumerable<double> enumerable)
        {
            var list = enumerable.ToList();
            if (list.Max() - list.Min() < tolForDoubleResolution)
            {
                return new Tuple<string, double?>(list.Average().ToString("F5", CultureInfo.InvariantCulture), list.Average());
            }
            else
            {
                var possibleReturn = string.Join(" or ", list.Select(b => b.ToString("F5", CultureInfo.InvariantCulture)));
                return (ExcelCompatible && possibleReturn.Length > 32000) ? new Tuple<string, double?>("(too many)", null) : new Tuple<string, double?>(possibleReturn, null);
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
                return (ExcelCompatible && possibleReturn.Length > 32000) ? new Tuple<string, int?>("(too many)", null) : new Tuple<string, int?>(possibleReturn, null);
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
                return (ExcelCompatible && possibleReturn.Length > 32000) ? new Tuple<string, string>("(too many)", null) : new Tuple<string, string>(possibleReturn, null);
            }
        }

        #endregion Private Methods
    }
}