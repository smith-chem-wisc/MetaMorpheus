using Chemistry;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class PeptideSpectralMatch
    {
        public const double ToleranceForScoreDifferentiation = 1e-9;

        private const double ToleranceForDoubleResolution = 1e-6;
        private Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>> _CompactPeptides = new Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>();

        public PeptideSpectralMatch(CompactPeptideBase peptide, int notch, double score, int scanIndex, IScan scan, DigestionParams digestionParams)
        {
            ScanIndex = scanIndex;
            FullFilePath = scan.FullFilePath;
            ScanNumber = scan.OneBasedScanNumber;
            PrecursorScanNumber = scan.OneBasedPrecursorScanNumber;
            ScanRetentionTime = scan.RetentionTime;
            ScanExperimentalPeaks = scan.NumPeaks;
            TotalIonCurrent = scan.TotalIonCurrent;
            ScanPrecursorCharge = scan.PrecursorCharge;
            ScanPrecursorMonoisotopicPeakMz = scan.PrecursorMonoisotopicPeakMz;
            ScanPrecursorMass = scan.PrecursorMass;
            AddOrReplace(peptide, score, notch, true);
            AllScores = new List<double>();
            DigestionParams = digestionParams;
            MatchedIonSeriesDict = new Dictionary<ProductType, int[]>();
            MatchedIonMassToChargeRatioDict = new Dictionary<ProductType, double[]>();
            MatchedIonIntensitiesDict = new Dictionary<ProductType, double[]>();
            ProductMassErrorDa = new Dictionary<ProductType, double[]>();
            ProductMassErrorPpm = new Dictionary<ProductType, double[]>();
            MatchedFragmentIons = new List<MatchedFragmentIon>();
        }

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
        public IEnumerable<KeyValuePair<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>> CompactPeptides { get { return _CompactPeptides.AsEnumerable(); } }
        public int NumDifferentCompactPeptides { get { return _CompactPeptides.Count; } }
        public FdrInfo FdrInfo { get; private set; }
        public double Score { get; private set; }
        public double DeltaScore { get; private set; }
        public double RunnerUpScore { get; set; }
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
        public Dictionary<ProductType, int[]> MatchedIonSeriesDict { get; internal set; }
        public Dictionary<ProductType, double[]> MatchedIonMassToChargeRatioDict { get; internal set; }
        public Dictionary<ProductType, double[]> MatchedIonIntensitiesDict { get; internal set; }
        public string ProteinAccesion { get; private set; }
        public string Organism { get; private set; }
        public Dictionary<string, int> ModsIdentified { get; private set; }
        public Dictionary<ProductType, double[]> ProductMassErrorDa { get; internal set; }
        public Dictionary<ProductType, double[]> ProductMassErrorPpm { get; internal set; }
        public readonly DigestionParams DigestionParams;
        public List<double> AllScores { get; set; }
        public List<MatchedFragmentIon> MatchedFragmentIons { get; private set; }

        public double[] Features
        {
            get
            {
                return new[] { Math.Round(Score), Score - Math.Round(Score) };
            }
        }

        public static string GetTabSeparatedHeader()
        {
            return String.Join("\t", DataDictionary(null, null).Keys);
        }

        public void SetMatchedFragments(List<MatchedFragmentIon> matchedFragmentIons)
        {
            MatchedFragmentIons = matchedFragmentIons;
        }

        public void AddOrReplace(CompactPeptideBase compactPeptide, double score, int notch, bool reportAllAmbiguity)
        {
            if (score - Score > ToleranceForScoreDifferentiation) //if new score beat the old score, overwrite it
            {
                _CompactPeptides = new Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>
                {
                    { compactPeptide, new  Tuple<int, HashSet<PeptideWithSetModifications>>(notch,null)}
                };
                if (Score - RunnerUpScore > ToleranceForScoreDifferentiation)
                {
                    RunnerUpScore = Score;
                }
                Score = score;
            }
            else if (score - Score > -ToleranceForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
            {
                _CompactPeptides[compactPeptide] = new Tuple<int, HashSet<PeptideWithSetModifications>>(notch, null);
            }
            else if (Score - RunnerUpScore > ToleranceForScoreDifferentiation)
            {
                RunnerUpScore = score;
            }
        }

        public void CompactCompactPeptides()
        {
            List<Tuple<CompactPeptideBase, int>> cps = new List<Tuple<CompactPeptideBase, int>>();
            foreach (KeyValuePair<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>> kvp in CompactPeptides)
            {
                //Change CPWM to reflect actual CP
                Tuple<CompactPeptideBase, int> tempTuple = new Tuple<CompactPeptideBase, int>(kvp.Key, kvp.Value.Item1);
                if (!cps.Contains(tempTuple))
                {
                    cps.Add(tempTuple);
                }
            }
            _CompactPeptides = new Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>();
            foreach (Tuple<CompactPeptideBase, int> cp in cps)
            {
                _CompactPeptides[cp.Item1] = new Tuple<int, HashSet<PeptideWithSetModifications>>(cp.Item2, null);
            }
        }

        public void MatchToProteinLinkedPeptides(Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching)
        {
            foreach (var cpKey in _CompactPeptides.Keys.ToList())
            {
                _CompactPeptides[cpKey] = new Tuple<int, HashSet<PeptideWithSetModifications>>(_CompactPeptides[cpKey].Item1, matching[cpKey]);
            }
            var pepsWithMods = CompactPeptides.SelectMany(b => b.Value.Item2);
            IsDecoy = CompactPeptides.Any(b => b.Value.Item2.All(c => c.Protein.IsDecoy));
            FullSequence = Resolve(pepsWithMods.Select(b => b.Sequence)).Item2;
            BaseSequence = Resolve(pepsWithMods.Select(b => b.BaseSequence)).Item2;
            PeptideLength = Resolve(pepsWithMods.Select(b => b.Length)).Item2;
            OneBasedStartResidueInProtein = Resolve(pepsWithMods.Select(b => b.OneBasedStartResidueInProtein)).Item2;
            OneBasedEndResidueInProtein = Resolve(pepsWithMods.Select(b => b.OneBasedEndResidueInProtein)).Item2;
            ProteinLength = Resolve(pepsWithMods.Select(b => b.Protein.Length)).Item2;
            PeptideMonisotopicMass = Resolve(pepsWithMods.Select(b => b.MonoisotopicMass)).Item2;
            ProteinAccesion = Resolve(pepsWithMods.Select(b => b.Protein.Accession)).Item2;
            Organism = Resolve(pepsWithMods.Select(b => b.Protein.Organism)).Item2;
            ModsIdentified = Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus)).Item2;
            ModsChemicalFormula = Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus.Select(c => (c.Value as ModificationWithMassAndCf)))).Item2;
            Notch = Resolve(CompactPeptides.Select(b => b.Value.Item1)).Item2;
        }

        public override string ToString()
        {
            return ToString(new Dictionary<string, int>());
        }

        public string ToString(IReadOnlyDictionary<string, int> ModstoWritePruned)
        {
            return String.Join("\t", DataDictionary(this, ModstoWritePruned).Values);
        }

        public static Dictionary<string, string> DataDictionary(PeptideSpectralMatch peptide, IReadOnlyDictionary<string, int> ModsToWritePruned)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            AddBasicMatchData(s, peptide);
            AddPeptideSequenceData(s, peptide, ModsToWritePruned);
            AddMatchedIonsData(s, peptide);
            AddMatchScoreData(s, peptide);
            return s;
        }

        public void CalculateDeltaScore(double scoreCutoff)
        {
            DeltaScore = Score - Math.Max(RunnerUpScore, scoreCutoff);
        }

        public void SetFdrValues(int cumulativeTarget, int cumulativeDecoy, double tempQValue, int cumulativeTargetNotch, int cumulativeDecoyNotch, double tempQValueNotch, double maximumLikelihood, double eValue, double eScore, bool calculateEValue)
        {
            FdrInfo = new FdrInfo
            {
                CumulativeTarget = cumulativeTarget,
                CumulativeDecoy = cumulativeDecoy,
                QValue = tempQValue,
                CumulativeTargetNotch = cumulativeTargetNotch,
                CumulativeDecoyNotch = cumulativeDecoyNotch,
                QValueNotch = tempQValueNotch,
                MaximumLikelihood = maximumLikelihood,
                EScore = eScore,
                EValue = eValue,
                CalculateEValue = calculateEValue
            };
        }

        private static (string, ChemicalFormula) Resolve(IEnumerable<IEnumerable<ModificationWithMassAndCf>> enumerable)
        {
            ChemicalFormula f = new ChemicalFormula();
            {
                var firstEnum = enumerable.First();
                foreach (var mod in firstEnum)
                {
                    if (mod == null)
                    {
                        return ("unknown", null);
                    }
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
                    {
                        return ("unknown", null);
                    }
                    fhere.Add(mod.chemicalFormula);
                }
                if (!f.Equals(fhere))
                {
                    equals = false;
                }
                formulas.Add(fhere);
            }
            if (!equals)
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", formulas.Select(b => b.Formula)));
                return (returnString, null);
            }
            else
            {
                return (f.Formula, f);
            }
        }

        private static Tuple<string, Dictionary<string, int>> Resolve(IEnumerable<Dictionary<int, ModificationWithMass>> enumerable)
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
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", enumerable.Select(b => string.Join(" ", b.Values.Select(c => c.id).OrderBy(c => c)))));
                return new Tuple<string, Dictionary<string, int>>(returnString, null);
            }
            else
            {
                return new Tuple<string, Dictionary<string, int>>(string.Join(" ", enumerable.First().Values.Select(c => c.id).OrderBy(c => c)), ok);
            }
        }

        private static Tuple<string, double?> ResolveF2(IEnumerable<double> enumerable)
        {
            var list = enumerable.ToList();
            if (list.Max() - list.Min() < ToleranceForDoubleResolution)
            {
                return new Tuple<string, double?>(list.Average().ToString("F2", CultureInfo.InvariantCulture), list.Average());
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F2", CultureInfo.InvariantCulture))));
                return new Tuple<string, double?>(returnString, null);
            }
        }

        private static Tuple<string, double?> Resolve(IEnumerable<double> enumerable)
        {
            var list = enumerable.ToList();
            if (list.Max() - list.Min() < ToleranceForDoubleResolution)
            {
                return new Tuple<string, double?>(list.Average().ToString("F5", CultureInfo.InvariantCulture), list.Average());
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))));
                return new Tuple<string, double?>(returnString, null);
            }
        }

        private static Tuple<string, int?> Resolve(IEnumerable<int> enumerable)
        {
            var list = enumerable.ToList();
            var first = list[0];
            if (list.All(b => first.Equals(b)))
            {
                return new Tuple<string, int?>(first.ToString(CultureInfo.InvariantCulture), first);
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString(CultureInfo.InvariantCulture))));
                return new Tuple<string, int?>(returnString, null);
            }
        }

        private static Tuple<string, string> Resolve(IEnumerable<string> enumerable)
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
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list));
                return new Tuple<string, string>(returnString, null);
            }
        }

        private static void AddBasicMatchData(Dictionary<string, string> s, PeptideSpectralMatch peptide)
        {
            s["File Name"] = peptide == null ? " " : Path.GetFileNameWithoutExtension(peptide.FullFilePath);
            s["Scan Number"] = peptide == null ? " " : peptide.ScanNumber.ToString(CultureInfo.InvariantCulture);
            s["Scan Retention Time"] = peptide == null ? " " : peptide.ScanRetentionTime.ToString("F5", CultureInfo.InvariantCulture);
            s["Num Experimental Peaks"] = peptide == null ? " " : peptide.ScanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture);
            s["Total Ion Current"] = peptide == null ? " " : peptide.TotalIonCurrent.ToString("F5", CultureInfo.InvariantCulture);
            s["Precursor Scan Number"] = peptide == null ? " " : peptide.PrecursorScanNumber.HasValue ? peptide.PrecursorScanNumber.Value.ToString(CultureInfo.InvariantCulture) : "unknown";
            s["Precursor Charge"] = peptide == null ? " " : peptide.ScanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture);
            s["Precursor MZ"] = peptide == null ? " " : peptide.ScanPrecursorMonoisotopicPeakMz.ToString("F5", CultureInfo.InvariantCulture);
            s["Precursor Mass"] = peptide == null ? " " : peptide.ScanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture);
            s["Score"] = peptide == null ? " " : peptide.Score.ToString("F3", CultureInfo.InvariantCulture);
            s["Delta Score"] = peptide == null ? " " : peptide.DeltaScore.ToString("F3", CultureInfo.InvariantCulture);
            s["Notch"] = peptide == null ? " " : Resolve(peptide.CompactPeptides.Select(b => b.Value.Item1)).Item1; // Notch
            s["Different Peak Matches"] = peptide == null ? " " : peptide.NumDifferentCompactPeptides.ToString("F5", CultureInfo.InvariantCulture);
        }

        private static void AddPeptideSequenceData(Dictionary<string, string> s, PeptideSpectralMatch peptide, IReadOnlyDictionary<string, int> ModsToWritePruned)
        {
            bool pepWithModsIsNull = peptide == null || peptide.CompactPeptides.First().Value.Item2 == null;

            var pepsWithMods = pepWithModsIsNull ? null : peptide.CompactPeptides.SelectMany(b => b.Value.Item2)
                .OrderBy(p => p.Sequence)
                .ThenBy(p => p.Protein.Accession)
                .ThenBy(p => p.OneBasedStartResidueInProtein).ToList();

            s["Peptides Sharing Same Peaks"] = pepWithModsIsNull ? " " :
                GlobalVariables.CheckLengthOfOutput(string.Join("|", peptide.CompactPeptides.Select(b => b.Value.Item2.Count.ToString(CultureInfo.InvariantCulture))));
            s["Base Sequence"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.BaseSequence)).Item1;
            s["Full Sequence"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Sequence)).Item1;
            s["Essential Sequence"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.EssentialSequence(ModsToWritePruned))).Item1;
            s["Mods"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus)).Item1;
            s["Mods Chemical Formulas"] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(b => string.Join("|", b.AllModsOneIsNterminus.OrderBy(c => c.Key)
                    .Where(c => c.Value is ModificationWithMassAndCf).Select(c => (c.Value as ModificationWithMassAndCf).chemicalFormula.Formula)))).Item1;
            s["Mods Combined Chemical Formula"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus.Select(c => (c.Value as ModificationWithMassAndCf)))).Item1;
            s["Num Variable Mods"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NumVariableMods)).Item1;
            s["Missed Cleavages"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MissedCleavages.ToString(CultureInfo.InvariantCulture))).Item1;
            s["Peptide Monoisotopic Mass"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MonoisotopicMass)).Item1;
            s["Mass Diff (Da)"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => peptide.ScanPrecursorMass - b.MonoisotopicMass)).Item1;
            s["Mass Diff (ppm)"] = pepWithModsIsNull ? " " : ResolveF2(pepsWithMods.Select(b => ((peptide.ScanPrecursorMass - b.MonoisotopicMass) / b.MonoisotopicMass * 1e6))).Item1;
            s["Protein Accession"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.Accession)).Item1;
            s["Protein Name"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.FullName)).Item1;
            s["Gene Name"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => string.Join(", ", b.Protein.GeneNames.Select(d => d.Item1 + ":" + d.Item2)))).Item1;
            s["Sequence Variations"] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(b => string.Join(", ", b.Protein.SequenceVariations
                    .Where(d => peptide.OneBasedStartResidueInProtein <= d.OneBasedBeginPosition && d.OneBasedBeginPosition <= peptide.OneBasedEndResidueInProtein)
                    .Select(d => d.OriginalSequence + d.OneBasedBeginPosition.ToString() + d.VariantSequence)))).Item1;
            s["Organism Name"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.Organism)).Item1;
            s["Contaminant"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.IsContaminant ? "Y" : "N")).Item1;
            s["Decoy"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.IsDecoy ? "Y" : "N")).Item1;
            s["Peptide Description"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.PeptideDescription)).Item1;
            s["Start and End Residues In Protein"] =
                pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => ("[" + b.OneBasedStartResidueInProtein.ToString(CultureInfo.InvariantCulture) + " to " +
                    b.OneBasedEndResidueInProtein.ToString(CultureInfo.InvariantCulture) + "]"))).Item1;
            s["Previous Amino Acid"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.PreviousAminoAcid.ToString())).Item1;
            s["Next Amino Acid"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NextAminoAcid.ToString())).Item1;

            string allScores = " ";
            string theoreticalsSearched = " ";
            if (!pepWithModsIsNull && peptide.FdrInfo != null && peptide.FdrInfo.CalculateEValue)
            {
                allScores = string.Join(";", peptide.AllScores.Select(p => p.ToString("F2", CultureInfo.InvariantCulture)));
                theoreticalsSearched = peptide.AllScores.Count.ToString();
            }

            s["All Scores"] = allScores;
            s["Theoreticals Searched"] = theoreticalsSearched;
            s["Decoy/Contaminant/Target"] = pepWithModsIsNull ? " " : peptide.IsDecoy ? "D" : pepsWithMods.Any(c => c.Protein.IsContaminant) ? "C" : "T";
        }

        private static void AddMatchedIonsData(Dictionary<string, string> s, PeptideSpectralMatch peptide)
        {
            //sb for writing, format for double.ToString, header for input dictionary key
            (StringBuilder sb, string format, string header)[] matchedIonInfo = new(StringBuilder, string, string)[]
            {
                (new StringBuilder(), "F5", "Matched Ion Mass-To-Charge Ratios"),
                (new StringBuilder(), "F5", "Matched Ion Mass Diff (Da)"),
                (new StringBuilder(), "F2", "Matched Ion Mass Diff (Ppm)"),
                (new StringBuilder(), "F0", "Matched Ion Intensities")
            };

            const string blankEntry = " "; //the space prevents sorting issues in excel
            string ionCounts = blankEntry;
            string ionSeries = blankEntry;

            if (peptide != null && peptide.MatchedIonMassToChargeRatioDict.Any())
            {
                StringBuilder seriesStringBuilder = new StringBuilder(); //this one is used for series, all others are done in the dictionary

                //Dictionary allows for easy iteration over the series, where the key is the interesting data and the value is a tuple of the string to build and the string used for double.ToString rounding
                Dictionary<Dictionary<ProductType, double[]>, (StringBuilder sb, string format, string header)> peptideInfoToStringBuilderDict
                    = new Dictionary<Dictionary<ProductType, double[]>, (StringBuilder sb, string format, string header)>
                {
                    { peptide.MatchedIonMassToChargeRatioDict,  matchedIonInfo[0]},
                    { peptide.ProductMassErrorDa, matchedIonInfo[1] },
                    { peptide.ProductMassErrorPpm, matchedIonInfo[2] },
                    { peptide.MatchedIonIntensitiesDict, matchedIonInfo[3] }
                };

                const string delimiter = ", "; //using ", " instead of "," improves human readability

                foreach (ProductType productType in peptide.MatchedIonSeriesDict.Keys)
                {
                    string ionType = productType.ToString()[0].ToString().ToLower(); //gets the first char of the type (ie b, y, c, z)

                    //Ion series
                    string[] seriesToWrite = peptide.MatchedIonSeriesDict[productType].Select(x => ionType + x.ToString() + "+1").ToArray(); //assumes charge of +1
                    seriesStringBuilder.Append("[" + string.Join(delimiter, seriesToWrite) + "];");

                    //All other data present in dictionary
                    //add bracket to the beginning of each stringbuilder in the dictionary for this product type
                    foreach (var kvp in peptideInfoToStringBuilderDict)
                    {
                        kvp.Value.sb.Append("[");
                    }

                    //add all the data to the stringbuilder
                    for (int i = 0; i < seriesToWrite.Length - 1; i++) //-1 so no delimiter for the last entry
                    {
                        string currentSeriesToWrite = seriesToWrite[i];
                        foreach (var kvp in peptideInfoToStringBuilderDict)
                        {
                            kvp.Value.sb.Append(currentSeriesToWrite + ":" + kvp.Key[productType][i].ToString(kvp.Value.format, CultureInfo.InvariantCulture) + delimiter);
                        }
                    }
                    //add the last data from each array without the delimiter
                    if (seriesToWrite.Length != 0)//check it's not empty
                    {
                        string currentSeriesToWrite = seriesToWrite.Last();
                        foreach (var kvp in peptideInfoToStringBuilderDict)
                        {
                            kvp.Value.sb.Append(currentSeriesToWrite + ":" + kvp.Key[productType].Last().ToString(kvp.Value.format, CultureInfo.InvariantCulture)); //no delimiter here!
                        }
                    }

                    foreach (var kvp in peptideInfoToStringBuilderDict) //add bracket to the end of each stringbuilder in the dictionary
                    {
                        kvp.Value.sb.Append("];");
                    }
                }
                ionCounts = string.Join(";", peptide.MatchedIonMassToChargeRatioDict.Select(b => b.Value.Count(c => c > 0)));
                ionSeries = "[" + GlobalVariables.CheckLengthOfOutput(seriesStringBuilder.ToString()) + "]";
            }
            else
            {
                foreach (var info in matchedIonInfo)
                {
                    info.sb.Append(blankEntry);
                }
            }

            //write into input dictionary
            s["Matched Ion Counts"] = ionCounts;
            s["Matched Ion Series"] = ionSeries;
            foreach (var info in matchedIonInfo)
            {
                s[info.header] = info.sb.ToString();
            }
        }

        private static void AddMatchScoreData(Dictionary<string, string> s, PeptideSpectralMatch peptide)
        {
            string localizedScores = " ";
            string improvementPossible = " ";
            if (peptide != null && peptide.LocalizedScores != null)
            {
                localizedScores = GlobalVariables.CheckLengthOfOutput(("[" + string.Join(",", peptide.LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]"));
                improvementPossible = (peptide.LocalizedScores.Max() - peptide.Score).ToString("F3", CultureInfo.InvariantCulture);
            }
            s["Localized Scores"] = localizedScores;
            s["Improvement Possible"] = improvementPossible;

            string cumulativeTarget = " ";
            string cumulativeDecoy = " ";
            string qValue = " ";
            string cumulativeTargetNotch = " ";
            string cumulativeDecoyNotch = " ";
            string qValueNotch = " ";
            string eValue = " ";
            string eScore = " ";
            if (peptide != null && peptide.FdrInfo != null)
            {
                cumulativeTarget = peptide.FdrInfo.CumulativeTarget.ToString(CultureInfo.InvariantCulture);
                cumulativeDecoy = peptide.FdrInfo.CumulativeDecoy.ToString(CultureInfo.InvariantCulture);
                qValue = peptide.FdrInfo.QValue.ToString("F6", CultureInfo.InvariantCulture);
                cumulativeTargetNotch = peptide.FdrInfo.CumulativeTargetNotch.ToString(CultureInfo.InvariantCulture);
                cumulativeDecoyNotch = peptide.FdrInfo.CumulativeDecoyNotch.ToString(CultureInfo.InvariantCulture);
                qValueNotch = peptide.FdrInfo.QValueNotch.ToString("F6", CultureInfo.InvariantCulture);
                if (peptide.FdrInfo.CalculateEValue)
                {
                    eValue = peptide.FdrInfo.EValue.ToString("F6", CultureInfo.InvariantCulture);
                    eScore = peptide.FdrInfo.EScore.ToString("F6", CultureInfo.InvariantCulture);
                }
            }
            s["Cumulative Target"] = cumulativeTarget;
            s["Cumulative Decoy"] = cumulativeDecoy;
            s["QValue"] = qValue;
            s["Cumulative Target Notch"] = cumulativeTargetNotch;
            s["Cumulative Decoy Notch"] = cumulativeDecoyNotch;
            s["QValue Notch"] = qValueNotch;
            s["eValue"] = eValue;
            s["eScore"] = eScore;
        }
    }
}