using Chemistry;
using EngineLayer.FdrAnalysis;
using Proteomics;
using Proteomics.Fragmentation;
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
        private const double ToleranceForDoubleResolution = 1e-6;

        private List<(int Notch, PeptideWithSetModifications Pwsm)> _bestMatchingPeptides;

        public const double ToleranceForScoreDifferentiation = 1e-9;

        public PeptideSpectralMatch(PeptideWithSetModifications peptide, int notch, double score, int scanIndex, IScan scan, DigestionParams digestionParams, List<MatchedFragmentIon> matchedFragmentIons)
        {
            _bestMatchingPeptides = new List<(int, PeptideWithSetModifications)>();
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
            AllScores = new List<double>();
            DigestionParams = digestionParams;
            PeptidesToMatchingFragments = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();

            AddOrReplace(peptide, score, notch, true, matchedFragmentIons);
        }

        // these fields will be null if they are ambiguous
        public ChemicalFormula ModsChemicalFormula { get; private set; }
        public string FullSequence { get; private set; }
        public int? Notch { get; private set; }
        public string BaseSequence { get; private set; }
        public int? PeptideLength { get; private set; }
        public int? OneBasedStartResidueInProtein { get; private set; }
        public int? OneBasedEndResidueInProtein { get; private set; }
        public double? PeptideMonisotopicMass { get; private set; }
        public int? ProteinLength { get; private set; }
        public string ProteinAccession { get; private set; }
        public string Organism { get; private set; }
        public List<MatchedFragmentIon> MatchedFragmentIons { get; protected set; }

        // these should never be null under normal circumstances
        public Dictionary<string, int> ModsIdentified { get; private set; }
        public List<double> LocalizedScores { get; internal set; }
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
        public int NumDifferentMatchingPeptides { get { return _bestMatchingPeptides.Count; } }
        public FdrInfo FdrInfo { get; private set; }
        public double Score { get; private set; }
        public double DeltaScore { get; private set; }
        public double RunnerUpScore { get; set; }
        public bool IsDecoy { get; private set; }
        public bool IsContaminant { get; private set; }
        public readonly DigestionParams DigestionParams;
        public List<double> AllScores { get; set; }
        public Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>> PeptidesToMatchingFragments { get; private set; }

        public IEnumerable<(int Notch, PeptideWithSetModifications Peptide)> BestMatchingPeptides
        {
            get
            {
                return _bestMatchingPeptides.OrderBy(p => p.Item2.FullSequence)
                    .ThenBy(p => p.Item2.Protein.Accession)
                    .ThenBy(p => p.Item2.OneBasedStartResidueInProtein);
            }
        }

        /// <summary>
        /// Used for Percolator output
        /// </summary>
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

        public void AddOrReplace(PeptideWithSetModifications pwsm, double newScore, int notch, bool reportAllAmbiguity, List<MatchedFragmentIon> matchedFragmentIons)
        {
            if (newScore - Score > ToleranceForScoreDifferentiation) //if new score beat the old score, overwrite it
            {
                _bestMatchingPeptides.Clear();
                _bestMatchingPeptides.Add((notch, pwsm));

                if (Score - RunnerUpScore > ToleranceForScoreDifferentiation)
                {
                    RunnerUpScore = Score;
                }

                Score = newScore;

                PeptidesToMatchingFragments.Clear();
                PeptidesToMatchingFragments.Add(pwsm, matchedFragmentIons);
            }
            else if (newScore - Score > -ToleranceForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
            {
                _bestMatchingPeptides.Add((notch, pwsm));

                if (!PeptidesToMatchingFragments.ContainsKey(pwsm))
                {
                    PeptidesToMatchingFragments.Add(pwsm, matchedFragmentIons);
                }
            }
            else if (Score - RunnerUpScore > ToleranceForScoreDifferentiation)
            {
                RunnerUpScore = newScore;
            }
        }

        public override string ToString()
        {
            return ToString(new Dictionary<string, int>());
        }

        public string ToString(IReadOnlyDictionary<string, int> ModstoWritePruned)
        {
            return String.Join("\t", DataDictionary(this, ModstoWritePruned).Values);
        }

        public static Dictionary<string, string> DataDictionary(PeptideSpectralMatch psm, IReadOnlyDictionary<string, int> ModsToWritePruned)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            AddBasicMatchData(s, psm);
            AddPeptideSequenceData(s, psm, ModsToWritePruned);
            AddMatchedIonsData(s, psm);
            AddMatchScoreData(s, psm);
            return s;
        }

        public void CalculateDeltaScore(double scoreCutoff)
        {
            DeltaScore = Score - Math.Max(RunnerUpScore, scoreCutoff);
        }

        public void SetFdrValues(double cumulativeTarget, double cumulativeDecoy, double qValue, double cumulativeTargetNotch, double cumulativeDecoyNotch, double qValueNotch, double maximumLikelihood, double eValue, double eScore, bool calculateEValue)
        {
            FdrInfo = new FdrInfo
            {
                CumulativeTarget = cumulativeTarget,
                CumulativeDecoy = cumulativeDecoy,
                QValue = qValue,
                CumulativeTargetNotch = cumulativeTargetNotch,
                CumulativeDecoyNotch = cumulativeDecoyNotch,
                QValueNotch = qValueNotch,
                MaximumLikelihood = maximumLikelihood,
                EScore = eScore,
                EValue = eValue,
                CalculateEValue = calculateEValue
            };
        }

        /// <summary>
        /// This method saves properties of this PSM for internal use. It is NOT used for any output.
        /// These resolved fields are (usually) null if there is more than one option.
        /// e.g., if this PSM can be explained by more than one base sequence, the BaseSequence property will be null
        /// </summary>
        public void ResolveAllAmbiguities()
        {
            IsDecoy = _bestMatchingPeptides.Any(p => p.Pwsm.Protein.IsDecoy);
            IsContaminant = _bestMatchingPeptides.Any(p => p.Pwsm.Protein.IsContaminant);

            FullSequence = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.FullSequence)).ResolvedValue;
            BaseSequence = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.BaseSequence)).ResolvedValue;
            PeptideLength = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.Length)).ResolvedValue;
            OneBasedStartResidueInProtein = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.OneBasedStartResidueInProtein)).ResolvedValue;
            OneBasedEndResidueInProtein = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.OneBasedEndResidueInProtein)).ResolvedValue;
            ProteinLength = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.Protein.Length)).ResolvedValue;
            PeptideMonisotopicMass = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.MonoisotopicMass)).ResolvedValue;
            ProteinAccession = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.Protein.Accession)).ResolvedValue;
            Organism = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.Protein.Organism)).ResolvedValue;
            ModsIdentified = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.AllModsOneIsNterminus)).ResolvedValue;
            ModsChemicalFormula = Resolve(_bestMatchingPeptides.Select(b => b.Pwsm.AllModsOneIsNterminus.Select(c => (c.Value)))).ResolvedValue;
            Notch = Resolve(_bestMatchingPeptides.Select(b => b.Notch)).ResolvedValue;

            // if the PSM matches a target and a decoy and they are the SAME SEQUENCE, remove the decoy
            if (IsDecoy)
            {
                bool removedPeptides = false;
                var hits = _bestMatchingPeptides.GroupBy(p => p.Pwsm.FullSequence);

                foreach (var hit in hits)
                {
                    if (hit.Any(p => p.Pwsm.Protein.IsDecoy) && hit.Any(p => !p.Pwsm.Protein.IsDecoy))
                    {
                        // at least one peptide with this sequence is a target and at least one is a decoy
                        // remove the decoys with this sequence
                        _bestMatchingPeptides.RemoveAll(p => p.Pwsm.FullSequence == hit.Key && p.Pwsm.Protein.IsDecoy);
                        removedPeptides = true;
                    }
                }

                if (removedPeptides)
                {
                    ResolveAllAmbiguities();
                }
            }

            // TODO: technically, different peptide options for this PSM can have different matched ions
            // we can write a Resolve method for this if we want...
            MatchedFragmentIons = PeptidesToMatchingFragments.First().Value;
        }

        /// <summary>
        /// This method is used by protein parsimony to remove PeptideWithSetModifications objects that have non-parsimonious protein associations
        /// </summary>
        public void TrimProteinMatches(HashSet<Protein> parsimoniousProteins)
        {
            if (IsDecoy)
            {
                if (_bestMatchingPeptides.Any(p => parsimoniousProteins.Contains(p.Pwsm.Protein) && p.Pwsm.Protein.IsDecoy))
                {
                    _bestMatchingPeptides.RemoveAll(p => !parsimoniousProteins.Contains(p.Item2.Protein));
                }
                // else do nothing
            }
            else
            {
                _bestMatchingPeptides.RemoveAll(p => !parsimoniousProteins.Contains(p.Item2.Protein));
            }

            ResolveAllAmbiguities();
        }

        /// <summary>
        /// This method is used by protein parsimony to add PeptideWithSetModifications objects for modification-agnostic parsimony
        /// </summary>
        public void AddProteinMatch((int, PeptideWithSetModifications) peptideWithNotch)
        {
            _bestMatchingPeptides.Add(peptideWithNotch);
            ResolveAllAmbiguities();
        }

        /// <summary>
        /// Resolve Methods()
        /// if all 'values' are the same this returns the one value otherwise you get a separated list of all values in their original order.
        /// for example:
        /// Notches 1,1,1,1 returns as 1
        /// Notches 1,0,1,0 returns as 1|0|1|0
        /// </summary>
        private static (string ResolvedString, ChemicalFormula ResolvedValue) Resolve(IEnumerable<IEnumerable<Modification>> enumerable)
        {
            ChemicalFormula f = new ChemicalFormula();
            {
                var firstEnum = enumerable.First();
                foreach (var mod in firstEnum)
                {
                    if (mod == null || mod.ChemicalFormula == null)
                    {
                        return ("unknown", null);
                    }
                    f.Add(mod.ChemicalFormula);
                }
            }
            bool equals = true;
            List<ChemicalFormula> formulas = new List<ChemicalFormula>();
            foreach (var anEnum in enumerable)
            {
                ChemicalFormula fhere = new ChemicalFormula();
                foreach (var mod in anEnum)
                {
                    if (mod == null || mod.ChemicalFormula == null)
                    {
                        return ("unknown", null);
                    }
                    fhere.Add(mod.ChemicalFormula);
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

        private static (string ResolvedString, Dictionary<string, int> ResolvedValue) Resolve(IEnumerable<Dictionary<int, Modification>> enumerable)
        {
            Dictionary<string, int> ok = enumerable.First().Values.OrderBy(b => b.IdWithMotif).GroupBy(b => b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());
            bool notEqual = false;
            foreach (var ha in enumerable)
            {
                Dictionary<string, int> okTest = ha.Values.OrderBy(b => b.IdWithMotif).GroupBy(b => b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());
                if (!ok.SequenceEqual(okTest))
                {
                    notEqual = true;
                    break;
                }
            }
            if (notEqual)
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", enumerable.Select(b => string.Join(" ", b.Values.Select(c => c.IdWithMotif).OrderBy(c => c)))));
                return (returnString, null);
            }
            else
            {
                return (string.Join(" ", enumerable.First().Values.Select(c => c.IdWithMotif).OrderBy(c => c)), ok);
            }
        }

        private static (string ResolvedString, double? ResolvedValue) ResolveF2(IEnumerable<double> enumerable)
        {
            var list = enumerable.ToList();
            if (list.Max() - list.Min() < ToleranceForDoubleResolution)
            {
                return (list.Average().ToString("F2", CultureInfo.InvariantCulture), list.Average());
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F2", CultureInfo.InvariantCulture))));
                return (returnString, null);
            }
        }

        private static (string ResolvedString, double? ResolvedValue) Resolve(IEnumerable<double> enumerable)
        {
            var list = enumerable.ToList();
            if (list.Max() - list.Min() < ToleranceForDoubleResolution)
            {
                return (list.Average().ToString("F5", CultureInfo.InvariantCulture), list.Average());
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))));
                return (returnString, null);
            }
        }

        private static (string ResolvedString, int? ResolvedValue) Resolve(IEnumerable<int> enumerable)
        {
            var list = enumerable.ToList();
            var first = list[0];
            if (list.All(b => first.Equals(b)))
            {
                return (first.ToString(CultureInfo.InvariantCulture), first);
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString(CultureInfo.InvariantCulture))));
                return (returnString, null);
            }
        }

        private static (string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable)
        {
            var list = enumerable.ToList();
            string first = list.FirstOrDefault(b => b != null);
            // Only first if list is either all null or all equal to the first
            if (list.All(b => b == null) || list.All(b => first.Equals(b)))
            {
                return (first, first);
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list));
                return (returnString, null);
            }
        }

        private static void AddBasicMatchData(Dictionary<string, string> s, PeptideSpectralMatch psm)
        {
            s["File Name"] = psm == null ? " " : Path.GetFileNameWithoutExtension(psm.FullFilePath);
            s["Scan Number"] = psm == null ? " " : psm.ScanNumber.ToString(CultureInfo.InvariantCulture);
            s["Scan Retention Time"] = psm == null ? " " : psm.ScanRetentionTime.ToString("F5", CultureInfo.InvariantCulture);
            s["Num Experimental Peaks"] = psm == null ? " " : psm.ScanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture);
            s["Total Ion Current"] = psm == null ? " " : psm.TotalIonCurrent.ToString("F5", CultureInfo.InvariantCulture);
            s["Precursor Scan Number"] = psm == null ? " " : psm.PrecursorScanNumber.HasValue ? psm.PrecursorScanNumber.Value.ToString(CultureInfo.InvariantCulture) : "unknown";
            s["Precursor Charge"] = psm == null ? " " : psm.ScanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture);
            s["Precursor MZ"] = psm == null ? " " : psm.ScanPrecursorMonoisotopicPeakMz.ToString("F5", CultureInfo.InvariantCulture);
            s["Precursor Mass"] = psm == null ? " " : psm.ScanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture);
            s["Score"] = psm == null ? " " : psm.Score.ToString("F3", CultureInfo.InvariantCulture);
            s["Delta Score"] = psm == null ? " " : psm.DeltaScore.ToString("F3", CultureInfo.InvariantCulture);
            s["Notch"] = psm == null ? " " : Resolve(psm.BestMatchingPeptides.Select(p => p.Notch)).ResolvedString;
            s["Different Peak Matches"] = psm == null ? " " : psm.NumDifferentMatchingPeptides.ToString("F5", CultureInfo.InvariantCulture);
        }

        private static void AddPeptideSequenceData(Dictionary<string, string> s, PeptideSpectralMatch psm, IReadOnlyDictionary<string, int> ModsToWritePruned)
        {
            bool pepWithModsIsNull = psm == null || psm.BestMatchingPeptides == null || !psm.BestMatchingPeptides.Any();

            var pepsWithMods = pepWithModsIsNull ? null : psm.BestMatchingPeptides.Select(p => p.Peptide).ToList();

            s["Base Sequence"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.BaseSequence)).ResolvedString;
            s["Full Sequence"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.FullSequence)).ResolvedString;
            s["Essential Sequence"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.EssentialSequence(ModsToWritePruned))).ResolvedString;
            s["Mods"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus)).ResolvedString;
            s["Mods Chemical Formulas"] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(p => p.AllModsOneIsNterminus.Select(v => v.Value))).ResolvedString;
            s["Mods Combined Chemical Formula"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus.Select(c => (c.Value as Modification)))).ResolvedString;
            s["Num Variable Mods"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NumVariableMods)).Item1;
            s["Missed Cleavages"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MissedCleavages.ToString(CultureInfo.InvariantCulture))).ResolvedString;
            s["Peptide Monoisotopic Mass"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MonoisotopicMass)).ResolvedString;
            s["Mass Diff (Da)"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => psm.ScanPrecursorMass - b.MonoisotopicMass)).ResolvedString;
            s["Mass Diff (ppm)"] = pepWithModsIsNull ? " " : ResolveF2(pepsWithMods.Select(b => ((psm.ScanPrecursorMass - b.MonoisotopicMass) / b.MonoisotopicMass * 1e6))).ResolvedString;
            s["Protein Accession"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.Accession)).ResolvedString;
            s["Protein Name"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.FullName)).ResolvedString;
            s["Gene Name"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => string.Join(", ", b.Protein.GeneNames.Select(d => d.Item1 + ":" + d.Item2)))).ResolvedString;
            s["Sequence Variations"] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(b => string.Join(", ", b.Protein.SequenceVariations
                    .Where(d => psm.OneBasedStartResidueInProtein <= d.OneBasedBeginPosition && d.OneBasedBeginPosition <= psm.OneBasedEndResidueInProtein)
                    .Select(d => d.OriginalSequence + d.OneBasedBeginPosition.ToString() + d.VariantSequence)))).ResolvedString;
            s["Organism Name"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.Organism)).ResolvedString;
            s["Contaminant"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.IsContaminant ? "Y" : "N")).ResolvedString;
            s["Decoy"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.IsDecoy ? "Y" : "N")).ResolvedString;
            s["Peptide Description"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.PeptideDescription)).ResolvedString;
            s["Start and End Residues In Protein"] =
                pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => ("[" + b.OneBasedStartResidueInProtein.ToString(CultureInfo.InvariantCulture) + " to " +
                    b.OneBasedEndResidueInProtein.ToString(CultureInfo.InvariantCulture) + "]"))).ResolvedString;
            s["Previous Amino Acid"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.PreviousAminoAcid.ToString())).ResolvedString;
            s["Next Amino Acid"] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NextAminoAcid.ToString())).ResolvedString;

            string allScores = " ";
            string theoreticalsSearched = " ";
            if (!pepWithModsIsNull && psm.FdrInfo != null && psm.FdrInfo.CalculateEValue)
            {
                allScores = string.Join(";", psm.AllScores.Select(p => p.ToString("F2", CultureInfo.InvariantCulture)));
                theoreticalsSearched = psm.AllScores.Count.ToString();
            }

            s["All Scores"] = allScores;
            s["Theoreticals Searched"] = theoreticalsSearched;
            s["Decoy/Contaminant/Target"] = pepWithModsIsNull ? " " : psm.IsDecoy ? "D" : psm.IsContaminant ? "C" : "T";
        }

        /// <summary>
        ///
        /// </summary>
        private static void AddMatchedIonsData(Dictionary<string, string> s, PeptideSpectralMatch psm)
        {
            bool nullPsm = (psm == null);

            StringBuilder seriesStringBuilder = new StringBuilder();
            StringBuilder mzStringBuilder = new StringBuilder();
            StringBuilder fragmentDaErrorStringBuilder = new StringBuilder();
            StringBuilder fragmentPpmErrorStringBuilder = new StringBuilder();
            StringBuilder fragmentIntensityStringBuilder = new StringBuilder();
            List<StringBuilder> stringBuilders = new List<StringBuilder> { seriesStringBuilder, mzStringBuilder, fragmentDaErrorStringBuilder, fragmentPpmErrorStringBuilder, fragmentIntensityStringBuilder };

            if (!nullPsm)
            {
                var matchedIons = psm.MatchedFragmentIons;
                if (matchedIons == null)
                {
                    matchedIons = psm.PeptidesToMatchingFragments.First().Value;
                }

                // using ", " instead of "," improves human readability
                const string delimiter = ", ";

                var matchedIonsGroupedByProductType = matchedIons.GroupBy(i => i.NeutralTheoreticalProduct.ProductType).OrderBy(i => i.Key).ToList();

                foreach (var productType in matchedIonsGroupedByProductType)
                {
                    var products = productType.OrderBy(p => p.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber)
                        .ToList();

                    stringBuilders.ForEach(p => p.Append("["));

                    for (int i = 0; i < products.Count; i++)
                    {
                        MatchedFragmentIon ion = products[i];
                        string ionLabel;

                        double massError = ion.Mz.ToMass(ion.Charge) - ion.NeutralTheoreticalProduct.NeutralMass;
                        double ppmMassError = massError / ion.NeutralTheoreticalProduct.NeutralMass * 1e6;

                        if (ion.NeutralTheoreticalProduct.NeutralLoss == 0)
                        {
                            // no neutral loss
                            ionLabel = ion.NeutralTheoreticalProduct.ProductType + "" + ion.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber + "+" + ion.Charge;
                        }
                        else
                        {
                            // ion label with neutral loss
                            ionLabel = "(" + ion.NeutralTheoreticalProduct.ProductType + "" + ion.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber
                                + "-" + ion.NeutralTheoreticalProduct.NeutralLoss.ToString("F2") + ")" + "+" + ion.Charge;
                        }

                        // append ion label
                        seriesStringBuilder.Append(ionLabel);

                        // append experimental m/z
                        mzStringBuilder.Append(ionLabel + ":" + ion.Mz.ToString("F5"));

                        // append absolute mass error
                        fragmentDaErrorStringBuilder.Append(ionLabel + ":" + massError.ToString("F5"));

                        // append ppm mass error
                        fragmentPpmErrorStringBuilder.Append(ionLabel + ":" + ppmMassError.ToString("F2"));

                        // append fragment ion intensity
                        fragmentIntensityStringBuilder.Append(ionLabel + ":" + ion.Intensity.ToString("F0"));

                        // append delimiter ", "
                        if (i < products.Count - 1)
                        {
                            stringBuilders.ForEach(p => p.Append(delimiter));
                        }
                    }

                    // append product type delimiter
                    stringBuilders.ForEach(p => p.Append("];"));
                }
            }

            // save ion series strings to output dictionary
            s["Matched Ion Series"] = nullPsm ? " " : seriesStringBuilder.ToString().TrimEnd(';');
            s["Matched Ion Mass-To-Charge Ratios"] = nullPsm ? " " : mzStringBuilder.ToString().TrimEnd(';');
            s["Matched Ion Mass Diff (Da)"] = nullPsm ? " " : fragmentDaErrorStringBuilder.ToString().TrimEnd(';');
            s["Matched Ion Mass Diff (Ppm)"] = nullPsm ? " " : fragmentPpmErrorStringBuilder.ToString().TrimEnd(';');
            s["Matched Ion Intensities"] = nullPsm ? " " : fragmentIntensityStringBuilder.ToString().TrimEnd(';');

            // number of matched ions
            s["Matched Ion Counts"] = nullPsm ? " " : psm.MatchedFragmentIons.Count.ToString();
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