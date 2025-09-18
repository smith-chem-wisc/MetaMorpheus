using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Chemistry;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;

namespace EngineLayer
{
    public static class PsmTsvWriter
    {
        internal const double ToleranceForDoubleResolutionF5 = 1e-6;
        internal const double ToleranceForDoubleResolutionF2 = 1e-3;

        /// <summary>
        /// Resolve Methods()
        /// if all 'values' are the same this returns the one value otherwise you get a separated list of all values in their original order.
        /// for example:
        /// Notches 1,1,1,1 returns as 1
        /// Notches 1,0,1,0 returns as 1|0|1|0
        /// </summary>
        internal static (string ResolvedString, ChemicalFormula ResolvedValue) Resolve(IEnumerable<IEnumerable<Modification>> enumerable)
        {
            var list = enumerable.ToList();
            ChemicalFormula firstChemFormula = new ChemicalFormula();
            foreach (var firstMods in list[0])
            {
                if (firstMods == null || firstMods.ChemicalFormula == null)
                {
                    return ("unknown", null);
                }
                firstChemFormula.Add(firstMods.ChemicalFormula);
            }

            bool equals = true;
            List<ChemicalFormula> formulas = new List<ChemicalFormula>();
            foreach (var anEnum in list)
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
                if (!firstChemFormula.Equals(fhere))
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
                return (firstChemFormula.Formula, firstChemFormula);
            }
        }

        internal static (string ResolvedString, Dictionary<string, int> ResolvedValue) Resolve(IEnumerable<Dictionary<int, Modification>> enumerable)
        {
            var list = enumerable.ToList();
            Dictionary<string, int> firstDict = list[0].Values.OrderBy(b => b.IdWithMotif).GroupBy(b => b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());

            bool equals = true;
            foreach (var dict in list)
            {
                Dictionary<string, int> okTest = dict.Values.OrderBy(b => b.IdWithMotif).GroupBy(b => b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());
                if (!firstDict.SequenceEqual(okTest))
                {
                    equals = false;
                    break;
                }
            }
            if (!equals)
            {
                var returnString = string.Join("|", list.Select(b => string.Join(" ", b.Values.Select(c => c.IdWithMotif).OrderBy(c => c))));
                returnString = GlobalVariables.CheckLengthOfOutput(returnString);
                return (returnString, null);
            }
            else
            {
                return (string.Join(" ", list[0].Values.Select(c => c.IdWithMotif).OrderBy(c => c)), firstDict);
            }
        }

        internal static (string ResolvedString, double? ResolvedValue) Resolve(IEnumerable<double> enumerable)
        {
            var list = enumerable.ToList();
            double min = list.Min();
            double max = list.Max();
            double avg = list.Average();

            // If all values are within tolerance, check if all are integers
            if (max - min < ToleranceForDoubleResolutionF5)
            {
                bool allInts = list.All(x => Math.Abs(x - Math.Round(x)) < ToleranceForDoubleResolutionF2);
                if (!allInts)
                    return (avg.ToString("F5", CultureInfo.InvariantCulture), avg);

                int intVal = (int)Math.Round(avg);
                return (intVal.ToString(CultureInfo.InvariantCulture), intVal);
            }

            var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))));
            return (returnString, null);
        }

        internal static (string ResolvedString, int? ResolvedValue) Resolve(IEnumerable<int> enumerable)
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

        internal static (string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable)
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

        internal static (string ResolvedString, string ResolvedValue) Resolve(IEnumerable<string> enumerable, string ambiguousIfNull)
        {
            var list = enumerable.ToList();
            string first = list.FirstOrDefault(b => b != null);
            // Only first if list is either all null or all equal to the first
            if (list.All(b => b == null) || list.All(b => first.Equals(b)))
            {
                return (first, first);
            }
            // use only distinct names if all of the base sequences are the same
            else if (ambiguousIfNull != null)
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Distinct()));
                return (returnString, null);
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list));
                return (returnString, null);
            }
        }

        internal static void AddBasicMatchData(Dictionary<string, string> s, SpectralMatch psm, bool includeOneOverK0Column = false)
        {
            s[SpectrumMatchFromTsvHeader.FileName] = psm == null ? " " : Path.GetFileNameWithoutExtension(psm.FullFilePath);
            s[SpectrumMatchFromTsvHeader.Ms2ScanNumber] = psm == null ? " " : psm.ScanNumber.ToString(CultureInfo.InvariantCulture);
            s[SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime] = psm == null ? " " : psm.ScanRetentionTime.ToString("F5", CultureInfo.InvariantCulture);
            s[SpectrumMatchFromTsvHeader.NumExperimentalPeaks] = psm == null ? " " : psm.ScanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture);
            s[SpectrumMatchFromTsvHeader.TotalIonCurrent] = psm == null ? " " : psm.TotalIonCurrent.ToString("F5", CultureInfo.InvariantCulture);
            s[SpectrumMatchFromTsvHeader.PrecursorScanNum] = psm == null ? " " : psm.PrecursorScanNumber.HasValue ? psm.PrecursorScanNumber.Value.ToString(CultureInfo.InvariantCulture) : "unknown";
            s[SpectrumMatchFromTsvHeader.PrecursorCharge] = psm == null ? " " : psm.ScanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture);
            s[SpectrumMatchFromTsvHeader.PrecursorIntensity] = psm == null ? " " : psm.PrecursorScanIntensity.ToString("F5", CultureInfo.InvariantCulture);
            s[SpectrumMatchFromTsvHeader.PrecursorMz] = psm == null ? " " : psm.ScanPrecursorMonoisotopicPeakMz.ToString("F5", CultureInfo.InvariantCulture);
            s[SpectrumMatchFromTsvHeader.PrecursorMass] = psm == null ? " " : psm.ScanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture);
            if ( includeOneOverK0Column) // This information is only written if one or more spectra have a K0 value, otherwise it is not included in the output
                s[SpectrumMatchFromTsvHeader.OneOverK0] = psm == null ? " " : psm.ScanOneOverK0.HasValue ? psm.ScanOneOverK0.Value.ToString("F5", CultureInfo.InvariantCulture) : "N/A";
            s[SpectrumMatchFromTsvHeader.Score] = psm == null ? " " : psm.Score.ToString("F3", CultureInfo.InvariantCulture);
            s[SpectrumMatchFromTsvHeader.DeltaScore] = psm == null ? " " : psm.DeltaScore.ToString("F3", CultureInfo.InvariantCulture);
            s[SpectrumMatchFromTsvHeader.Notch] = psm == null ? " " : Resolve(psm.BestMatchingBioPolymersWithSetMods.Select(p => p.Notch / MassDiffAcceptor.NotchScalar)).ResolvedString;
        }

        internal static void AddPeptideSequenceData(Dictionary<string, string> s, SpectralMatch sm, IReadOnlyDictionary<string, int> ModsToWritePruned)
        {
            bool pepWithModsIsNull = sm == null || sm.BestMatchingBioPolymersWithSetMods == null || !sm.BestMatchingBioPolymersWithSetMods.Any();

            List<IBioPolymerWithSetMods> pepsWithMods = pepWithModsIsNull ? null : sm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer).ToList();

            s[SpectrumMatchFromTsvHeader.BaseSequence] = pepWithModsIsNull ? " " : sm.BaseSequence ?? Resolve(pepWithModsIsNull ? null : pepsWithMods.Select(b => b.BaseSequence)).ResolvedString;
            s[SpectrumMatchFromTsvHeader.FullSequence] = pepWithModsIsNull ? " " : sm.FullSequence != null ? sm.FullSequence : Resolve(pepWithModsIsNull ? null : pepsWithMods.Select(b => b.FullSequence)).ResolvedString;
            s[SpectrumMatchFromTsvHeader.EssentialSequence] = pepWithModsIsNull ? " " : sm.EssentialSequence != null ? sm.EssentialSequence : Resolve(pepWithModsIsNull ? null : pepsWithMods.Select(b => b.EssentialSequence(ModsToWritePruned))).ResolvedString;
            string geneString = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => string.Join(", ", b.Parent.GeneNames.Select(d => $"{d.Item1}:{d.Item2}"))), sm.FullSequence).ResolvedString;
            s[SpectrumMatchFromTsvHeader.AmbiguityLevel] = ProteoformLevelClassifier.ClassifyPrSM(s[SpectrumMatchFromTsvHeader.FullSequence], geneString);
            s[SpectrumMatchFromTsvHeader.SpectrumMatchCount] = pepWithModsIsNull ? " " : sm.PsmCount.ToString();
            s[SpectrumMatchFromTsvHeader.Mods] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus)).ResolvedString;
            s[SpectrumMatchFromTsvHeader.ModsChemicalFormulas] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(p => p.AllModsOneIsNterminus.Select(v => v.Value))).ResolvedString;
            s[SpectrumMatchFromTsvHeader.ModsCombinedChemicalFormula] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus.Select(c => c.Value as Modification))).ResolvedString;
            s[SpectrumMatchFromTsvHeader.NumVariableMods] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NumVariableMods)).ResolvedString;
            s[SpectrumMatchFromTsvHeader.MissedCleavages] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MissedCleavages.ToString(CultureInfo.InvariantCulture))).ResolvedString;
            s[SpectrumMatchFromTsvHeader.MonoisotopicMass] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MonoisotopicMass)).ResolvedString;
            s[SpectrumMatchFromTsvHeader.MassDiffDa] = pepWithModsIsNull ? " " : Resolve(sm.PrecursorMassErrorDa).ResolvedString;
            s[SpectrumMatchFromTsvHeader.MassDiffPpm] = pepWithModsIsNull ? " " : Resolve(sm.PrecursorMassErrorPpm).ResolvedString;
            s[SpectrumMatchFromTsvHeader.Accession] = pepWithModsIsNull ? " " : sm.Accession != null ? sm.Accession : Resolve(pepsWithMods.Select(b => b.Parent.Accession), sm.FullSequence).ResolvedString;
            s[SpectrumMatchFromTsvHeader.Name] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.FullName), sm.FullSequence).ResolvedString;
            s[SpectrumMatchFromTsvHeader.GeneName] = geneString;
            s[SpectrumMatchFromTsvHeader.OrganismName] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.Organism)).ResolvedString;

            if (sm is PeptideSpectralMatch psm || sm is null && GlobalVariables.AnalyteType != AnalyteType.Oligo)
            {
                s[SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations] = pepWithModsIsNull ? " " :
               Resolve(pepsWithMods.Select(p => p as PeptideWithSetModifications)
                   .Select(b => string.Join(", ", b.Protein.AppliedSequenceVariations
                   .Where(av => b.IntersectsAndIdentifiesVariation(av).identifies)
                   .Select(av => b.SequenceVariantString(av, b.IntersectsAndIdentifiesVariation(av).intersects))))).ResolvedString;
                s[SpectrumMatchFromTsvHeader.SpliceSites] = pepWithModsIsNull ? " " :
                    Resolve(pepsWithMods.Select(p => p as PeptideWithSetModifications)
                        .Select(b => string.Join(", ", b.Protein.SpliceSites
                        .Where(d => Includes(b, d))
                        .Select(d => $"{d.OneBasedBeginPosition.ToString()}-{d.OneBasedEndPosition.ToString()}")))).ResolvedString;
            }

            s[SpectrumMatchFromTsvHeader.Contaminant] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.IsContaminant ? "Y" : "N")).ResolvedString;
            s[SpectrumMatchFromTsvHeader.Decoy] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.IsDecoy ? "Y" : "N")).ResolvedString;
            s[SpectrumMatchFromTsvHeader.Description] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Description)).ResolvedString;
            s[SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(b => $"[{b.OneBasedStartResidue.ToString(CultureInfo.InvariantCulture)} to {b.OneBasedEndResidue.ToString(CultureInfo.InvariantCulture)}]"), sm.FullSequence).ResolvedString;
            s[SpectrumMatchFromTsvHeader.PreviousResidue] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.PreviousResidue.ToString())).ResolvedString;
            s[SpectrumMatchFromTsvHeader.NextResidue] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NextResidue.ToString())).ResolvedString;

            string theoreticalsSearched = " ";
            s[SpectrumMatchFromTsvHeader.TheoreticalsSearched] = theoreticalsSearched;
            s[SpectrumMatchFromTsvHeader.DecoyContaminantTarget] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.IsDecoy ? "D" : b.Parent.IsContaminant ? "C" : "T")).ResolvedString;
        }

        /// <summary>
        /// Determines whether a peptide includes a splice site
        /// </summary>
        /// <param name="pep"></param>
        /// <param name="site"></param>
        /// <returns></returns>
        private static bool Includes(PeptideWithSetModifications pep, SpliceSite site)
        {
            return pep.OneBasedStartResidue <= site.OneBasedBeginPosition && pep.OneBasedEndResidue >= site.OneBasedEndPosition;
        }

        internal static void AddMatchedIonsData(Dictionary<string, string> s, List<MatchedFragmentIon> matchedIons)
        {
            bool nullPsm = matchedIons == null;

            StringBuilder seriesStringBuilder = new StringBuilder();
            StringBuilder mzStringBuilder = new StringBuilder();
            StringBuilder fragmentDaErrorStringBuilder = new StringBuilder();
            StringBuilder fragmentPpmErrorStringBuilder = new StringBuilder();
            StringBuilder fragmentIntensityStringBuilder = new StringBuilder();
            List<StringBuilder> stringBuilders = new List<StringBuilder> { seriesStringBuilder, mzStringBuilder, fragmentDaErrorStringBuilder, fragmentPpmErrorStringBuilder, fragmentIntensityStringBuilder };

            if (!nullPsm)
            {
                // using ", " instead of "," improves human readability
                const string delimiter = ", ";

                var matchedIonsGroupedByProductType = matchedIons.GroupBy(x => new { x.NeutralTheoreticalProduct.ProductType, x.NeutralTheoreticalProduct.SecondaryProductType }).ToList();

                foreach (var productType in matchedIonsGroupedByProductType)
                {
                    var products = productType.OrderBy(p => p.NeutralTheoreticalProduct.FragmentNumber)
                        .ToList();

                    stringBuilders.ForEach(p => p.Append("["));

                    for (int i = 0; i < products.Count; i++)
                    {
                        MatchedFragmentIon ion = products[i];
                        string ionLabel;

                        double massError = ion.Mz.ToMass(ion.Charge) - ion.NeutralTheoreticalProduct.NeutralMass;
                        double ppmMassError = massError / ion.NeutralTheoreticalProduct.NeutralMass * 1e6;

                        ionLabel = ion.Annotation;

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
            s[SpectrumMatchFromTsvHeader.MatchedIonSeries] = nullPsm ? " " : seriesStringBuilder.ToString().TrimEnd(';');
            s[SpectrumMatchFromTsvHeader.MatchedIonMzRatios] = nullPsm ? " " : mzStringBuilder.ToString().TrimEnd(';');
            s[SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa] = nullPsm ? " " : fragmentDaErrorStringBuilder.ToString().TrimEnd(';');
            s[SpectrumMatchFromTsvHeader.MatchedIonMassDiffPpm] = nullPsm ? " " : fragmentPpmErrorStringBuilder.ToString().TrimEnd(';');
            s[SpectrumMatchFromTsvHeader.MatchedIonIntensities] = nullPsm ? " " : fragmentIntensityStringBuilder.ToString().TrimEnd(';');

            // number of matched ions
            s[SpectrumMatchFromTsvHeader.MatchedIonCounts] = nullPsm ? " " : matchedIons.Count.ToString();
        }

        internal static void AddMatchScoreData(Dictionary<string, string> s, SpectralMatch peptide, bool writePeptideLevelFdr = false)
        {
            string spectralAngle = peptide == null ? " " : peptide.SpectralAngle.ToString("F4");
            string localizedScores = " ";
            string improvementPossible = " ";
            if (peptide != null && peptide.LocalizedScores != null)
            {
                localizedScores = GlobalVariables.CheckLengthOfOutput("[" + string.Join(",", peptide.LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]");
                improvementPossible = (peptide.LocalizedScores.Max() - peptide.Score).ToString("F3", CultureInfo.InvariantCulture);
            }
            s[SpectrumMatchFromTsvHeader.SpectralAngle] = spectralAngle;
            s[SpectrumMatchFromTsvHeader.LocalizedScores] = localizedScores;
            s[SpectrumMatchFromTsvHeader.ImprovementPossible] = improvementPossible;

            string cumulativeTarget = " ";
            string cumulativeDecoy = " ";
            string qValue = " ";
            string cumulativeTargetNotch = " ";
            string cumulativeDecoyNotch = " ";
            string qValueNotch = " ";
            string PEP = " ";
            string PEP_Qvalue = " ";

            if (peptide != null && peptide.GetFdrInfo(writePeptideLevelFdr) != null)
            {
                cumulativeTarget = peptide.GetFdrInfo(writePeptideLevelFdr).CumulativeTarget.ToString(CultureInfo.InvariantCulture);
                cumulativeDecoy = peptide.GetFdrInfo(writePeptideLevelFdr).CumulativeDecoy.ToString(CultureInfo.InvariantCulture);
                qValue = peptide.GetFdrInfo(writePeptideLevelFdr).QValue.ToString("F6", CultureInfo.InvariantCulture);
                PEP = peptide.GetFdrInfo(writePeptideLevelFdr).PEP.ToString();
                PEP_Qvalue = peptide.GetFdrInfo(writePeptideLevelFdr).PEP_QValue.ToString();

                // ambiguous notch, has never been resolved by our disambiguation, so take the best of the notches for the fdr columns. 
                if (peptide.Notch == null && peptide.GetFdrInfo(writePeptideLevelFdr).QValueNotch > 1)
                {
                    var min = peptide.BestMatchingBioPolymersWithSetMods.MinBy(b => writePeptideLevelFdr ? b.PeptideQValueNotch : b.QValueNotch);
                    if (min != null && (writePeptideLevelFdr ? min.PeptideQValueNotch : min.QValueNotch).HasValue)
                    {
                        cumulativeTargetNotch = (writePeptideLevelFdr ? min.PeptideCumulativeTargetNotch : min.CumulativeTargetNotch).Value.ToString(CultureInfo.InvariantCulture);
                        cumulativeDecoyNotch = (writePeptideLevelFdr ? min.PeptideCumulativeDecoyNotch : min.CumulativeDecoyNotch).Value.ToString(CultureInfo.InvariantCulture);
                        qValueNotch = (writePeptideLevelFdr ? min.PeptideQValueNotch : min.QValueNotch).Value.ToString("F6", CultureInfo.InvariantCulture);
                    }
                    else
                    {
                        throw new Exception("If the psm is notch ambiguous, the q value notch in the hypothesis should never be null");
                    }
                }
                else
                {
                    cumulativeTargetNotch = peptide.GetFdrInfo(writePeptideLevelFdr).CumulativeTargetNotch.ToString("F6", CultureInfo.InvariantCulture);
                    cumulativeDecoyNotch = peptide.GetFdrInfo(writePeptideLevelFdr).CumulativeDecoyNotch.ToString("F6", CultureInfo.InvariantCulture);
                    qValueNotch = peptide.GetFdrInfo(writePeptideLevelFdr).QValueNotch.ToString("F6", CultureInfo.InvariantCulture);
                }
            }

            s[SpectrumMatchFromTsvHeader.CumulativeTarget] = cumulativeTarget;
            s[SpectrumMatchFromTsvHeader.CumulativeDecoy] = cumulativeDecoy;
            s[SpectrumMatchFromTsvHeader.QValue] = qValue;
            s[SpectrumMatchFromTsvHeader.CumulativeTargetNotch] = cumulativeTargetNotch;
            s[SpectrumMatchFromTsvHeader.CumulativeDecoyNotch] = cumulativeDecoyNotch;
            s[SpectrumMatchFromTsvHeader.QValueNotch] = qValueNotch;
            s[SpectrumMatchFromTsvHeader.PEP] = PEP;
            s[SpectrumMatchFromTsvHeader.PEP_QValue] = PEP_Qvalue;
        }
    }
}