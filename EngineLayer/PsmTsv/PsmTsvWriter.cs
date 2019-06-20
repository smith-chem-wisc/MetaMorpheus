using Chemistry;
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

        internal static (string ResolvedString, double? ResolvedValue) ResolveF2(IEnumerable<double> enumerable)
        {
            var list = enumerable.ToList();
            if (list.Max() - list.Min() < ToleranceForDoubleResolutionF2)
            {
                return (list.Average().ToString("F2", CultureInfo.InvariantCulture), list.Average());
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F2", CultureInfo.InvariantCulture))));
                return (returnString, null);
            }
        }

        internal static (string ResolvedString, double? ResolvedValue) Resolve(IEnumerable<double> enumerable)
        {
            var list = enumerable.ToList();
            if (list.Max() - list.Min() < ToleranceForDoubleResolutionF5)
            {
                return (list.Average().ToString("F5", CultureInfo.InvariantCulture), list.Average());
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))));
                return (returnString, null);
            }
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

        internal static void AddBasicMatchData(Dictionary<string, string> s, PeptideSpectralMatch psm)
        {
            s[PsmTsvHeader.FileName] = psm == null ? " " : Path.GetFileNameWithoutExtension(psm.FullFilePath);
            s[PsmTsvHeader.Ms2ScanNumber] = psm == null ? " " : psm.ScanNumber.ToString(CultureInfo.InvariantCulture);
            s[PsmTsvHeader.Ms2ScanRetentionTime] = psm == null ? " " : psm.ScanRetentionTime.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.NumExperimentalPeaks] = psm == null ? " " : psm.ScanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.TotalIonCurrent] = psm == null ? " " : psm.TotalIonCurrent.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.PrecursorScanNum] = psm == null ? " " : psm.PrecursorScanNumber.HasValue ? psm.PrecursorScanNumber.Value.ToString(CultureInfo.InvariantCulture) : "unknown";
            s[PsmTsvHeader.PrecursorCharge] = psm == null ? " " : psm.ScanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.PrecursorMz] = psm == null ? " " : psm.ScanPrecursorMonoisotopicPeakMz.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.PrecursorMass] = psm == null ? " " : psm.ScanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.Score] = psm == null ? " " : psm.Score.ToString("F3", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.DeltaScore] = psm == null ? " " : psm.DeltaScore.ToString("F3", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.Notch] = psm == null ? " " : Resolve(psm.BestMatchingPeptides.Select(p => p.Notch)).ResolvedString;
        }

        internal static void AddPeptideSequenceData(Dictionary<string, string> s, PeptideSpectralMatch psm, IReadOnlyDictionary<string, int> ModsToWritePruned)
        {
            bool pepWithModsIsNull = psm == null || psm.BestMatchingPeptides == null || !psm.BestMatchingPeptides.Any();

            List<PeptideWithSetModifications> pepsWithMods = pepWithModsIsNull ? null : psm.BestMatchingPeptides.Select(p => p.Peptide).ToList();

            s[PsmTsvHeader.BaseSequence] = pepWithModsIsNull ? " " : (psm.BaseSequence ?? Resolve(pepWithModsIsNull ? null : pepsWithMods.Select(b => b.BaseSequence)).ResolvedString);
            s[PsmTsvHeader.FullSequence] = pepWithModsIsNull ? " " : (psm.FullSequence != null ? psm.FullSequence : Resolve(pepWithModsIsNull ? null : pepsWithMods.Select(b => b.FullSequence)).ResolvedString);
            s[PsmTsvHeader.EssentialSequence] = pepWithModsIsNull ? " " : (psm.EssentialSequence != null ? psm.EssentialSequence : Resolve(pepWithModsIsNull ? null : pepsWithMods.Select(b => b.EssentialSequence(ModsToWritePruned))).ResolvedString);
            s[PsmTsvHeader.PsmCount] = pepWithModsIsNull ? " " : psm.PsmCount.ToString();
            s[PsmTsvHeader.Mods] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus)).ResolvedString;
            s[PsmTsvHeader.ModsChemicalFormulas] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(p => p.AllModsOneIsNterminus.Select(v => v.Value))).ResolvedString;
            s[PsmTsvHeader.ModsCombinedChemicalFormula] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus.Select(c => (c.Value as Modification)))).ResolvedString;
            s[PsmTsvHeader.NumVariableMods] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NumVariableMods)).ResolvedString;
            s[PsmTsvHeader.MissedCleavages] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MissedCleavages.ToString(CultureInfo.InvariantCulture))).ResolvedString;
            s[PsmTsvHeader.PeptideMonoMass] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MonoisotopicMass)).ResolvedString;
            s[PsmTsvHeader.MassDiffDa] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => psm.ScanPrecursorMass - b.MonoisotopicMass)).ResolvedString;
            s[PsmTsvHeader.MassDiffPpm] = pepWithModsIsNull ? " " : ResolveF2(pepsWithMods.Select(b => ((psm.ScanPrecursorMass - b.MonoisotopicMass) / b.MonoisotopicMass * 1e6))).ResolvedString;
            s[PsmTsvHeader.ProteinAccession] = pepWithModsIsNull ? " " : (psm.ProteinAccession != null ? psm.ProteinAccession : Resolve(pepsWithMods.Select(b => b.Protein.Accession), psm.FullSequence).ResolvedString);
            s[PsmTsvHeader.ProteinName] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.FullName), psm.FullSequence).ResolvedString;
            s[PsmTsvHeader.GeneName] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => string.Join(", ", b.Protein.GeneNames.Select(d => $"{d.Item1}:{d.Item2}"))), psm.FullSequence).ResolvedString;
            s[PsmTsvHeader.OrganismName] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.Organism)).ResolvedString;
            s[PsmTsvHeader.IntersectingSequenceVariations] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(b => string.Join(", ", b.Protein.AppliedSequenceVariations
                    .Where(av => IntersectsWithVariation(b, av, false))
                    .Select(av => SequenceVariantString(b, av))))).ResolvedString;
            s[PsmTsvHeader.IdentifiedSequenceVariations] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(b => string.Join(", ", b.Protein.AppliedSequenceVariations
                    .Where(av => IntersectsWithVariation(b, av, true))
                    .Select(av => SequenceVariantString(b, av))))).ResolvedString;
            s[PsmTsvHeader.SpliceSites] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(b => string.Join(", ", b.Protein.SpliceSites
                    .Where(d => Includes(b, d))
                    .Select(d => $"{d.OneBasedBeginPosition.ToString()}-{d.OneBasedEndPosition.ToString()}")))).ResolvedString;
            s[PsmTsvHeader.Contaminant] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.IsContaminant ? "Y" : "N")).ResolvedString;
            s[PsmTsvHeader.Decoy] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Protein.IsDecoy ? "Y" : "N")).ResolvedString;
            s[PsmTsvHeader.PeptideDesicription] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.PeptideDescription)).ResolvedString;
            s[PsmTsvHeader.StartAndEndResiduesInProtein] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(b => ($"[{b.OneBasedStartResidueInProtein.ToString(CultureInfo.InvariantCulture)} to {b.OneBasedEndResidueInProtein.ToString(CultureInfo.InvariantCulture)}]")), psm.FullSequence).ResolvedString;
            s[PsmTsvHeader.PreviousAminoAcid] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.PreviousAminoAcid.ToString())).ResolvedString;
            s[PsmTsvHeader.NextAminoAcid] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NextAminoAcid.ToString())).ResolvedString;

            string theoreticalsSearched = " ";
            s[PsmTsvHeader.TheoreticalsSearched] = theoreticalsSearched;
            s[PsmTsvHeader.DecoyContaminantTarget] = pepWithModsIsNull ? " " : psm.IsDecoy ? "D" : psm.IsContaminant ? "C" : "T";
        }

        /// <summary>
        /// Determines whether a peptide includes a splice site
        /// </summary>
        /// <param name="pep"></param>
        /// <param name="site"></param>
        /// <returns></returns>
        private static bool Includes(PeptideWithSetModifications pep, SpliceSite site)
        {
            return pep.OneBasedStartResidueInProtein <= site.OneBasedBeginPosition && pep.OneBasedEndResidueInProtein >= site.OneBasedEndPosition;
        }

        /// <summary>
        /// Checks for an intersection between a peptide and applied variant that shows a sequence change.
        /// </summary>
        /// <param name="pep"></param>
        /// <param name="appliedVariation"></param>
        /// <returns></returns>
        private static bool IntersectsWithVariation(PeptideWithSetModifications pep, SequenceVariation appliedVariation, bool checkUnique)
        {
            // does it intersect?
            int intersectOneBasedStart = Math.Max(pep.OneBasedStartResidueInProtein, appliedVariation.OneBasedBeginPosition);
            int intersectOneBasedEnd = Math.Min(pep.OneBasedEndResidueInProtein, appliedVariation.OneBasedEndPosition);
            if (intersectOneBasedEnd < intersectOneBasedStart)
            {
                return false;
            }
            else if (!checkUnique)
            {
                return true;
            }
            else
            {
                // if the original sequence is too short or long, the intersect of the peptide and variant is unique
                int intersectSize = intersectOneBasedEnd - intersectOneBasedStart + 1;
                int variantZeroBasedStart = intersectOneBasedStart - appliedVariation.OneBasedBeginPosition;
                bool origSeqIsShort = appliedVariation.OriginalSequence.Length - variantZeroBasedStart < intersectSize;
                bool origSeqIsLong = appliedVariation.OriginalSequence.Length > intersectSize && pep.OneBasedEndResidueInProtein > intersectOneBasedEnd;
                if (origSeqIsShort || origSeqIsLong)
                {
                    return true;
                }

                // is the variant sequence intersecting the peptide different than the original sequence?
                string originalAtIntersect = appliedVariation.OriginalSequence.Substring(intersectOneBasedStart - appliedVariation.OneBasedBeginPosition, intersectSize);
                string variantAtIntersect = appliedVariation.VariantSequence.Substring(intersectOneBasedStart - appliedVariation.OneBasedBeginPosition, intersectSize);
                return originalAtIntersect != variantAtIntersect;
            }
        }

        /// <summary>
        /// Makes the string representing a detected sequence variation, including any modifications on a variant amino acid
        /// </summary>
        /// <param name="p"></param>
        /// <param name="d"></param>
        /// <returns></returns>
        private static string SequenceVariantString(PeptideWithSetModifications p, SequenceVariation applied)
        {
            var modsOnVariantOneIsNTerm = p.AllModsOneIsNterminus
                .Where(kv => kv.Key == 1 && applied.OneBasedBeginPosition == 1 || applied.OneBasedBeginPosition <= kv.Key - 2 + p.OneBasedStartResidueInProtein && kv.Key - 2 + p.OneBasedStartResidueInProtein <= applied.OneBasedEndPosition)
                .ToDictionary(kv => kv.Key - applied.OneBasedBeginPosition + 1, kv => kv.Value);
            PeptideWithSetModifications variantWithAnyMods = new PeptideWithSetModifications(p.Protein, p.DigestionParams, applied.OneBasedBeginPosition, applied.OneBasedEndPosition, p.CleavageSpecificityForFdrCategory, p.PeptideDescription, p.MissedCleavages, modsOnVariantOneIsNTerm, p.NumFixedMods);
            return $"{applied.OriginalSequence}{applied.OneBasedBeginPosition}{variantWithAnyMods.FullSequence}";
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
            s[PsmTsvHeader.MatchedIonSeries] = nullPsm ? " " : seriesStringBuilder.ToString().TrimEnd(';');
            s[PsmTsvHeader.MatchedIonMzRatios] = nullPsm ? " " : mzStringBuilder.ToString().TrimEnd(';');
            s[PsmTsvHeader.MatchedIonMassDiffDa] = nullPsm ? " " : fragmentDaErrorStringBuilder.ToString().TrimEnd(';');
            s[PsmTsvHeader.MatchedIonMassDiffPpm] = nullPsm ? " " : fragmentPpmErrorStringBuilder.ToString().TrimEnd(';');
            s[PsmTsvHeader.MatchedIonIntensities] = nullPsm ? " " : fragmentIntensityStringBuilder.ToString().TrimEnd(';');

            // number of matched ions
            s[PsmTsvHeader.MatchedIonCounts] = nullPsm ? " " : matchedIons.Count.ToString();
        }

        internal static void AddMatchScoreData(Dictionary<string, string> s, PeptideSpectralMatch peptide)
        {
            string localizedScores = " ";
            string improvementPossible = " ";
            if (peptide != null && peptide.LocalizedScores != null)
            {
                localizedScores = GlobalVariables.CheckLengthOfOutput(("[" + string.Join(",", peptide.LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]"));
                improvementPossible = (peptide.LocalizedScores.Max() - peptide.Score).ToString("F3", CultureInfo.InvariantCulture);
            }
            s[PsmTsvHeader.LocalizedScores] = localizedScores;
            s[PsmTsvHeader.ImprovementPossible] = improvementPossible;

            string cumulativeTarget = " ";
            string cumulativeDecoy = " ";
            string qValue = " ";
            string cumulativeTargetNotch = " ";
            string cumulativeDecoyNotch = " ";
            string qValueNotch = " ";
            string PEP = " ";
            string PEP_Qvalue = " ";
            if (peptide != null && peptide.FdrInfo != null)
            {
                cumulativeTarget = peptide.FdrInfo.CumulativeTarget.ToString(CultureInfo.InvariantCulture);
                cumulativeDecoy = peptide.FdrInfo.CumulativeDecoy.ToString(CultureInfo.InvariantCulture);
                qValue = peptide.FdrInfo.QValue.ToString("F6", CultureInfo.InvariantCulture);
                cumulativeTargetNotch = peptide.FdrInfo.CumulativeTargetNotch.ToString(CultureInfo.InvariantCulture);
                cumulativeDecoyNotch = peptide.FdrInfo.CumulativeDecoyNotch.ToString(CultureInfo.InvariantCulture);
                qValueNotch = peptide.FdrInfo.QValueNotch.ToString("F6", CultureInfo.InvariantCulture);
                PEP = peptide.FdrInfo.PEP.ToString();
                PEP_Qvalue = peptide.FdrInfo.PEP_QValue.ToString();
            }
            s[PsmTsvHeader.CumulativeTarget] = cumulativeTarget;
            s[PsmTsvHeader.CumulativeDecoy] = cumulativeDecoy;
            s[PsmTsvHeader.QValue] = qValue;
            s[PsmTsvHeader.CumulativeTargetNotch] = cumulativeTargetNotch;
            s[PsmTsvHeader.CumulativeDecoyNotch] = cumulativeDecoyNotch;
            s[PsmTsvHeader.QValueNotch] = qValueNotch;
            s[PsmTsvHeader.PEP] = PEP;
            s[PsmTsvHeader.PEP_QValue] = PEP_Qvalue;
        }
    }
}