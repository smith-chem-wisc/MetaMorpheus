﻿using Chemistry;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Omics;
using Omics.Modifications;


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

        internal static void AddBasicMatchData(Dictionary<string, string> s, SpectralMatch psm)
        {
            s[PsmTsvHeader.FileName] = psm == null ? " " : Path.GetFileNameWithoutExtension(psm.FullFilePath);
            s[PsmTsvHeader.Ms2ScanNumber] = psm == null ? " " : psm.ScanNumber.ToString(CultureInfo.InvariantCulture);
            s[PsmTsvHeader.Ms2ScanRetentionTime] = psm == null ? " " : psm.ScanRetentionTime.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.NumExperimentalPeaks] = psm == null ? " " : psm.ScanExperimentalPeaks.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.TotalIonCurrent] = psm == null ? " " : psm.TotalIonCurrent.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.PrecursorScanNum] = psm == null ? " " : psm.PrecursorScanNumber.HasValue ? psm.PrecursorScanNumber.Value.ToString(CultureInfo.InvariantCulture) : "unknown";
            s[PsmTsvHeader.PrecursorCharge] = psm == null ? " " : psm.ScanPrecursorCharge.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.PrecursorIntensity] = psm == null ? " " : psm.PrecursorScanIntensity.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.PrecursorMz] = psm == null ? " " : psm.ScanPrecursorMonoisotopicPeakMz.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.PrecursorMass] = psm == null ? " " : psm.ScanPrecursorMass.ToString("F5", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.Score] = psm == null ? " " : psm.Score.ToString("F3", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.DeltaScore] = psm == null ? " " : psm.DeltaScore.ToString("F3", CultureInfo.InvariantCulture);
            s[PsmTsvHeader.Notch] = psm == null ? " " : Resolve(psm.BestMatchingBioPolymersWithSetMods.Select(p => p.Notch)).ResolvedString;
        }

        internal static void AddPeptideSequenceData(Dictionary<string, string> s, SpectralMatch sm, IReadOnlyDictionary<string, int> ModsToWritePruned)
        {
            bool pepWithModsIsNull = sm == null || sm.BestMatchingBioPolymersWithSetMods == null || !sm.BestMatchingBioPolymersWithSetMods.Any();

            List<IBioPolymerWithSetMods> pepsWithMods = pepWithModsIsNull ? null : sm.BestMatchingBioPolymersWithSetMods.Select(p => p.Peptide).ToList();

            s[PsmTsvHeader.BaseSequence] = pepWithModsIsNull ? " " : (sm.BaseSequence ?? Resolve(pepWithModsIsNull ? null : pepsWithMods.Select(b => b.BaseSequence)).ResolvedString);
            s[PsmTsvHeader.FullSequence] = pepWithModsIsNull ? " " : (sm.FullSequence != null ? sm.FullSequence : Resolve(pepWithModsIsNull ? null : pepsWithMods.Select(b => b.FullSequence)).ResolvedString);
            s[PsmTsvHeader.EssentialSequence] = pepWithModsIsNull ? " " : (sm.EssentialSequence != null ? sm.EssentialSequence : Resolve(pepWithModsIsNull ? null : pepsWithMods.Select(b => b.EssentialSequence(ModsToWritePruned))).ResolvedString);
            string geneString = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => string.Join(", ", b.Parent.GeneNames.Select(d => $"{d.Item1}:{d.Item2}"))), sm.FullSequence).ResolvedString;
            s[PsmTsvHeader.AmbiguityLevel] = ProteoformLevelClassifier.ClassifyPrSM(s[PsmTsvHeader.FullSequence], geneString);
            s[PsmTsvHeader.PsmCount] = pepWithModsIsNull ? " " : sm.PsmCount.ToString();
            s[PsmTsvHeader.Mods] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus)).ResolvedString;
            s[PsmTsvHeader.ModsChemicalFormulas] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(p => p.AllModsOneIsNterminus.Select(v => v.Value))).ResolvedString;
            s[PsmTsvHeader.ModsCombinedChemicalFormula] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.AllModsOneIsNterminus.Select(c => (c.Value as Modification)))).ResolvedString;
            s[PsmTsvHeader.NumVariableMods] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NumVariableMods)).ResolvedString;
            s[PsmTsvHeader.MissedCleavages] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MissedCleavages.ToString(CultureInfo.InvariantCulture))).ResolvedString;
            s[PsmTsvHeader.PeptideMonoMass] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.MonoisotopicMass)).ResolvedString;
            s[PsmTsvHeader.MassDiffDa] = pepWithModsIsNull ? " " : Resolve(sm.PrecursorMassErrorDa).ResolvedString;
            s[PsmTsvHeader.MassDiffPpm] = pepWithModsIsNull ? " " : Resolve(sm.PrecursorMassErrorPpm).ResolvedString;
            s[PsmTsvHeader.ProteinAccession] = pepWithModsIsNull ? " " : (sm.Accession != null ? sm.Accession : Resolve(pepsWithMods.Select(b => b.Parent.Accession), sm.FullSequence).ResolvedString);
            s[PsmTsvHeader.ProteinName] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.FullName), sm.FullSequence).ResolvedString;
            s[PsmTsvHeader.GeneName] = geneString;
            s[PsmTsvHeader.OrganismName] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.Organism)).ResolvedString;

            if (sm is PeptideSpectralMatch psm || sm is null)
            {
                s[PsmTsvHeader.IdentifiedSequenceVariations] = pepWithModsIsNull ? " " :
               Resolve(pepsWithMods.Select(p => p as PeptideWithSetModifications)
                   .Select(b => string.Join(", ", b.Protein.AppliedSequenceVariations
                   .Where(av => b.IntersectsAndIdentifiesVariation(av).identifies)
                   .Select(av => b.SequenceVariantString(av, b.IntersectsAndIdentifiesVariation(av).intersects))))).ResolvedString;
                s[PsmTsvHeader.SpliceSites] = pepWithModsIsNull ? " " :
                    Resolve(pepsWithMods.Select(p => p as PeptideWithSetModifications)
                        .Select(b => string.Join(", ", b.Protein.SpliceSites
                        .Where(d => Includes(b, d))
                        .Select(d => $"{d.OneBasedBeginPosition.ToString()}-{d.OneBasedEndPosition.ToString()}")))).ResolvedString;
            }

            s[PsmTsvHeader.Contaminant] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.IsContaminant ? "Y" : "N")).ResolvedString;
            s[PsmTsvHeader.Decoy] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.IsDecoy ? "Y" : "N")).ResolvedString;
            s[PsmTsvHeader.PeptideDesicription] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Description)).ResolvedString;
            s[PsmTsvHeader.StartAndEndResiduesInProtein] = pepWithModsIsNull ? " " :
                Resolve(pepsWithMods.Select(b => ($"[{b.OneBasedStartResidue.ToString(CultureInfo.InvariantCulture)} to {b.OneBasedEndResidue.ToString(CultureInfo.InvariantCulture)}]")), sm.FullSequence).ResolvedString;
            s[PsmTsvHeader.PreviousAminoAcid] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.PreviousResidue.ToString())).ResolvedString;
            s[PsmTsvHeader.NextAminoAcid] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.NextResidue.ToString())).ResolvedString;

            string theoreticalsSearched = " ";
            s[PsmTsvHeader.TheoreticalsSearched] = theoreticalsSearched;
            s[PsmTsvHeader.DecoyContaminantTarget] = pepWithModsIsNull ? " " : Resolve(pepsWithMods.Select(b => b.Parent.IsDecoy ? "D" : b.Parent.IsContaminant ? "C" : "T")).ResolvedString;
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
            s[PsmTsvHeader.MatchedIonSeries] = nullPsm ? " " : seriesStringBuilder.ToString().TrimEnd(';');
            s[PsmTsvHeader.MatchedIonMzRatios] = nullPsm ? " " : mzStringBuilder.ToString().TrimEnd(';');
            s[PsmTsvHeader.MatchedIonMassDiffDa] = nullPsm ? " " : fragmentDaErrorStringBuilder.ToString().TrimEnd(';');
            s[PsmTsvHeader.MatchedIonMassDiffPpm] = nullPsm ? " " : fragmentPpmErrorStringBuilder.ToString().TrimEnd(';');
            s[PsmTsvHeader.MatchedIonIntensities] = nullPsm ? " " : fragmentIntensityStringBuilder.ToString().TrimEnd(';');

            // number of matched ions
            s[PsmTsvHeader.MatchedIonCounts] = nullPsm ? " " : matchedIons.Count.ToString();
        }

        internal static void AddMatchScoreData(Dictionary<string, string> s, SpectralMatch peptide, bool writePeptideLevelFdr = false)
        {
            string spectralAngle = peptide == null ? " " : peptide.SpectralAngle.ToString("F4");
            string localizedScores = " ";
            string improvementPossible = " ";
            if (peptide != null && peptide.LocalizedScores != null)
            {
                localizedScores = GlobalVariables.CheckLengthOfOutput(("[" + string.Join(",", peptide.LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]"));
                improvementPossible = (peptide.LocalizedScores.Max() - peptide.Score).ToString("F3", CultureInfo.InvariantCulture);
            }
            s[PsmTsvHeader.SpectralAngle] = spectralAngle;
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

            if (peptide != null && peptide.GetFdrInfo(writePeptideLevelFdr) != null)
            {
                cumulativeTarget = peptide.GetFdrInfo(writePeptideLevelFdr).CumulativeTarget.ToString(CultureInfo.InvariantCulture);
                cumulativeDecoy = peptide.GetFdrInfo(writePeptideLevelFdr).CumulativeDecoy.ToString(CultureInfo.InvariantCulture);
                qValue = peptide.GetFdrInfo(writePeptideLevelFdr).QValue.ToString("F6", CultureInfo.InvariantCulture);
                cumulativeTargetNotch = peptide.GetFdrInfo(writePeptideLevelFdr).CumulativeTargetNotch.ToString(CultureInfo.InvariantCulture);
                cumulativeDecoyNotch = peptide.GetFdrInfo(writePeptideLevelFdr).CumulativeDecoyNotch.ToString(CultureInfo.InvariantCulture);
                qValueNotch = peptide.GetFdrInfo(writePeptideLevelFdr).QValueNotch.ToString("F6", CultureInfo.InvariantCulture);
                PEP = peptide.GetFdrInfo(writePeptideLevelFdr).PEP.ToString();
                PEP_Qvalue = peptide.GetFdrInfo(writePeptideLevelFdr).PEP_QValue.ToString();
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