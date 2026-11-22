using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using GuiFunctions.Util;
using EngineLayer.ModificationAnalysis;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// <see cref="CommonParameters.IdentificationFilter"/>: the user picks the value that decides identification, the PEP
    /// q-value, the q-value notch or the pooled q-value, and every step of a search uses it at the one
    /// <see cref="CommonParameters.FilterThreshold"/>. Where PEP was not trained for a match, the PEP q-value filter judges
    /// that match by its q-value notch instead, and the search says so. Unnamed, the filter is the PEP q-value bottom-up and
    /// the q-value notch for top-down and RNA; a task file saved before FilterThreshold existed keeps what it meant.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // the whole-search tests set FdrAnalysisEngine.QvalueThresholdOverride, a process-wide static
    public static class IdentificationFilterTests
    {
        private static CommonParameters Chosen(FilterType filter, double threshold = 0.01) =>
            new CommonParameters(identificationFilter: filter, filterThreshold: threshold);

        /// <summary>
        /// <paramref name="count"/> target matches to <paramref name="sequence"/>, digested with <paramref name="protease"/>.
        /// The first <paramref name="confident"/> have q-value, notch and PEP q-value 0.001, the rest 0.5. When PEP was not
        /// trained, every PEP q-value is FdrInfo's initial 2, as FdrAnalysisEngine leaves it.
        /// </summary>
        private static List<SpectralMatch> Matches(string sequence, string protease, int count, int confident, bool pepTrained,
            CleavageSpecificity specificity = CleavageSpecificity.Full)
        {
            var peptide = new Protein(sequence, "acc_" + sequence)
                .Digest(new DigestionParams(protease: protease, minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .First();
            peptide.CleavageSpecificityForFdrCategory = specificity;
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());

            // a match's digestion agent comes from the CommonParameters it is made with
            var searchParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease, minPeptideLength: 1));
            var psms = new List<SpectralMatch>();
            for (int i = 0; i < count; i++)
            {
                var psm = new PeptideSpectralMatch(peptide, 0, 10, i, scan, searchParameters, new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                double q = i < confident ? 0.001 : 0.5;
                psm.SetFdrValues(i + 1, 0, q, i + 1, 0, q, pep: 0, pepQValue: pepTrained ? q : 2);
                psm.PsmFdrInfo.PepTrained = pepTrained;
                psms.Add(psm);
            }
            return psms;
        }

        private static SpectralMatch NotchAmbiguousPsm(double psmHypothesisQValue, double peptideHypothesisQValue)
        {
            var psm = Matches("PEPTIDEK", "trypsin", 1, 0, pepTrained: false)[0];
            var alternate = new Protein("PEPTIDER", "acc_alternate")
                .Digest(new DigestionParams(protease: "trypsin", minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .First();
            psm.AddOrReplace(alternate, psm.Score, 1, reportAllAmbiguity: true, matchedFragmentIons: new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            Set(psm, q: 0.5, notch: 2, pepQ: 2, pepTrained: false);
            psm.PeptideFdrInfo = new FdrInfo { QValue = 0.5, QValueNotch = 2, PEP_QValue = 2 };

            foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
            {
                hypothesis.QValueNotch = psmHypothesisQValue;
                hypothesis.PeptideQValueNotch = peptideHypothesisQValue;
            }

            return psm;
        }

        private static void Set(SpectralMatch psm, double q, double notch, double pepQ, bool pepTrained = true, double pep = 0)
        {
            psm.SetFdrValues(1, 0, q, 1, 0, notch, pep, pepQ);
            psm.PsmFdrInfo.PepTrained = pepTrained;
        }

        private static int CountFor(FilteredPsms filtered, string protease) =>
            filtered.FilteredPsmsList.Count(p => p.DigestionParams.DigestionAgent.Name == protease);

        // Spec 1: the user picks the filter, and holds it to the one threshold.
        [Test]
        public static void TheSearchNamesItsFilter_AndHoldsItToTheOneThreshold()
        {
            var trained = Matches("PEPTIDEK", "trypsin", 300, 150, pepTrained: true);

            var pep = FilteredPsms.Filter(trained, Chosen(FilterType.PepQValue));
            Assert.That((pep.FilterType, pep.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));

            var notch = FilteredPsms.Filter(trained, Chosen(FilterType.QValueNotch, 1.0));
            Assert.That((notch.FilterType, notch.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 1.0)));

            var pooled = FilteredPsms.Filter(trained, Chosen(FilterType.QValue, 0.05));
            Assert.That((pooled.FilterType, pooled.FilterThreshold), Is.EqualTo((FilterType.QValue, 0.05)));

            var combined = FilteredPsms.Filter(trained, Chosen(FilterType.QValueAndQValueNotch, 0.05));
            Assert.That((combined.FilterType, combined.FilterThreshold), Is.EqualTo((FilterType.QValueAndQValueNotch, 0.05)));
        }

        // The maintainer's request: unnamed, the filter is the PEP q-value bottom-up and the q-value notch for top-down
        // (and RNA), each at 0.01.
        [Test]
        [TestCase(false, false, FilterType.PepQValue)]
        [TestCase(false, true, FilterType.QValueNotch)]
        [TestCase(true, false, FilterType.QValueNotch)]
        [TestCase(true, true, FilterType.QValueNotch)]
        public static void DefaultIdentificationFilterBySearchType(bool isRna, bool isTopDown, FilterType expected)
        {
            Assert.That(CommonParameters.DefaultIdentificationFilter(isRna, isTopDown), Is.EqualTo(expected));
        }

        [Test]
        public static void TheDefaultFilter_IsThePepQValueBottomUp_AndTheNotchForTopDownAndRna()
        {
            Assert.That((new CommonParameters().IdentificationFilter, new CommonParameters().FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
            var topDown = new CommonParameters(digestionParams: new DigestionParams(protease: "top-down"));
            Assert.That((topDown.IdentificationFilter, topDown.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.01)));
            var rna = new CommonParameters(digestionParams: new global::Transcriptomics.Digestion.RnaDigestionParams("RNase T1"));
            Assert.That((rna.IdentificationFilter, rna.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.01)));

            // A named filter always wins, whatever the analyte
            Assert.That(new CommonParameters(digestionParams: new DigestionParams(protease: "top-down"), identificationFilter: FilterType.PepQValue).IdentificationFilter,
                Is.EqualTo(FilterType.PepQValue));
            Assert.That(new CommonParameters(identificationFilter: FilterType.QValueNotch).IdentificationFilter, Is.EqualTo(FilterType.QValueNotch));

            // The bottom-up default falls back to the notch where PEP was not trained
            var psms = Matches("PEPTIDEK", "trypsin", 300, 0, pepTrained: false);
            foreach (var psm in psms.Take(120)) Set(psm, q: 0.5, notch: 0.001, pepQ: 2, pepTrained: false);
            var filtered = FilteredPsms.Filter(psms, new CommonParameters());
            Assert.That((filtered.FilterType, filtered.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(120), "untrained: the notch alone decides");

            // The notch is the notch alone, not the q-value and the notch together
            Assert.That(FilteredPsms.Filter(psms, Chosen(FilterType.QValueNotch)).FilteredPsmsList, Has.Count.EqualTo(120));
        }

        // Spec 1: the filter and its threshold are saved in the TOML, and read back; the old thresholds are not written.
        [Test]
        [TestCase(FilterType.PepQValue, 0.01)]
        [TestCase(FilterType.QValueNotch, 0.05)]
        [TestCase(FilterType.QValue, 0.005)]
        public static void TheChoiceIsSavedInTheToml_WithItsOneThreshold(FilterType filter, double threshold)
        {
            string saved = TomlOf(Chosen(filter, threshold));
            string[] lines = saved.Split(Environment.NewLine);
            Assert.That(lines, Has.One.Match($@"^\s*IdentificationFilter = ""{filter}""$"));
            Assert.That(lines, Has.One.Match($@"^\s*FilterThreshold = {threshold.ToString(System.Globalization.CultureInfo.InvariantCulture)}$"));
            Assert.That(lines, Has.None.Match(@"^\s*(Pep)?QValueThreshold ="), "the old thresholds are not written");

            var reloadedTask = LoadTask(saved);
            var reloaded = reloadedTask.CommonParameters;
            Assert.That((reloaded.IdentificationFilter, reloaded.FilterThreshold), Is.EqualTo((filter, threshold)));
        }

        // Spec 2: every filter is held to the one threshold, and a PEP fallback to the notch at that threshold.
        [Test]
        public static void EachFilterIsHeldToTheOneThreshold_AndAPepFallbackUsesIt()
        {
            var psm = Matches("PEPTIDEK", "trypsin", 1, 0, pepTrained: true)[0];

            Set(psm, q: 0.5, notch: 0.03, pepQ: 0.5);
            Assert.That(FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.QValueNotch, 0.05)).FilteredPsmsList,
                Has.Count.EqualTo(1), "the notch at the threshold");

            Set(psm, q: 0.5, notch: 0.5, pepQ: 0.03);
            Assert.That(FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.PepQValue, 0.05)).FilteredPsmsList,
                Has.Count.EqualTo(1), "the PEP q-value at the threshold");

            Set(psm, q: 0.5, notch: 0.03, pepQ: 2, pepTrained: false);
            var fellBack = FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.PepQValue, 0.05));
            Assert.That(fellBack.FilteredPsmsList, Has.Count.EqualTo(1), "untrained: the notch, at the threshold");
            Assert.That(fellBack.GetFilterValue(psm), Is.EqualTo(0.03), "the value handed to quantification is the one it was judged by");

            Set(psm, q: 0.001, notch: 0.07, pepQ: 2, pepTrained: false);
            Assert.That(FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.PepQValue, 0.05)).FilteredPsmsList,
                Is.Empty, "the fallback is the notch alone, not the q-value");
        }

        // The pooled q-value is a choice of its own: the pooled q-value alone, not the notch and not PEP.
        [Test]
        public static void ThePooledQValue_JudgesByThePooledQValueAlone()
        {
            var psm = Matches("PEPTIDEK", "trypsin", 1, 0, pepTrained: true)[0];

            Set(psm, q: 0.03, notch: 0.5, pepQ: 0.5);
            var pooled = FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.QValue, 0.05));
            Assert.That(pooled.FilteredPsmsList, Has.Count.EqualTo(1), "the pooled q-value at the threshold");
            Assert.That(pooled.GetFilterValue(psm), Is.EqualTo(0.03));
            Assert.That(pooled.GetFilterTypeString(), Is.EqualTo("q-value"));
            Assert.That(pooled.PepFilteringNotDone, Is.Null, "no PEP, so no PEP fallback");

            Set(psm, q: 0.5, notch: 0.001, pepQ: 0.001);
            Assert.That(FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.QValue, 0.05)).FilteredPsmsList, Is.Empty,
                "not the notch, not the PEP q-value");

            Set(psm, q: 0.001, notch: 0.5, pepQ: 2, pepTrained: false);
            Assert.That(FilteredPsms.Filter(new[] { psm }, Chosen(FilterType.QValue, 0.05)).FilteredPsmsList, Has.Count.EqualTo(1),
                "PEP training does not matter to it");
        }

        [Test]
        public static void CombinedQValueAndQValueNotchFilterRequiresBothValuesAtTheThreshold()
        {
            var psms = Matches("PEPTIDEK", "trypsin", 3, 0, pepTrained: true);
            Set(psms[0], q: 0.01, notch: 0.01, pepQ: 0.01);
            Set(psms[1], q: 0.011, notch: 0.001, pepQ: 0.01);
            Set(psms[2], q: 0.001, notch: 0.011, pepQ: 0.01);

            var filtered = FilteredPsms.Filter(psms, Chosen(FilterType.QValueNotch, 0.01),
                includeAmbiguousMods: false,
                filterType: FilterType.QValueAndQValueNotch);

            Assert.That(filtered.FilteredPsmsList, Is.EqualTo(new[] { psms[0] }));
            Assert.That(filtered.GetFilterValue(psms[0]), Is.EqualTo(0.01));
        }

        [Test]
        public static void NotchAmbiguousFiltering_DefaultMatchesThePreviousPsmPredicate_AndOptInKeepsPassingHypotheses()
        {
            var ambiguous = NotchAmbiguousPsm(psmHypothesisQValue: 0.001, peptideHypothesisQValue: 0.001);
            var ordinaryPass = Matches("ORDINARYK", "trypsin", 1, 0, pepTrained: false)[0];
            Set(ordinaryPass, q: 0.5, notch: 0.01, pepQ: 2, pepTrained: false);
            var ordinaryFail = Matches("FAILINGK", "trypsin", 1, 0, pepTrained: false)[0];
            Set(ordinaryFail, q: 0.5, notch: 0.2, pepQ: 2, pepTrained: false);
            var ambiguousFail = NotchAmbiguousPsm(psmHypothesisQValue: 0.2, peptideHypothesisQValue: 0.2);
            const double threshold = 0.05;
            var input = new[] { ordinaryPass, ordinaryFail, ambiguous, ambiguousFail };

            bool PreviousPredicate(SpectralMatch psm) => psm.PsmFdrInfo.QValueNotch <= threshold
                || (Math.Abs(psm.PsmFdrInfo.QValueNotch - 2) < 0.001
                    && psm.BestMatchingBioPolymersWithSetMods.Any(h => h.QValueNotch.HasValue && h.QValueNotch <= threshold));

            var defaultFiltered = FilteredPsms.Filter(input, Chosen(FilterType.QValueNotch, threshold),
                includeAmbiguous: true, includeAmbiguousMods: true);
            var optInFiltered = FilteredPsms.Filter(input, Chosen(FilterType.QValueNotch, threshold),
                includeAmbiguous: true, includeAmbiguousMods: true, includeNotchAmbiguous: true);

            Assert.That(input.Where(PreviousPredicate), Is.EqualTo(new[] { ordinaryPass, ambiguous }));
            Assert.That(defaultFiltered.FilteredPsmsList, Is.EqualTo(new[] { ordinaryPass }), "the new option must remain opt-in");
            Assert.That(optInFiltered.FilteredPsmsList, Is.EqualTo(new[] { ordinaryPass, ambiguous }));
            Assert.That(ambiguous.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(2),
                "filtering must not perform notch disambiguation or mutate hypotheses");
        }

        [Test]
        public static void NotchAmbiguousFiltering_DropsWhenEveryHypothesisFails()
        {
            var ambiguous = NotchAmbiguousPsm(psmHypothesisQValue: 0.2, peptideHypothesisQValue: 0.2);

            var filtered = FilteredPsms.Filter(new[] { ambiguous }, Chosen(FilterType.QValueNotch, 0.05),
                includeAmbiguous: true, includeAmbiguousMods: true, includeNotchAmbiguous: true);

            Assert.That(filtered.FilteredPsmsList, Is.Empty);
        }

        [Test]
        public static void NotchAmbiguousFiltering_UsesPeptideLevelHypothesisNotchValues()
        {
            var ambiguous = NotchAmbiguousPsm(psmHypothesisQValue: 0.2, peptideHypothesisQValue: 0.001);

            var filtered = FilteredPsms.Filter(new[] { ambiguous }, Chosen(FilterType.QValueNotch, 0.05),
                includeAmbiguous: true, includeAmbiguousMods: true, filterAtPeptideLevel: true, includeNotchAmbiguous: true);

            Assert.That(filtered.FilteredPsmsList, Is.EqualTo(new[] { ambiguous }));
        }

        // Spec 5 (README C): one protease trained PEP, the other did not.
        [Test]
        public static void PepNotTrainedForOneProtease_ItsMatchesAreJudgedByTheNotch_AndTheSentenceNamesIt()
        {
            var both = Matches("PEPTIDEK", "trypsin", 300, 150, pepTrained: true)
                .Concat(Matches("PEPTIDER", "Arg-C", 200, 100, pepTrained: false)).ToList();

            var filtered = FilteredPsms.Filter(both, Chosen(FilterType.PepQValue));

            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(CountFor(filtered, "trypsin"), Is.EqualTo(150));
            Assert.That(CountFor(filtered, "Arg-C"), Is.EqualTo(100), "judged by the notch, not dropped on a PEP q-value of 2");
            Assert.That(filtered.PepFilteringNotDone, Is.EqualTo(
                "PEP filtering was requested, but PEP was not trained for 200 target PSMs (Arg-C: 200), so those were filtered by q-value notch <= 0.01."));
        }

        // Spec 5: a semi- or non-specific search runs FDR, and trains PEP, once per cleavage-specificity category.
        [Test]
        public static void PepNotTrainedForOneCategory_IsReportedForThatCategory()
        {
            var both = Matches("PEPTIDEK", "trypsin", 300, 150, pepTrained: true, specificity: CleavageSpecificity.Full)
                .Concat(Matches("PEPTIDER", "trypsin", 200, 100, pepTrained: false, specificity: CleavageSpecificity.Semi)).ToList();

            var filtered = FilteredPsms.Filter(both, Chosen(FilterType.PepQValue));

            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(250));
            Assert.That(filtered.PepFilteringNotDone, Does.Contain("(trypsin (Semi): 200)"));
        }

        // Spec 5 (README E): fewer than 100 matches are filtered, not kept whole.
        [Test]
        public static void FewerThan100Matches_AreFiltered()
        {
            var filtered = FilteredPsms.Filter(Matches("PEPTIDER", "Arg-C", 50, 20, pepTrained: false), Chosen(FilterType.PepQValue));

            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(20));
            Assert.That(filtered.PepFilteringNotDone, Does.StartWith("PEP filtering was requested, but PEP was not trained for 50 target PSMs"));
        }

        // Spec 3 (README F): the Spectrum Match Count column follows the chosen filter.
        [Test]
        public static void SpectrumMatchCount_FollowsTheChosenFilter()
        {
            var psms = Matches("PEPTIDEK", "trypsin", 300, 0, pepTrained: true);
            for (int i = 0; i < 150; i++) Set(psms[i], q: i < 100 ? 0.001 : 0.5, notch: i < 100 ? 0.001 : 0.5, pepQ: 0.001);

            FdrAnalysisEngine.CountPsm(psms, Chosen(FilterType.PepQValue));

            Assert.That(psms.Select(p => p.PsmCount).Distinct(), Is.EqualTo(new[] { 150 }), "the 150 that pass the PEP q-value");
        }

        /// <summary>One protein, one peptide, one PSM at the given score, notch q-value and PEP, PEP trained.</summary>
        private static (ProteinGroup Group, SpectralMatch Psm) OneGroup(string accession, string sequence, bool decoy, double score, double notch, double pep)
        {
            var protein = new Protein(sequence, accession, isDecoy: decoy);
            var peptide = protein.Digest(new DigestionParams(protease: "trypsin", minPeptideLength: 1), new List<Modification>(), new List<Modification>()).First();
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());
            var psm = new PeptideSpectralMatch(peptide, 0, score, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            Set(psm, q: 0.001, notch: notch, pepQ: 0.001, pep: pep);
            var group = new ProteinGroup(new HashSet<IBioPolymer> { protein }, new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });
            return (group, psm);
        }

        private static Dictionary<string, double> ProteinQValues(CommonParameters commonParameters)
        {
            // X: best score 30, notch 0.0001, PEP 0.5. Y: best score 50, notch 0.005, PEP 0.0001. A decoy between them on
            // both. Ranked by the notch: X, decoy, Y. Ranked by PEP: Y, decoy, X.
            var all = new[]
            {
                OneGroup("X", "ACDEFK", decoy: false, score: 30, notch: 0.0001, pep: 0.5),
                OneGroup("Y", "GHILMK", decoy: false, score: 50, notch: 0.005, pep: 0.0001),
                OneGroup("DECOY_Z", "NQSTVK", decoy: true, score: 40, notch: 0.003, pep: 0.001),
            };
            var results = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(all.Select(a => a.Group).ToList(), all.Select(a => a.Psm).ToList(),
                noOneHitWonders: false, treatModPeptidesAsDifferentPeptides: true, mergeIndistinguishableProteinGroups: false, commonParameters,
                new List<(string, CommonParameters)>(), new List<string>()).Run();
            return results.SortedAndScoredProteinGroups.ToDictionary(g => g.ProteinGroupName, g => g.QValue);
        }

        // Spec 3: protein scoring and protein FDR use the chosen filter and stop inferring their own.
        [Test]
        public static void ProteinFdr_UsesTheChosenFilter()
        {
            var byPep = ProteinQValues(Chosen(FilterType.PepQValue));
            Assert.That((byPep["Y"], byPep["X"]), Is.EqualTo((0.0, 0.5)), "ranked by PEP");

            var byNotch = ProteinQValues(Chosen(FilterType.QValueNotch));
            Assert.That((byNotch["X"], byNotch["Y"]), Is.EqualTo((0.0, 0.5)), "ranked by the notch");

            // The pooled q-value ranks protein groups by their best q-value notch, as every filter but a trained PEP does
            var byPooled = ProteinQValues(Chosen(FilterType.QValue));
            Assert.That((byPooled["X"], byPooled["Y"]), Is.EqualTo((0.0, 0.5)), "ranked by the notch");
        }

        /// <summary>A match to the first peptide of <paramref name="protein"/>, on its own scan.</summary>
        private static SpectralMatch Match(Protein protein, int notch, double score, int scanNumber)
        {
            var peptide = protein.Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>()).First();
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                scanNumber, 2, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=" + scanNumber, double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());
            var psm = new PeptideSpectralMatch(peptide, notch, score, scanNumber, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            return psm;
        }

        // Spec 5b: the engines' counts in results.txt follow the chosen filter at its threshold and name it.
        [Test]
        public static void TheEngineCounts_FollowTheChosenFilter_AndNameIt()
        {
            // Through FdrAnalysisEngine.Run, which computes the q-values itself: 10 decoys in notch 0 outscore 100 targets
            // in notch 1, so every target fails the pooled q-value and passes the notch.
            var target = new Protein("PEPTIDEK", "target", databaseFilePath: "temp");
            var decoy = new Protein("EDITPEPK", "DECOY_target", isDecoy: true, databaseFilePath: "temp");
            var psms = Enumerable.Range(1, 10).Select(i => Match(decoy, notch: 0, score: 100, scanNumber: i))
                .Concat(Enumerable.Range(11, 100).Select(i => Match(target, notch: 1, score: 50, scanNumber: i)))
                .ToList();
            var commonParameters = Chosen(FilterType.QValueNotch);
            var fdr = (FdrAnalysisResults)new FdrAnalysisEngine(psms, 2, commonParameters,
                new List<(string, CommonParameters)> { ("File", commonParameters) }, new List<string>(), doPEP: false).Run();
            Assert.That(psms.Where(p => !p.IsDecoy).All(p => p.FdrInfo.QValue > 0.01 && p.FdrInfo.QValueNotch <= 0.01));
            Assert.That(fdr.PsmsWithin1PercentFdr, Is.EqualTo(100));
            Assert.That(fdr.ToString(), Does.Contain("PSMs with q-value notch <= 0.01: 100"));

            // The same matches under the pooled q-value: every target fails it
            var pooledParameters = Chosen(FilterType.QValue);
            var pooled = (FdrAnalysisResults)new FdrAnalysisEngine(psms, 2, pooledParameters,
                new List<(string, CommonParameters)> { ("File", pooledParameters) }, new List<string>(), doPEP: false).Run();
            Assert.That(pooled.PsmsWithin1PercentFdr, Is.EqualTo(0));
            Assert.That(pooled.ToString(), Does.Contain("PSMs with q-value <= 0.01: 0"));

            // The localized-modification counts.
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif);
            var mod = new Modification(_originalId: "mod1", _modificationType: "mt", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
            var protein = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: new Dictionary<int, List<Modification>> { { 2, new List<Modification> { mod } } });
            var peptide = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0,
                new Dictionary<int, Modification> { { 2, mod } }, 0);
            var scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), 0, 0, "", new CommonParameters());
            var modPsm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            modPsm.ResolveAllAmbiguities();
            Set(modPsm, q: 0.5, notch: 0.001, pepQ: 2, pepTrained: false);
            var mods = (ModificationAnalysisResults)new ModificationAnalysisEngine(new List<SpectralMatch> { modPsm }, commonParameters,
                new List<(string, CommonParameters)> { ("", commonParameters) }, new List<string>()).Run();
            Assert.That(mods.CountOfEachModSeenOnProteins.Values.Sum() + mods.CountOfModsSeenAndLocalized.Values.Sum(), Is.GreaterThan(0));
            Assert.That(mods.ToString(), Does.Contain("Localized mods seen with q-value notch <= 0.01:"));

            // The protein groups, at the chosen threshold.
            var groups = ProteinQValues(Chosen(FilterType.QValueNotch, 0.05));
            Assert.That(groups["X"], Is.EqualTo(0));
            var proteinResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(new List<ProteinGroup>(), new List<SpectralMatch>(),
                false, true, false, Chosen(FilterType.QValueNotch, 0.05), new List<(string, CommonParameters)>(), new List<string>()).Run();
            Assert.That(proteinResults.ToString(), Does.Contain("Number of proteins with q-value <= 0.05: 0"));
        }

        // Spec 5b: mzIdentML's pass flags follow the chosen filter, for matches and protein groups.
        [Test]
        public static void MzIdentMl_PassThreshold_FollowsTheChosenFilter()
        {
            var protein = new Protein("PEPTIDEK", "ACC1", databaseFilePath: "temp");
            var notchOnly = Match(protein, notch: 0, score: 10, scanNumber: 1);
            Set(notchOnly, q: 0.5, notch: 0.001, pepQ: 2, pepTrained: false);
            var qValueOnly = Match(protein, notch: 0, score: 10, scanNumber: 2);
            Set(qValueOnly, q: 0.001, notch: 0.5, pepQ: 2, pepTrained: false);
            var psms = new List<SpectralMatch> { notchOnly, qValueOnly };
            var group = new ProteinGroup(new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods>(psms.SelectMany(p => p.BestMatchingBioPolymersWithSetMods).Select(p => p.SpecificBioPolymer)),
                new HashSet<IBioPolymerWithSetMods>())
            {
                QValue = 0.03
            };
            group.SequenceCoverageFraction.Add(0.5);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "IdentificationFilterTests_passThreshold.mzID");
            try
            {
                var filter = FilteredPsms.Filter(psms, Chosen(FilterType.QValueNotch, 0.05), includeHighQValuePsms: true);
                MzIdentMLWriter.WriteMzIdentMl(psms, new List<ProteinGroup> { group }, new List<Modification>(), new List<Modification>(),
                    null, new List<DigestionAgent>(), new MzLibUtil.PpmTolerance(20), new MzLibUtil.PpmTolerance(20), 0, path, true,
                    identificationFilter: filter);
                var elements = XDocument.Load(path).Descendants().ToList();
                string[] Of(string name) => elements.Where(e => e.Name.LocalName == name).Select(e => (string)e.Attribute("passThreshold")).ToArray();
                Assert.That(Of("SpectrumIdentificationItem"), Is.EqualTo(new[] { "true", "false" }), "the notch decides each match");
                Assert.That(Of("ProteinDetectionHypothesis"), Is.EqualTo(new[] { "true" }), "a group at 0.03 passes a 0.05 filter");
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        /// <summary>
        /// The FDR engine records, on each match and at each level, whether PEP was trained for it; a later pass that
        /// trains nothing keeps the record. The HeLa subset, searched as the PEP tests search it.
        /// </summary>
        [Test]
        public static void TheFdrEngineRecordsWhetherItTrainedPep_AndALaterPassKeepsTheRecord()
        {
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams());
            var fsp = new List<(string, CommonParameters)> { ("TaGe_SA_HeLa_04_subset_longestSeq.mzML", commonParameters) };
            string dataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            var msDataFile = new MyFileManager(true).LoadFile(dataFile, commonParameters);
            List<Protein> proteins = UsefulProteomicsDatabases.ProteinDbLoader.LoadProteinFasta(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"), true, UsefulProteomicsDatabases.DecoyType.Reverse, false, out _,
                UsefulProteomicsDatabases.ProteinDbLoader.UniprotAccessionRegex, UsefulProteomicsDatabases.ProteinDbLoader.UniprotFullNameRegex, UsefulProteomicsDatabases.ProteinDbLoader.UniprotFullNameRegex,
                UsefulProteomicsDatabases.ProteinDbLoader.UniprotGeneNameRegex, UsefulProteomicsDatabases.ProteinDbLoader.UniprotOrganismRegex, -1);
            var scans = MetaMorpheusTask.GetMs2Scans(msDataFile, dataFile, commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] psmArray = new PeptideSpectralMatch[scans.Length];
            new EngineLayer.ClassicSearch.ClassicSearchEngine(psmArray, scans, new List<Modification>(), new List<Modification>(), null, null, null,
                proteins, new SinglePpmAroundZeroSearchMode(5), commonParameters, fsp, null, new List<string>(), false).Run();
            var psms = psmArray.Where(p => p != null).ToList();

            var qValueOverride = typeof(FdrAnalysisEngine).GetProperty("QvalueThresholdOverride")!;
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                try
                {
                    qValueOverride.SetValue(null, false);
                    new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>()).Run();
                    Assert.That(psms.Count, Is.LessThanOrEqualTo(1000));
                    Assert.That(psms.Select(p => (p.PsmFdrInfo.PepTrained, p.PeptideFdrInfo.PepTrained)).Distinct(), Is.EqualTo(new[] { (false, false) }),
                        "1000 matches or fewer: PEP is not trained");

                    // The override's (decoys + 1) / targets q-values leave no positive example below the default 0.005.
                    qValueOverride.SetValue(null, true);
                    commonParameters.QValueCutoffForPepCalculation = 0.1;
                    new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>()).Run();
                    Assert.That(psms.Select(p => p.PsmFdrInfo.PEP).Distinct().Count(), Is.GreaterThan(1), "a trained model, not PEPs left at 0");
                    Assert.That(psms.Select(p => (p.PsmFdrInfo.PepTrained, p.PeptideFdrInfo.PepTrained)).Distinct(), Is.EqualTo(new[] { (true, true) }));

                    // A later pass that reuses the PEPs (the individual-file and variant tables) trains nothing and keeps the record.
                    new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>(), doPEP: false).Run();
                    Assert.That(psms.All(p => p.PsmFdrInfo.PepTrained && p.PeptideFdrInfo.PepTrained));
                }
                finally
                {
                    qValueOverride.SetValue(null, false);
                }
            }
        }

        /// <summary>
        /// A whole SearchTask on the two TaGe snips: 258 target PSMs and no decoys. PEP is not trained; with
        /// <paramref name="trainPep"/> training is attempted and fails.
        /// </summary>
        private static (string Folder, string[] Results, List<string> Warnings) Search(string name, CommonParameters commonParameters, bool trainPep,
            bool writeIndividualFiles = false, bool writeIdentificationFiles = false) =>
            Search(name, new SearchTask { CommonParameters = commonParameters }, trainPep, writeIndividualFiles, writeIdentificationFiles);

        private static (string Folder, string[] Results, List<string> Warnings) Search(string name, SearchTask task, bool trainPep,
            bool writeIndividualFiles = false, bool writeIdentificationFiles = false)
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;
            string outputFolder = Path.Combine(testDirectory, "IdentificationFilterTests_" + name);
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            task.SearchParameters = new SearchParameters
            {
                WriteIndividualFiles = writeIndividualFiles,
                WriteMzId = writeIdentificationFiles,
                WritePepXml = writeIdentificationFiles
            };
            var spectraFiles = new List<string>
            {
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML"),
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip_2.mzML")
            };
            var databases = new List<DbForTask> { new DbForTask(Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta"), false) };
            var warnings = new List<string>();
            EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);

            var qValueOverride = typeof(FdrAnalysisEngine).GetProperty("QvalueThresholdOverride")!;
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                MetaMorpheusTask.WarnHandler += onWarn;
                try
                {
                    qValueOverride.SetValue(null, trainPep);
                    Directory.CreateDirectory(outputFolder);
                    task.RunTask(outputFolder, databases, spectraFiles, name);
                }
                finally
                {
                    qValueOverride.SetValue(null, false);
                    MetaMorpheusTask.WarnHandler -= onWarn;
                }
            }
            return (outputFolder, File.ReadAllLines(Path.Combine(outputFolder, "results.txt")), warnings);
        }

        // Spec 5 (README H, J): a failed training is not a trained PEP.
        [Test]
        public static void SearchTask_PepTrainingFailed_FallsBackToTheNotch_AndSaysSo()
        {
            var (folder, results, warnings) = Search("trainingFailed", Chosen(FilterType.PepQValue), trainPep: true);
            foreach (string line in results.Where(l => l.StartsWith("PEP filtering") || l.StartsWith("All target") || l.StartsWith("Posterior")))
                TestContext.WriteLine("results.txt: " + (line.Length > 200 ? line.Substring(0, 200) + "..." : line));

            Assert.That(results, Has.Some.StartsWith("Posterior error probability analysis failed"));
            int psmLine = Array.FindIndex(results, l => l.StartsWith("All target PSMs with "));
            Assert.That(results[psmLine], Does.StartWith("All target PSMs with pep q-value <= 0.01: "));
            Assert.That(results[psmLine - 1], Is.EqualTo(
                "PEP filtering was requested, but PEP was not trained for 258 target PSMs (trypsin: 258), so those were filtered by q-value notch <= 0.01."));
            Assert.That(warnings, Has.One.EqualTo(results[psmLine - 1]), "one warning per search");

            Directory.Delete(folder, true);
        }

        // Spec 4 (README I): the per-file lines follow the same rule as the all-files lines.
        [Test]
        public static void SearchTask_PerFileLines_FollowTheSameRule_AndCarryTheSentence()
        {
            var (folder, results, warnings) = Search("perFile", Chosen(FilterType.PepQValue), trainPep: false, writeIndividualFiles: true);
            foreach (string line in results.Where(l => l.StartsWith("PEP filtering") || l.Contains(" - Target ") || l.StartsWith("All target")))
                TestContext.WriteLine("results.txt: " + line);

            var countLines = Enumerable.Range(0, results.Length)
                .Where(i => results[i].Contains(" - Target PSMs with ") || results[i].Contains(" - Target peptides with ")
                    || results[i].StartsWith("All target PSMs with ") || results[i].StartsWith("All target peptides with "))
                .ToList();
            Assert.That(countLines, Has.Count.EqualTo(6), "all files and two files, PSMs and peptides");
            foreach (int i in countLines)
            {
                Assert.That(results[i], Does.Contain("pep q-value <= 0.01: "), "the same filter and threshold on every line");
                Assert.That(results[i - 1], Does.StartWith("PEP filtering was requested, but "), "the fallback sentence above every count it qualifies");
            }
            Assert.That(warnings.Count(w => w.StartsWith("PEP filtering was requested")), Is.EqualTo(1), "one warning per search");

            Directory.Delete(folder, true);
        }

        // Spec 5b (README I): mzIdentML and pepXML state the chosen filter and threshold.
        [Test]
        public static void SearchTask_TheChosenFilter_ReachesMzIdentMlAndPepXml()
        {
            // At 0.05, so a file that still said the hard-coded 0.01 would fail.
            var (folder, _, _) = Search("writers", Chosen(FilterType.QValueNotch, 0.05),
                trainPep: false, writeIdentificationFiles: true);
            string individual = Path.Combine(folder, "Individual File Results");

            XDocument mzId = XDocument.Load(Path.Combine(individual, "TaGe_SA_A549_3_snip.mzID"));
            var thresholds = mzId.Descendants().Where(e => e.Name.LocalName == "cvParam" && ((string)e.Attribute("name"))?.EndsWith("FDR threshold") == true).ToList();
            Assert.That(thresholds.Select(e => (string)e.Attribute("value")), Is.EqualTo(new[] { "0.05", "0.05" }));

            string pepXml = File.ReadAllText(Path.Combine(individual, "TaGe_SA_A549_3_snip.pep.XML"));
            Assert.That(pepXml, Does.Contain("name=\"Identification Filter\" value=\"q-value notch &lt;= 0.05\""));

            Directory.Delete(folder, true);
        }

        // Spec 8: a chosen filter at 1 is allowed, with a warning at the start and a line in results.txt.
        [Test]
        public static void SearchTask_AFilterAtOne_RunsAndWarnsThatItKeepsEveryMatch()
        {
            var (folder, results, warnings) = Search("atOne", Chosen(FilterType.QValueNotch, 1.0), trainPep: false);

            const string sentence = "This search filters by q-value notch <= 1, which keeps every match.";
            Assert.That(warnings, Has.One.EqualTo(sentence));
            Assert.That(results, Has.One.EqualTo(sentence));

            Directory.Delete(folder, true);
        }

        // The pooled q-value is a choice: the constructor takes it, and a task TOML that names it loads and searches with it.
        [Test]
        public static void ThePooledQValue_IsAChoice_AndATomlThatNamesItLoads()
        {
            Assert.That(new CommonParameters(identificationFilter: FilterType.QValue).IdentificationFilter, Is.EqualTo(FilterType.QValue));

            string saved = TomlOf(Chosen(FilterType.QValueNotch, 0.02));
            var reloaded = LoadCommon(string.Join(Environment.NewLine, saved.Split(Environment.NewLine)
                .Select(l => l.TrimStart().StartsWith("IdentificationFilter") ? "IdentificationFilter = \"QValue\"" : l)));
            Assert.That((reloaded.IdentificationFilter, reloaded.FilterThreshold), Is.EqualTo((FilterType.QValue, 0.02)));
            Assert.That(FilteredPsms.Describe(reloaded), Is.EqualTo("q-value <= 0.02"));
        }

        private static string TestData(string file) => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", file);

        private static SearchTask LoadTask(string tomlText)
        {
            string toml = Path.Combine(TestContext.CurrentContext.TestDirectory, $"IdentificationFilterTests_{Guid.NewGuid():N}.toml");
            try
            {
                File.WriteAllText(toml, tomlText);
                return Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig);
            }
            finally
            {
                File.Delete(toml);
            }
        }

        private static CommonParameters LoadCommon(string tomlText) => LoadTask(tomlText).CommonParameters;

        private static string TomlOf(CommonParameters commonParameters)
        {
            string toml = Path.Combine(TestContext.CurrentContext.TestDirectory, $"IdentificationFilterTests_{Guid.NewGuid():N}.toml");
            try
            {
                Nett.Toml.WriteFile(new SearchTask { CommonParameters = commonParameters }, toml, MetaMorpheusTask.tomlConfig);
                return string.Join(Environment.NewLine, File.ReadAllLines(toml));
            }
            finally
            {
                File.Delete(toml);
            }
        }

        /// <summary>
        /// A search task TOML as an older MetaMorpheus wrote it: today's file for <paramref name="digestionParams"/> without
        /// IdentificationFilter and FilterThreshold, and with <paramref name="oldKeys"/> (whole lines, such as
        /// "QValueThreshold = 1.0") in the [CommonParameters] table, before its first key when <paramref name="first"/> and
        /// after its last otherwise. A filter is named only when <paramref name="filter"/> is given.
        /// </summary>
        private static string OldToml(string[] oldKeys, string filter = null, bool first = false, IDigestionParams digestionParams = null)
        {
            var lines = TomlOf(new CommonParameters(digestionParams: digestionParams)).Split(Environment.NewLine)
                .Where(l => !l.TrimStart().StartsWith("IdentificationFilter") && !l.TrimStart().StartsWith("FilterThreshold"))
                .ToList();
            var keys = oldKeys.ToList();
            if (filter != null)
            {
                keys.Insert(first ? keys.Count : 0, $"IdentificationFilter = \"{filter}\"");
            }

            int table = lines.FindIndex(l => l.Trim() == "[CommonParameters]");
            Assert.That(table, Is.GreaterThanOrEqualTo(0), "the fixture writes a [CommonParameters] table");
            // A table's keys end at the next table header
            int end = lines.FindIndex(table + 1, l => l.TrimStart().StartsWith("["));
            int at = first ? table + 1 : end < 0 ? lines.Count : end;
            while (!first && at > table + 1 && string.IsNullOrWhiteSpace(lines[at - 1])) at--;
            lines.InsertRange(at, keys);
            return string.Join(Environment.NewLine, lines);
        }

        private static readonly IDigestionParams TopDown = new DigestionParams(protease: "top-down");

        // Review of #2982 (nbollis), loading old TOMLs. A task saved with the old GUI's PEP box ticked names no filter and holds
        // QValueThreshold 1.0 and PepQValueThreshold 0.01 (TestData/Task2-SearchTaskconfig.toml, as on master). It loads as the PEP
        // q-value filter at 0.01, and is not reported as a filter at 1.
        [Test]
        public static void AnOldGuiPepTask_LoadsAsThePepFilter()
        {
            var commonTask = Toml.ReadFile<SearchTask>(TestData("Task2-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            var common = commonTask.CommonParameters;

            Assert.That((common.IdentificationFilter, common.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
            Assert.That(SearchTask.FilterKeepsEveryMatch(common), Is.Null, "a migrated PEP task filters at 0.01, not at 1");
            Assert.That((common.Clone().IdentificationFilter, common.Clone().FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)), "a copy keeps the filter");
        }

        // Review of #2982 (Alexander-Sol), his reproduction: an unmodified old PEP TOML (master's Task2-SearchTaskconfig.toml,
        // byte for byte) admitted ten trained matches with PEP q-value and notch 0.5, as a notch filter at its saved
        // QValueThreshold of 1.0. Read as the PEP filter at its saved 0.01, it rejects all ten.
        [Test]
        public static void AnUnmodifiedOldPepToml_RejectsTenTrainedMatchesAtOneHalf()
        {
            var common = Toml.ReadFile<SearchTask>(TestData("Task2-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig).CommonParameters;
            var psms = Matches("PEPTIDEK", "trypsin", 10, 0, pepTrained: true);
            Assume.That(psms.Select(p => (p.FdrInfo.PEP_QValue, p.FdrInfo.QValueNotch, p.FdrInfo.PepTrained)).Distinct(), Is.EqualTo(new[] { (0.5, 0.5, true) }));

            var filtered = FilteredPsms.Filter(psms, common);

            Assert.That((filtered.FilterType, filtered.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
            Assert.That(filtered.FilteredPsmsList, Is.Empty, "none of the ten passes the saved 0.01");
            Assert.That(psms.Where(p => FilteredPsms.Passes(p, false, common.IdentificationFilter, common.FilterThreshold)), Is.Empty);
        }

        // Review of #2982 (Alexander-Sol): opening and saving an old PEP task in the Search task window keeps it a PEP task.
        // The window's load and save rule is TaskWindowIdentificationFilter: the filter choice and one threshold box.
        [Test]
        public static void TheSearchWindow_ShowsAnOldPepTaskAsThePepFilter_AndSavesItsFilterByName()
        {
            var old = Toml.ReadFile<SearchTask>(TestData("Task2-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig).CommonParameters;

            var (shown, thresholdText) = TaskWindowIdentificationFilter.Load(old, isRna: false, isTopDown: false);
            Assert.That((shown, thresholdText), Is.EqualTo((FilterType.PepQValue, "0.01")), "the PEP q-value, at the saved 0.01");

            var (filter, threshold) = TaskWindowIdentificationFilter.ForSave(shown, thresholdText);
            Assert.That((filter, threshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
            string saved = TomlOf(new CommonParameters(identificationFilter: filter, filterThreshold: threshold));
            Assert.That(saved.Split(Environment.NewLine), Has.One.Match(@"^\s*IdentificationFilter = ""PepQValue""$"));
            Assert.That(saved.Split(Environment.NewLine), Has.One.Match(@"^\s*FilterThreshold = 0.01$"));
            var reopenedTask = LoadTask(saved);
            var reopened = reopenedTask.CommonParameters;
            Assert.That(TaskWindowIdentificationFilter.Load(reopened, false, false), Is.EqualTo((FilterType.PepQValue, "0.01")));

            // A new task shows the new bottom-up default, and a top-down or RNA one the q-value notch
            var fresh = new SearchTask().CommonParameters;
            Assert.That(TaskWindowIdentificationFilter.Load(fresh, isRna: false, isTopDown: false), Is.EqualTo((FilterType.PepQValue, "0.01")));
            Assert.That(TaskWindowIdentificationFilter.Load(fresh, isRna: false, isTopDown: true), Is.EqualTo((FilterType.QValueNotch, "0.01")));
            Assert.That(TaskWindowIdentificationFilter.Load(fresh, isRna: true, isTopDown: false), Is.EqualTo((FilterType.QValueNotch, "0.01")));
        }

        // Row 5: an old default task (QValueThreshold 0.01, PepQValueThreshold 1.0, no filter named) keeps the notch at 0.01,
        // not the new PEP q-value default, wherever the old keys sit in its table.
        [Test]
        [TestCase(false)]
        [TestCase(true)]
        public static void AnOldDefaultTask_LoadsAsTheNotch(bool keysFirst)
        {
            var commonTask = LoadTask(OldToml(new[] { "QValueThreshold = 0.01", "PepQValueThreshold = 1.0" }, first: keysFirst));
            var common = commonTask.CommonParameters;

            Assert.That((common.IdentificationFilter, common.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.01)));
        }

        // Row 4, with the old keys in either order and before or after the rest of the table
        [Test]
        [TestCase(false)]
        [TestCase(true)]
        public static void AnOldPepTask_ReadsAsThePepFilter_WhereverItsKeysSit(bool keysFirst)
        {
            var commonTask = LoadTask(OldToml(new[] { "PepQValueThreshold = 0.005", "QValueThreshold = 1.0" }, first: keysFirst));
            var common = commonTask.CommonParameters;

            Assert.That((common.IdentificationFilter, common.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.005)));
        }

        // Row 5 for a top-down task: the old keys keep it the notch at its QValueThreshold.
        [Test]
        public static void AnOldTopDownTask_WithoutAFilter_StaysTheNotchAtItsQValueThreshold()
        {
            var common = LoadCommon(OldToml(new[] { "QValueThreshold = 0.02", "PepQValueThreshold = 1.0" }, digestionParams: TopDown));

            Assume.That(common.DigestionParams.DigestionAgent.Name, Is.EqualTo("top-down"));
            Assert.That((common.IdentificationFilter, common.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.02)));
        }

        // A key the old file lacks means its old default: QValueThreshold 0.01, PepQValueThreshold 1.0.
        [Test]
        public static void AnOldTask_WithOneOldKey_ReadsTheOtherAsItsOldDefault()
        {
            var onlyQ = LoadCommon(OldToml(new[] { "QValueThreshold = 0.05" }));
            Assert.That((onlyQ.IdentificationFilter, onlyQ.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.05)));

            var onlyPepTask = LoadTask(OldToml(new[] { "PepQValueThreshold = 0.005" }));
            var onlyPep = onlyPepTask.CommonParameters;
            Assert.That((onlyPep.IdentificationFilter, onlyPep.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.005)), "0.005 is below the old default 0.01");
        }

        // The twelve test TOMLs that name both old keys are kept as they are, as old task files. Each keeps what it meant.
        [Test]
        [TestCase(@"DatabaseTests\Task1-SearchTaskconfig.toml", FilterType.QValueNotch, 0.01)]
        [TestCase(@"TestData\Task1-SearchTaskconfig.toml", FilterType.QValueNotch, 0.01)]
        [TestCase(@"TestData\Task2-SearchTaskconfig.toml", FilterType.PepQValue, 0.01)]
        [TestCase(@"TestData\ModificationAnalysis\SearchTaskconfig.toml", FilterType.QValueNotch, 0.01)]
        [TestCase(@"TestData\NonSpecificSearchToml.toml", FilterType.QValueNotch, 0.01)]
        [TestCase(@"TestData\RnaSearchTask_MultipleDecon.toml", FilterType.QValueNotch, 0.05)]
        [TestCase(@"TMT_test\TMT11_LowCID_SearchTask.toml", FilterType.QValueNotch, 0.01)]
        [TestCase(@"TopDownTestData\TopDownSearchToml.toml", FilterType.QValueNotch, 0.01)]
        [TestCase(@"GlycoTestData\TaskForOglycoPruneDb.toml", FilterType.QValueNotch, 0.01)]
        public static void TheOldSearchTaskFixtures_KeepTheirFilter(string file, FilterType filter, double threshold)
        {
            var common = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, file), MetaMorpheusTask.tomlConfig).CommonParameters;
            Assert.That((common.IdentificationFilter, common.FilterThreshold), Is.EqualTo((filter, threshold)));
        }

        [Test]
        [TestCase(@"GlycoTestData\GlycoSnip.toml")]
        [TestCase(@"GlycoTestData\Glyco_Isobaric_testing.toml")]
        public static void TheOldGlycoTaskFixtures_KeepTheNotch(string file)
        {
            var common = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, file), MetaMorpheusTask.tomlConfig).CommonParameters;
            Assert.That((common.IdentificationFilter, common.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.01)));
        }

        [Test]
        public static void TheOldGptmdTaskFixture_KeepsTheNotch()
        {
            var common = Toml.ReadFile<GptmdTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\ModificationAnalysis\GPTMDTaskconfig.toml"), MetaMorpheusTask.tomlConfig).CommonParameters;
            Assert.That((common.IdentificationFilter, common.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.01)));
        }

        [Test]
        public static void AnOldGuiPepTask_SavedAgain_NamesItsFilterAndThreshold()
        {
            var old = Toml.ReadFile<SearchTask>(TestData("Task2-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            string saved = TomlOf(old.CommonParameters);

            Assert.That(saved.Split(Environment.NewLine), Has.One.Match(@"^\s*IdentificationFilter = ""PepQValue""$"));
            Assert.That(saved.Split(Environment.NewLine), Has.One.Match(@"^\s*FilterThreshold = 0.01$"));
            Assert.That(saved.Split(Environment.NewLine), Has.None.Match(@"^\s*(Pep)?QValueThreshold ="));
            var reloadedTask = LoadTask(saved);
            var reloaded = reloadedTask.CommonParameters;
            Assert.That((reloaded.IdentificationFilter, reloaded.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
        }

        // A copy keeps the record of an old file, and resolves as the original does
        [Test]
        public static void ACopy_OfAnOldTask_ResolvesAsTheOriginal()
        {
            var old = LoadCommon(OldToml(new[] { "QValueThreshold = 0.02", "PepQValueThreshold = 1.0" }, digestionParams: TopDown));
            foreach (var copy in new[] { old.Clone(), old.CloneWithNewDissociationType(MassSpectrometry.DissociationType.ETD),
                old.CloneWithNewTerminus(), old.CloneWithNewTotalPartitions(2) })
            {
                Assert.That((copy.IdentificationFilter, copy.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.02)));
            }
            var file = MetaMorpheusTask.SetAllFileSpecificCommonParams(old, new FileSpecificParameters(), "file");
            Assert.That((file.IdentificationFilter, file.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.02)));
        }

        [Test]
        public static void SearchTask_AnOldGuiPepTask_UsesThePepFilter()
        {
            var commonTask = LoadTask(OldToml(new[] { "QValueThreshold = 1.0", "PepQValueThreshold = 0.01" }));

            var (folder, results, warnings) = Search("legacyToml", commonTask, trainPep: false);
            Assert.That(results, Has.Some.StartWith("All target PSMs with pep q-value <= 0.01: "));
            Assert.That(warnings, Has.None.Contain("keeps every match"));

            Directory.Delete(folder, true);
        }

        // Review of #2982 (nbollis), the peptide/PSM PepTrained mismatch. With more than 1000 matches but no more than 1000
        // peptides, the FDR engine trained PEP but computed no peptide-level PEP q-values, so the PSM level was stamped trained
        // and the peptide level not: under the PEP filter the PSM tables used PEP and the peptide tables the notch. Both levels
        // now follow the training.
        [Test]
        public static void PepTrainedWithFewPeptides_BothLevelsFollowTheTraining()
        {
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams());
            var fsp = new List<(string, CommonParameters)> { ("TaGe_SA_HeLa_04_subset_longestSeq.mzML", commonParameters) };
            string dataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            var msDataFile = new MyFileManager(true).LoadFile(dataFile, commonParameters);
            List<Protein> proteins = UsefulProteomicsDatabases.ProteinDbLoader.LoadProteinFasta(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"), true, UsefulProteomicsDatabases.DecoyType.Reverse, false, out _,
                UsefulProteomicsDatabases.ProteinDbLoader.UniprotAccessionRegex, UsefulProteomicsDatabases.ProteinDbLoader.UniprotFullNameRegex, UsefulProteomicsDatabases.ProteinDbLoader.UniprotFullNameRegex,
                UsefulProteomicsDatabases.ProteinDbLoader.UniprotGeneNameRegex, UsefulProteomicsDatabases.ProteinDbLoader.UniprotOrganismRegex, -1);
            var scans = MetaMorpheusTask.GetMs2Scans(msDataFile, dataFile, commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] psmArray = new PeptideSpectralMatch[scans.Length];
            new EngineLayer.ClassicSearch.ClassicSearchEngine(psmArray, scans, new List<Modification>(), new List<Modification>(), null, null, null,
                proteins, new SinglePpmAroundZeroSearchMode(5), commonParameters, fsp, null, new List<string>(), false).Run();
            var psms = psmArray.Where(p => p != null).ToList();
            psms.ForEach(p => p.ResolveAllAmbiguities());
            int peptides = psms.Select(p => p.FullSequence).Distinct().Count();
            Assume.That(psms.Count, Is.GreaterThan(peptides), "the fixture needs more matches than peptides");

            var threshold = typeof(FdrAnalysisEngine).GetProperty(nameof(FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue))!;
            int original = FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue;
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                try
                {
                    // More matches than the threshold, but no more peptides: the engine's second branch.
                    threshold.SetValue(null, peptides);
                    commonParameters.QValueCutoffForPepCalculation = 0.1;
                    new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>()).Run();
                }
                finally
                {
                    threshold.SetValue(null, original);
                }
            }

            Assert.That(psms.Select(p => p.PsmFdrInfo.PepTrained).Distinct(), Is.EqualTo(new[] { true }), "the fixture trains PEP");
            Assert.That(psms.Select(p => p.PeptideFdrInfo.PepTrained).Distinct(), Is.EqualTo(new[] { true }), "the peptide level follows the training");
            var peptideStandIns = FilteredPsms.Filter(psms, commonParameters, includeDecoys: true, includeContaminants: true, includeAmbiguous: true,
                includeHighQValuePsms: true, filterAtPeptideLevel: true).FilteredPsmsList;
            Assert.That(peptideStandIns.Select(p => p.PeptideFdrInfo.PEP_QValue), Has.All.LessThan(2), "peptide-level PEP q-values were computed");
        }

        // Review of #2982 (nbollis), SetFdrValues: replacing a match's q-values keeps the record of whether its PEP was trained.
        [Test]
        public static void SetFdrValues_KeepsTheRecordOfWhetherPepWasTrained()
        {
            var psm = Matches("PEPTIDEK", "trypsin", 1, 1, pepTrained: true)[0];
            Assume.That(psm.PsmFdrInfo.PepTrained, Is.True);

            psm.SetFdrValues(1, 0, 0.002, 1, 0, 0.003, 0.01, 0.004);

            Assert.That(psm.PsmFdrInfo.PepTrained, Is.True);
            Assert.That((psm.FdrInfo.QValue, psm.FdrInfo.QValueNotch, psm.FdrInfo.PEP_QValue), Is.EqualTo((0.002, 0.003, 0.004)));
        }

        // Review of #2982 (nbollis), rounding: results.txt states the threshold exactly; 0.005 is not 0.01.
        [Test]
        public static void SearchTask_ResultsTxt_StatesTheThresholdExactly()
        {
            var (folder, results, _) = Search("exactThreshold", Chosen(FilterType.QValueNotch, 0.005), trainPep: false);

            var countLines = results.Where(l => l.StartsWith("All target ") && l.Contains(" with q-value notch <= ")).ToList();
            Assert.That(countLines, Has.Count.EqualTo(3), "PSMs, peptides and protein groups");
            Assert.That(countLines, Has.All.Contain("q-value notch <= 0.005"));

            Directory.Delete(folder, true);
        }

        // The pooled q-value can be asked for by name, as the other two can.
        [Test]
        public static void ThePooledQValue_CanBeAskedForByName()
        {
            var psms = Matches("PEPTIDEK", "trypsin", 10, 4, pepTrained: true);

            var filtered = FilteredPsms.Filter(psms, new CommonParameters(), filterType: FilterType.QValue);
            Assert.That((filtered.FilterType, filtered.FilterThreshold, filtered.FilteredPsmsList.Count), Is.EqualTo((FilterType.QValue, 0.01, 4)));
            Assert.That(FilteredPsms.Passes(psms[0], false, FilterType.QValue, 0.01), Is.True);
            Assert.That(FilteredPsms.Passes(psms[9], false, FilterType.QValue, 0.01), Is.False);
        }
    }
}
