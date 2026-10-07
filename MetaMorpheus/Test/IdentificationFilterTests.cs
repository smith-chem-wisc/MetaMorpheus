using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// The tiered identification filter (<see cref="IdentificationFilter"/>): PEP q-value where PEP was trained, else the
    /// q-value notch, else the q-value; strict &lt;; opt-in; decided by what was computed, never by counts.
    /// </summary>
    [TestFixture]
    public static class IdentificationFilterTests
    {
        private const string Residues = "ACDEFGHILMNQSTVWY"; // no K, R or P, so each sequence digests to itself

        // The fallback reasons name "PSMs"/"peptides" from the global analyte type, which other tests change.
        private static AnalyteType _analyteTypeBefore;

        [SetUp]
        public static void UsePeptideLabels()
        {
            _analyteTypeBefore = GlobalVariables.AnalyteType;
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
        }

        [TearDown]
        public static void RestoreAnalyteType()
        {
            GlobalVariables.AnalyteType = _analyteTypeBefore;
        }

        /// <summary>A unique, undigestable sequence for index i (up to 289 * 17 of them).</summary>
        private static string Sequence(int i) =>
            "PEPT" + Residues[i % 17] + Residues[(i / 17) % 17] + Residues[(i / 289) % 17];

        private static SpectralMatch MakeMatch(int i, bool decoy, double q, double notch, double pep, double pepQ,
            double? peptideQ = null, double? peptideNotch = null, double? peptidePepQ = null)
        {
            var protein = new Protein(Sequence(i), (decoy ? "DECOY_" : "") + "ACC" + i, isDecoy: decoy);
            var peptide = protein.Digest(new DigestionParams(protease: "trypsin", minPeptideLength: 1),
                new List<Modification>(), new List<Modification>()).First();

            var fakeScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                i + 1, 2, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=" + (i + 1), double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            var scan = new Ms2ScanWithSpecificMass(fakeScan, 2, 0, "File", new CommonParameters());

            SpectralMatch psm = new PeptideSpectralMatch(peptide, 0, 10, i, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.PsmFdrInfo = new FdrInfo { QValue = q, QValueNotch = notch, PEP = pep, PEP_QValue = pepQ };
            psm.PeptideFdrInfo = new FdrInfo
            {
                QValue = peptideQ ?? q,
                QValueNotch = peptideNotch ?? notch,
                PEP = pep,
                PEP_QValue = peptidePepQ ?? pepQ
            };
            return psm;
        }

        /// <summary>n target matches whose values are all well below 0.01, with PEP varying (trained) or not.</summary>
        private static List<SpectralMatch> Matches(int n, bool pepTrained, double pepWhenUntrained = 0, double pepQWhenUntrained = 2)
        {
            return Enumerable.Range(0, n)
                .Select(i => pepTrained
                    ? MakeMatch(i, false, q: 0.001, notch: 0.001, pep: 0.0001 * (i + 1), pepQ: 0.001)
                    : MakeMatch(i, false, q: 0.001, notch: 0.001, pep: pepWhenUntrained, pepQ: pepQWhenUntrained))
                .ToList();
        }

        [Test]
        [TestCase(0.01, 1.0, false, TestName = "Defaults (q 0.01, PEP q 1.0): tiered mode off")]
        [TestCase(0.01, 0.05, false, TestName = "PEP q above q: tiered mode off")]
        [TestCase(0.01, 0.01, true, TestName = "PEP q equal to q: tiered mode on (equal no longer does nothing)")]
        [TestCase(0.05, 0.01, true, TestName = "PEP q below q: tiered mode on")]
        [TestCase(1.0, 0.01, true, TestName = "GUI PEP box (q 1.0, PEP q 0.01): tiered mode on")]
        [TestCase(2.0, 1.0, false, TestName = "PEP q of 1: tiered mode off")]
        public static void TieredModeIsOptIn(double q, double pepQ, bool expected)
        {
            Assert.That(IdentificationFilter.IsTieredMode(q, pepQ), Is.EqualTo(expected));
            Assert.That(IdentificationFilter.IsTieredMode(new CommonParameters(qValueThreshold: q, pepQValueThreshold: pepQ)), Is.EqualTo(expected));
        }

        [Test]
        public static void EveryTierUsesTheSmallerThreshold()
        {
            // The GUI sends q = 1.0 when its PEP box is ticked: a fallback held to 1.0 would filter nothing.
            Assert.That(IdentificationFilter.Threshold(1.0, 0.01), Is.EqualTo(0.01));
            var untrained = Matches(10, pepTrained: false);
            var tier = IdentificationFilter.Resolve(untrained, false, new CommonParameters(qValueThreshold: 1.0, pepQValueThreshold: 0.01));
            Assert.That(tier.FilterType, Is.EqualTo(FilterType.QValueNotch));
            Assert.That(tier.Threshold, Is.EqualTo(0.01));
        }

        [Test]
        [TestCase(false)]
        [TestCase(true)]
        public static void TrainedPepIsUsedAtEitherLevel(bool peptideLevel)
        {
            var tier = IdentificationFilter.Resolve(Matches(20, pepTrained: true), peptideLevel, 0.01);
            Assert.That(tier.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(tier.FallbackReason, Is.Null);
            Assert.That(tier.PeptideLevel, Is.EqualTo(peptideLevel));
            Assert.That(tier.Describe(), Is.EqualTo("pep q-value < 0.01"));
        }

        [Test]
        public static void UntrainedPepFallsBackToTheNotchAndSaysWhy()
        {
            // 640 PSMs: PEP is trained only above 1000, so every PEP q-value is still the initial 2.
            var tier = IdentificationFilter.Resolve(Matches(640, pepTrained: false), false, 0.01);
            Assert.That(tier.FilterType, Is.EqualTo(FilterType.QValueNotch));
            Assert.That(tier.FallbackReason, Is.EqualTo("PEP not trained: 640 PSMs"));
            Assert.That(tier.DescribeWithReason(), Is.EqualTo("q-value notch < 0.01 (PEP not trained: 640 PSMs)"));
        }

        [Test]
        public static void FailedPepTrainingFallsBackToTheNotch()
        {
            // PEPAnalysisEngine returned early: every PEP stayed 0, yet the FDR engine still computed PEP q-values
            // (a score ranking). Those must not be used as PEP q-values.
            var matches = Matches(50, pepTrained: false, pepWhenUntrained: 0, pepQWhenUntrained: 0.001);
            var tier = IdentificationFilter.Resolve(matches, false, 0.01);
            Assert.That(tier.FilterType, Is.EqualTo(FilterType.QValueNotch));
            Assert.That(tier.FallbackReason, Is.EqualTo("PEP training failed: every PEP is 0"));
        }

        [Test]
        public static void PepQValuesOnDecoysAloneDoNotCountAsTrained()
        {
            var matches = Matches(10, pepTrained: false);
            matches.Add(MakeMatch(100, true, q: 0.5, notch: 0.5, pep: 0.9, pepQ: 0.5));
            Assert.That(IdentificationFilter.Resolve(matches, false, 0.01).FilterType, Is.EqualTo(FilterType.QValueNotch));
        }

        [Test]
        public static void NoNotchFallsBackToTheQValueAlone()
        {
            var matches = Enumerable.Range(0, 10).Select(i => MakeMatch(i, false, q: 0.001, notch: 2, pep: 0, pepQ: 2)).ToList();
            var tier = IdentificationFilter.Resolve(matches, false, 0.01);
            Assert.That(tier.FilterType, Is.EqualTo(FilterType.QValue));
            Assert.That(tier.FallbackReason, Is.EqualTo("PEP not trained: 10 PSMs; no q-value notch was computed"));
            // q only: the notch (2) is not read
            Assert.That(matches.All(tier.Passes), Is.True);
        }

        [Test]
        public static void PeptideLevelFallsBackOnItsOwn()
        {
            // PEP trained on PSMs (more than 1000) but not on peptides (1000 or fewer): the peptide-level PEP q-value
            // stays 2. The PSM level keeps PEP; the peptide level falls back to the notch instead of losing everything.
            var matches = Enumerable.Range(0, 30)
                .Select(i => MakeMatch(i, false, q: 0.001, notch: 0.001, pep: 0.0001 * (i + 1), pepQ: 0.001, peptidePepQ: 2))
                .ToList();
            var psmTier = IdentificationFilter.Resolve(matches, false, 0.01);
            var peptideTier = IdentificationFilter.Resolve(matches, true, 0.01);
            Assert.That(psmTier.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(peptideTier.FilterType, Is.EqualTo(FilterType.QValueNotch));
            Assert.That(peptideTier.FallbackReason, Is.EqualTo("PEP not trained: 30 peptides"));
            Assert.That(matches.All(peptideTier.Passes), Is.True);
        }

        [Test]
        [TestCase(FilterType.PepQValue)]
        [TestCase(FilterType.QValueNotch)]
        [TestCase(FilterType.QValue)]
        public static void ExactlyAtTheThresholdFails(FilterType filterType)
        {
            var atThreshold = new FdrInfo { QValue = 0.01, QValueNotch = 0.01, PEP_QValue = 0.01 };
            var below = new FdrInfo { QValue = 0.0099, QValueNotch = 0.0099, PEP_QValue = 0.0099 };
            Assert.That(IdentificationFilter.Passes(atThreshold, filterType, 0.01), Is.False);
            Assert.That(IdentificationFilter.Passes(below, filterType, 0.01), Is.True);
        }

        [Test]
        public static void EachTierReadsOnlyItsOwnValue()
        {
            var fdr = new FdrInfo { QValue = 0.3, QValueNotch = 0.2, PEP_QValue = 0.1 };
            Assert.That(IdentificationFilter.GetValue(fdr, FilterType.PepQValue), Is.EqualTo(0.1));
            Assert.That(IdentificationFilter.GetValue(fdr, FilterType.QValueNotch), Is.EqualTo(0.2));
            Assert.That(IdentificationFilter.GetValue(fdr, FilterType.QValue), Is.EqualTo(0.3));
        }

        [Test]
        public static void FewerThan100MatchesAreFilteredNotPassedThrough()
        {
            // Before: PEP mode with fewer than 100 matches set the threshold to 1 and filtered nothing.
            var matches = Matches(50, pepTrained: false);
            matches.Add(MakeMatch(60, false, q: 0.5, notch: 0.5, pep: 0, pepQ: 2));
            var filtered = FilteredPsms.Filter(matches, new CommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.01));

            Assert.That(filtered.FilteringNotPerformed, Is.False);
            Assert.That(filtered.Tier, Is.Not.Null);
            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.QValueNotch));
            Assert.That(filtered.FilterThreshold, Is.EqualTo(0.01));
            Assert.That(filtered.Count(), Is.EqualTo(50));
            Assert.That(filtered.TargetPsmsAboveThreshold, Is.EqualTo(50));
            Assert.That(filtered.GetThresholdString(), Is.EqualTo("< 0.01"));
        }

        [Test]
        public static void BetweenOneHundredAndOneThousandUntrainedMatchesAreNotAllDropped()
        {
            // Before: PEP mode with 100-1000 matches filtered on PEP q-values that were never computed (all 2): nothing passed.
            var matches = Matches(640, pepTrained: false);
            var filtered = FilteredPsms.Filter(matches, new CommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.01));
            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.QValueNotch));
            Assert.That(filtered.Count(), Is.EqualTo(640));

            var peptides = FilteredPsms.Filter(matches, new CommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.01), filterAtPeptideLevel: true);
            Assert.That(peptides.FilterType, Is.EqualTo(FilterType.QValueNotch));
            Assert.That(peptides.Count(), Is.EqualTo(640));
        }

        [Test]
        public static void TieredFilterIsStrictAndReadsTheChosenValue()
        {
            var commonParameters = new CommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.01);
            var matches = Matches(20, pepTrained: true);
            var atThreshold = MakeMatch(30, false, q: 0.001, notch: 0.001, pep: 0.5, pepQ: 0.01);
            var poorQGoodPep = MakeMatch(31, false, q: 0.5, notch: 0.5, pep: 0.001, pepQ: 0.002);
            matches.Add(atThreshold);
            matches.Add(poorQGoodPep);

            var filtered = FilteredPsms.Filter(matches, commonParameters);
            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(filtered.Contains(atThreshold), Is.False, "exactly at the threshold must fail");
            Assert.That(filtered.Contains(poorQGoodPep), Is.True, "the PEP tier reads only the PEP q-value");
            Assert.That(filtered.GetFilterValue(poorQGoodPep), Is.EqualTo(0.002));
        }

        [Test]
        public static void DefaultSettingsKeepTheLegacyFilter()
        {
            // Defaults (q 0.01, PEP q 1.0): q-value AND notch, inclusive, no tier.
            var commonParameters = new CommonParameters();
            var atThreshold = MakeMatch(0, false, q: 0.01, notch: 0.01, pep: 0, pepQ: 2);
            var badNotch = MakeMatch(1, false, q: 0.001, notch: 0.02, pep: 0, pepQ: 2);
            var good = MakeMatch(2, false, q: 0.001, notch: 0.001, pep: 0, pepQ: 2);
            var filtered = FilteredPsms.Filter(new List<SpectralMatch> { atThreshold, badNotch, good }, commonParameters);

            Assert.That(filtered.Tier, Is.Null);
            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.QValue));
            Assert.That(filtered.Contains(atThreshold), Is.True, "legacy comparison stays inclusive");
            Assert.That(filtered.Contains(badNotch), Is.False, "legacy filter still needs the notch too");
            Assert.That(filtered.Contains(good), Is.True);
            Assert.That(filtered.GetThresholdString(), Is.EqualTo("<= 0.01"));
            Assert.That(filtered.GetFilterValue(good), Is.EqualTo(0.001));
            Assert.That(filtered.PassesThreshold(0.01), Is.True);
        }

        [Test]
        public static void CountPsmUsesTheTierUnderTieredMode()
        {
            // q fails, notch passes: counted by the tiered filter's notch tier, not by the legacy q AND notch rule.
            var matches = Enumerable.Range(0, 5).Select(i => MakeMatch(i, false, q: 0.5, notch: 0.001, pep: 0, pepQ: 2)).ToList();
            FdrAnalysisEngine.CountPsm(matches, new CommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.01));
            Assert.That(matches.All(m => m.PsmCount == 1), Is.True);

            var legacy = Enumerable.Range(0, 5).Select(i => MakeMatch(i, false, q: 0.5, notch: 0.001, pep: 0, pepQ: 2)).ToList();
            FdrAnalysisEngine.CountPsm(legacy);
            Assert.That(legacy.All(m => m.PsmCount == 0), Is.True);
        }

        [Test]
        public static void ProteinScoringUsesTheSameTier()
        {
            // Untrained PEP: the protein engine must use the notch tier, as the PSM tables do; before, it filtered on
            // the never-computed PEP q-values only when every one was 2, and never decided per level.
            var commonParameters = new CommonParameters(qValueThreshold: 0.01, pepQValueThreshold: 0.01);
            var matches = Enumerable.Range(0, 5).Select(i => MakeMatch(i, false, q: 0.5, notch: 0.001, pep: 0, pepQ: 2)).ToList();
            var filtered = FilteredPsms.Filter(matches, commonParameters);
            Assert.That(filtered.Count(), Is.EqualTo(5), "the notch tier keeps matches whose q-value alone would fail");

            var parsimony = (ProteinParsimonyResults)new ProteinParsimonyEngine(
                matches, false, commonParameters, null, new List<string>()).Run();
            var scored = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(parsimony.ProteinGroups, matches,
                false, false, true, commonParameters, null, new List<string>()).Run();
            Assert.That(scored.SortedAndScoredProteinGroups.Count(p => !p.IsDecoy), Is.EqualTo(5));
            Assert.That(scored.SortedAndScoredProteinGroups.All(p => p.AllPsmsBelowOnePercentFDR.Count == 1), Is.True);
        }
    }
}
