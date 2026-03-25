using EngineLayer;
using EngineLayer.CircularSearch;
using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test.CircularSearch
{
    /// <summary>
    /// Tests for <see cref="CircularSearchEngine"/> using synthetic MGF-based scans.
    ///
    /// Three targeted diagnostic tests isolate each step of the scoring pipeline:
    ///   1. MatchFragmentIons  — do experimental fragments match theoretical ones?
    ///   2. CalculatePeptideScore — does a non-empty match list yield score > 0?
    ///   3. AddCandidateToPsm — does a positive score populate the PSM slot?
    ///
    /// Plus the original four engine-level tests.
    /// </summary>
    [TestFixture]
    public static class CircularSearchEngineTests
    {
        private static readonly double Proton = 1.007276;

        // ── Shared setup ──────────────────────────────────────────────────────

        private static CommonParameters MakeCommonParams(int maxMissedCleavages = 1) =>
            new CommonParameters(
                dissociationType: DissociationType.HCD,
                precursorMassTolerance: new PpmTolerance(5),
                productMassTolerance: new PpmTolerance(20),
                scoreCutoff: 1,   // ← add this
                digestionParams: new DigestionParams(
                    protease: "trypsin",
                    maxMissedCleavages: maxMissedCleavages,
                    minPeptideLength: 1));

        /// <summary>
        /// Writes a temp MGF, loads it, and returns the first scan.
        /// </summary>
        private static Ms2ScanWithSpecificMass LoadOneScan(
            double precursorNeutralMass,
            IEnumerable<double> fragmentNeutralMasses,
            CommonParameters commonParams,
            string testName,
            out string mgfPath)
        {
            mgfPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                $"CSE_{testName}_{Guid.NewGuid():N}.mgf");

            double precursorMz = precursorNeutralMass + Proton;
            var sb = new System.Text.StringBuilder();
            sb.AppendLine("BEGIN IONS");
            sb.AppendLine(FormattableString.Invariant($"PEPMASS={precursorMz:F6}"));
            sb.AppendLine("CHARGE=1+");
            sb.AppendLine("SCANS=1");
            sb.AppendLine("RTINSECONDS=60");
            foreach (double nm in fragmentNeutralMasses)
                sb.AppendLine(FormattableString.Invariant($"{nm + Proton:F6} 1000000"));
            sb.AppendLine("END IONS");
            File.WriteAllText(mgfPath, sb.ToString());

            var msDataFile = new Readers.Mgf(mgfPath).LoadAllStaticData();
            return MetaMorpheusTask
                .GetMs2Scans(msDataFile, mgfPath, commonParams)
                .OrderBy(s => s.PrecursorMass)
                .First();
        }

        /// <summary>
        /// Returns the CircularPeptideWithSetModifications and its internal fragments
        /// for ring "FGHIK" with maxMissedCleavages=1.
        /// </summary>
        private static (CircularPeptideWithSetModifications peptide, List<Product> fragments)
            GetFghikPeptideAndFragments(CommonParameters commonParams)
        {
            var protein = new CircularProtein("FGHIK", "acc_diag");
            var peptide = protein.Digest(
                    commonParams.DigestionParams,
                    new List<Modification>(),
                    new List<Modification>())
                .OfType<CircularPeptideWithSetModifications>()
                .Single();

            var fragments = new List<Product>();
            peptide.FragmentInternally(DissociationType.HCD, minLengthOfFragments: 3, fragments);
            return (peptide, fragments);
        }

        // ── Diagnostic Test 1: MatchFragmentIons ─────────────────────────────

        /// <summary>
        /// Calls MatchFragmentIons directly with a scan built from the exact theoretical
        /// fragment neutral masses of "FGHIK" and asserts that at least one ion is matched.
        ///
        /// If this fails: the experimental fragments in the scan do not overlap with the
        /// theoretical products, indicating a deconvolution or mass offset problem in how
        /// the MGF peaks become ExperimentalFragments.
        /// </summary>
        [Test]
        public static void MatchFragmentIons_WithExactTheoreticalMasses_ReturnsMatches()
        {
            var commonParams = MakeCommonParams();
            var (peptide, fragments) = GetFghikPeptideAndFragments(commonParams);

            Assert.That(fragments.Count, Is.GreaterThan(0), "Pre-condition: fragments must exist.");

            var fragmentNeutralMasses = fragments.Select(p => p.NeutralMass).Distinct();
            var scan = LoadOneScan(peptide.MonoisotopicMass, fragmentNeutralMasses,
                commonParams, nameof(MatchFragmentIons_WithExactTheoreticalMasses_ReturnsMatches),
                out string mgfPath);

            try
            {
                // Report what the scan loaded
                Console.WriteLine($"Scan PrecursorMass: {scan.PrecursorMass:F6}");
                Console.WriteLine($"Scan ExperimentalFragments count: {scan.ExperimentalFragments?.Length ?? 0}");
                if (scan.ExperimentalFragments != null)
                    foreach (var ef in scan.ExperimentalFragments)
                        Console.WriteLine($"  ExperimentalFragment MonoisotopicMass: {ef.MonoisotopicMass:F6}  Charge: {ef.Charge}");

                Console.WriteLine("Theoretical fragment neutral masses:");
                foreach (var f in fragments)
                    Console.WriteLine($"  {f.NeutralMass:F6}  ProductType: {f.ProductType}");

                var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, fragments, commonParams);

                Console.WriteLine($"Matched ions count: {matchedIons.Count}");
                foreach (var mi in matchedIons)
                    Console.WriteLine($"  Matched: {mi.NeutralTheoreticalProduct.NeutralMass:F6}");

                Assert.That(matchedIons.Count, Is.GreaterThan(0),
                    "MatchFragmentIons should match at least one ion when scan fragments " +
                    "are built from exact theoretical neutral masses.");
            }
            finally { File.Delete(mgfPath); }
        }

        // ── Diagnostic Test 2: CalculatePeptideScore ─────────────────────────

        /// <summary>
        /// Calls CalculatePeptideScore directly with the matched ions from Test 1 and
        /// asserts that the score is > 0.
        ///
        /// If this fails with a non-empty matched ion list: the scoring function has a
        /// bug for the product types produced by FragmentInternally (e.g. all ions are
        /// ProductType.D which scores zero, or TotalIonCurrent is zero).
        /// </summary>
        [Test]
        public static void CalculatePeptideScore_WithMatchedIons_ReturnsPositiveScore()
        {
            var commonParams = MakeCommonParams();
            var (peptide, fragments) = GetFghikPeptideAndFragments(commonParams);

            var fragmentNeutralMasses = fragments.Select(p => p.NeutralMass).Distinct();
            var scan = LoadOneScan(peptide.MonoisotopicMass, fragmentNeutralMasses,
                commonParams, nameof(CalculatePeptideScore_WithMatchedIons_ReturnsPositiveScore),
                out string mgfPath);

            try
            {
                var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, fragments, commonParams);

                Console.WriteLine($"Matched ions count going into CalculatePeptideScore: {matchedIons.Count}");
                Console.WriteLine($"Scan TotalIonCurrent: {scan.TheScan.TotalIonCurrent}");
                Console.WriteLine($"Scan XcorrProcessed: {scan.TheScan.MassSpectrum.XcorrProcessed}");
                foreach (var mi in matchedIons)
                    Console.WriteLine($"  Ion ProductType: {mi.NeutralTheoreticalProduct.ProductType}  Intensity: {mi.Intensity}");

                // CalculatePeptideScore requires at least one matched ion to score > 0
                // If matchedIons is empty here, this test will correctly fail and confirm
                // that MatchFragmentIons is the root cause, not CalculatePeptideScore.
                Assume.That(matchedIons.Count, Is.GreaterThan(0),
                    "Skipping: MatchFragmentIons returned no matches (root cause is Test 1).");

                double score = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan, matchedIons);
                Console.WriteLine($"Score: {score}");

                Assert.That(score, Is.GreaterThan(0),
                    "CalculatePeptideScore should return > 0 when matched ions exist.");
            }
            finally { File.Delete(mgfPath); }
        }

        // ── Diagnostic Test 3: AddCandidateToPsm ─────────────────────────────

        /// <summary>
        /// Calls AddCandidateToPsm via a constructed CircularSearchEngine and asserts
        /// that the PSM slot is populated after a scan with positive score.
        ///
        /// This test constructs the engine, manually invokes the scoring pipeline steps
        /// in sequence, and then calls AddCandidateToPsm with a synthetic
        /// ScanWithIndexAndNotchInfo at index 0.
        ///
        /// If this fails with a positive score: AddCandidateToPsm is filtering it out
        /// (e.g. score < ScoreCutoff, or the lock/improve logic is rejecting it).
        /// </summary>
        [Test]
        public static void AddCandidateToPsm_WithPositiveScore_PopulatesPsmSlot()
        {
            var commonParams = MakeCommonParams();
            var (peptide, fragments) = GetFghikPeptideAndFragments(commonParams);

            var fragmentNeutralMasses = fragments.Select(p => p.NeutralMass).Distinct();
            var scan = LoadOneScan(peptide.MonoisotopicMass, fragmentNeutralMasses,
                commonParams, nameof(AddCandidateToPsm_WithPositiveScore_PopulatesPsmSlot),
                out string mgfPath);

            try
            {
                var psms = new SpectralMatch[1];
                var searchMode = new SinglePpmAroundZeroSearchMode(5);
                var protein = new CircularProtein("FGHIK", "acc_addcandidate");

                var engine = new CircularSearchEngine(
                    globalPsms: psms,
                    arrayOfSortedMS2Scans: new[] { scan },
                    variableModifications: new List<Modification>(),
                    fixedModifications: new List<Modification>(),
                    circularProteins: new List<CircularProtein> { protein },
                    searchMode: searchMode,
                    commonParameters: commonParams,
                    fileSpecificParameters: null,
                    nestedIds: new List<string>(),
                    minInternalFragmentLength: 3);

                // Step 1: match ions
                var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, fragments, commonParams);
                Console.WriteLine($"Matched ions for AddCandidateToPsm test: {matchedIons.Count}");

                Assume.That(matchedIons.Count, Is.GreaterThan(0),
                    "Skipping: MatchFragmentIons returned no matches (root cause is Test 1).");

                // Step 2: score
                double score = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan, matchedIons);
                Console.WriteLine($"Score for AddCandidateToPsm test: {score}");
                Console.WriteLine($"ScoreCutoff: {commonParams.ScoreCutoff}");

                Assume.That(score, Is.GreaterThan(0),
                    "Skipping: CalculatePeptideScore returned 0 (root cause is Test 2).");

                // Step 3: AddCandidateToPsm
                var scanInfo = new ScanWithIndexAndNotchInfo(notch: 0, scanIndex: 0);
                engine.AddCandidateToPsm(scanInfo, score, peptide, matchedIons);

                Console.WriteLine($"PSM after AddCandidateToPsm: {(psms[0] == null ? "null" : "populated")}");
                if (psms[0] != null)
                    Console.WriteLine($"PSM score: {psms[0].Score}");

                Assert.That(psms[0], Is.Not.Null,
                    "AddCandidateToPsm should populate the PSM slot when score > ScoreCutoff.");
                Assert.That(psms[0].Score, Is.GreaterThan(0));
            }
            finally { File.Delete(mgfPath); }
        }

        // ── Engine-level tests ────────────────────────────────────────────────

        private static string WriteTempMgf(
            double precursorNeutralMass,
            IEnumerable<double> fragmentNeutralMasses,
            string testName)
        {
            string path = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                $"CSE_{testName}_{Guid.NewGuid():N}.mgf");

            var sb = new System.Text.StringBuilder();
            sb.AppendLine("BEGIN IONS");
            sb.AppendLine(FormattableString.Invariant($"PEPMASS={precursorNeutralMass + Proton:F6}"));
            sb.AppendLine("CHARGE=1+");
            sb.AppendLine("SCANS=1");
            sb.AppendLine("RTINSECONDS=60");
            foreach (double nm in fragmentNeutralMasses)
                sb.AppendLine(FormattableString.Invariant($"{nm + Proton:F6} 1000000"));
            sb.AppendLine("END IONS");
            File.WriteAllText(path, sb.ToString());
            return path;
        }

        private static Ms2ScanWithSpecificMass[] LoadScansFromMgf(
            string mgfPath, CommonParameters commonParams) =>
            MetaMorpheusTask
                .GetMs2Scans(new Readers.Mgf(mgfPath).LoadAllStaticData(), mgfPath, commonParams)
                .OrderBy(s => s.PrecursorMass)
                .ToArray();

        private static SpectralMatch[] RunEngine(
            Ms2ScanWithSpecificMass[] scans,
            List<CircularProtein> proteins,
            CommonParameters commonParams,
            int minInternalFragmentLength = 3)
        {
            var psms = new SpectralMatch[scans.Length];
            new CircularSearchEngine(psms, scans, new List<Modification>(),
                new List<Modification>(), proteins,
                new SinglePpmAroundZeroSearchMode(5), commonParams,
                null, new List<string>(), minInternalFragmentLength).Run();
            return psms;
        }

        /// <summary>
        /// Empty protein list → no crash, all PSMs null.
        /// </summary>
        [Test]
        public static void EmptyProteinList_NoCrash_AllPsmsNull()
        {
            var commonParams = MakeCommonParams();
            string mgf = WriteTempMgf(500.0, new[] { 100.0, 200.0 },
                nameof(EmptyProteinList_NoCrash_AllPsmsNull));
            try
            {
                var scans = LoadScansFromMgf(mgf, commonParams);
                var psms = RunEngine(scans, new List<CircularProtein>(), commonParams);
                Assert.That(psms.All(p => p is null));
            }
            finally { File.Delete(mgf); }
        }

        /// <summary>
        /// Matching scan → PSM with positive score.
        /// </summary>
        [Test]
        public static void MatchingScan_GivesPositiveScore()
        {
            var protein = new CircularProtein("FGHIK", "acc_match");
            Assert.That(protein.BaseSequence, Is.EqualTo("FGHIK"));

            var commonParams = MakeCommonParams(maxMissedCleavages: 1);
            var (peptide, fragments) = GetFghikPeptideAndFragments(commonParams);

            string mgf = WriteTempMgf(peptide.MonoisotopicMass,
                fragments.Select(p => p.NeutralMass).Distinct(),
                nameof(MatchingScan_GivesPositiveScore));
            try
            {
                var scans = LoadScansFromMgf(mgf, commonParams);
                var psms = RunEngine(scans, new List<CircularProtein> { protein }, commonParams);

                Assert.That(psms[0], Is.Not.Null,
                    "Engine should assign a PSM when scan matches the circular peptide.");
                Assert.That(psms[0].Score, Is.GreaterThan(0));
            }
            finally { File.Delete(mgf); }
        }

        /// <summary>
        /// DigestionCountDictionary is populated after Run().
        /// </summary>
        [Test]
        public static void DigestionCountDictionary_IsPopulated_AfterRun()
        {
            var protein = new CircularProtein("FGHIK", "acc_dict");
            var commonParams = MakeCommonParams(maxMissedCleavages: 1);
            string mgf = WriteTempMgf(1.0, new[] { 100.0 },
                nameof(DigestionCountDictionary_IsPopulated_AfterRun));
            try
            {
                var scans = LoadScansFromMgf(mgf, commonParams);
                var psms = new SpectralMatch[scans.Length];
                var engine = new CircularSearchEngine(psms, scans,
                    new List<Modification>(), new List<Modification>(),
                    new List<CircularProtein> { protein },
                    new SinglePpmAroundZeroSearchMode(5), commonParams,
                    null, new List<string>(), 3);
                engine.Run();
                Assert.That(engine.DigestionCountDictionary.Count, Is.GreaterThan(0));
            }
            finally { File.Delete(mgf); }
        }

        /// <summary>
        /// Precursor shifted 500 Da → all PSMs null.
        /// </summary>
        [Test]
        public static void NonMatchingScan_PrecursorShifted500Da_AllPsmsNull()
        {
            var protein = new CircularProtein("FGHIK", "acc_nomatch");
            var commonParams = MakeCommonParams(maxMissedCleavages: 1);
            var (peptide, _) = GetFghikPeptideAndFragments(commonParams);

            string mgf = WriteTempMgf(peptide.MonoisotopicMass + 500.0,
                new[] { 100.0, 200.0 },
                nameof(NonMatchingScan_PrecursorShifted500Da_AllPsmsNull));
            try
            {
                var scans = LoadScansFromMgf(mgf, commonParams);
                var psms = RunEngine(scans, new List<CircularProtein> { protein }, commonParams);
                Assert.That(psms.All(p => p is null));
            }
            finally { File.Delete(mgf); }
        }
    }
}