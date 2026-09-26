using System.Collections.Generic;
using System.Linq;
using System.Text;
using Chemistry;
using EngineLayer;
using EngineLayer.Truncation;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Phase 1 tests: the TruncationAcceptor precursor filter (#7) and the Pass 2 dual single-series
    /// scoring engine (#4-#8, #4a, #6). Uses a fully in-memory synthetic fixture — three top-down
    /// parent proteoforms and eight hand-built scans whose expected outcomes are encoded per test.
    /// </summary>
    [TestFixture]
    public class TruncationSearchEngineTests
    {
        // ---------- TruncationAcceptor (decision #7): 0 < M_obs <= M_theo + tolerance ----------

        [Test]
        public void Acceptor_ParentHeavierThanObserved_Accepted()
        {
            var acceptor = new TruncationAcceptor(new AbsoluteTolerance(0.5));
            Assert.That(acceptor.Accepts(1000.0, 1010.0), Is.EqualTo(0)); // truncation candidate
        }

        [Test]
        public void Acceptor_IntactEqualMass_Accepted()
        {
            var acceptor = new TruncationAcceptor(new AbsoluteTolerance(0.5));
            Assert.That(acceptor.Accepts(1000.0, 1000.0), Is.EqualTo(0));
        }

        [Test]
        public void Acceptor_ObservedExceedsTheoWithinTolerance_Accepted()
        {
            var acceptor = new TruncationAcceptor(new AbsoluteTolerance(0.5));
            Assert.That(acceptor.Accepts(1000.3, 1000.0), Is.EqualTo(0));
        }

        [Test]
        public void Acceptor_ObservedExceedsTheoBeyondTolerance_Rejected()
        {
            var acceptor = new TruncationAcceptor(new AbsoluteTolerance(0.5));
            Assert.That(acceptor.Accepts(1000.6, 1000.0), Is.EqualTo(-1));
        }

        [Test]
        public void Acceptor_NonPositiveObserved_Rejected()
        {
            var acceptor = new TruncationAcceptor(new AbsoluteTolerance(0.5));
            Assert.That(acceptor.Accepts(0.0, 1000.0), Is.EqualTo(-1));
            Assert.That(acceptor.Accepts(-5.0, 1000.0), Is.EqualTo(-1));
        }

        [Test]
        public void Acceptor_SingleNotch()
        {
            var acceptor = new TruncationAcceptor(new AbsoluteTolerance(0.5));
            Assert.That(acceptor.NumNotches, Is.EqualTo(1));
        }

        [Test]
        public void Acceptor_AllowedIntervalFromObserved_IsLowerBoundedAndUnboundedAbove()
        {
            var acceptor = new TruncationAcceptor(new AbsoluteTolerance(0.5));
            var interval = acceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(1000.0).Single();
            Assert.That(interval.Minimum, Is.EqualTo(999.5).Within(1e-9)); // M_obs - tol
            Assert.That(double.IsPositiveInfinity(interval.Maximum), Is.True);
            Assert.That(interval.Notch, Is.EqualTo(0));
        }

        [Test]
        public void Acceptor_AllowedIntervalFromTheoretical_IsZeroToTheoPlusTolerance()
        {
            var acceptor = new TruncationAcceptor(new AbsoluteTolerance(0.5));
            var interval = acceptor.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(1000.0).Single();
            Assert.That(interval.Minimum, Is.EqualTo(0));
            Assert.That(interval.Maximum, Is.EqualTo(1000.5).Within(1e-9)); // M_theo + tol
            Assert.That(interval.Notch, Is.EqualTo(0));
        }

        [Test]
        public void Acceptor_ProseString_DescribesRule()
        {
            var prose = new TruncationAcceptor(new AbsoluteTolerance(0.5)).ToProseString();
            // The prose must actually convey the acceptance rule (0 < observed <= theoretical + tolerance),
            // not merely mention "observed". Assert the structural pieces - both mass terms and both
            // comparison operators - while staying agnostic to how the tolerance itself renders.
            Assert.That(prose, Does.Contain("observed"));
            Assert.That(prose, Does.Contain("theoretical"));
            Assert.That(prose, Does.Contain("0 <"));
            Assert.That(prose, Does.Contain("<="));
        }

        // ---------- Pass 2 engine: synthetic fixture (3 parents, 8 scans) ----------

        // 19 standard residues, no Met (avoids initiator-methionine handling during top-down digestion).
        private const string AlphabetP1 = "ACDEFGHIKLNPQRSTVWY";
        private const string AlphabetP2 = "YWVTSRQPNLKIHGFEDCA";
        private const string AlphabetP3 = "DEKRHSTNQGAVLIFYWPC";

        private CommonParameters _cp;
        private PeptideWithSetModifications _p1; // 50 residues
        private PeptideWithSetModifications _p2; // 80 residues
        private PeptideWithSetModifications _p3; // 30 residues
        private List<TruncationParentSelection> _results;

        [OneTimeSetUp]
        public void BuildFixtureAndRun()
        {
            _cp = new CommonParameters(
                dissociationType: DissociationType.HCD,
                scoreCutoff: 3,
                precursorMassTolerance: new PpmTolerance(10),
                productMassTolerance: new PpmTolerance(20));

            string p1Seq = RepeatTo(AlphabetP1, 50);
            string p2Seq = RepeatTo(AlphabetP2, 80);
            string p3Seq = RepeatTo(AlphabetP3, 30);

            _p1 = BuildTopDownProteoform(p1Seq, "P1");
            _p2 = BuildTopDownProteoform(p2Seq, "P2");
            _p3 = BuildTopDownProteoform(p3Seq, "P3");

            // Truncated forms used to source each scan's single-series ions.
            var p1ToResidue40 = BuildTopDownProteoform(p1Seq.Substring(0, 40), "P1_1-40"); // C-terminal truncation of P1
            var p2FromResidue6 = BuildTopDownProteoform(p2Seq.Substring(5), "P2_6-80");    // N-terminal truncation of P2
            var p3FromResidue2 = BuildTopDownProteoform(p3Seq.Substring(1), "P3_2-30");    // N-terminal truncation of P3

            // Two real C-series ions of full P1 used as opposing-series noise in S7.
            List<double> p1OpposingCIons = SeriesMasses(_p1, FragmentationTerminus.C).Take(2).ToList();

            var scans = new[]
            {
                BuildScan(1, _p1.MonoisotopicMass, new List<double> { 500.0 }),                                 // S1 intact P1
                BuildScan(2, p1ToResidue40.MonoisotopicMass, SeriesMasses(p1ToResidue40, FragmentationTerminus.N)), // S2 C-trunc -> N series
                BuildScan(3, p2FromResidue6.MonoisotopicMass, SeriesMasses(p2FromResidue6, FragmentationTerminus.C)), // S3 N-trunc -> C series
                BuildScan(4, p3FromResidue2.MonoisotopicMass, SeriesMasses(p3FromResidue2, FragmentationTerminus.C)), // S4 N-trunc -> C series
                BuildScan(5, 2000.0, new List<double> { 123.4567, 456.7891, 789.0123 }),                        // S5 no plausible parent
                BuildScan(6, _p2.MonoisotopicMass, new List<double> { 500.0 }),                                 // S6 intact P2
                BuildScan(7, p1ToResidue40.MonoisotopicMass,                                                    // S7 like S2 + opposing-series noise
                    SeriesMasses(p1ToResidue40, FragmentationTerminus.N).Concat(p1OpposingCIons).ToList()),
                BuildScan(8, _p2.MonoisotopicMass + 1000.0, SeriesMasses(_p2, FragmentationTerminus.C).Take(5).ToList()) // S8 heavier than every parent
            };

            // Pass 1 hits: S1 -> P1 intact, S6 -> P2 intact (keyed by one-based scan number).
            var pass1 = new Dictionary<int, IReadOnlyList<SpectralMatch>>
            {
                { 1, new[] { Pass1Match(_p1, scans[0], notch: 0) } },
                { 6, new[] { Pass1Match(_p2, scans[5], notch: 0) } }
            };

            var parents = new List<TruncationParent>
            {
                new(_p1, "P1", "P1", false),
                new(_p2, "P2", "P2", false),
                new(_p3, "P3", "P3", false)
            };

            var engine = new TruncationSearchEngine(parents, scans, _cp,
                new TruncationAcceptor(_cp.PrecursorMassTolerance), pass1);
            _results = engine.Run();
        }

        [Test]
        public void ScansS1andS6_AreIntactSkipped()
        {
            Assert.That(_results[0].Outcome, Is.EqualTo(TruncationScanOutcome.IntactInherited));
            Assert.That(_results[5].Outcome, Is.EqualTo(TruncationScanOutcome.IntactInherited));
        }

        [Test]
        public void ScanS2_PicksP1_AndNTermSeries()
        {
            var s2 = _results[1];
            Assert.That(s2.Outcome, Is.EqualTo(TruncationScanOutcome.Winner));
            Assert.That(s2.WinningParent.ProteinAccession, Is.EqualTo("P1"));
            Assert.That(s2.WinningSeries, Is.EqualTo(FragmentationTerminus.N));
        }

        [Test]
        public void ScanS3_PicksP2_AndCTermSeries()
        {
            var s3 = _results[2];
            Assert.That(s3.Outcome, Is.EqualTo(TruncationScanOutcome.Winner));
            Assert.That(s3.WinningParent.ProteinAccession, Is.EqualTo("P2"));
            Assert.That(s3.WinningSeries, Is.EqualTo(FragmentationTerminus.C));
        }

        [Test]
        public void ScanS4_PicksP3_AndCTermSeries()
        {
            var s4 = _results[3];
            Assert.That(s4.Outcome, Is.EqualTo(TruncationScanOutcome.Winner));
            Assert.That(s4.WinningParent.ProteinAccession, Is.EqualTo("P3"));
            Assert.That(s4.WinningSeries, Is.EqualTo(FragmentationTerminus.C));
        }

        [Test]
        public void ScanS5_NoWinner()
        {
            Assert.That(_results[4].Outcome, Is.EqualTo(TruncationScanOutcome.NoWinner));
        }

        [Test]
        public void ScanS7_OpposingIonsIgnored()
        {
            // Same winner as S2: the two wrong-series peaks neither disqualify P1 nor flip the winner.
            var s7 = _results[6];
            Assert.That(s7.Outcome, Is.EqualTo(TruncationScanOutcome.Winner));
            Assert.That(s7.WinningParent.ProteinAccession, Is.EqualTo("P1"));
            Assert.That(s7.WinningSeries, Is.EqualTo(FragmentationTerminus.N));

            // The opposing-series peaks must not be MATCHED: they can only dilute the intensity-normalized
            // score (they enlarge the total-intensity denominator), never raise it, and must not change the
            // integer matched-ion count -- i.e. the score moves by less than one ion. (Verified empirically:
            // exact equality fails because the two extra peaks shift the normalization by ~0.05; a matched
            // opposing ion would instead RAISE the score by ~1.) (#17)
            var s2 = _results[1];
            Assert.That(s7.Score, Is.LessThanOrEqualTo(s2.Score));
            Assert.That(s2.Score - s7.Score, Is.LessThan(1.0));
        }

        [Test]
        public void ScanS8_AllParentsTooLight()
        {
            Assert.That(_results[7].Outcome, Is.EqualTo(TruncationScanOutcome.NoWinner));
        }

        /// <summary>
        /// Top-down Pass 1 runs with ThreeMM, so an intact match can sit at notch 1: the observed precursor is one
        /// C13 heavier than the theoretical mass. It must still be inherited, at its own notch, not searched for
        /// truncations (#4a).
        /// </summary>
        [Test]
        public void IntactMatchAtNonzeroNotch_IsInheritedAtThatNotch()
        {
            Ms2ScanWithSpecificMass scan = BuildScan(1, _p1.MonoisotopicMass + Constants.C13MinusC12,
                SeriesMasses(_p1, FragmentationTerminus.Both));
            SpectralMatch pass1Match = Pass1Match(_p1, scan, notch: 1);
            var engine = new TruncationSearchEngine(new List<TruncationParent> { new(_p1, "P1", "P1", false) },
                new[] { scan }, _cp, new TruncationAcceptor(_cp.PrecursorMassTolerance),
                new Dictionary<int, IReadOnlyList<SpectralMatch>> { { 1, new[] { pass1Match } } });

            TruncationParentSelection selection = engine.Run().Single();

            Assert.That(selection.Outcome, Is.EqualTo(TruncationScanOutcome.IntactInherited));
            Assert.That(selection.IntactMatch, Is.SameAs(pass1Match));
            SpectralMatch inherited = TruncationPass3.InheritAsFullLength(selection.IntactMatch, scan, 0, _cp);
            inherited.ResolveAllAmbiguities();
            Assert.That(inherited.Notch, Is.EqualTo(1));
            Assert.That(inherited.BaseSequence, Is.EqualTo(_p1.BaseSequence));
        }

        /// <summary>
        /// A chimeric scan carries one entry per precursor, and Pass 1 can match both. Each entry must be
        /// inherited with its own Pass 1 match, not only the best match on the scan (#4a).
        /// </summary>
        [Test]
        public void ChimericScan_EachPrecursorInheritsItsOwnPass1Match()
        {
            var fragments = SeriesMasses(_p1, FragmentationTerminus.Both).Concat(SeriesMasses(_p2, FragmentationTerminus.Both)).ToList();
            Ms2ScanWithSpecificMass p1Entry = BuildScan(9, _p1.MonoisotopicMass, fragments);
            Ms2ScanWithSpecificMass p2Entry = BuildScan(9, _p2.MonoisotopicMass, fragments);
            SpectralMatch p1Match = Pass1Match(_p1, p1Entry, notch: 0, score: 20);
            SpectralMatch p2Match = Pass1Match(_p2, p2Entry, notch: 0, score: 30);
            var engine = new TruncationSearchEngine(
                new List<TruncationParent> { new(_p1, "P1", "P1", false), new(_p2, "P2", "P2", false) },
                new[] { p1Entry, p2Entry }, _cp, new TruncationAcceptor(_cp.PrecursorMassTolerance),
                new Dictionary<int, IReadOnlyList<SpectralMatch>> { { 9, new[] { p2Match, p1Match } } }); // best first

            List<TruncationParentSelection> selections = engine.Run();

            Assert.That(selections.Select(s => s.Outcome), Is.All.EqualTo(TruncationScanOutcome.IntactInherited));
            Assert.That(selections[0].IntactMatch, Is.SameAs(p1Match));
            Assert.That(selections[1].IntactMatch, Is.SameAs(p2Match));
        }

        /// <summary>
        /// Two parents with the same sequence under different accessions tie exactly in Pass 2. Both must reach
        /// Pass 3 so the truncation is reported once with both accessions, as decision #10 describes, instead of
        /// under whichever parent the candidate order put first.
        /// </summary>
        [Test]
        public void TiedParents_AreCarriedToPass3_AndCollapseToOneAmbiguousPsm()
        {
            string p1Seq = RepeatTo(AlphabetP1, 50);
            var isoform = BuildTopDownProteoform(p1Seq, "P1b");
            var p1ToResidue40 = BuildTopDownProteoform(p1Seq.Substring(0, 40), "P1_1-40");
            Ms2ScanWithSpecificMass scan = BuildScan(2, p1ToResidue40.MonoisotopicMass, SeriesMasses(p1ToResidue40, FragmentationTerminus.N));
            var engine = new TruncationSearchEngine(
                new List<TruncationParent> { new(_p1, "P1", "P1", false), new(isoform, "P1b", "P1b", false) },
                new[] { scan }, _cp, new TruncationAcceptor(_cp.PrecursorMassTolerance));

            TruncationParentSelection selection = engine.Run().Single();

            Assert.That(selection.Outcome, Is.EqualTo(TruncationScanOutcome.Winner));
            Assert.That(selection.TiedWinners.Select(t => t.Parent.ProteinAccession),
                Is.EquivalentTo(new[] { "P1", "P1b" }.Except(new[] { selection.WinningParent.ProteinAccession })));

            List<TruncationPsm> psms = TruncationPass3.ScoreTruncations(selection, _cp,
                SearchTask.GetMassDiffAcceptor(_cp.PrecursorMassTolerance, MassDiffAcceptorType.Exact, null));
            Assert.That(psms, Has.Count.EqualTo(2));

            SpectralMatch collapsed = TruncationPass3.CollapseDuplicateTruncations(psms).Single();
            collapsed.ResolveAllAmbiguities();
            Assert.That(collapsed.BestMatchingBioPolymersWithSetMods.Select(b => b.SpecificBioPolymer.Parent.Accession).Distinct(),
                Is.EquivalentTo(new[] { "P1", "P1b" }));
        }

        /// <summary>
        /// An inherited full-length PSM keeps every Pass 1 hypothesis, so a Pass 1 match that was ambiguous
        /// between two accessions stays ambiguous in the truncation output (#4a).
        /// </summary>
        [Test]
        public void InheritAsFullLength_KeepsPass1ProteinAmbiguity()
        {
            var isoform = BuildTopDownProteoform(RepeatTo(AlphabetP1, 50), "P1b");
            Ms2ScanWithSpecificMass scan = BuildScan(1, _p1.MonoisotopicMass, SeriesMasses(_p1, FragmentationTerminus.Both));
            SpectralMatch pass1Match = Pass1Match(_p1, scan, notch: 0);
            pass1Match.AddOrReplace(isoform, pass1Match.Score, 0, reportAllAmbiguity: true, new List<MatchedFragmentIon>());

            SpectralMatch inherited = TruncationPass3.InheritAsFullLength(pass1Match, scan, 0, _cp);
            inherited.ResolveAllAmbiguities();

            Assert.That(inherited.BestMatchingBioPolymersWithSetMods.Select(b => b.SpecificBioPolymer.Parent.Accession),
                Is.EquivalentTo(new[] { "P1", "P1b" }));
            Assert.That(inherited.BestMatchingBioPolymersWithSetMods.Select(b => ((PeptideWithSetModifications)b.SpecificBioPolymer).Description),
                Is.All.EqualTo(TruncationPass3.FullLength));
        }

        /// <summary>
        /// When every parent is excluded as oversized, the index is empty and every scan is a no-winner rather
        /// than an index-out-of-range on the missing heaviest parent.
        /// </summary>
        [Test]
        public void NoIndexableParents_EveryScanIsNoWinner()
        {
            Ms2ScanWithSpecificMass scan = BuildScan(2, _p1.MonoisotopicMass - 500, SeriesMasses(_p1, FragmentationTerminus.N));
            var engine = new TruncationSearchEngine(new List<TruncationParent> { new(_p1, "P1", "P1", false) },
                new[] { scan }, _cp, new TruncationAcceptor(_cp.PrecursorMassTolerance), maxFragmentSize: 1000);

            List<TruncationParentSelection> selections = engine.Run();

            Assert.That(engine.IndexedParentCount, Is.EqualTo(0));
            Assert.That(selections.Single().Outcome, Is.EqualTo(TruncationScanOutcome.NoWinner));
        }

        [Test]
        public void OversizedParent_ExcludedAndWarned()
        {
            var oversized = BuildTopDownProteoform(RepeatTo("WFYRLKHEQND", 400), "BIG"); // ~35+ kDa
            Assert.That(oversized.MonoisotopicMass, Is.GreaterThan(30000));

            var parents = new List<TruncationParent>
            {
                new(_p1, "P1", "P1", false),
                new(_p2, "P2", "P2", false),
                new(_p3, "P3", "P3", false),
                new(oversized, "BIG", "BIG", false)
            };

            var engine = new TruncationSearchEngine(parents, System.Array.Empty<Ms2ScanWithSpecificMass>(),
                _cp, new TruncationAcceptor(_cp.PrecursorMassTolerance));

            Assert.That(engine.ExcludedOversizedParentCount, Is.EqualTo(1));
            Assert.That(engine.Warnings, Has.Exactly(1).EqualTo("1 parent proteoforms exceeded MaxFragmentSize and were excluded."));
        }

        // ---------- parent loading / filtering (#2, #3) ----------

        [Test]
        public void PipeAmbiguousParent_BecomesMultipleParents()
        {
            var row = new Pass1ProteoformRow
            {
                FullSequence = "PEPTIDEPEPTIDEK|PEPTIDEPEPTIDER|PEPTIDEPEPTIDEM",
                ProteinAccession = "ACC1",
                PepQValue = 0.0,
                IsDecoy = false
            };

            var parents = TruncationParentBuilder.BuildParents(new[] { row }, 0.10);

            Assert.That(parents.Count, Is.EqualTo(3));
            Assert.That(parents.Select(p => p.Proteoform.BaseSequence).Distinct().Count(), Is.EqualTo(3));
            // All three alternatives tie back to the one originating row (#2).
            Assert.That(parents.All(p => ReferenceEquals(p.OriginatingId, row)), Is.True);
        }

        [Test]
        public void ParentFilter_UsesPepQWhenComputedElseNotchQ()
        {
            // PEP computed -> judged on PEP q-value, ignoring notch q-value.
            Assert.That(TruncationParentBuilder.PassesParentFilter(new Pass1ProteoformRow { PepQValue = 0.05, NotchQValue = 0.99 }, 0.10), Is.True);
            Assert.That(TruncationParentBuilder.PassesParentFilter(new Pass1ProteoformRow { PepQValue = 0.50, NotchQValue = 0.00 }, 0.10), Is.False);
            // PEP not computed -> fall back to notch q-value.
            Assert.That(TruncationParentBuilder.PassesParentFilter(new Pass1ProteoformRow { PepQValue = null, NotchQValue = 0.05 }, 0.10), Is.True);
            Assert.That(TruncationParentBuilder.PassesParentFilter(new Pass1ProteoformRow { PepQValue = null, NotchQValue = 0.50 }, 0.10), Is.False);
        }

        [Test]
        public void BuildParents_ExcludesRowsFailingFilter()
        {
            var rows = new[]
            {
                new Pass1ProteoformRow { FullSequence = "AAAAAAAAAA", ProteinAccession = "A", PepQValue = 0.01 },
                new Pass1ProteoformRow { FullSequence = "CCCCCCCCCC", ProteinAccession = "C", PepQValue = 0.50 }
            };

            var parents = TruncationParentBuilder.BuildParents(rows, 0.10);

            Assert.That(parents.Count, Is.EqualTo(1));
            Assert.That(parents[0].ProteinAccession, Is.EqualTo("A"));
        }

        [Test]
        public void PerScanRestriction_ExcludesParentFromDisallowedProtein()
        {
            // A C-terminal truncation of P1 (N-series ions) that would normally win for P1.
            string p1Seq = RepeatTo(AlphabetP1, 50);
            var p1 = BuildTopDownProteoform(p1Seq, "P1");
            var p1ToResidue40 = BuildTopDownProteoform(p1Seq.Substring(0, 40), "x");
            var scan = BuildScan(1, p1ToResidue40.MonoisotopicMass, SeriesMasses(p1ToResidue40, FragmentationTerminus.N));
            var scans = new[] { scan };
            var parents = new List<TruncationParent> { new(p1, "P1", "P1", false) };

            // Allowed set for this scan EXCLUDES P1 -> the only candidate is filtered out -> no winner.
            var disallow = new Dictionary<Ms2ScanWithSpecificMass, HashSet<string>> { { scan, new HashSet<string>() } };
            var excluded = new TruncationSearchEngine(parents, scans, _cp,
                new TruncationAcceptor(_cp.PrecursorMassTolerance), null,
                TruncationSearchEngine.DefaultMaxFragmentSize, disallow).Run();
            Assert.That(excluded[0].Outcome, Is.EqualTo(TruncationScanOutcome.NoWinner));

            // Allowed set INCLUDES P1 -> P1 wins (N-terminal-ion series, C-terminus chopped).
            var allow = new Dictionary<Ms2ScanWithSpecificMass, HashSet<string>> { { scan, new HashSet<string> { "P1" } } };
            var allowed = new TruncationSearchEngine(parents, scans, _cp,
                new TruncationAcceptor(_cp.PrecursorMassTolerance), null,
                TruncationSearchEngine.DefaultMaxFragmentSize, allow).Run();
            Assert.That(allowed[0].Outcome, Is.EqualTo(TruncationScanOutcome.Winner));
            Assert.That(allowed[0].WinningParent.ProteinAccession, Is.EqualTo("P1"));
        }

        // ---------- fixture helpers ----------

        private static string RepeatTo(string alphabet, int length)
        {
            var sb = new StringBuilder(length);
            while (sb.Length < length)
            {
                sb.Append(alphabet);
            }
            return sb.ToString().Substring(0, length);
        }

        private static PeptideWithSetModifications BuildTopDownProteoform(string sequence, string accession)
        {
            var protein = new Protein(sequence, accession);
            var digestionParams = new DigestionParams(protease: "top-down", minPeptideLength: 1, maxPeptideLength: 100000);
            return protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
        }

        private SpectralMatch Pass1Match(PeptideWithSetModifications form, Ms2ScanWithSpecificMass scan, int notch, double score = 10) =>
            new PeptideSpectralMatch(form, notch, score, 0, scan, _cp, new List<MatchedFragmentIon>());

        private List<double> SeriesMasses(PeptideWithSetModifications form, FragmentationTerminus terminus)
        {
            var products = new List<Product>();
            form.Fragment(_cp.DissociationType, terminus, products, _cp.FragmentationParameters);
            return products.Where(p => !double.IsNaN(p.NeutralMass)).Select(p => p.NeutralMass).ToList();
        }

        private Ms2ScanWithSpecificMass BuildScan(int scanNumber, double precursorMass, IReadOnlyList<double> fragmentNeutralMasses)
        {
            const double intensity = 1000.0;
            var ordered = fragmentNeutralMasses.OrderBy(m => m).ToList();

            double[] mz = ordered.Select(m => m.ToMz(1)).ToArray();
            double[] intensities = ordered.Select(_ => intensity).ToArray();
            var spectrum = new MzSpectrum(mz, intensities, false);

            double tic = intensities.Sum();
            var msDataScan = new MsDataScan(
                massSpectrum: spectrum,
                oneBasedScanNumber: scanNumber,
                msnOrder: 2,
                isCentroid: true,
                polarity: Polarity.Positive,
                retentionTime: scanNumber,
                scanWindowRange: new MzRange(0, 1_000_000),
                scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap,
                totalIonCurrent: tic <= 0 ? 1 : tic,
                injectionTime: 1.0,
                noiseData: null,
                nativeId: $"scan={scanNumber}");

            var envelopes = ordered
                .Select(m => new IsotopicEnvelope(new List<(double mz, double intensity)> { (m.ToMz(1), intensity) }, m, 1, intensity, 0))
                .ToArray();

            return new Ms2ScanWithSpecificMass(msDataScan, precursorMass.ToMz(1), 1, "synthetic", _cp,
                neutralExperimentalFragments: envelopes);
        }
    }
}
