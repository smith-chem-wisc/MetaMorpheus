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
            Assert.That(prose, Does.Contain("observed"));
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
            var pass1 = new Dictionary<int, double> { { 1, _p1.MonoisotopicMass }, { 6, _p2.MonoisotopicMass } };

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
        }

        [Test]
        public void ScanS8_AllParentsTooLight()
        {
            Assert.That(_results[7].Outcome, Is.EqualTo(TruncationScanOutcome.NoWinner));
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
