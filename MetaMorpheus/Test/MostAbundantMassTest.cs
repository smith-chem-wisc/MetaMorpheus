using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Tests for "most abundant mass" precursor selection (Strategy B) on the MetaMorpheus side:
    /// the <see cref="MostAbundantMassDiffAcceptor"/>, the <see cref="PrecursorMassMatchMode"/>
    /// parameter plumbing, and <see cref="Ms2ScanWithSpecificMass.PrecursorMassToMatch"/>.
    /// </summary>
    [TestFixture]
    public static class MostAbundantMassTest
    {
        private static readonly AverageResidue Averagine = new Averagine();

        [Test]
        public static void Acceptor_PinsExactMonoisotopic_AndRejectsOffByOne()
        {
            // A ~15 kDa proteoform: the observed precursor is the most abundant isotopic peak.
            const double peptideMono = 15000.0;
            double observedMostAbundant = peptideMono + Averagine.GetMostAbundantOffset(peptideMono);

            var acceptor = new MostAbundantMassDiffAcceptor("mostAbundant", new PpmTolerance(5), Averagine);

            // The correct monoisotopic candidate is accepted (notch 0)...
            Assert.That(acceptor.Accepts(observedMostAbundant, peptideMono), Is.EqualTo(0));

            // ...but candidates that are off by ±1 neutron in the monoisotopic mass are rejected,
            // because their theoretical most-abundant peak no longer lines up with the observed one.
            Assert.That(acceptor.Accepts(observedMostAbundant, peptideMono + Constants.C13MinusC12), Is.EqualTo(-1));
            Assert.That(acceptor.Accepts(observedMostAbundant, peptideMono - Constants.C13MinusC12), Is.EqualTo(-1));
        }

        [Test]
        public static void Acceptor_ContrastsWithOneMmNotch_WhichAcceptsOffByOne()
        {
            // Illustrates what most-abundant mode replaces: the OneMM notch acceptor (which operates
            // in monoisotopic space) deliberately accepts a candidate one neutron lighter than the
            // observed monoisotopic mass — exactly the off-by-one ambiguity this work removes.
            const double observedMono = 15000.0;
            var oneMm = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.OneMM, null);

            Assert.That(oneMm.Accepts(observedMono, observedMono), Is.GreaterThanOrEqualTo(0));
            Assert.That(oneMm.Accepts(observedMono, observedMono - 1.0029), Is.GreaterThanOrEqualTo(0)); // off-by-one accepted
        }

        [Test]
        public static void Acceptor_TheoreticalInterval_IsCenteredOnMostAbundantMass()
        {
            const double peptideMono = 20000.0;
            var tol = new PpmTolerance(5);
            var acceptor = new MostAbundantMassDiffAcceptor("mostAbundant", tol, Averagine);

            double expectedCenter = peptideMono + Averagine.GetMostAbundantOffset(peptideMono);

            var intervals = System.Linq.Enumerable.ToList(acceptor.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMono));
            Assert.That(intervals.Count, Is.EqualTo(1));
            Assert.That(intervals[0].Notch, Is.EqualTo(0));
            Assert.That(intervals[0].Minimum, Is.EqualTo(tol.GetMinimumValue(expectedCenter)).Within(1e-6));
            Assert.That(intervals[0].Maximum, Is.EqualTo(tol.GetMaximumValue(expectedCenter)).Within(1e-6));
        }

        [Test]
        public static void SearchTask_GetMassDiffAcceptor_ReturnsMostAbundantWhenModeSet()
        {
            var acc = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.OneMM, null,
                PrecursorMassMatchMode.MostAbundant, Averagine);
            Assert.That(acc, Is.TypeOf<MostAbundantMassDiffAcceptor>());

            // Default mode is unaffected.
            var def = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.OneMM, null);
            Assert.That(def, Is.Not.TypeOf<MostAbundantMassDiffAcceptor>());
        }

        [Test]
        public static void Ms2Scan_PrecursorMassToMatch_DefaultsToPrecursorMass()
        {
            var scan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap,
                double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null,
                DissociationType.AnyActivationType, 1, null);

            // No precursorMassToMatch supplied → falls back to the monoisotopic PrecursorMass.
            var defaultScan = new Ms2ScanWithSpecificMass(scan, 1500.0.ToMz(2), 2, "", new CommonParameters());
            Assert.That(defaultScan.PrecursorMassToMatch, Is.EqualTo(defaultScan.PrecursorMass).Within(1e-9));

            // Explicit precursorMassToMatch is carried through without altering PrecursorMass.
            double matchMass = defaultScan.PrecursorMass + 5.0;
            var mostAbundantScan = new Ms2ScanWithSpecificMass(scan, 1500.0.ToMz(2), 2, "", new CommonParameters(),
                precursorMassToMatch: matchMass);
            Assert.That(mostAbundantScan.PrecursorMassToMatch, Is.EqualTo(matchMass).Within(1e-9));
            Assert.That(mostAbundantScan.PrecursorMass, Is.EqualTo(defaultScan.PrecursorMass).Within(1e-9));
        }

        [Test]
        public static void CommonParameters_PrecursorMassMatchMode_RoundTripsAndClones()
        {
            var cp = new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);
            Assert.That(cp.PrecursorMassMatchMode, Is.EqualTo(PrecursorMassMatchMode.MostAbundant));

            // Clone preserves it.
            Assert.That(cp.Clone().PrecursorMassMatchMode, Is.EqualTo(PrecursorMassMatchMode.MostAbundant));

            // Toml round-trip preserves it.
            string toml = Toml.WriteString(cp, MetaMorpheusTask.tomlConfig);
            var restored = Toml.ReadString<CommonParameters>(toml, MetaMorpheusTask.tomlConfig);
            Assert.That(restored.PrecursorMassMatchMode, Is.EqualTo(PrecursorMassMatchMode.MostAbundant));

            // Default stays Monoisotopic.
            Assert.That(new CommonParameters().PrecursorMassMatchMode, Is.EqualTo(PrecursorMassMatchMode.Monoisotopic));
        }
    }
}
