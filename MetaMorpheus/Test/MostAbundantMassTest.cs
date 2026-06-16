using System;
using System.Collections.Generic;
using System.Linq;
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
        public static void Acceptor_OnApexCandidate_GetsNotchZero()
        {
            // A ~15 kDa proteoform: the observed precursor is the most abundant isotopic peak.
            const double peptideMono = 15000.0;
            double observedMostAbundant = peptideMono + Averagine.GetMostAbundantOffset(peptideMono);

            var acceptor = new MostAbundantMassDiffAcceptor("mostAbundant", new PpmTolerance(5), Averagine);

            // The candidate whose predicted apex lands exactly on the observed peak is the on-apex
            // match → notch 0.
            Assert.That(acceptor.Accepts(observedMostAbundant, peptideMono), Is.EqualTo(0));
        }

        [Test]
        public static void Acceptor_ToleratesApexMisprediction_WithinNotchSet_RejectsBeyond()
        {
            // The averagine apex can miss the experimental apex by ±1–2 neutrons; the notch set
            // (default ±2) deliberately accepts those candidates (with a nonzero notch), while
            // candidates beyond the window are rejected.
            const double peptideMono = 15000.0;
            double observedMostAbundant = peptideMono + Averagine.GetMostAbundantOffset(peptideMono);
            var acceptor = new MostAbundantMassDiffAcceptor("mostAbundant", new PpmTolerance(5), Averagine, maxApexOffsetNeutrons: 2);

            // ±1 and ±2 neutron apex offsets are accepted (nonzero notch, never the -1 sentinel)...
            foreach (int k in new[] { -2, -1, 1, 2 })
            {
                int notch = acceptor.Accepts(observedMostAbundant, peptideMono - k * Constants.C13MinusC12);
                Assert.That(notch, Is.GreaterThanOrEqualTo(0).Or.LessThan(-1), $"k={k} should be accepted with a notch != -1");
            }

            // ...but a +3 neutron offset (beyond the ±2 window) is rejected.
            Assert.That(acceptor.Accepts(observedMostAbundant, peptideMono - 3 * Constants.C13MinusC12), Is.EqualTo(-1));

            Assert.That(acceptor.NumNotches, Is.EqualTo(5));
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
        public static void Acceptor_TheoreticalIntervals_SpanApexPlusMinusTwoNeutrons()
        {
            const double peptideMono = 20000.0;
            var tol = new PpmTolerance(5);
            var acceptor = new MostAbundantMassDiffAcceptor("mostAbundant", tol, Averagine, maxApexOffsetNeutrons: 2);

            double apex = peptideMono + Averagine.GetMostAbundantOffset(peptideMono);

            var intervals = System.Linq.Enumerable.ToList(acceptor.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMono));
            Assert.That(intervals.Count, Is.EqualTo(5));

            // The on-apex interval (notch 0) is centered on the predicted most-abundant mass.
            var onApex = intervals.Find(i => i.Notch == 0);
            Assert.That(onApex, Is.Not.Null);
            Assert.That(onApex.Minimum, Is.EqualTo(tol.GetMinimumValue(apex)).Within(1e-6));
            Assert.That(onApex.Maximum, Is.EqualTo(tol.GetMaximumValue(apex)).Within(1e-6));

            // The extreme intervals sit at apex ± 2 neutrons.
            Assert.That(intervals.Exists(i => i.Contains(apex + 2 * Constants.C13MinusC12)), Is.True);
            Assert.That(intervals.Exists(i => i.Contains(apex - 2 * Constants.C13MinusC12)), Is.True);
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
        public static void DotAcceptor_AcceptsShiftedCandidateAtApex_AndToleratesApex()
        {
            const double peptideMono = 10000.0;
            const double phospho = 79.96633;
            int phosphoNotch = (int)Math.Round(phospho * MassDiffAcceptor.NotchScalar);
            var acc = new MostAbundantDotMassDiffAcceptor("g", new[] { 0.0, phospho }, new PpmTolerance(5), Averagine, 2);

            Assert.That(acc.NumNotches, Is.EqualTo(2));

            // Unmodified candidate at its apex → notch 0.
            double obsUnmod = peptideMono + Averagine.GetMostAbundantOffset(peptideMono);
            Assert.That(acc.Accepts(obsUnmod, peptideMono), Is.EqualTo(0));

            // Phospho-shifted candidate at the apex of (peptide + phospho) → the phospho shift's notch.
            double shiftedMono = peptideMono + phospho;
            double obsPhospho = shiftedMono + Averagine.GetMostAbundantOffset(shiftedMono);
            Assert.That(acc.Accepts(obsPhospho, peptideMono), Is.EqualTo(phosphoNotch));

            // A +1 neutron apex misprediction on the phospho candidate is still accepted under the same notch.
            Assert.That(acc.Accepts(obsPhospho + Constants.C13MinusC12, peptideMono), Is.EqualTo(phosphoNotch));
        }

        [Test]
        public static void GptmdFilterTypes_ResolveToActiveFilters_TomlSettable()
        {
            // Names (toml-serializable) resolve to filter instances; unknown names are ignored.
            var p = new GptmdParameters
            {
                GptmdFilters = new List<IGptmdFilter>(),
                GptmdFilterTypes = new List<string> { nameof(ImprovedScoreFilter), "bogus" }
            };
            var active = p.GetActiveFilters();
            Assert.That(active.Any(f => f is ImprovedScoreFilter), Is.True);
            Assert.That(active.Count, Is.EqualTo(1)); // "bogus" ignored

            // In-memory (GUI) filter + named filter of same type dedupe to one.
            p.GptmdFilters = new List<IGptmdFilter> { new ImprovedScoreFilter() };
            Assert.That(p.GetActiveFilters().Count, Is.EqualTo(1));

            // Default (no names, no in-memory filters) → no filtering, every mod added.
            Assert.That(new GptmdParameters { GptmdFilters = new(), GptmdFilterTypes = new() }.GetActiveFilters(), Is.Empty);
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
