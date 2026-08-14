using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.DatabaseLoading;
using EngineLayer.Gptmd;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
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

        // Averagine most-abundant offset for a monoisotopic mass (diff-to-monoisotopic at the nearest bin).
        private static double ApexOffset(double mono) => Averagine.GetDiffToMonoisotopic(Averagine.GetMostIntenseMassIndex(mono));

        [Test]
        public static void Acceptor_OnApexCandidate_GetsNotchZero()
        {
            // A ~15 kDa proteoform: the observed precursor is the most abundant isotopic peak.
            const double peptideMono = 15000.0;
            double observedMostAbundant = peptideMono + ApexOffset(peptideMono);

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
            double observedMostAbundant = peptideMono + ApexOffset(peptideMono);
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

            double apex = peptideMono + ApexOffset(peptideMono);

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
        public static void SearchTask_GetMassDiffAcceptor_WithMostAbundantType_ReturnsMostAbundantMassDiffAcceptor()
        {
            var acc = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.MostAbundant_Exact, null,
                averagineModel: Averagine);
            Assert.That(acc, Is.TypeOf<MostAbundantMassDiffAcceptor>());
            Assert.That(acc.FileNameAddition, Is.EqualTo("mostAbundant"));

            acc = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.MostAbundant_PlusMinusOne, null,
                averagineModel: Averagine);
            Assert.That(acc, Is.TypeOf<MostAbundantMassDiffAcceptor>());
            Assert.That(acc.FileNameAddition, Is.EqualTo("mostAbundant_1"));

            acc = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.MostAbundant_PlusMinusTwo, null,
                averagineModel: Averagine);
            Assert.That(acc, Is.TypeOf<MostAbundantMassDiffAcceptor>());
            Assert.That(acc.FileNameAddition, Is.EqualTo("mostAbundant_2"));
        }

        /// <summary>
        /// The scan stores observations, not search decisions: PrecursorMass is always the monoisotopic
        /// mass, and PrecursorMostAbundantMass is always the observed apex (0 when none was observed).
        /// The search mode — not the scan — decides which of the two is selected on.
        /// </summary>
        [Test]
        public static void Ms2Scan_PrecursorMasses_AreObservations_SearchModeSelectsBetweenThem()
        {
            var scan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap,
                double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null,
                DissociationType.AnyActivationType, 1, null);

            var monoParams = new CommonParameters();
            var apexParams = new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);

            // No apex observed (e.g. a scan-header precursor): the observation is 0, and BOTH modes select
            // the monoisotopic mass — a most-abundant search cannot match on a peak that was never seen.
            var noEnvelopeScan = new Ms2ScanWithSpecificMass(scan, 1500.0.ToMz(2), 2, "", monoParams);
            Assert.That(noEnvelopeScan.PrecursorMostAbundantMass, Is.EqualTo(0));
            Assert.That(noEnvelopeScan.GetPrecursorMassForSearch(monoParams), Is.EqualTo(noEnvelopeScan.PrecursorMass).Within(1e-9));
            Assert.That(noEnvelopeScan.GetPrecursorMassForSearch(apexParams), Is.EqualTo(noEnvelopeScan.PrecursorMass).Within(1e-9));

            // Apex observed: it is recorded without touching PrecursorMass, which keeps its monoisotopic
            // meaning. Only a most-abundant search selects the apex.
            double apexMass = noEnvelopeScan.PrecursorMass + 5.0;
            var apexScan = new Ms2ScanWithSpecificMass(scan, 1500.0.ToMz(2), 2, "", monoParams,
                precursorMostAbundantMass: apexMass);
            Assert.That(apexScan.PrecursorMass, Is.EqualTo(noEnvelopeScan.PrecursorMass).Within(1e-9));
            Assert.That(apexScan.PrecursorMostAbundantMass, Is.EqualTo(apexMass).Within(1e-9));
            Assert.That(apexScan.GetPrecursorMassForSearch(monoParams), Is.EqualTo(apexScan.PrecursorMass).Within(1e-9));
            Assert.That(apexScan.GetPrecursorMassForSearch(apexParams), Is.EqualTo(apexMass).Within(1e-9));
        }

        /// <summary>
        /// The regression behind nbollis's review: a most-abundant search accepts candidates whose
        /// deconvoluted monoisotopic peak is off by whole isotopologues, so (ScanPrecursorMass -
        /// peptideMonoisotopic) is NOT a chemical mass difference. Any consumer that reads it as one — the
        /// localization engine, the mass-shift histogram, calibration — must go through
        /// GetObservedMonoisotopicMass, or it interprets a ~1-2 Da isotope-assignment offset as a
        /// modification.
        /// </summary>
        [Test]
        public static void GetObservedMonoisotopicMass_RemovesApexAndNeutronOffsets([Values(-2, -1, 0, 1, 2)] int neutronOffset)
        {
            var protein = new Protein("PEPTIDESEQUENCEK", "prot1");
            var peptide = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 16,
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0, null);
            double mono = peptide.MonoisotopicMass;

            var msDataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap,
                double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null,
                DissociationType.AnyActivationType, 1, null);

            // What the detector sees: the apex, mispredicted by k neutrons, plus a little instrument drift.
            // The deconvoluted monoisotopic peak is off by the same k — this is exactly the PSM population
            // most-abundant mode exists to rescue.
            const double driftPpm = 2.0;
            double observedApex = (mono + ApexOffset(mono)) * (1 + driftPpm / 1e6) + neutronOffset * Constants.C13MinusC12;
            double deconvolutedMono = mono * (1 + driftPpm / 1e6) + neutronOffset * Constants.C13MinusC12;

            var apexParams = new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);
            var apexScan = new Ms2ScanWithSpecificMass(msDataScan, deconvolutedMono.ToMz(1), 1, "", apexParams,
                precursorMostAbundantMass: observedApex);
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, apexScan, apexParams, new List<MatchedFragmentIon>());

            // The match is one the apex acceptor really accepts.
            var acceptor = new MostAbundantMassDiffAcceptor("mostAbundant", new PpmTolerance(5), Averagine);
            Assert.That(acceptor.Accepts(observedApex, mono), Is.GreaterThanOrEqualTo(0), "should be an accepted match");

            // Undoing the apex offset and the k neutrons leaves the monoisotopic mass plus only the drift, so
            // the mass difference a downstream engine sees is ~0 ppm — not k * 1.00335 Da.
            double observedMono = psm.GetObservedMonoisotopicMass(mono, apexParams);
            Assert.That((observedMono - mono) / mono * 1e6, Is.EqualTo(driftPpm).Within(0.5));

            // Read raw, that same PSM would show a whole-neutron mass difference — the phantom modification.
            double rawDifference = psm.ScanPrecursorMass - mono;
            Assert.That(rawDifference, Is.EqualTo(neutronOffset * Constants.C13MinusC12).Within(0.01));

            // Monoisotopic mode is the identity — baseline behaviour is untouched.
            Assert.That(psm.GetObservedMonoisotopicMass(mono, new CommonParameters()), Is.EqualTo(psm.ScanPrecursorMass));
        }

        /// <summary>
        /// G-PTM-D cannot use the correction above: with an unknown modification mass in the residual, the
        /// k-neutron apex misprediction can no longer be recovered by rounding — the modification would be
        /// folded into k. Each candidate mass is therefore tested on its own, apex-to-apex, so a real
        /// modification is still found when the apex was mispredicted, and an unrelated one is still rejected.
        /// </summary>
        [Test]
        public static void MatchesCandidateMass_FindsModification_DespiteApexMisprediction()
        {
            const double phospho = 79.96633;
            const double oxidation = 15.99491;

            var protein = new Protein("PEPTIDESEQUENCEK", "prot1");
            var peptide = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 16,
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0, null);
            double mono = peptide.MonoisotopicMass;

            // A phosphorylated proteoform observed at its apex, mispredicted by +1 neutron.
            double phosphoMono = mono + phospho;
            double observedApex = phosphoMono + ApexOffset(phosphoMono) + Constants.C13MinusC12;

            var msDataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap,
                double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null,
                DissociationType.AnyActivationType, 1, null);

            var apexParams = new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);
            var scan = new Ms2ScanWithSpecificMass(msDataScan, phosphoMono.ToMz(1), 1, "", apexParams,
                precursorMostAbundantMass: observedApex);
            var psm = new PeptideSpectralMatch(peptide, 0, 10, 0, scan, apexParams, new List<MatchedFragmentIon>());

            var tolerance = new PpmTolerance(5);

            // The real modification is supported, despite the apex being off by a neutron.
            Assert.That(psm.MatchesCandidateMass(mono + phospho, tolerance, apexParams), Is.True);

            // An unrelated modification is not.
            Assert.That(psm.MatchesCandidateMass(mono + oxidation, tolerance, apexParams), Is.False);

            // In monoisotopic mode the test is the plain monoisotopic comparison.
            var monoParams = new CommonParameters();
            Assert.That(psm.MatchesCandidateMass(psm.ScanPrecursorMass, tolerance, monoParams), Is.True);
            Assert.That(psm.MatchesCandidateMass(psm.ScanPrecursorMass + 10, tolerance, monoParams), Is.False);
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
            double obsUnmod = peptideMono + ApexOffset(peptideMono);
            Assert.That(acc.Accepts(obsUnmod, peptideMono), Is.EqualTo(0));

            // Phospho-shifted candidate at the apex of (peptide + phospho) → the phospho shift's notch.
            double shiftedMono = peptideMono + phospho;
            double obsPhospho = shiftedMono + ApexOffset(shiftedMono);
            Assert.That(acc.Accepts(obsPhospho, peptideMono), Is.EqualTo(phosphoNotch));

            // A +1 neutron apex misprediction on the phospho candidate is still accepted under the same notch.
            Assert.That(acc.Accepts(obsPhospho + Constants.C13MinusC12, peptideMono), Is.EqualTo(phosphoNotch));
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

        [Test]
        public static void MostAbundantMassErrorPpm_TightOnApex_NullInMonoisotopicMode()
        {
            var protein = new Protein("PEPTIDESEQUENCEK", "prot1");
            var peptide = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 16,
                CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0, null);
            double mono = peptide.MonoisotopicMass;
            double apex = mono + ApexOffset(mono); // observed most-abundant lands exactly on the predicted apex

            var msDataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap,
                double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null,
                DissociationType.AnyActivationType, 1, null);
            var noIons = new List<MatchedFragmentIon>();

            // Most-abundant mode: the scan observes the apex; the monoisotopic PrecursorMz is left as mono.
            // Reported most-abundant error is ~0 (on apex).
            var maParams = new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);
            var maScan = new Ms2ScanWithSpecificMass(msDataScan, mono.ToMz(1), 1, "", maParams, precursorMostAbundantMass: apex);
            var maPsm = new PeptideSpectralMatch(peptide, 0, 10, 0, maScan, maParams, noIons);

            Assert.That(maPsm.MostAbundantMassErrorPpm, Is.Not.Null);
            Assert.That(maPsm.MostAbundantMassErrorPpm[0], Is.EqualTo(0).Within(0.5)); // tight, in ppm

            // Monoisotopic mode: no most-abundant error is reported (column omitted).
            var monoParams = new CommonParameters();
            var monoScan = new Ms2ScanWithSpecificMass(msDataScan, mono.ToMz(1), 1, "", monoParams);
            var monoPsm = new PeptideSpectralMatch(peptide, 0, 10, 0, monoScan, monoParams, noIons);
            Assert.That(monoPsm.MostAbundantMassErrorPpm, Is.Null);

            // The optional columns appear in the header only when requested.
            Assert.That(SpectralMatch.GetTabSeparatedHeader(includeMostAbundantColumn: true),
                Does.Contain("Most Abundant Mass Diff (ppm)"));
            Assert.That(SpectralMatch.GetTabSeparatedHeader(includeMostAbundantColumn: true),
                Does.Contain("Precursor Most Abundant Mass"));
            Assert.That(SpectralMatch.GetTabSeparatedHeader(), Does.Not.Contain("Most Abundant Mass Diff (ppm)"));
            Assert.That(SpectralMatch.GetTabSeparatedHeader(), Does.Not.Contain("Precursor Most Abundant Mass"));

            // End-to-end writer path: the emitted row cells align under the new header columns and
            // carry the reported values (row and header stay column-aligned).
            var header = SpectralMatch.GetTabSeparatedHeader(includeMostAbundantColumn: true).Split('\t');
            var row = maPsm.ToString(new Dictionary<string, int>(), false, false, false, includeMostAbundantColumn: true).Split('\t');
            int col = Array.IndexOf(header, "Most Abundant Mass Diff (ppm)");
            Assert.That(col, Is.GreaterThanOrEqualTo(0));
            Assert.That(double.Parse(row[col], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(0).Within(0.5));

            // The absolute observed-apex column carries the same mass the search matched on.
            int massCol = Array.IndexOf(header, "Precursor Most Abundant Mass");
            Assert.That(massCol, Is.GreaterThanOrEqualTo(0));
            Assert.That(double.Parse(row[massCol], System.Globalization.CultureInfo.InvariantCulture),
                Is.EqualTo(maScan.PrecursorMostAbundantMass).Within(1e-4));
        }

        [Test]
        public static void Acceptor_ObservedIntervals_RecoverCandidateMonoWindows()
        {
            // The indexed (ModernSearch) path: GetAllowedPrecursorMassIntervalsFromObservedMass converts
            // an observed most-abundant mass back to candidate-monoisotopic windows, one per apex notch.
            const double peptideMono = 20000.0;
            var tol = new PpmTolerance(5);
            var acceptor = new MostAbundantMassDiffAcceptor("mostAbundant", tol, Averagine, maxApexOffsetNeutrons: 2);

            double observedApex = peptideMono + ApexOffset(peptideMono);
            var intervals = acceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(observedApex).ToList();
            Assert.That(intervals.Count, Is.EqualTo(5));

            // The on-apex window (notch 0) recovers the candidate's true monoisotopic mass.
            var onApex = intervals.Find(i => i.Notch == 0);
            Assert.That(onApex, Is.Not.Null);
            Assert.That(onApex.Contains(peptideMono), Is.True);
        }

        [Test]
        public static void DotAcceptor_Intervals_CoverShifts_AndRejectsFarCandidate()
        {
            const double peptideMono = 10000.0;
            const double phospho = 79.96633;
            var tol = new PpmTolerance(5);
            var acc = new MostAbundantDotMassDiffAcceptor("g", new[] { 0.0, phospho }, tol, Averagine, maxApexOffsetNeutrons: 1);
            int phosphoNotch = (int)Math.Round(phospho * MassDiffAcceptor.NotchScalar);

            // 2 shifts × (2·1+1) apex offsets = 6 theoretical windows; the phospho apex is covered.
            var theo = acc.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMono).ToList();
            Assert.That(theo.Count, Is.EqualTo(6));
            double shiftedApex = (peptideMono + phospho) + ApexOffset(peptideMono + phospho);
            Assert.That(theo.Exists(i => i.Notch == phosphoNotch && i.Contains(shiftedApex)), Is.True);

            // The observed-side conversion recovers the unmodified candidate's mono under notch 0.
            var obs = acc.GetAllowedPrecursorMassIntervalsFromObservedMass(peptideMono + ApexOffset(peptideMono)).ToList();
            Assert.That(obs.Count, Is.EqualTo(6));
            Assert.That(obs.Exists(i => i.Notch == 0 && i.Contains(peptideMono)), Is.True);

            // A candidate far from every shifted apex is rejected.
            Assert.That(acc.Accepts(peptideMono + 500.0, peptideMono), Is.EqualTo(-1));
        }

        [Test]
        public static void Acceptors_ExposeMaxApexOffset_ToString_AndProseString()
        {
            var mass = new MostAbundantMassDiffAcceptor("mostAbundant", new PpmTolerance(5), Averagine, maxApexOffsetNeutrons: 2);
            Assert.That(mass.MaxApexOffsetNeutrons, Is.EqualTo(2));
            Assert.That(mass.ToString(), Does.Contain("mostAbundant"));
            Assert.That(mass.ToProseString(), Does.Contain("most abundant"));

            var dot = new MostAbundantDotMassDiffAcceptor("g", new[] { 0.0, 79.96633 }, new PpmTolerance(5), Averagine, maxApexOffsetNeutrons: 2);
            Assert.That(dot.MaxApexOffsetNeutrons, Is.EqualTo(2));
            Assert.That(dot.ToString(), Does.Contain("mostAbundantDot"));
            Assert.That(dot.ToProseString(), Does.Contain("most abundant"));
        }

        [Test]
        public static void Acceptors_NonPositiveCandidateMass_UseZeroApexOffset_DoNotThrow()
        {
            // A candidate monoisotopic mass <= 0 must take the guarded ApexOffset branch (return 0,
            // not index the averagine model out of range). Routed through Accepts(scan, peptideMass).
            var mass = new MostAbundantMassDiffAcceptor("mostAbundant", new PpmTolerance(5), Averagine);
            Assert.That(mass.Accepts(0.0, 0.0), Is.EqualTo(0).Or.EqualTo(-1));

            var dot = new MostAbundantDotMassDiffAcceptor("g", new[] { 0.0 }, new PpmTolerance(5), Averagine);
            Assert.That(dot.Accepts(0.0, 0.0), Is.EqualTo(0).Or.EqualTo(-1));
        }

        [Test]
        public static void GptmdRunsInMostAbundantMode()
        {
            // GPTMD end to end in most-abundant mode: the notch/prose acceptors switch to the apex
            // (MostAbundantDotMassDiffAcceptor) variant. Deconvolution is on by default for GPTMD, so
            // only the match-mode needs flipping.
            var gptmdTask = new GptmdTask();
            gptmdTask.CommonParameters.PrecursorMassMatchMode = PrecursorMassMatchMode.MostAbundant;

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestGptmdMostAbundant");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            Directory.CreateDirectory(outputFolder);

            Assert.DoesNotThrow(() =>
                gptmdTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) },
                    new List<string> { myFile }, "test"));
            Assert.That(Directory.GetFiles(outputFolder, "*", SearchOption.AllDirectories).Length, Is.GreaterThan(0));

            Directory.Delete(outputFolder, true);
            string taskSettings = Path.Combine(TestContext.CurrentContext.TestDirectory, "Task Settings");
            if (Directory.Exists(taskSettings)) { Directory.Delete(taskSettings, true); }
        }

        [Test]
        public static void GetMs2Scans_RecordsApexObservation_RegardlessOfSearchMode()
        {
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            var fileManager = new MyFileManager(true);

            var monoParams = new CommonParameters(maxThreadsToUsePerFile: 1, doPrecursorDeconvolution: true,
                precursorMassMatchMode: PrecursorMassMatchMode.Monoisotopic);
            var apexParams = new CommonParameters(maxThreadsToUsePerFile: 1, doPrecursorDeconvolution: true,
                precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);

            var msDataFile = fileManager.LoadFile(origDataFile, monoParams);

            var monoScans = MetaMorpheusTask.GetMs2Scans(msDataFile, origDataFile, monoParams).ToArray();
            var apexScans = MetaMorpheusTask.GetMs2Scans(msDataFile, origDataFile, apexParams).ToArray();

            Assert.That(apexScans.Length, Is.GreaterThan(0));

            // The apex is an observation of the deconvoluted envelope, so it is recorded in BOTH modes —
            // the scan does not change shape depending on the search being run.
            Assert.That(monoScans.Any(s => s.PrecursorMostAbundantMass > 0), Is.True);
            Assert.That(apexScans.Any(s => s.PrecursorMostAbundantMass > 0), Is.True);
            // And in both modes PrecursorMass keeps its monoisotopic meaning.
            Assert.That(monoScans.Select(s => s.PrecursorMass), Is.EqualTo(apexScans.Select(s => s.PrecursorMass)).Within(1e-6));

            // Only the search mode changes what gets matched: a most-abundant search selects a mass that
            // departs from the monoisotopic one, a monoisotopic search never does.
            Assert.That(apexScans.All(s => Math.Abs(s.GetPrecursorMassForSearch(monoParams) - s.PrecursorMass) < 1e-6), Is.True);
            Assert.That(apexScans.Any(s => Math.Abs(s.GetPrecursorMassForSearch(apexParams) - s.PrecursorMass) > 1e-6), Is.True);
        }

        /// <summary>
        /// End-to-end proof that ClassicSearchEngine selects candidates by the APEX in most-abundant mode.
        /// <para>
        /// ClassicSearchEngine:149 passes <c>specificBioPolymer.MonoisotopicMass</c> into GetAcceptableScans,
        /// which reads as "still monoisotopic". It is not: that value is the INPUT to the theory-side apex
        /// conversion, not a bypass of it. GetAcceptableScans (:253) hands it to the acceptor's
        /// GetAllowedPrecursorMassIntervalsFromTheoreticalMass, which returns windows centred on
        /// mono + averagineApexOffset (MostAbundantMassDiffAcceptor:120-128); those windows are then compared
        /// against MyScanPrecursorMasses (:258, :268, :276), which was built from GetPrecursorMassForSearch
        /// (:46) — the observed APEX. So the comparison is apex-to-apex at both ends.
        /// </para>
        /// <para>
        /// This test forces the distinction: the scan is given a deliberately WRONG monoisotopic mass
        /// (+3.5 Da, far outside tolerance and not a whole number of neutrons) and a CORRECT apex. A
        /// monoisotopic search must therefore find nothing; a most-abundant search must find the peptide.
        /// If selection were really driven by the monoisotopic mass, the most-abundant search would return
        /// null too.
        /// </para>
        /// </summary>
        [Test]
        public static void ClassicSearch_SelectsScansByApex_NotMonoisotopic()
        {
            // ~9.6 kDa tryptic peptide: big enough that the averagine apex is several neutrons above the
            // monoisotopic mass, so "apex" and "monoisotopic" are unambiguously different numbers.
            const string sequence = "PEPTIDEASYGLVNFQWMCH" + "PEPTIDEASYGLVNFQWMCH"
                                  + "PEPTIDEASYGLVNFQWMCH" + "PEPTIDEASYGLVNFQWMCH"
                                  + "PEPTIDE" + "K";
            var protein = new Protein(sequence, "accession");
            var digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5, maxPeptideLength: 200);
            var peptide = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>(), null, null).First();

            double mono = peptide.MonoisotopicMass;
            double apex = mono + ApexOffset(mono);
            Assert.That(apex - mono, Is.GreaterThan(1.0), "test is only meaningful if the apex differs from the monoisotopic mass");

            // Theoretical fragments, handed to the scan directly as neutral-mass envelopes so that fragment
            // matching is guaranteed and the test isolates PRECURSOR selection.
            var products = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            var fragmentEnvelopes = products.Where(p => !double.IsNaN(p.NeutralMass))
                .Select(p => new IsotopicEnvelope(p.NeutralMass, 1.0, 1)).ToArray();

            var spectrum = new MzSpectrum(
                products.Where(p => !double.IsNaN(p.NeutralMass)).Select(p => p.NeutralMass.ToMz(1)).OrderBy(m => m).ToArray(),
                products.Where(p => !double.IsNaN(p.NeutralMass)).Select(p => 1.0).ToArray(), false);
            var dataScan = new MsDataScan(spectrum, 2, 2, true, Polarity.Positive, 1.0, new MzRange(100, 12000), "f",
                MZAnalyzerType.Orbitrap, spectrum.SumOfAllY, null, null, "scan=2", 500.0, null, null, 500.0, null,
                DissociationType.HCD, 1, null);

            // The observation: monoisotopic mass deliberately wrong by 3.5 Da, apex exactly right.
            const double monoError = 3.5;
            var scan = new Ms2ScanWithSpecificMass(dataScan, (mono + monoError).ToMz(1), 1, "f.mzML",
                new CommonParameters(), neutralExperimentalFragments: fragmentEnvelopes,
                precursorMostAbundantMass: apex);

            Assert.That(scan.PrecursorMass, Is.EqualTo(mono + monoError).Within(1e-6));
            Assert.That(scan.PrecursorMostAbundantMass, Is.EqualTo(apex).Within(1e-6));

            var scans = new[] { scan };
            var proteinList = new List<Protein> { protein };
            var tol = new PpmTolerance(5);

            var monoParams = new CommonParameters(digestionParams: digestionParams, scoreCutoff: 1,
                precursorMassTolerance: tol, dissociationType: DissociationType.HCD,
                precursorMassMatchMode: PrecursorMassMatchMode.Monoisotopic);
            var apexParams = new CommonParameters(digestionParams: digestionParams, scoreCutoff: 1,
                precursorMassTolerance: tol, dissociationType: DissociationType.HCD,
                precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);

            var fsp = new List<(string FileName, CommonParameters Parameters)> { ("f.mzML", monoParams) };

            SpectralMatch[] monoPsms = new PeptideSpectralMatch[1];
            new ClassicSearchEngine(monoPsms, scans, new List<Modification>(), new List<Modification>(), null, null, null,
                proteinList, new SinglePpmAroundZeroSearchMode(5), monoParams, fsp, null, new List<string>(), false).Run();

            SpectralMatch[] apexPsms = new PeptideSpectralMatch[1];
            new ClassicSearchEngine(apexPsms, scans, new List<Modification>(), new List<Modification>(), null, null, null,
                proteinList, new MostAbundantMassDiffAcceptor("mostAbundant", tol, Averagine), apexParams, fsp, null,
                new List<string>(), false).Run();

            // The monoisotopic search cannot reach this scan: its monoisotopic mass is 3.5 Da off.
            Assert.That(monoPsms[0], Is.Null,
                "monoisotopic search selected a scan whose monoisotopic mass is 3.5 Da from the candidate");

            // The most-abundant search finds it, because selection used the observed apex on the scan side
            // and mono + averagineApexOffset on the theory side.
            Assert.That(apexPsms[0], Is.Not.Null,
                "most-abundant search failed to select the scan by its apex — selection is still monoisotopic");
            Assert.That(apexPsms[0].Notch, Is.EqualTo(0), "on-apex match should carry notch 0");
            apexPsms[0].ResolveAllAmbiguities();
            Assert.That(apexPsms[0].BaseSequence, Is.EqualTo(peptide.BaseSequence));

            // CONTROL: the monoisotopic search is not broken — give it a scan whose monoisotopic mass IS
            // correct and it finds the same peptide. So the null above is caused by the mass it matched on,
            // not by anything else in this setup (protease, score cutoff, fragment matching).
            var correctMonoScan = new Ms2ScanWithSpecificMass(dataScan, mono.ToMz(1), 1, "f.mzML",
                new CommonParameters(), neutralExperimentalFragments: fragmentEnvelopes,
                precursorMostAbundantMass: apex);
            SpectralMatch[] controlPsms = new PeptideSpectralMatch[1];
            new ClassicSearchEngine(controlPsms, new[] { correctMonoScan }, new List<Modification>(), new List<Modification>(),
                null, null, null, proteinList, new SinglePpmAroundZeroSearchMode(5), monoParams, fsp, null,
                new List<string>(), false).Run();
            Assert.That(controlPsms[0], Is.Not.Null, "control: monoisotopic search should find a correct monoisotopic mass");
            controlPsms[0].ResolveAllAmbiguities();
            Assert.That(controlPsms[0].BaseSequence, Is.EqualTo(peptide.BaseSequence));

            // NEGATIVE CONTROL: the most-abundant search is not simply permissive. Move the observed apex
            // 7.5 Da away — beyond the ±2-neutron notch set (±2.007 Da) and not a whole number of neutrons,
            // so it cannot be rescued by any k — and the same search now finds nothing.
            var wrongApexScan = new Ms2ScanWithSpecificMass(dataScan, (mono + monoError).ToMz(1), 1, "f.mzML",
                new CommonParameters(), neutralExperimentalFragments: fragmentEnvelopes,
                precursorMostAbundantMass: apex + 7.5);
            SpectralMatch[] negativePsms = new PeptideSpectralMatch[1];
            new ClassicSearchEngine(negativePsms, new[] { wrongApexScan }, new List<Modification>(), new List<Modification>(),
                null, null, null, proteinList, new MostAbundantMassDiffAcceptor("mostAbundant", tol, Averagine),
                apexParams, fsp, null, new List<string>(), false).Run();
            Assert.That(negativePsms[0], Is.Null,
                "negative control: most-abundant search accepted an apex 7.5 Da outside every notch");
        }

        /// <summary>
        /// Q1 (review question): does TRADITIONAL G-PTM-D use missed-monoisotopic notches?
        /// <para>
        /// <b>No.</b> GptmdTask.GetAcceptableMassShifts (:282-291) builds the shift set from exactly four
        /// sources — G-PTM-D mod masses, mod-replacement differences, combos, and a single zero notch — and
        /// none of them is a neutron offset. The resulting <see cref="DotMassDiffAcceptor"/> therefore places
        /// one tight-ppm window per shift, with no isotope-level tolerance anywhere, including on the zero
        /// notch. A deconvolution that misassigns the monoisotopic peak by one neutron is simply lost.
        /// </para>
        /// </summary>
        [Test]
        public static void Q1_TraditionalGptmd_HasNoMissedMonoisotopicNotches()
        {
            const double peptideMono = 1200.0;
            const double oxidation = 15.994915;
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            var gptmdMods = new List<Modification>
            {
                new Modification(_originalId: "ox", _modificationType: "mt", _target: motifM,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: oxidation)
            };

            var shifts = GptmdTask.GetAcceptableMassShifts(new List<Modification>(), new List<Modification>(),
                gptmdMods, new List<Tuple<double, double>>()).ToList();

            // The shift set is chemistry only: {0, +15.994915}. Nothing near a neutron.
            Assert.That(shifts.Count, Is.EqualTo(2));
            Assert.That(shifts.Any(s => Math.Abs(s) < 1e-9), Is.True, "shift set should contain the zero notch");
            Assert.That(shifts.Any(s => Math.Abs(s - oxidation) < 1e-9), Is.True, "shift set should contain the mod mass");
            Assert.That(shifts.Any(s => Math.Abs(Math.Abs(s) - Constants.C13MinusC12) < 0.05), Is.False,
                "the traditional G-PTM-D shift set contains a neutron-sized offset");

            var acceptor = new DotMassDiffAcceptor("gptmd", shifts, new PpmTolerance(10));

            // POSITIVE CONTROLS — the two real shifts are accepted at tight ppm.
            Assert.That(acceptor.Accepts(peptideMono, peptideMono), Is.GreaterThanOrEqualTo(0),
                "positive control: unmodified peptide should match the zero notch");
            Assert.That(acceptor.Accepts(peptideMono + oxidation, peptideMono), Is.GreaterThanOrEqualTo(0),
                "positive control: oxidised peptide should match the oxidation notch");

            // NEGATIVE CONTROLS — one neutron off either shift is rejected outright. This is the answer to
            // Q1: there is no missed-monoisotopic tolerance, on the mod notch or on the zero notch.
            Assert.That(acceptor.Accepts(peptideMono + Constants.C13MinusC12, peptideMono), Is.EqualTo(-1),
                "traditional G-PTM-D accepted a +1-neutron precursor on the zero notch");
            Assert.That(acceptor.Accepts(peptideMono + oxidation + Constants.C13MinusC12, peptideMono), Is.EqualTo(-1),
                "traditional G-PTM-D accepted a +1-neutron precursor on the oxidation notch");
            Assert.That(acceptor.Accepts(peptideMono - Constants.C13MinusC12, peptideMono), Is.EqualTo(-1));
        }

        /// <summary>
        /// Q2 (review question): does MOST-ABUNDANT G-PTM-D use the analogous notches?
        /// <para>
        /// <b>Yes, but they are not the same thing.</b> <see cref="MostAbundantDotMassDiffAcceptor"/> emits
        /// k ∈ {0, ±1, ±2} apex windows around <i>every</i> shift, so five windows per shift where the
        /// traditional acceptor has one. They compensate for a different error: a missed-monoisotopic notch
        /// covers the DECONVOLUTION picking the wrong isotopologue as monoisotopic (an observation-side
        /// error), whereas the apex notches cover AVERAGINE mispredicting which isotopologue is tallest (a
        /// theory-side error). All five windows for a shift share that shift's notch, so per-shift FDR
        /// grouping is preserved.
        /// </para>
        /// </summary>
        [Test]
        public static void Q2_MostAbundantGptmd_HasApexNotchesAroundEveryShift()
        {
            const double peptideMono = 9500.0;
            const double oxidation = 15.994915;
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            var gptmdMods = new List<Modification>
            {
                new Modification(_originalId: "ox", _modificationType: "mt", _target: motifM,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: oxidation)
            };
            var shifts = GptmdTask.GetAcceptableMassShifts(new List<Modification>(), new List<Modification>(),
                gptmdMods, new List<Tuple<double, double>>()).ToList();

            var tolerance = new PpmTolerance(10);
            var apexAcceptor = new MostAbundantDotMassDiffAcceptor("gptmd", shifts, tolerance, Averagine);
            var traditional = new DotMassDiffAcceptor("gptmd", shifts, tolerance);
            double spacing = Constants.C13MinusC12;

            foreach (double shift in new[] { 0.0, oxidation })
            {
                double shiftedMono = peptideMono + shift;
                double apex = shiftedMono + ApexOffset(shiftedMono);
                int expectedNotch = apexAcceptor.Accepts(apex, peptideMono);

                // POSITIVE CONTROLS — on-apex and ±1, ±2 neutrons are all accepted, and all report the SAME
                // notch, because they are one modification hypothesis differing only by apex misprediction.
                Assert.That(expectedNotch, Is.GreaterThanOrEqualTo(0),
                    $"positive control: on-apex candidate for shift {shift} should be accepted");
                foreach (int k in new[] { -2, -1, 1, 2 })
                {
                    Assert.That(apexAcceptor.Accepts(apex + k * spacing, peptideMono), Is.EqualTo(expectedNotch),
                        $"shift {shift}: k={k} should be accepted under the same notch as k=0");
                }

                // NEGATIVE CONTROL — ±3 neutrons is outside the set and must be rejected.
                foreach (int k in new[] { -3, 3 })
                {
                    Assert.That(apexAcceptor.Accepts(apex + k * spacing, peptideMono), Is.EqualTo(-1),
                        $"shift {shift}: k={k} is outside the ±2 apex set and must be rejected");
                }

                // THE ASYMMETRY, stated as an assertion: the traditional acceptor rejects exactly the
                // ±1-neutron case that the most-abundant acceptor accepts.
                Assert.That(traditional.Accepts(shiftedMono + spacing, peptideMono), Is.EqualTo(-1));
                Assert.That(apexAcceptor.Accepts(apex + spacing, peptideMono), Is.GreaterThanOrEqualTo(0));
            }

            // NumNotches counts SHIFTS, not windows — the ±k windows do not inflate the FDR grouping.
            Assert.That(apexAcceptor.NumNotches, Is.EqualTo(shifts.Count));
            Assert.That(apexAcceptor.NumNotches, Is.EqualTo(traditional.NumNotches));
        }

        /// <summary>
        /// The error source the Q1/Q2 asymmetry introduces, characterised rather than asserted away.
        /// <para>
        /// Because the zero notch also carries ±2 apex windows, and <see cref="MostAbundantDotMassDiffAcceptor.Accepts"/>
        /// returns the FIRST matching shift in |shift|-ascending order (so the zero notch is tested first), a
        /// real modification whose mass is close to a whole number of neutrons can be absorbed by the zero
        /// notch instead of being reported as that modification. Deamidation (+0.98402) sits 0.01934 Da from
        /// one neutron (1.00336), so whether the two are separable depends entirely on the peptide mass:
        /// </para>
        /// <list type="bullet">
        /// <item>~1200 Da at 10 ppm (±0.012 Da): 0.01934 Da is OUTSIDE tolerance → correctly reported as deamidation.</item>
        /// <item>~9500 Da at 10 ppm (±0.095 Da): 0.01934 Da is INSIDE tolerance → absorbed by the zero notch.</item>
        /// </list>
        /// <para>
        /// Crucially, the TRADITIONAL acceptor is not safe here either — it fails in the opposite direction.
        /// At 9500 Da its deamidation window is ±0.095 Da wide, so a precursor off by exactly one neutron (a
        /// pure deconvolution isotope error, no chemistry) lands inside it and is reported AS deamidation.
        /// So at high mass: traditional G-PTM-D invents a modification that is not there (false positive),
        /// most-abundant G-PTM-D absorbs one that is (false negative). Neither is introduced by most-abundant
        /// mode — at this mass and tolerance the two masses are not separable at all, and the mode only
        /// decides which way the ambiguity falls. Documented here so the boundary is known rather than
        /// discovered.
        /// </para>
        /// </summary>
        [Test]
        public static void ApexNotches_CanAbsorbANearNeutronMod_AtHighMassOnly()
        {
            const double deamidation = 0.984016;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            var gptmdMods = new List<Modification>
            {
                new Modification(_originalId: "deam", _modificationType: "mt", _target: motifN,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: deamidation)
            };
            var shifts = GptmdTask.GetAcceptableMassShifts(new List<Modification>(), new List<Modification>(),
                gptmdMods, new List<Tuple<double, double>>()).ToList();
            var tolerance = new PpmTolerance(10);
            var acceptor = new MostAbundantDotMassDiffAcceptor("gptmd", shifts, tolerance, Averagine);

            // The zero notch sorts first, so it is tested before the deamidation notch.
            Assert.That(shifts.First(), Is.EqualTo(0).Within(1e-9));
            int zeroNotch = (int)Math.Round(0.0 * MassDiffAcceptor.NotchScalar);
            int deamidationNotch = (int)Math.Round(deamidation * MassDiffAcceptor.NotchScalar);
            Assert.That(deamidationNotch, Is.Not.EqualTo(zeroNotch));

            // The gap that decides it.
            double gap = Constants.C13MinusC12 - deamidation;
            Assert.That(gap, Is.EqualTo(0.01934).Within(1e-4));

            // POSITIVE CONTROL — small peptide: the gap exceeds tolerance, so the zero notch's k=+1 window
            // does NOT reach the deamidated apex, and deamidation is correctly identified.
            const double smallMono = 1200.0;
            double smallDeamidatedApex = (smallMono + deamidation) + ApexOffset(smallMono + deamidation);
            Assert.That(tolerance.GetMaximumValue(smallMono) - smallMono, Is.LessThan(gap),
                "premise: at 1200 Da the tolerance must be tighter than the deamidation/neutron gap");
            Assert.That(acceptor.Accepts(smallDeamidatedApex, smallMono), Is.EqualTo(deamidationNotch),
                "at low mass a deamidated precursor should be reported as deamidation, not absorbed by notch 0");

            // NEGATIVE CONTROL (the hazard) — large proteoform: the same gap is now inside tolerance, the
            // zero notch matches first, and the modification is never attributed.
            const double largeMono = 9500.0;
            double largeDeamidatedApex = (largeMono + deamidation) + ApexOffset(largeMono + deamidation);
            Assert.That(tolerance.GetMaximumValue(largeMono) - largeMono, Is.GreaterThan(gap),
                "premise: at 9500 Da the tolerance must be looser than the deamidation/neutron gap");
            Assert.That(acceptor.Accepts(largeDeamidatedApex, largeMono), Is.EqualTo(zeroNotch),
                "expected the zero notch to absorb the near-neutron mod at high mass");

            // AND THE TRADITIONAL ACCEPTOR IS NOT SAFE EITHER — it fails in the opposite direction.
            // At 9500 Da the deamidation window is ±0.095 Da wide, so a precursor that is off by exactly one
            // neutron (a pure deconvolution isotope error, no chemistry at all) lands inside it and is
            // reported AS deamidation. Traditional G-PTM-D therefore invents a modification that is not
            // there (false positive); most-abundant G-PTM-D loses one that is (false negative). Neither is
            // introduced by most-abundant mode: at this mass and tolerance the two masses are simply not
            // separable, and the mode only decides which way the ambiguity falls.
            var traditional = new DotMassDiffAcceptor("gptmd", shifts, tolerance);
            Assert.That(traditional.Accepts(largeMono + deamidation, largeMono), Is.EqualTo(deamidationNotch),
                "positive control: traditional G-PTM-D should identify a genuinely deamidated precursor");
            Assert.That(traditional.Accepts(largeMono + Constants.C13MinusC12, largeMono), Is.EqualTo(deamidationNotch),
                "at high mass traditional G-PTM-D reports a pure +1-neutron isotope error as deamidation");

            // NEGATIVE CONTROL for that claim — at low mass the traditional acceptor correctly rejects the
            // same +1-neutron precursor, confirming the failure above is a mass-scale effect and not a
            // property of the acceptor.
            Assert.That(traditional.Accepts(smallMono + Constants.C13MinusC12, smallMono), Is.EqualTo(-1),
                "at low mass a +1-neutron precursor should not be mistaken for deamidation");
        }

        /// <summary>
        /// The sharpest form of the apex-notch error source: the ±k windows can make two DIFFERENT
        /// modifications collide, not just a modification and "unmodified".
        /// <para>
        /// Acetylation (+42.010565) and carbamylation (+43.005814) are 0.995 Da apart — a full Dalton, and
        /// completely unambiguous under the traditional acceptor, which places one tight window per shift.
        /// But acetyl's <c>k = +1</c> apex window sits at 43.014, only <b>0.0081 Da</b> from carbamyl. Since
        /// shifts are tested in |shift|-ascending order, acetyl is reached first, so above ~811 Da (at
        /// 10 ppm) a genuinely carbamylated precursor is reported as acetylated.
        /// </para>
        /// <para>
        /// This collision is created by the notch set: it does not exist without the ±k windows. It is the
        /// concrete reason NOT to add missed-monoisotopic notches to traditional G-PTM-D, and a known
        /// limitation of the most-abundant G-PTM-D acceptor, which already carries them.
        /// </para>
        /// </summary>
        [Test]
        public static void ApexNotches_CanConfuseTwoModifications_AcetylVersusCarbamyl()
        {
            const double acetyl = 42.010565;
            const double carbamyl = 43.005814;
            ModificationMotif.TryGetMotif("K", out ModificationMotif motifK);
            var gptmdMods = new List<Modification>
            {
                new Modification(_originalId: "acetyl", _modificationType: "mt", _target: motifK,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: acetyl),
                new Modification(_originalId: "carbamyl", _modificationType: "mt", _target: motifK,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: carbamyl)
            };
            var shifts = GptmdTask.GetAcceptableMassShifts(new List<Modification>(), new List<Modification>(),
                gptmdMods, new List<Tuple<double, double>>()).ToList();
            var tolerance = new PpmTolerance(10);

            int acetylNotch = (int)Math.Round(acetyl * MassDiffAcceptor.NotchScalar);
            int carbamylNotch = (int)Math.Round(carbamyl * MassDiffAcceptor.NotchScalar);

            // The two mods are a full Dalton apart; the collision is manufactured by the +1 neutron window.
            Assert.That(carbamyl - acetyl, Is.EqualTo(0.99525).Within(1e-5));
            double collisionGap = Math.Abs((acetyl + Constants.C13MinusC12) - carbamyl);
            Assert.That(collisionGap, Is.EqualTo(0.00811).Within(1e-4));

            // TRADITIONAL — no ±k windows, so the two mods never collide, at any mass.
            var traditional = new DotMassDiffAcceptor("gptmd", shifts, tolerance);
            foreach (double mono in new[] { 1200.0, 9500.0 })
            {
                Assert.That(traditional.Accepts(mono + acetyl, mono), Is.EqualTo(acetylNotch),
                    "positive control: traditional should identify acetylation");
                Assert.That(traditional.Accepts(mono + carbamyl, mono), Is.EqualTo(carbamylNotch),
                    $"traditional confused carbamyl with acetyl at {mono} Da");
            }

            var apexAcceptor = new MostAbundantDotMassDiffAcceptor("gptmd", shifts, tolerance, Averagine);

            // The crossover is set by the tolerance evaluated at the WINDOW mass (mono + ~43 Da), not at the
            // peptide mass, so it lands near 770 Da at 10 ppm — below essentially every real peptide.
            const double smallMono = 500.0;
            Assert.That(tolerance.GetMaximumValue(smallMono + carbamyl) - (smallMono + carbamyl),
                Is.LessThan(collisionGap),
                "premise: at 500 Da the tolerance must be tighter than the acetyl+1n / carbamyl gap");

            // POSITIVE CONTROL — a genuinely acetylated precursor is identified as acetylation.
            double smallAcetylApex = (smallMono + acetyl) + ApexOffset(smallMono + acetyl);
            Assert.That(apexAcceptor.Accepts(smallAcetylApex, smallMono), Is.EqualTo(acetylNotch));

            // NEGATIVE CONTROL — below the crossover the two mods are still resolved correctly.
            double smallCarbamylApex = (smallMono + carbamyl) + ApexOffset(smallMono + carbamyl);
            Assert.That(apexAcceptor.Accepts(smallCarbamylApex, smallMono), Is.EqualTo(carbamylNotch),
                "below the crossover carbamylation should still be identified correctly");

            // THE HAZARD — above the crossover, a genuinely CARBAMYLATED precursor is reported as ACETYL,
            // because acetyl sorts first and its k=+1 window reaches it. Demonstrated here at a routine
            // TRYPTIC peptide mass, not only at proteoform scale.
            const double largeMono = 1500.0;
            Assert.That(tolerance.GetMaximumValue(largeMono + carbamyl) - (largeMono + carbamyl),
                Is.GreaterThan(collisionGap),
                "premise: at 1500 Da the tolerance must be looser than the collision gap");
            double largeCarbamylApex = (largeMono + carbamyl) + ApexOffset(largeMono + carbamyl);
            Assert.That(apexAcceptor.Accepts(largeCarbamylApex, largeMono), Is.EqualTo(acetylNotch),
                "expected acetyl's +1-neutron window to capture the carbamylated precursor");
            Assert.That(apexAcceptor.Accepts(largeCarbamylApex, largeMono), Is.Not.EqualTo(carbamylNotch));
        }

        /// <summary>
        /// PROTOTYPE MEASUREMENT — does resolving apex collisions in favour of the on-apex (k = 0)
        /// hypothesis improve G-PTM-D, and what does it cost in depth?
        /// <para>
        /// Today <see cref="MostAbundantDotMassDiffAcceptor.Accepts"/> loops shift-major (for each shift, for
        /// each k) and returns the first match, so the smallest |shift| wins — in practice the zero notch.
        /// The alternative is k-major (for each k, for each shift), which resolves every collision toward the
        /// hypothesis whose apex was predicted correctly. This test simulates both orderings over the same
        /// acceptor inputs and reports where they disagree.
        /// </para>
        /// <para>
        /// <b>Depth is identical by construction:</b> both orderings enumerate the same window set, so the
        /// ACCEPTED/REJECTED decision is unchanged for every input. Only the notch attributed differs. The
        /// question is therefore purely which error we prefer, not how many identifications we keep.
        /// </para>
        /// </summary>
        [Test]
        public static void Prototype_KMajorOrdering_ChangesAttributionNotDepth()
        {
            const double deamidation = 0.984016;
            const double acetyl = 42.010565;
            const double carbamyl = 43.005814;
            var shifts = new[] { 0.0, deamidation, acetyl, carbamyl }.OrderBy(Math.Abs).ToArray();
            var tolerance = new PpmTolerance(10);
            double spacing = Constants.C13MinusC12;
            int[] ks = { 0, -1, 1, -2, 2 };
            const double mono = 9500.0;

            int NotchOf(double shift) => (int)Math.Round(shift * MassDiffAcceptor.NotchScalar);
            double ApexOf(double shift) => (mono + shift) + ApexOffset(mono + shift);

            // Current production ordering, mirrored from MostAbundantDotMassDiffAcceptor:78-93.
            int ShiftMajor(double observed)
            {
                foreach (double s in shifts)
                    foreach (int k in ks)
                        if (tolerance.Within(observed, ApexOf(s) + k * spacing)) return NotchOf(s);
                return -1;
            }
            // Proposed alternative: every shift is tried at k = 0 before any shift is tried at k = ±1.
            int KMajor(double observed)
            {
                foreach (int k in ks)
                    foreach (double s in shifts)
                        if (tolerance.Within(observed, ApexOf(s) + k * spacing)) return NotchOf(s);
                return -1;
            }

            var scenarios = new (string Name, double Observed, int Truth)[]
            {
                ("unmodified, on apex",            ApexOf(0.0),                    NotchOf(0.0)),
                ("unmodified, apex off by +1",     ApexOf(0.0) + spacing,          NotchOf(0.0)),
                ("deamidated, on apex",            ApexOf(deamidation),            NotchOf(deamidation)),
                ("acetylated, on apex",            ApexOf(acetyl),                 NotchOf(acetyl)),
                ("acetylated, apex off by +1",     ApexOf(acetyl) + spacing,       NotchOf(acetyl)),
                ("carbamylated, on apex",          ApexOf(carbamyl),               NotchOf(carbamyl)),
            };

            int shiftMajorCorrect = 0, kMajorCorrect = 0;
            foreach (var sc in scenarios)
            {
                int a = ShiftMajor(sc.Observed), b = KMajor(sc.Observed);

                // DEPTH: both orderings accept exactly the same inputs. This is the cost measurement.
                Assert.That(a >= 0, Is.EqualTo(b >= 0), $"depth changed for '{sc.Name}'");

                if (a == sc.Truth) shiftMajorCorrect++;
                if (b == sc.Truth) kMajorCorrect++;
                TestContext.WriteLine($"{sc.Name,-30} truth={sc.Truth,8}  shift-major={a,8}{(a == sc.Truth ? " ok " : " WRONG")}  k-major={b,8}{(b == sc.Truth ? " ok" : " WRONG")}");
            }
            TestContext.WriteLine($"correct attributions: shift-major {shiftMajorCorrect}/{scenarios.Length}, k-major {kMajorCorrect}/{scenarios.Length}");

            // Every accepted input stays accepted under either ordering: reordering costs ZERO depth.
            Assert.That(scenarios.All(sc => ShiftMajor(sc.Observed) >= 0 && KMajor(sc.Observed) >= 0), Is.True);

            // The orderings genuinely disagree — this is not a no-op.
            Assert.That(scenarios.Any(sc => ShiftMajor(sc.Observed) != KMajor(sc.Observed)), Is.True,
                "the two orderings should disagree on at least one collision scenario");
        }

        /// <summary>
        /// FDR consequences of the apex notch set — the question that matters most for top-down, where
        /// depth must be bought without losing FDR control.
        /// <para>
        /// Two separate properties, and they land differently:
        /// </para>
        /// <list type="number">
        /// <item><b>Validity is preserved.</b> Targets and decoys traverse the SAME acceptor with the SAME
        /// five windows per shift, so the null distribution widens symmetrically and the per-notch q-value
        /// stays calibrated. GPTMDTask:182 passes <c>NumNotches</c> to FdrAnalysisEngine, and that count is
        /// identical in both modes — the ±k windows do not silently change the FDR grouping.</item>
        /// <item><b>Stratification is LOST relative to the plain most-abundant search.</b> The search
        /// acceptor makes k itself the notch, so a confident on-apex match and a ±2 apex-misprediction match
        /// are FDR-controlled separately. The G-PTM-D acceptor collapses all k onto the shift's notch, so
        /// those two populations are pooled into one q-value. Valid, but less discriminating.</item>
        /// </list>
        /// </summary>
        [Test]
        public static void FdrStratification_SearchSeparatesApexNotches_GptmdPoolsThem()
        {
            const double peptideMono = 9500.0;
            const double oxidation = 15.994915;
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            var gptmdMods = new List<Modification>
            {
                new Modification(_originalId: "ox", _modificationType: "mt", _target: motifM,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: oxidation)
            };
            var shifts = GptmdTask.GetAcceptableMassShifts(new List<Modification>(), new List<Modification>(),
                gptmdMods, new List<Tuple<double, double>>()).ToList();
            var tolerance = new PpmTolerance(10);
            double spacing = Constants.C13MinusC12;

            var searchAcceptor = new MostAbundantMassDiffAcceptor("search", tolerance, Averagine);
            var gptmdAcceptor = new MostAbundantDotMassDiffAcceptor("gptmd", shifts, tolerance, Averagine);
            var traditional = new DotMassDiffAcceptor("gptmd", shifts, tolerance);

            double apex = peptideMono + ApexOffset(peptideMono);

            // SEARCH acceptor: k is the notch, so on-apex and ±2 land in DIFFERENT FDR groups.
            int searchK0 = searchAcceptor.Accepts(apex, peptideMono);
            int searchK2 = searchAcceptor.Accepts(apex + 2 * spacing, peptideMono);
            Assert.That(searchK0, Is.GreaterThanOrEqualTo(0));
            Assert.That(searchK2, Is.GreaterThanOrEqualTo(0));
            Assert.That(searchK2, Is.Not.EqualTo(searchK0),
                "the plain most-abundant search should FDR-control each apex offset separately");
            Assert.That(searchAcceptor.NumNotches, Is.EqualTo(5), "one notch per k");

            // G-PTM-D acceptor: the shift is the notch, so on-apex and ±2 are POOLED into one FDR group.
            int gptmdK0 = gptmdAcceptor.Accepts(apex, peptideMono);
            int gptmdK2 = gptmdAcceptor.Accepts(apex + 2 * spacing, peptideMono);
            Assert.That(gptmdK0, Is.GreaterThanOrEqualTo(0));
            Assert.That(gptmdK2, Is.EqualTo(gptmdK0),
                "G-PTM-D pools every apex offset for a shift into that shift's notch");

            // Validity plumbing: the notch COUNT handed to FdrAnalysisEngine (GPTMDTask:182) is unchanged by
            // the apex windows, so most-abundant mode does not silently redefine the FDR grouping.
            Assert.That(gptmdAcceptor.NumNotches, Is.EqualTo(traditional.NumNotches));
            Assert.That(gptmdAcceptor.NumNotches, Is.EqualTo(shifts.Count));

            // And the pooling is what keeps the notch groups populated: stratifying by k as well would
            // multiply the group count by 5, which GPTMDTask:181 already warns is the failure mode
            // ("it takes multiple files to get enough PSMs for all the different notches").
            Assert.That(searchAcceptor.NumNotches * shifts.Count, Is.EqualTo(gptmdAcceptor.NumNotches * 5),
                "stratifying G-PTM-D by shift AND k would give 5x the FDR groups");
        }

        /// <summary>
        /// Candidate inflation — the historical reason missed-monoisotopic notches were kept out of
        /// traditional G-PTM-D, measured against the REAL default G-PTM-D modification list.
        /// <para>
        /// This is a different failure from misassignment. <see cref="MassDiffAcceptor.Accepts"/> returns the
        /// FIRST matching shift, so at most one notch is wrong. But GetPossibleMods (GptmdEngine:253-282) has
        /// no early exit — it yields EVERY modification whose candidate mass matches. Widening each test from
        /// a ppm window to ±2 neutrons widens the accepted mass span to ~4 Da, and since modification masses
        /// cluster roughly a Dalton apart, the shortlist inflates. Every shortlisted mod is then placed at
        /// every legal site, re-fragmented and re-scored, and any survivor is written into the database.
        /// </para>
        /// <para>
        /// The unmodified case is the sharp one: a peptide carrying NO modification should propose nothing.
        /// </para>
        /// </summary>
        [Test]
        public static void ApexNotches_InflateTheGptmdCandidateShortlist()
        {
            // The actual default G-PTM-D list: Common Artifact + Common Biological + Metal.
            var defaultIds = new GptmdParameters().ListOfModsGptmd.Select(t => t.Item2).ToHashSet();
            var gptmdMods = GlobalVariables.AllModsKnown
                .Where(m => m.ValidModification && m.MonoisotopicMass.HasValue && defaultIds.Contains(m.IdWithMotif))
                .ToList();
            Assert.That(gptmdMods.Count, Is.GreaterThan(20), "expected the real default G-PTM-D mod list");

            const double peptideMono = 9500.0;
            var tolerance = new PpmTolerance(10);
            var monoParams = new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.Monoisotopic);
            var apexParams = new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);

            // One observation carrying BOTH masses; only the interpretation changes between modes.
            SpectralMatch PsmFor(double trueModMass)
            {
                double trueMono = peptideMono + trueModMass;
                double trueApex = trueMono + ApexOffset(trueMono);
                var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1,
                    true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null,
                    "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
                var scan = new Ms2ScanWithSpecificMass(dataScan, trueMono.ToMz(1), 1, "f", new CommonParameters(),
                    precursorMostAbundantMass: trueApex);
                ModificationMotif.TryGetMotif("X", out ModificationMotif motifX);
                var pep = new PeptideWithSetModifications("PEPTIDEK", new Dictionary<string, Modification>());
                return new PeptideSpectralMatch(pep, 0, 0, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            }

            // Count exactly what GetPossibleMods counts: one MatchesCandidateMass call per mod (line 261).
            int Shortlist(SpectralMatch psm, CommonParameters cp, out int distinctMasses)
            {
                var hits = gptmdMods
                    .Where(m => psm.MatchesCandidateMass(peptideMono + m.MonoisotopicMass.Value, tolerance, cp))
                    .ToList();
                distinctMasses = hits.Select(m => Math.Round(m.MonoisotopicMass.Value, 4)).Distinct().Count();
                return hits.Count;
            }

            // --- Case 1: a genuinely OXIDISED proteoform -------------------------------------------------
            const double oxidation = 15.994915;
            var oxPsm = PsmFor(oxidation);
            int monoOx = Shortlist(oxPsm, monoParams, out int monoOxMasses);
            int apexOx = Shortlist(oxPsm, apexParams, out int apexOxMasses);
            TestContext.WriteLine($"oxidised  : monoisotopic {monoOx} mods ({monoOxMasses} masses) -> most-abundant {apexOx} mods ({apexOxMasses} masses)");

            Assert.That(monoOx, Is.GreaterThan(0), "positive control: the true modification must be shortlisted");
            Assert.That(apexOx, Is.GreaterThan(0), "positive control: the true modification must be shortlisted");
            Assert.That(apexOx, Is.GreaterThanOrEqualTo(monoOx));
            Assert.That(apexOxMasses, Is.GreaterThan(monoOxMasses),
                "the apex windows should admit modification masses the monoisotopic test rejects");

            // --- Case 2: an UNMODIFIED peptide ------------------------------------------------------------
            var cleanPsm = PsmFor(0.0);
            int monoClean = Shortlist(cleanPsm, monoParams, out int monoCleanMasses);
            int apexClean = Shortlist(cleanPsm, apexParams, out int apexCleanMasses);
            TestContext.WriteLine($"unmodified: monoisotopic {monoClean} mods ({monoCleanMasses} masses) -> most-abundant {apexClean} mods ({apexCleanMasses} masses)");

            // NEGATIVE CONTROL — a monoisotopic search proposes nothing for an unmodified peptide, because no
            // real modification has a mass within 0.095 Da of zero.
            Assert.That(monoClean, Is.EqualTo(0),
                "monoisotopic G-PTM-D should propose no modification for an unmodified peptide");

            // THE HAZARD — most-abundant proposes modifications for a peptide that carries none, because the
            // ±k windows reach the ~1 and ~2 Da region where many modification masses live.
            Assert.That(apexClean, Is.GreaterThan(0),
                "expected the apex windows to admit near-neutron modifications for an unmodified peptide");
            TestContext.WriteLine("  spurious on unmodified peptide: " + string.Join(", ",
                gptmdMods.Where(m => cleanPsm.MatchesCandidateMass(peptideMono + m.MonoisotopicMass.Value, tolerance, apexParams))
                         .Select(m => $"{m.IdWithMotif} ({m.MonoisotopicMass.Value:F4})").Distinct().Take(12)));
        }

        /// <summary>
        /// The monoisotopic fallback in GetPrecursorMassForSearch is deliberate — a most-abundant search
        /// cannot match on a peak that was never observed — but it is silent, and a run where most scans take
        /// it is not the search the user asked for. GetMonoisotopicFallbackWarning reports it. Warn only when
        /// there is something to warn about: never in monoisotopic mode (where the fallback does not exist),
        /// and never when every precursor has an observed apex.
        /// </summary>
        [Test]
        public static void MonoisotopicFallbackWarning_ReportsOnlyWhenScansLackAnApex()
        {
            var monoParams = new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.Monoisotopic);
            var apexParams = new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);

            Ms2ScanWithSpecificMass MakeScan(double? apex)
            {
                var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1,
                    true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null,
                    "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
                return new Ms2ScanWithSpecificMass(dataScan, 1000.0.ToMz(1), 1, "f.mzML", new CommonParameters(),
                    precursorMostAbundantMass: apex);
            }

            var withApex = MakeScan(1005.0);
            var withoutApex = MakeScan(null);   // no envelope deconvoluted -> PrecursorMostAbundantMass == 0
            Assert.That(withoutApex.PrecursorMostAbundantMass, Is.EqualTo(0));

            var mixed = new[] { withApex, withoutApex, withoutApex, withApex };

            // Monoisotopic mode: the fallback does not exist, so there is nothing to warn about.
            Assert.That(mixed.GetMonoisotopicFallbackWarning(monoParams), Is.Null);

            // Most-abundant mode with a complete scan set: no fallback taken, no warning.
            Assert.That(new[] { withApex, withApex }.GetMonoisotopicFallbackWarning(apexParams), Is.Null);

            // Most-abundant mode with scans that lack an apex: warn, and say how many.
            string warning = mixed.GetMonoisotopicFallbackWarning(apexParams, "myfile.mzML");
            Assert.That(warning, Is.Not.Null, "most-abundant search silently fell back to monoisotopic matching");
            Assert.That(warning, Does.Contain("2 of 4"));
            Assert.That(warning, Does.Contain("50.0%"));
            Assert.That(warning, Does.Contain("myfile.mzML"));
            Assert.That(warning, Does.Contain("monoisotopic"));

            // The file name is optional and simply omitted when absent.
            Assert.That(mixed.GetMonoisotopicFallbackWarning(apexParams), Does.Contain("2 of 4"));
            Assert.That(mixed.GetMonoisotopicFallbackWarning(apexParams), Does.Not.Contain(" in "));

            // Degenerate inputs must not throw or invent a warning.
            Assert.That(((Ms2ScanWithSpecificMass[])null).GetMonoisotopicFallbackWarning(apexParams), Is.Null);
            Assert.That(Array.Empty<Ms2ScanWithSpecificMass>().GetMonoisotopicFallbackWarning(apexParams), Is.Null);
            Assert.That(new[] { withoutApex, null }.GetMonoisotopicFallbackWarning(apexParams), Does.Contain("1 of 1"));
        }

        /// <summary>
        /// The warning is actually emitted by a real task run, not merely computable. Runs a search whose theory
        /// (most-abundant mass diff acceptor) and observed (precursor mass match mode) sides are both most-abundant
        /// over a file whose precursors deconvolute, with precursor deconvolution turned OFF so every scan falls
        /// back, and asserts the task raised the warning.
        /// </summary>
        [Test]
        public static void SearchTask_EmitsMonoisotopicFallbackWarning_WhenNoApexObserved()
        {
            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => { lock (warnings) { warnings.Add(e.S); } };
            MetaMorpheusTask.WarnHandler += handler;
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestMostAbundantFallbackWarning");
            try
            {
                var searchTask = new SearchTask();
                searchTask.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.MostAbundant_Exact;
                searchTask.CommonParameters = new CommonParameters(
                    precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant,
                    doPrecursorDeconvolution: false);   // no envelope -> no observed apex on any scan

                string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
                string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
                Directory.CreateDirectory(outputFolder);

                searchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) },
                    new List<string> { myFile }, "test");

                Assert.That(warnings.Any(w => w != null && w.Contains("no observed most-abundant peak")), Is.True,
                    "a most-abundant search fell back to monoisotopic matching without warning");
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
                if (Directory.Exists(outputFolder)) { Directory.Delete(outputFolder, true); }
                string taskSettings = Path.Combine(TestContext.CurrentContext.TestDirectory, "Task Settings");
                if (Directory.Exists(taskSettings)) { Directory.Delete(taskSettings, true); }
            }
        }

        /// <summary>
        /// End-to-end proof that GptmdEngine discovers modifications through the APEX in most-abundant mode.
        /// <para>
        /// The G-PTM-D analogue of <see cref="ClassicSearch_SelectsScansByApex_NotMonoisotopic"/>. Mod
        /// discovery cannot reuse GetObservedMonoisotopicMass, because there the residual holds
        /// (k neutrons + the UNKNOWN modification mass) and rounding would destroy the modification. Instead
        /// GetPossibleMods (GptmdEngine:253-282) tests each candidate mass — peptide + modification — on its
        /// own via PrecursorMassExtensions.MatchesCandidateMass, apex-to-apex, never adjusting the
        /// observation.
        /// </para>
        /// <para>
        /// The scan here carries a deliberately WRONG monoisotopic mass and a CORRECT apex for the MODIFIED
        /// candidate, so a monoisotopic G-PTM-D run must discover nothing while a most-abundant run must
        /// discover the modification. Fragment scoring is deliberately uninformative (single-peak spectrum,
        /// no filters) so the test isolates precursor-driven mod selection, which is the only thing
        /// most-abundant mode changes.
        /// </para>
        /// </summary>
        [Test]
        public static void Gptmd_DiscoversModByApex_NotMonoisotopic()
        {
            // ~9.6 kDa peptide with four N residues, so the apex sits several neutrons above the
            // monoisotopic mass and the mod has legal sites.
            const string sequence = "PEPTIDEASYGLVNFQWMCH" + "PEPTIDEASYGLVNFQWMCH"
                                  + "PEPTIDEASYGLVNFQWMCH" + "PEPTIDEASYGLVNFQWMCH"
                                  + "PEPTIDE" + "K";
            const int nCount = 4;
            var protein = new Protein(sequence, "accession");
            var digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5, maxPeptideLength: 200);
            var peptide = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>(), null, null).First();

            // Deamidation-like mass on N: far from any whole number of neutrons, so an accepted match
            // cannot be an isotope artifact in disguise.
            const double modMass = 21.981943;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            var gptmdModifications = new List<Modification>
            {
                new Modification(_originalId: "testmod", _modificationType: "mt", _target: motifN,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: modMass)
            };

            double modifiedMono = peptide.MonoisotopicMass + modMass;
            double modifiedApex = modifiedMono + ApexOffset(modifiedMono);
            Assert.That(modifiedApex - modifiedMono, Is.GreaterThan(1.0),
                "test is only meaningful if the apex differs from the monoisotopic mass");

            var tolerance = new PpmTolerance(10);
            var monoParams = new CommonParameters(digestionParams: digestionParams,
                precursorMassMatchMode: PrecursorMassMatchMode.Monoisotopic);
            var apexParams = new CommonParameters(digestionParams: digestionParams,
                precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant);
            var fsp = new List<(string FileName, CommonParameters Parameters)> { ("filepath", monoParams) };
            var combos = new List<Tuple<double, double>>();

            // Uninformative single-peak spectrum: every legal site scores the same, so nothing but the
            // precursor decides which mods are shortlisted.
            MsDataScan MakeDataScan() => new MsDataScan(
                new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive,
                double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN,
                null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);

            SpectralMatch MakePsm(double observedMono, double observedApex, CommonParameters cp)
            {
                var dataScan = MakeDataScan();
                var scan = new Ms2ScanWithSpecificMass(dataScan, observedMono.ToMz(1), 1, "filepath", cp,
                    precursorMostAbundantMass: observedApex);
                SpectralMatch psm = new PeptideSpectralMatch(peptide, 0, 0, 0, scan, cp, new List<MatchedFragmentIon>());
                psm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
                psm.SetMs2Scan(scan.TheScan);
                return psm;
            }

            GptmdResults RunGptmd(SpectralMatch psm, CommonParameters cp) =>
                (GptmdResults)new GptmdEngine(new List<SpectralMatch> { psm }, gptmdModifications, combos,
                    new Dictionary<string, Tolerance> { { "filepath", tolerance } }, cp, fsp,
                    new List<string>(), null).Run();

            // Monoisotopic mass deliberately wrong by 3.5 Da; apex exactly right for peptide + mod.
            const double monoError = 3.5;

            // A monoisotopic G-PTM-D run cannot reach the modification: the mass it reads is 3.5 Da off.
            var monoResults = RunGptmd(MakePsm(modifiedMono + monoError, modifiedApex, monoParams), monoParams);
            Assert.That(monoResults.Mods.Count, Is.EqualTo(0),
                "monoisotopic G-PTM-D discovered a mod from a precursor mass that is 3.5 Da off");

            // A most-abundant run finds it, because GetPossibleMods tested (peptide + mod) apex-to-apex.
            var apexResults = RunGptmd(MakePsm(modifiedMono + monoError, modifiedApex, apexParams), apexParams);
            Assert.That(apexResults.Mods.Count, Is.EqualTo(1),
                "most-abundant G-PTM-D failed to discover the mod by its apex — discovery is still monoisotopic");
            Assert.That(apexResults.Mods["accession"].Count, Is.EqualTo(nCount));
            Assert.That(apexResults.Mods["accession"].Select(m => m.Item2.MonoisotopicMass.Value),
                Is.All.EqualTo(modMass).Within(1e-6));

            // CONTROL: the monoisotopic path is not broken — give it a correct monoisotopic mass and it
            // discovers the same mod at the same sites. So the empty result above is caused by the mass it
            // matched on, not by the mod list, the motif, or the FDR/score plumbing.
            var controlResults = RunGptmd(MakePsm(modifiedMono, modifiedApex, monoParams), monoParams);
            Assert.That(controlResults.Mods.Count, Is.EqualTo(1),
                "control: monoisotopic G-PTM-D should discover the mod from a correct monoisotopic mass");
            Assert.That(controlResults.Mods["accession"].Count, Is.EqualTo(nCount));

            // NEGATIVE CONTROL: most-abundant discovery is not simply permissive. Move the observed apex
            // 7.5 Da away — beyond the ±2-neutron apex tolerance (±2.007 Da) and not a whole number of
            // neutrons, so no k can rescue it — and no mod is proposed for any site.
            var negativeResults = RunGptmd(MakePsm(modifiedMono + monoError, modifiedApex + 7.5, apexParams), apexParams);
            Assert.That(negativeResults.Mods.Count, Is.EqualTo(0),
                "negative control: most-abundant G-PTM-D proposed a mod for an apex 7.5 Da outside every notch");

            // NEGATIVE CONTROL 2: a mod whose mass the observation does not support is not proposed, even
            // though the apex is otherwise correct. Guards against the ±k windows swallowing the mod mass.
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN2);
            var wrongMassMods = new List<Modification>
            {
                new Modification(_originalId: "wrongmass", _modificationType: "mt", _target: motifN2,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: modMass + 7.5)
            };
            var wrongModResults = (GptmdResults)new GptmdEngine(
                new List<SpectralMatch> { MakePsm(modifiedMono + monoError, modifiedApex, apexParams) },
                wrongMassMods, combos, new Dictionary<string, Tolerance> { { "filepath", tolerance } },
                apexParams, fsp, new List<string>(), null).Run();
            Assert.That(wrongModResults.Mods.Count, Is.EqualTo(0),
                "negative control: most-abundant G-PTM-D proposed a mod whose mass the apex does not support");
        }
    }
}
