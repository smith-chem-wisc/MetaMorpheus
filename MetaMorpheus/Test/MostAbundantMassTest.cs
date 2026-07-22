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
        public static void SearchTask_GetMassDiffAcceptor_ReturnsMostAbundantWhenModeSet()
        {
            var acc = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.OneMM, null,
                PrecursorMassMatchMode.MostAbundant, Averagine);
            Assert.That(acc, Is.TypeOf<MostAbundantMassDiffAcceptor>());

            // Default mode is unaffected.
            var def = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.OneMM, null);
            Assert.That(def, Is.Not.TypeOf<MostAbundantMassDiffAcceptor>());
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
