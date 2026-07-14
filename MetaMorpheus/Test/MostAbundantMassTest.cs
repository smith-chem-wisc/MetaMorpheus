using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
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
        /// The acceptor — not the scan — decides which of the two a given search selects on.
        /// </summary>
        [Test]
        public static void Ms2Scan_PrecursorMasses_AreObservations_AcceptorSelectsBetweenThem()
        {
            var scan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap,
                double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null,
                DissociationType.AnyActivationType, 1, null);

            var monoAcceptor = new SinglePpmAroundZeroSearchMode(5);
            var apexAcceptor = new MostAbundantMassDiffAcceptor("mostAbundant", new PpmTolerance(5), Averagine);

            // No apex observed (e.g. a scan-header precursor): the observation is 0, and BOTH acceptors
            // select the monoisotopic mass — a most-abundant search cannot invent a peak that wasn't seen.
            var noEnvelopeScan = new Ms2ScanWithSpecificMass(scan, 1500.0.ToMz(2), 2, "", new CommonParameters());
            Assert.That(noEnvelopeScan.PrecursorMostAbundantMass, Is.EqualTo(0));
            Assert.That(noEnvelopeScan.GetPrecursorMassForSearch(monoAcceptor), Is.EqualTo(noEnvelopeScan.PrecursorMass).Within(1e-9));
            Assert.That(noEnvelopeScan.GetPrecursorMassForSearch(apexAcceptor), Is.EqualTo(noEnvelopeScan.PrecursorMass).Within(1e-9));

            // Apex observed: it is recorded without touching PrecursorMass, which keeps its monoisotopic
            // meaning. Only the most-abundant acceptor selects the apex.
            double apexMass = noEnvelopeScan.PrecursorMass + 5.0;
            var apexScan = new Ms2ScanWithSpecificMass(scan, 1500.0.ToMz(2), 2, "", new CommonParameters(),
                precursorMostAbundantMass: apexMass);
            Assert.That(apexScan.PrecursorMass, Is.EqualTo(noEnvelopeScan.PrecursorMass).Within(1e-9));
            Assert.That(apexScan.PrecursorMostAbundantMass, Is.EqualTo(apexMass).Within(1e-9));
            Assert.That(apexScan.GetPrecursorMassForSearch(monoAcceptor), Is.EqualTo(apexScan.PrecursorMass).Within(1e-9));
            Assert.That(apexScan.GetPrecursorMassForSearch(apexAcceptor), Is.EqualTo(apexMass).Within(1e-9));
        }

        /// <summary>
        /// The regression behind nbollis's review: a most-abundant search accepts candidates whose
        /// deconvoluted monoisotopic peak is off by whole isotopologues, so (ScanPrecursorMass -
        /// peptideMonoisotopic) is NOT a chemical mass difference. Any consumer reading it as one — the
        /// localization engine, G-PTM-D mod assignment, calibration — must go through the acceptor first,
        /// or it interprets a ~1-2 Da isotope-assignment offset as a modification.
        /// </summary>
        [Test]
        public static void ToObservedMonoisotopicMass_RemovesApexAndNeutronOffsets([Values(-2, -1, 0, 1, 2)] int neutronOffset)
        {
            const double peptideMono = 15000.0;
            var acceptor = new MostAbundantMassDiffAcceptor("mostAbundant", new PpmTolerance(5), Averagine);

            // What the detector sees: the apex, mispredicted by k neutrons, plus a little instrument drift.
            const double driftPpm = 2.0;
            double observedApex = (peptideMono + ApexOffset(peptideMono)) * (1 + driftPpm / 1e6)
                                  + neutronOffset * Constants.C13MinusC12;
            Assert.That(acceptor.Accepts(observedApex, peptideMono), Is.GreaterThanOrEqualTo(0), "should be an accepted match");

            // Undoing the apex offset and the k neutrons leaves the monoisotopic mass plus only the drift,
            // so the mass difference a downstream engine sees is ~0 — not k * 1.00335 Da.
            double observedMono = acceptor.ToObservedMonoisotopicMass(observedApex, peptideMono);
            double massDifferencePpm = (observedMono - peptideMono) / peptideMono * 1e6;
            Assert.That(massDifferencePpm, Is.EqualTo(driftPpm).Within(0.5));

            // A monoisotopic acceptor is the identity — baseline behaviour is untouched.
            var monoAcceptor = new SinglePpmAroundZeroSearchMode(5);
            Assert.That(monoAcceptor.ToObservedMonoisotopicMass(observedApex, peptideMono), Is.EqualTo(observedApex));
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

            // The optional column appears in the header only when requested.
            Assert.That(SpectralMatch.GetTabSeparatedHeader(includeMostAbundantColumn: true),
                Does.Contain("Most Abundant Mass Diff (ppm)"));
            Assert.That(SpectralMatch.GetTabSeparatedHeader(), Does.Not.Contain("Most Abundant Mass Diff (ppm)"));

            // End-to-end writer path: the emitted row cell aligns under the new header column and
            // carries the reported value (row and header stay column-aligned).
            var header = SpectralMatch.GetTabSeparatedHeader(includeMostAbundantColumn: true).Split('\t');
            var row = maPsm.ToString(new Dictionary<string, int>(), false, false, false, includeMostAbundantColumn: true).Split('\t');
            int col = Array.IndexOf(header, "Most Abundant Mass Diff (ppm)");
            Assert.That(col, Is.GreaterThanOrEqualTo(0));
            Assert.That(double.Parse(row[col], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(0).Within(0.5));
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

            // Only the acceptor changes what gets matched: the apex acceptor selects a mass that departs
            // from the monoisotopic one, the monoisotopic acceptor never does.
            var monoAcceptor = new SinglePpmAroundZeroSearchMode(5);
            var apexAcceptor = new MostAbundantMassDiffAcceptor("mostAbundant", new PpmTolerance(5), Averagine);
            Assert.That(apexScans.All(s => Math.Abs(s.GetPrecursorMassForSearch(monoAcceptor) - s.PrecursorMass) < 1e-6), Is.True);
            Assert.That(apexScans.Any(s => Math.Abs(s.GetPrecursorMassForSearch(apexAcceptor) - s.PrecursorMass) > 1e-6), Is.True);
        }
    }
}
