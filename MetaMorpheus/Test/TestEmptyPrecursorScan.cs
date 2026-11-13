using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// A damaged raw file can keep a scan's header but lose its peaks, leaving an MS1 scan that exists yet holds no
    /// peaks. Reading MS2 scans whose precursor is such a scan used to leave a null slot in _GetMs2Scans, which
    /// crashed every task that reads MS2 scans with "Value cannot be null. (Parameter 'source')".
    /// </summary>
    [TestFixture]
    public static class TestEmptyPrecursorScan
    {
        private const double PrecursorMz = 500.25;
        private const int PrecursorCharge = 2;

        /// <summary>
        /// Scan 1: MS1 with peaks. Scan 2: MS2 of scan 1.
        /// Scan 3: MS1 with no peaks. Scan 4: MS2 of scan 3.
        /// Every MS2 carries a header precursor m/z and charge, as Thermo files do.
        /// </summary>
        private static MsDataFile FileWithAnEmptyMs1()
        {
            var ms1WithPeaks = new MzSpectrum(new[] { 499.0, PrecursorMz, 500.75, 501.25 }, new[] { 10.0, 1000.0, 600.0, 200.0 }, false);
            var emptyMs1 = new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false);
            var ms2Peaks = new MzSpectrum(new[] { 175.119, 276.155, 389.239, 502.323 }, new[] { 100.0, 200.0, 300.0, 400.0 }, false);

            var scans = new[]
            {
                Ms1(ms1WithPeaks, 1, 1.0),
                Ms2(new MzSpectrum(ms2Peaks.XArray, ms2Peaks.YArray, false), 2, 1.01, precursorScanNumber: 1),
                Ms1(emptyMs1, 3, 1.02),
                Ms2(new MzSpectrum(ms2Peaks.XArray, ms2Peaks.YArray, false), 4, 1.03, precursorScanNumber: 3),
            };
            return new GenericMsDataFile(scans, new SourceFile("no nativeID format", "mzML format", null, null, null));
        }

        private static MsDataScan Ms1(MzSpectrum spectrum, int scanNumber, double rt) =>
            new MsDataScan(spectrum, scanNumber, 1, true, Polarity.Positive, rt, new MzRange(300, 2000), "ms1",
                MZAnalyzerType.Orbitrap, spectrum.SumOfAllY, 1.0, null, "scan=" + scanNumber);

        private static MsDataScan Ms2(MzSpectrum spectrum, int scanNumber, double rt, int precursorScanNumber) =>
            new MsDataScan(spectrum, scanNumber, 2, true, Polarity.Positive, rt, new MzRange(100, 2000), "ms2",
                MZAnalyzerType.Orbitrap, spectrum.SumOfAllY, 1.0, null, "scan=" + scanNumber,
                selectedIonMz: PrecursorMz, selectedIonChargeStateGuess: PrecursorCharge, selectedIonIntensity: 1000,
                isolationMZ: PrecursorMz, isolationWidth: 2, dissociationType: DissociationType.HCD,
                oneBasedPrecursorScanNumber: precursorScanNumber, selectedIonMonoisotopicGuessMz: PrecursorMz);

        /// <summary>Runs the action and returns the warnings it raised that name this test's file.</summary>
        private static List<string> WarningsFor(string fileName, Action action)
        {
            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (_, e) =>
            {
                if (e.S != null && e.S.StartsWith(fileName)) lock (warnings) warnings.Add(e.S);
            };
            MetaMorpheusTask.WarnHandler += handler;
            try { action(); }
            finally { MetaMorpheusTask.WarnHandler -= handler; }
            return warnings;
        }

        [Test]
        public static void EmptyPrecursorMs1FallsBackOnTheScanHeader()
        {
            const string fileName = "EmptyMs1HeaderFallback.mzML";
            List<Ms2ScanWithSpecificMass> scans = null;

            var warnings = WarningsFor(fileName, () =>
                scans = MetaMorpheusTask.GetMs2Scans(FileWithAnEmptyMs1(), fileName, new CommonParameters(useProvidedPrecursorInfo: true)).ToList());

            // the scan under the intact MS1 is unaffected
            Assert.That(scans.Any(s => s.OneBasedScanNumber == 2));

            // the scan under the empty MS1 is still searched, on the precursor recorded in its header
            var fromHeader = scans.Where(s => s.OneBasedScanNumber == 4).ToList();
            Assert.That(fromHeader, Has.Count.EqualTo(1));
            Assert.That(fromHeader[0].PrecursorCharge, Is.EqualTo(PrecursorCharge));
            Assert.That(fromHeader[0].PrecursorMonoisotopicPeakMz, Is.EqualTo(PrecursorMz).Within(1e-6));

            // one summary warning for the file, not one per scan
            Assert.That(warnings, Has.Count.EqualTo(1));
            Assert.That(warnings[0], Does.Contain("1 MS2 scan(s) have a precursor MS1 scan with no peaks"));
        }

        [Test]
        public static void EmptyPrecursorMs1IsSkippedWithoutTheScanHeader()
        {
            const string fileName = "EmptyMs1NoHeader.mzML";
            List<Ms2ScanWithSpecificMass> scans = null;

            var warnings = WarningsFor(fileName, () =>
                scans = MetaMorpheusTask.GetMs2Scans(FileWithAnEmptyMs1(), fileName, new CommonParameters(useProvidedPrecursorInfo: false)).ToList());

            Assert.That(scans.Any(s => s.OneBasedScanNumber == 4), Is.False);
            Assert.That(warnings, Has.Count.EqualTo(1));
            Assert.That(warnings[0], Does.Contain("skipped"));
        }

        [Test]
        public static void EveryScanSlotIsFilled()
        {
            var slots = MetaMorpheusTask._GetMs2Scans(FileWithAnEmptyMs1(), "EmptyMs1Slots.mzML", new CommonParameters());

            Assert.That(slots, Has.Length.EqualTo(2));
            Assert.That(slots, Has.None.Null);
        }

        [Test]
        public static void CrosslinkScanReadingSurvivesAnEmptyPrecursorMs1()
        {
            List<Ms2ScanWithSpecificMass> scans = null;

            Assert.DoesNotThrow(() => scans = MetaMorpheusTask.GetMs2ScansWrapByScanNum(FileWithAnEmptyMs1(), "EmptyMs1Xl.mzML",
                new CommonParameters(), out _));
            Assert.That(scans.Select(s => s.OneBasedScanNumber), Is.EquivalentTo(new[] { 2, 4 }));
        }
    }
}
