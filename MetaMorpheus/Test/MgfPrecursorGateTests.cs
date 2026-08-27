using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Whether a search may deconvolute precursors out of an mgf. Third-party mgf files carry precursor
    /// information only in the scan header; one mzLib writes carries its MS1 scans too.
    /// </summary>
    [TestFixture]
    internal static class MgfPrecursorGateTests
    {
        private static MsDataScan Ms1(int scanNumber) =>
            new(new MzSpectrum(new[] { 570.0, 571.0, 572.0, 573.0 }, new[] { 100.0, 200.0, 300.0, 400.0 }, false),
                oneBasedScanNumber: scanNumber, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.0 + scanNumber, scanWindowRange: new MzRange(500, 600), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1000.0, injectionTime: null,
                noiseData: null, nativeId: "scan=" + scanNumber);

        private static MsDataScan Ms2(int scanNumber, int? precursorScanNumber, double? isolationWidth) =>
            new(new MzSpectrum(new[] { 110.05, 220.11 }, new[] { 500.0, 750.0 }, false),
                oneBasedScanNumber: scanNumber, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.0 + scanNumber, scanWindowRange: new MzRange(100, 250), scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1250.0, injectionTime: null,
                noiseData: null, nativeId: "scan=" + scanNumber, selectedIonMz: 571.8069,
                selectedIonChargeStateGuess: 2, selectedIonIntensity: 999999.0,
                isolationMZ: 572.0, isolationWidth: isolationWidth, dissociationType: DissociationType.HCD,
                oneBasedPrecursorScanNumber: precursorScanNumber, selectedIonMonoisotopicGuessMz: 571.8069);

        private static MsDataFile FileOf(params MsDataScan[] scans) => new GenericMsDataFile(scans, null);

        /// <summary>
        /// The premise the gate protects: precursor deconvolution needs an MS1 to read and an isolation
        /// range to read it over. Anything less and the search must fall back to the scan header, which is
        /// what every third-party mgf has always required.
        /// </summary>
        [Test]
        public static void GateRequiresAnMs1ScanAndAResolvableIsolationRange()
        {
            Assert.That(SearchTask.HasDeconvolutablePrecursorScans(FileOf(Ms1(1), Ms2(2, 1, 3.0))), Is.True,
                "a file with MS1 scans, precursor scan numbers and isolation widths is deconvolutable");

            Assert.That(SearchTask.HasDeconvolutablePrecursorScans(FileOf(Ms2(1, null, 3.0))), Is.False,
                "no MS1 scans at all -- the ordinary third-party mgf");

            Assert.That(SearchTask.HasDeconvolutablePrecursorScans(FileOf(Ms1(1), Ms2(2, null, 3.0))), Is.False,
                "MS1 present but nothing points at it");

            Assert.That(SearchTask.HasDeconvolutablePrecursorScans(FileOf(Ms1(1), Ms2(2, 1, null))), Is.False,
                "no isolation width, so IsolationRange is null and deconvolution would return nothing");

            Assert.That(SearchTask.HasDeconvolutablePrecursorScans(FileOf(Ms1(1), Ms2(2, 99, 3.0))), Is.False,
                "precursor scan number names a scan that is not in the file");

            Assert.That(SearchTask.HasDeconvolutablePrecursorScans(FileOf(Ms1(1))), Is.False,
                "no MS2 scans to search");
        }

        /// <summary>
        /// One partially annotated scan takes the whole file down the header-only path, rather than leaving
        /// some scans deconvoluted and others not for no visible reason.
        /// </summary>
        [Test]
        public static void GateIsAllOrNothingAcrossTheFile()
        {
            Assert.That(SearchTask.HasDeconvolutablePrecursorScans(
                    FileOf(Ms1(1), Ms2(2, 1, 3.0), Ms2(3, 1, 3.0))), Is.True);

            Assert.That(SearchTask.HasDeconvolutablePrecursorScans(
                    FileOf(Ms1(1), Ms2(2, 1, 3.0), Ms2(3, 1, null))), Is.False,
                "one unannotated scan out of two");
        }

        /// <summary>
        /// A round trip through the writer and reader has to satisfy the gate, or the mgf output option
        /// silently costs the next task its precursors.
        /// </summary>
        [Test]
        public static void AnMgfWrittenByMzLibSatisfiesTheGateAndAnMs2OnlyOneDoesNot()
        {
            var file = FileOf(Ms1(1), Ms2(2, 1, 3.0));
            string withMs1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "gateWithMs1.mgf");
            string ms2Only = Path.Combine(TestContext.CurrentContext.TestDirectory, "gateMs2Only.mgf");

            try
            {
                file.ExportAsMgf(withMs1);
                file.ExportAsMgf(ms2Only, includeMs1Scans: false);

                var readWithMs1 = MsDataFileReader.GetDataFile(withMs1);
                readWithMs1.LoadAllStaticData();
                Assert.That(SearchTask.HasDeconvolutablePrecursorScans(readWithMs1), Is.True);

                var readMs2Only = MsDataFileReader.GetDataFile(ms2Only);
                readMs2Only.LoadAllStaticData();
                Assert.That(SearchTask.HasDeconvolutablePrecursorScans(readMs2Only), Is.False);
            }
            finally
            {
                File.Delete(withMs1);
                File.Delete(ms2Only);
            }
        }

        /// <summary>
        /// The regression that motivated the content-based gate. SearchTask used to switch precursor
        /// deconvolution off for anything typed as <see cref="Mgf"/>, on the premise that an mgf has no MS1
        /// scans to deconvolute. An mgf mzLib wrote does, and the search silently fell back to one
        /// scan-header precursor per MS2 scan.
        /// </summary>
        [Test]
        [NonParallelizable]
        public static void SearchDeconvolutesPrecursorsOnlyFromAnMgfThatCarriesMs1Scans()
        {
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "MgfSearchGate");
            Directory.CreateDirectory(unitTestFolder);

            try
            {
                var source = MsDataFileReader.GetDataFile(
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML"));
                source.LoadAllStaticData();

                string withMs1 = Path.Combine(unitTestFolder, "withMs1.mgf");
                string ms2Only = Path.Combine(unitTestFolder, "ms2Only.mgf");
                source.ExportAsMgf(withMs1);
                source.ExportAsMgf(ms2Only, includeMs1Scans: false);

                string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta");

                (int precursors, int ms2Scans) RunSearch(string spectraPath, string name)
                {
                    string outputFolder = Path.Combine(unitTestFolder, name);
                    Directory.CreateDirectory(outputFolder);
                    new SearchTask().RunTask(outputFolder, new List<DbForTask> { new(database, false) },
                        new List<string> { spectraPath }, "test");

                    string[] results = File.ReadAllLines(Path.Combine(outputFolder, "results.txt"));
                    int Read(string label) => int.Parse(
                        results.First(l => l.StartsWith("All " + label + ": ")).Split(": ")[1]);
                    return (Read("Precursors"), Read("MS2 Scans"));
                }

                var deconvoluted = RunSearch(withMs1, "withMs1");
                var headerOnly = RunSearch(ms2Only, "ms2Only");

                Assert.That(deconvoluted.ms2Scans, Is.EqualTo(headerOnly.ms2Scans),
                    "both files hold the same MS2 scans");
                Assert.That(deconvoluted.precursors, Is.GreaterThan(deconvoluted.ms2Scans),
                    "MS1 scans are present, so deconvolution should find more precursors than there are scans");
                Assert.That(headerOnly.precursors, Is.EqualTo(headerOnly.ms2Scans),
                    "no MS1 scans, so the search must fall back to one scan-header precursor per scan");
            }
            finally
            {
                Directory.Delete(unitTestFolder, true);
            }
        }
    }
}
