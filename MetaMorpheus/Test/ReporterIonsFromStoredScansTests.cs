using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// TMT reporter ions have to be read from the scans as the file stores them. The search read them from its loaded
    /// copy, which MyFileManager has already trimmed: MS2 peaks when TrimMsMsPeaks is on (the default keeps the 200 most
    /// intense and drops anything under 1% of the base peak), and MS3 peaks always, because the filter's MSn switch is
    /// never turned off. A reporter ion the trimming removed was written as 0 and summed into peptide and protein values.
    /// </summary>
    [TestFixture]
    public static class ReporterIonsFromStoredScansTests
    {
        private static string TestFile(string name) => Path.Combine(TestContext.CurrentContext.TestDirectory, "TMT_test", name);

        private static string[] Formatted(double[] intensities) =>
            intensities.Select(i => i.ToString("F1", CultureInfo.InvariantCulture)).ToArray();

        /// <summary>
        /// The reporter values the search would extract from <paramref name="file"/> for the MS2 scan
        /// <paramref name="ms2ScanNumber"/>: from its most intense MS3 child when it has one, otherwise from itself.
        /// </summary>
        private static double[] Extract(MsDataFile file, int ms2ScanNumber, IsobaricMassTag tag)
        {
            var scans = file.GetAllScansList();
            var child = scans
                .Where(s => s.MsnOrder == 3 && s.OneBasedPrecursorScanNumber == ms2ScanNumber)
                .MaxBy(s => s.TotalIonCurrent);
            var source = child ?? scans.Single(s => s.OneBasedScanNumber == ms2ScanNumber);
            return tag.GetReporterIonIntensities(source.MassSpectrum);
        }

        [Test]
        public static void Ms2ReporterIons_AreReadFromTheStoredScan_NotTheTrimmedOne()
        {
            var searchTask = Toml.ReadFile<SearchTask>(TestFile("TMT-Task1-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            Assert.That(searchTask.CommonParameters.TrimMsMsPeaks, Is.True, "the fixture only means anything with trimming on");
            var tag = IsobaricMassTag.GetIsobaricMassTag(searchTask.SearchParameters.MultiplexModId);
            string mzml = TestFile("VA084TQ_6.mzML");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestReporterIonsFromStoredScansMs2");

            try
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("search", searchTask) },
                    new List<string> { mzml }, new List<DbForTask> { new DbForTask(TestFile("mouseTMT.fasta"), false) },
                    outputFolder).Run();

                string[] lines = File.ReadAllLines(Path.Combine(outputFolder, "search", "AllPSMs.psmtsv"));
                string[] header = lines[0].Split('\t');
                int scanColumn = System.Array.IndexOf(header, "Scan Number");
                int channels = tag.ReporterIonMzs.Length;
                var written = lines.Skip(1).Select(l => l.Split('\t'))
                    .Select(c => (scan: int.Parse(c[scanColumn], CultureInfo.InvariantCulture), reporters: c[^channels..]))
                    .ToList();

                var stored = MsDataFileReader.GetDataFile(mzml).LoadAllStaticData();
                var loaded = new MyFileManager(true).LoadFile(mzml, searchTask.CommonParameters);

                Assert.That(written, Is.Not.Empty);
                Assert.That(written.Any(w => !Formatted(Extract(stored, w.scan, tag)).SequenceEqual(Formatted(Extract(loaded, w.scan, tag)))),
                    Is.True, "the fixture only means anything if trimming at these settings removes a reporter from some PSM's scan");

                Assert.Multiple(() =>
                {
                    foreach (var (scan, reporters) in written)
                    {
                        Assert.That(reporters, Is.EqualTo(Formatted(Extract(stored, scan, tag))), $"scan {scan}");
                    }
                });
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// MS3 reporter scans are trimmed even with TrimMsMsPeaks off, because the load filter's MSn switch defaults to on
        /// and is never passed. A LowCID search loads without a filter, so this uses CID to reach the trimmed path. At the
        /// default 200 peaks no reporter in this file is lost (an MS3 spectrum is mostly reporters), so it keeps 10 peaks
        /// in each of 10 windows, a setting a search toml has used.
        /// </summary>
        [Test]
        public static void Ms3ReporterIons_AreReadFromTheStoredChildScan_NotTheTrimmedOne()
        {
            string mzml = TestFile("MS3_TMT11_Mouse_snip.mzML");
            var tag = IsobaricMassTag.GetIsobaricMassTag("TMT11");
            var commonParameters = new CommonParameters(dissociationType: DissociationType.CID,
                ms3childScanDissociationType: DissociationType.HCD, trimMsMsPeaks: false,
                numberOfPeaksToKeepPerWindow: 10, numberOfWindows: 10);

            var loaded = new MyFileManager(true).LoadFile(mzml, commonParameters);
            var stored = MsDataFileReader.GetDataFile(mzml).LoadAllStaticData();
            var storedByScanNumber = stored.GetAllScansList().ToDictionary(s => s.OneBasedScanNumber);

            var withChildren = MetaMorpheusTask.GetMs2Scans(loaded, mzml, commonParameters)
                .Where(s => s.ChildScans.Any(c => c.TheScan.MsnOrder == 3))
                .ToList();
            Assert.That(withChildren, Is.Not.Empty);
            Assert.That(withChildren.Any(s => !Formatted(Extract(stored, s.OneBasedScanNumber, tag))
                    .SequenceEqual(Formatted(Extract(loaded, s.OneBasedScanNumber, tag)))),
                Is.True, "the fixture only means anything if the loaded MS3 scans lost a reporter to trimming");

            Assert.Multiple(() =>
            {
                foreach (var scan in withChildren)
                {
                    scan.SetIsobaricMassTagReporterIonIntensities(tag, storedByScanNumber);
                    Assert.That(Formatted(scan.IsobaricMassTagReporterIonIntensities),
                        Is.EqualTo(Formatted(Extract(stored, scan.OneBasedScanNumber, tag))), $"scan {scan.OneBasedScanNumber}");
                }
            });
        }
    }
}
