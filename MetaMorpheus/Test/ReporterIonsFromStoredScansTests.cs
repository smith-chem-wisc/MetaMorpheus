using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Omics.Modifications;
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
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(searchTask.SearchParameters.MultiplexModId, GlobalVariables.AllModsKnown, out var tag), Is.True);
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
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag("TMT11", GlobalVariables.AllModsKnown, out var tag), Is.True);
            var commonParameters = new CommonParameters(dissociationType: DissociationType.CID,
                ms3childScanDissociationType: DissociationType.HCD, trimMsMsPeaks: false,
                numberOfPeaksToKeepPerWindow: 10, numberOfWindows: 10);

            var loaded = new MyFileManager(true).LoadFile(mzml, commonParameters);
            var stored = MsDataFileReader.GetDataFile(mzml).LoadAllStaticData();
            var storedScanFor = MyFileManager.PairWithStoredScans(loaded, stored);

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
                    scan.SetIsobaricMassTagReporterIonIntensities(tag, storedScanFor);
                    Assert.That(Formatted(scan.IsobaricMassTagReporterIonIntensities),
                        Is.EqualTo(Formatted(Extract(stored, scan.OneBasedScanNumber, tag))), $"scan {scan.OneBasedScanNumber}");
                }
            });
        }

        /// <summary>The scan of VA084TQ_6.mzML that is identified, and that the MGF fixture below writes twice.</summary>
        private const int RepeatedScan = 6;

        /// <summary>
        /// Writes <paramref name="mzml"/> as an MGF in which the block for <see cref="RepeatedScan"/> is followed by a copy
        /// with the same SCANS value, as an export that writes a scan once per candidate charge would. A static MGF load
        /// keeps both blocks. The copy's reporter peaks (below m/z 135) are three times as intense and its TIC is 1, so the
        /// two blocks can be told apart in the spectra and in AllPSMs.psmtsv.
        /// </summary>
        private static string WriteMgfWithRepeatedScan(string mzml, string mgf)
        {
            MgfMethods.WriteMgf(MsDataFileReader.GetDataFile(mzml).LoadAllStaticData(), mgf);
            string text = File.ReadAllText(mgf).Replace("\r\n", "\n");
            int start = text.LastIndexOf("BEGIN IONS", text.IndexOf($"\nSCANS={RepeatedScan}\n"), System.StringComparison.Ordinal);
            int end = text.IndexOf("END IONS", start, System.StringComparison.Ordinal) + "END IONS".Length;
            var copy = text[start..end].Split('\n').Select(line =>
            {
                if (line.StartsWith("TIC=")) return "TIC=1";
                string[] peak = line.Split(' ');
                return peak.Length == 2 && double.TryParse(peak[0], NumberStyles.Float, CultureInfo.InvariantCulture, out double mz) && mz < 135
                    ? peak[0] + " " + (3 * double.Parse(peak[1], CultureInfo.InvariantCulture)).ToString("R", CultureInfo.InvariantCulture)
                    : line;
            });
            File.WriteAllText(mgf, text.Insert(end, "\n" + string.Join("\n", copy)));
            return mgf;
        }

        /// <summary>
        /// Keying the stored scans by scan number threw "An item with the same key has already been added" on an MGF that
        /// repeats a SCANS value, so the whole search failed where master searched the file.
        /// </summary>
        [Test]
        public static void Mgf_WithARepeatedScanNumber_IsSearched_AndEachPsmReadsItsOwnBlock()
        {
            var searchTask = Toml.ReadFile<SearchTask>(TestFile("TMT-Task1-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(searchTask.SearchParameters.MultiplexModId, GlobalVariables.AllModsKnown, out var tag), Is.True);
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestReporterIonsRepeatedMgfScan");

            try
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
                Directory.CreateDirectory(folder);
                string mgf = WriteMgfWithRepeatedScan(TestFile("VA084TQ_6.mzML"), Path.Combine(folder, "VA084TQ_6.mgf"));
                string outputFolder = Path.Combine(folder, "output");
                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("search", searchTask) },
                    new List<string> { mgf }, new List<DbForTask> { new DbForTask(TestFile("mouseTMT.fasta"), false) },
                    outputFolder).Run();

                string[] lines = File.ReadAllLines(Path.Combine(outputFolder, "search", "AllPSMs.psmtsv"));
                string[] header = lines[0].Split('\t');
                int scanColumn = System.Array.IndexOf(header, "Scan Number");
                int ticColumn = System.Array.IndexOf(header, "Total Ion Current");
                int channels = tag.ReporterIonMzs.Length;
                var rows = lines.Skip(1).Select(l => l.Split('\t'))
                    .Where(c => int.Parse(c[scanColumn], CultureInfo.InvariantCulture) == RepeatedScan)
                    .ToList();
                var blocks = MsDataFileReader.GetDataFile(mgf).LoadAllStaticData().GetAllScansList()
                    .Where(s => s.OneBasedScanNumber == RepeatedScan)
                    .ToList();

                Assert.That(blocks, Has.Count.EqualTo(2));
                Assert.That(rows, Is.Not.Empty, "the fixture only means anything if a block with the repeated scan number is identified");
                Assert.Multiple(() =>
                {
                    foreach (var row in rows)
                    {
                        double tic = double.Parse(row[ticColumn], CultureInfo.InvariantCulture);
                        var block = blocks.MinBy(b => System.Math.Abs(b.TotalIonCurrent - tic));
                        Assert.That(row[^channels..], Is.EqualTo(Formatted(tag.GetReporterIonIntensities(block.MassSpectrum))), $"TIC {tic}");
                    }
                });
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// Two blocks with one scan number are two scans: each is read from its own stored block, never from the other.
        /// </summary>
        [Test]
        public static void RepeatedScanNumber_EachBlockIsReadFromItsOwnStoredBlock()
        {
            var searchTask = Toml.ReadFile<SearchTask>(TestFile("TMT-Task1-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(searchTask.SearchParameters.MultiplexModId, GlobalVariables.AllModsKnown, out var tag), Is.True);
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestReporterIonsRepeatedScanPairing");

            try
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
                Directory.CreateDirectory(folder);
                string mgf = WriteMgfWithRepeatedScan(TestFile("VA084TQ_6.mzML"), Path.Combine(folder, "VA084TQ_6.mgf"));
                var loaded = new MyFileManager(true).LoadFile(mgf, searchTask.CommonParameters);
                var stored = MsDataFileReader.GetDataFile(mgf).LoadAllStaticData();
                var storedScanFor = MyFileManager.PairWithStoredScans(loaded, stored);
                List<MsDataScan> loadedScans = loaded.GetAllScansList();
                List<MsDataScan> storedScans = stored.GetAllScansList();

                Assert.That(storedScanFor, Has.Count.EqualTo(loadedScans.Count));
                var repeated = Enumerable.Range(0, loadedScans.Count).Where(i => loadedScans[i].OneBasedScanNumber == RepeatedScan).ToList();
                Assert.That(repeated, Has.Count.EqualTo(2));
                var expected = repeated.Select(i => Formatted(tag.GetReporterIonIntensities(storedScans[i].MassSpectrum))).ToList();
                Assert.That(expected[0], Is.Not.EqualTo(expected[1]), "the fixture only means anything if the two blocks' reporters differ");

                var ms2Scans = MetaMorpheusTask.GetMs2Scans(loaded, mgf, searchTask.CommonParameters).ToList();
                Assert.Multiple(() =>
                {
                    for (int k = 0; k < repeated.Count; k++)
                    {
                        var fromBlock = ms2Scans.Where(s => ReferenceEquals(s.TheScan, loadedScans[repeated[k]])).ToList();
                        Assert.That(fromBlock, Is.Not.Empty, $"block {k}");
                        foreach (var scan in fromBlock)
                        {
                            scan.SetIsobaricMassTagReporterIonIntensities(tag, storedScanFor);
                            Assert.That(Formatted(scan.IsobaricMassTagReporterIonIntensities), Is.EqualTo(expected[k]), $"block {k}");
                        }
                    }
                });
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A LowCID search loads without a filter, so nothing is trimmed, but GetMs2Scans bins each MS2 spectrum in place for
        /// XCorr, and the search read its reporter ions from the binned spectrum: every reporter of every PSM was written as 0.
        /// They are read from the stored scans, here from the mzML and from an MGF written from it.
        /// </summary>
        [TestCase("mzML")]
        [TestCase("mgf")]
        public static void LowCidReporterIons_AreReadFromTheStoredScan_NotTheBinnedOne(string format)
        {
            var searchTask = Toml.ReadFile<SearchTask>(TestFile("TMT11_LowCID_SearchTask.toml"), MetaMorpheusTask.tomlConfig);
            Assert.That(searchTask.CommonParameters.DissociationType, Is.EqualTo(DissociationType.LowCID));
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(searchTask.SearchParameters.MultiplexModId, GlobalVariables.AllModsKnown, out var tag), Is.True);
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestReporterIonsLowCid_" + format);

            try
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
                Directory.CreateDirectory(folder);
                string spectra = TestFile("VA084TQ_6.mzML");
                if (format == "mgf")
                {
                    spectra = Path.Combine(folder, "VA084TQ_6.mgf");
                    MgfMethods.WriteMgf(MsDataFileReader.GetDataFile(TestFile("VA084TQ_6.mzML")).LoadAllStaticData(), spectra);
                }
                string outputFolder = Path.Combine(folder, "output");
                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("search", searchTask) },
                    new List<string> { spectra }, new List<DbForTask> { new DbForTask(TestFile("mouseTMT.fasta"), false) },
                    outputFolder).Run();

                string[] lines = File.ReadAllLines(Path.Combine(outputFolder, "search", "AllPSMs.psmtsv"));
                string[] header = lines[0].Split('\t');
                int scanColumn = System.Array.IndexOf(header, "Scan Number");
                int channels = tag.ReporterIonMzs.Length;
                var written = lines.Skip(1).Select(l => l.Split('\t'))
                    .Select(c => (scan: int.Parse(c[scanColumn], CultureInfo.InvariantCulture), reporters: c[^channels..]))
                    .ToList();

                var stored = MsDataFileReader.GetDataFile(spectra).LoadAllStaticData();
                var binned = new MyFileManager(true).LoadFile(spectra, searchTask.CommonParameters);
                _ = MetaMorpheusTask.GetMs2Scans(binned, spectra, searchTask.CommonParameters).ToList();

                Assert.That(written, Is.Not.Empty);
                Assert.That(written.Any(w => !Formatted(Extract(stored, w.scan, tag)).SequenceEqual(Formatted(Extract(binned, w.scan, tag)))),
                    Is.True, "the fixture only means anything if binning changes the reporters of some PSM's scan");

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
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// Pairing is by position, so it is checked: a position whose two scans have different scan numbers is left out, and
        /// a scan with no pair is read from its loaded spectrum, as master read every scan. Here the second read is the mzML
        /// the MGF was written from, which holds the repeated scan once, so every position after the repeat disagrees.
        /// </summary>
        [Test]
        public static void ScansWhoseTwoReadsDisagree_AreLeftUnpaired_AndReadFromTheLoadedSpectrum()
        {
            var searchTask = Toml.ReadFile<SearchTask>(TestFile("TMT-Task1-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(searchTask.SearchParameters.MultiplexModId, GlobalVariables.AllModsKnown, out var tag), Is.True);
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestReporterIonsUnpairedScans");

            try
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
                Directory.CreateDirectory(folder);
                string mzml = TestFile("VA084TQ_6.mzML");
                string mgf = WriteMgfWithRepeatedScan(mzml, Path.Combine(folder, "VA084TQ_6.mgf"));
                var loaded = new MyFileManager(true).LoadFile(mgf, searchTask.CommonParameters);
                List<MsDataScan> loadedScans = loaded.GetAllScansList();
                var storedScanFor = MyFileManager.PairWithStoredScans(loaded, MsDataFileReader.GetDataFile(mzml).LoadAllStaticData());

                int firstRepeat = loadedScans.FindIndex(s => s.OneBasedScanNumber == RepeatedScan);
                Assert.That(storedScanFor.Keys, Is.EquivalentTo(loadedScans.Take(firstRepeat + 1)));

                var unpaired = MetaMorpheusTask.GetMs2Scans(loaded, mgf, searchTask.CommonParameters)
                    .Where(s => s.OneBasedScanNumber > RepeatedScan)
                    .ToList();
                Assert.That(unpaired, Is.Not.Empty);
                Assert.Multiple(() =>
                {
                    foreach (var scan in unpaired)
                    {
                        scan.SetIsobaricMassTagReporterIonIntensities(tag);
                        string[] fromLoaded = Formatted(scan.IsobaricMassTagReporterIonIntensities);
                        scan.SetIsobaricMassTagReporterIonIntensities(tag, storedScanFor);
                        Assert.That(Formatted(scan.IsobaricMassTagReporterIonIntensities), Is.EqualTo(fromLoaded), $"scan {scan.OneBasedScanNumber}");
                    }
                });
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }
    }
}
