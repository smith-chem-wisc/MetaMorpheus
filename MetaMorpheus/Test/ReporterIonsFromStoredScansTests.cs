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

        /// <summary>
        /// Writes <paramref name="mzml"/> unindexed, with its native IDs renumbered to number its spectra from
        /// <paramref name="firstScan"/>. VA084TQ_6.mzML is a slice of a run (scan=31006 to scan=31014), so a read by scan
        /// number does not find the loaded scan in it; a whole run converted from .raw numbers its spectra from 1. Numbered
        /// from 0, a read by scan number finds a scan, but the next one.
        /// </summary>
        private static string WriteMzmlNumberedFrom(string mzml, string path, int firstScan)
        {
            string text = File.ReadAllText(mzml);
            int start = text.IndexOf("<mzML ", System.StringComparison.Ordinal);
            int end = text.IndexOf("</mzML>", System.StringComparison.Ordinal) + "</mzML>".Length;
            int first = int.Parse(System.Text.RegularExpressions.Regex.Match(text, @"scan=(\d+)").Groups[1].Value, CultureInfo.InvariantCulture);
            string renumbered = System.Text.RegularExpressions.Regex.Replace(text[start..end], @"scan=(\d+)",
                m => "scan=" + (int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture) - first + firstScan).ToString(CultureInfo.InvariantCulture));
            File.WriteAllText(path, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" + renumbered);
            return path;
        }

        /// <summary>
        /// One of the spectra files the tests below search: VA084TQ_6.mzML itself, or written from it in <paramref name="folder"/>.
        /// </summary>
        private static string Fixture(string fixture, string folder) => fixture switch
        {
            "mzML" => TestFile("VA084TQ_6.mzML"),
            "mzMLNumberedFromOne" => WriteMzmlNumberedFrom(TestFile("VA084TQ_6.mzML"), Path.Combine(folder, "VA084TQ_6.mzML"), 1),
            "mzMLNumberedFromZero" => WriteMzmlNumberedFrom(TestFile("VA084TQ_6.mzML"), Path.Combine(folder, "VA084TQ_6.mzML"), 0),
            "mzMLWithTwoIdsSwapped" => SwapSpectrumIds(WriteMzmlNumberedFrom(TestFile("VA084TQ_6.mzML"), Path.Combine(folder, "VA084TQ_6.mzML"), 1), 3, 4),
            "mgf" => WriteMgf(TestFile("VA084TQ_6.mzML"), Path.Combine(folder, "VA084TQ_6.mgf")),
            "mgfWithRepeatedScan" => WriteMgfWithRepeatedScan(TestFile("VA084TQ_6.mzML"), Path.Combine(folder, "VA084TQ_6.mgf")),
            "mzMLWithARepeatedScanNumber" => ReplaceSpectrumId(WriteMzmlNumberedFrom(TestFile("VA084TQ_6.mzML"), Path.Combine(folder, "VA084TQ_6.mzML"), 1), 9, "function=2 process=0 scan=8"),
            "mzMLWithANonNumericScanNumber" => ReplaceSpectrumId(WriteMzmlNumberedFrom(TestFile("VA084TQ_6.mzML"), Path.Combine(folder, "VA084TQ_6.mzML"), 1), 9, "controllerType=0 controllerNumber=1 scan=x9"),
            _ => throw new System.ArgumentException(fixture)
        };

        /// <summary>
        /// Gives the spectrum numbered <paramref name="scan"/> the native ID <paramref name="id"/>. The whole-file read numbers
        /// spectra by position and loads either file below; the dynamic connection's scan index reads the number from the
        /// native ID and throws: on a repeated scan= number, as a Waters conversion writes (function=1 ... scan=1 and
        /// function=2 ... scan=1), an ArgumentException; on a non-numeric one, a FormatException.
        /// </summary>
        private static string ReplaceSpectrumId(string mzml, int scan, string id)
        {
            File.WriteAllText(mzml, File.ReadAllText(mzml).Replace($"id=\"controllerType=0 controllerNumber=1 scan={scan}\"", $"id=\"{id}\""));
            return mzml;
        }

        /// <summary>
        /// Swaps the native IDs of two spectra, so every scan number is found by a read by scan number, two of them in the
        /// other's spectrum.
        /// </summary>
        private static string SwapSpectrumIds(string mzml, int a, int b)
        {
            string Id(int scan) => $"id=\"controllerType=0 controllerNumber=1 scan={scan}\"";
            File.WriteAllText(mzml, File.ReadAllText(mzml).Replace(Id(a), "SWAPPED").Replace(Id(b), Id(a)).Replace("SWAPPED", Id(b)));
            return mzml;
        }

        private static string WriteMgf(string mzml, string mgf)
        {
            MgfMethods.WriteMgf(MsDataFileReader.GetDataFile(mzml).LoadAllStaticData(), mgf);
            return mgf;
        }

        /// <summary>
        /// A LowCID search bins each MS2 spectrum in place for XCorr, so reporter ions read from the loaded spectrum are wrong.
        /// Read one scan at a time from the file, they are the stored scan's: here an mzML numbered from scan 1 and an MGF,
        /// both read that way (see the test below).
        /// </summary>
        [TestCase("mzMLNumberedFromOne")]
        [TestCase("mgf")]
        public static void LowCidReporterIons_ReadOneScanAtATime_AreTheStoredScans(string fixture)
        {
            var searchTask = Toml.ReadFile<SearchTask>(TestFile("TMT11_LowCID_SearchTask.toml"), MetaMorpheusTask.tomlConfig);
            Assert.That(searchTask.CommonParameters.DissociationType, Is.EqualTo(DissociationType.LowCID));
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(searchTask.SearchParameters.MultiplexModId, GlobalVariables.AllModsKnown, out var tag), Is.True);
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestReporterIonsLowCid_" + fixture);

            try
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
                Directory.CreateDirectory(folder);
                string spectra = Fixture(fixture, folder);
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
        /// The search read the whole file again to take reporter ions from it, a second copy beside the loaded one while the
        /// next file loaded. It now reads only the scans reporter ions come from, one at a time, where a read by scan number
        /// returns the loaded scan, and the file whole where it does not: VA084TQ_6.mzML (a slice, numbered from scan=31006),
        /// the same numbered from scan=0 or with two spectra's native IDs swapped (a read by scan number returns another
        /// spectrum), an MGF that repeats a SCANS= value, and an mzML whose scan index the dynamic connection cannot build (a
        /// repeated or non-numeric scan= number; the folder's deletion also shows the file was closed). Either way the values
        /// are the ones the whole-file read gives.
        /// </summary>
        [TestCase("mzMLNumberedFromOne", true)]
        [TestCase("mgf", true)]
        [TestCase("mzML", false)]
        [TestCase("mzMLNumberedFromZero", false)]
        [TestCase("mzMLWithTwoIdsSwapped", false)]
        [TestCase("mgfWithRepeatedScan", false)]
        [TestCase("mzMLWithARepeatedScanNumber", false)]
        [TestCase("mzMLWithANonNumericScanNumber", false)]
        public static void StoredScans_AreReadOneAtATime_WhereAReadByScanNumberFindsTheLoadedScan(string fixture, bool oneAtATime)
        {
            var searchTask = Toml.ReadFile<SearchTask>(TestFile("TMT-Task1-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(searchTask.SearchParameters.MultiplexModId, GlobalVariables.AllModsKnown, out var tag), Is.True);
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestReporterIonsOneAtATime_" + fixture);

            try
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
                Directory.CreateDirectory(folder);
                string spectra = Fixture(fixture, folder);
                var loaded = new MyFileManager(true).LoadFile(spectra, searchTask.CommonParameters);
                var ms2Scans = MetaMorpheusTask.GetMs2Scans(loaded, spectra, searchTask.CommonParameters).ToList();
                Assert.That(ms2Scans, Is.Not.Empty);

                Assert.That(MyFileManager.SetReporterIonIntensitiesFromStoredScans(spectra, loaded, ms2Scans, tag, 1), Is.EqualTo(oneAtATime));
                var read = ms2Scans.Select(s => Formatted(s.IsobaricMassTagReporterIonIntensities)).ToList();

                var storedScanFor = MyFileManager.PairWithStoredScans(loaded, MsDataFileReader.GetDataFile(spectra).LoadAllStaticData());
                var fromLoaded = new List<string[]>();
                var fromWholeFile = new List<string[]>();
                foreach (var scan in ms2Scans)
                {
                    scan.SetIsobaricMassTagReporterIonIntensities(tag);
                    fromLoaded.Add(Formatted(scan.IsobaricMassTagReporterIonIntensities));
                    scan.SetIsobaricMassTagReporterIonIntensities(tag, storedScanFor);
                    fromWholeFile.Add(Formatted(scan.IsobaricMassTagReporterIonIntensities));
                }

                Assert.That(fromWholeFile, Is.Not.EqualTo(fromLoaded), "the fixture only means anything if trimming removes a reporter");
                Assert.That(read, Is.EqualTo(fromWholeFile));
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// MS3 reporter scans read one at a time: MS3_TMT11_Mouse_snip.mzML is a slice (scan=28207 on), so it is read whole;
        /// renumbered from scan 1, a read by scan number finds each loaded MS3 child, and the values are the whole-file read's.
        /// The load keeps 10 peaks in each of 10 windows, as in the MS3 test above, so trimming removes reporters.
        /// </summary>
        [Test]
        public static void Ms3ReporterIons_AreReadOneAtATime_FromTheStoredChildScans()
        {
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag("TMT11", GlobalVariables.AllModsKnown, out var tag), Is.True);
            var commonParameters = new CommonParameters(dissociationType: DissociationType.CID,
                ms3childScanDissociationType: DissociationType.HCD, trimMsMsPeaks: false,
                numberOfPeaksToKeepPerWindow: 10, numberOfWindows: 10);
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestReporterIonsOneAtATime_Ms3");

            try
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
                Directory.CreateDirectory(folder);
                string mzml = WriteMzmlNumberedFrom(TestFile("MS3_TMT11_Mouse_snip.mzML"), Path.Combine(folder, "MS3_TMT11_Mouse_snip.mzML"), 1);
                var loaded = new MyFileManager(true).LoadFile(mzml, commonParameters);
                var ms2Scans = MetaMorpheusTask.GetMs2Scans(loaded, mzml, commonParameters).ToList();
                Assert.That(ms2Scans.Any(s => s.ReporterIonScan.MsnOrder == 3), Is.True);

                Assert.That(MyFileManager.SetReporterIonIntensitiesFromStoredScans(mzml, loaded, ms2Scans, tag, 1), Is.True);
                var read = ms2Scans.Select(s => Formatted(s.IsobaricMassTagReporterIonIntensities)).ToList();

                var storedScanFor = MyFileManager.PairWithStoredScans(loaded, MsDataFileReader.GetDataFile(mzml).LoadAllStaticData());
                var fromLoaded = new List<string[]>();
                var fromWholeFile = new List<string[]>();
                foreach (var scan in ms2Scans)
                {
                    scan.SetIsobaricMassTagReporterIonIntensities(tag);
                    fromLoaded.Add(Formatted(scan.IsobaricMassTagReporterIonIntensities));
                    scan.SetIsobaricMassTagReporterIonIntensities(tag, storedScanFor);
                    fromWholeFile.Add(Formatted(scan.IsobaricMassTagReporterIonIntensities));
                }

                Assert.That(fromWholeFile, Is.Not.EqualTo(fromLoaded), "the fixture only means anything if trimming removes a reporter");
                Assert.That(read, Is.EqualTo(fromWholeFile));
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A file with no scans to set has nothing to read: it counts as read one at a time, and nothing is read whole.
        /// </summary>
        [Test]
        public static void NoScans_NothingIsRead_AndItCountsAsOneAtATime()
        {
            string mzml = TestFile("VA084TQ_6.mzML");
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag("TMT11", GlobalVariables.AllModsKnown, out var tag), Is.True);
            var loaded = MsDataFileReader.GetDataFile(mzml).LoadAllStaticData();

            Assert.That(MyFileManager.SetReporterIonIntensitiesFromStoredScans(mzml, loaded, new List<Ms2ScanWithSpecificMass>(), tag, 1), Is.True);
        }

        /// <summary>
        /// An MGF's read by scan number throws for a number its SCANS= index does not hold. That is caught, and the file is
        /// read whole: here the highest-numbered loaded scan is renumbered past every SCANS= value, so it is left unpaired and
        /// read from its loaded spectrum, and every other scan from its stored block.
        /// </summary>
        [Test]
        public static void Mgf_AskedForAScanNumberItDoesNotHold_IsReadWhole()
        {
            var searchTask = Toml.ReadFile<SearchTask>(TestFile("TMT-Task1-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(searchTask.SearchParameters.MultiplexModId, GlobalVariables.AllModsKnown, out var tag), Is.True);
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestReporterIonsMgfMissingScan");

            try
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
                Directory.CreateDirectory(folder);
                string mgf = Fixture("mgf", folder);
                var loaded = new MyFileManager(true).LoadFile(mgf, searchTask.CommonParameters);
                var ms2Scans = MetaMorpheusTask.GetMs2Scans(loaded, mgf, searchTask.CommonParameters).ToList();
                Assert.That(ms2Scans, Is.Not.Empty);
                ms2Scans.Select(s => s.ReporterIonScan).MaxBy(s => s.OneBasedScanNumber).SetOneBasedScanNumber(999999);

                Assert.That(MyFileManager.SetReporterIonIntensitiesFromStoredScans(mgf, loaded, ms2Scans, tag, 1), Is.False);
                var read = ms2Scans.Select(s => Formatted(s.IsobaricMassTagReporterIonIntensities)).ToList();

                var storedScanFor = MyFileManager.PairWithStoredScans(loaded, MsDataFileReader.GetDataFile(mgf).LoadAllStaticData());
                var fromWholeFile = new List<string[]>();
                foreach (var scan in ms2Scans)
                {
                    scan.SetIsobaricMassTagReporterIonIntensities(tag, storedScanFor);
                    fromWholeFile.Add(Formatted(scan.IsobaricMassTagReporterIonIntensities));
                }
                Assert.That(read, Is.EqualTo(fromWholeFile));
            }
            finally
            {
                if (Directory.Exists(folder)) Directory.Delete(folder, true);
            }
        }

    }
}
