using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.Truncation;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Tests for the truncation search a <see cref="SearchTask"/> runs on its own matches when
    /// <see cref="CommonParameters.AddTruncations"/> is set (TaskLayer/SearchTask/TruncationSearch.cs).
    /// Engine-level tests are in the other Truncation*Tests files.
    /// </summary>
    [TestFixture]
    public class TruncationSearchTests
    {
        /// <summary>
        /// A pipe-ambiguous match has a null FullSequence. Parents are deduped on its sorted hypotheses instead, so two
        /// different ambiguous proteoforms stay separate (both seed parents, #2) while the same pair in either order
        /// is still one proteoform.
        /// </summary>
        [Test]
        public void ProteoformKey_KeepsDifferentPipeAmbiguousMatchesApart()
        {
            var cp = new CommonParameters();
            var scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new[] { 100.0 }, new[] { 1.0 }, false), 1, 2, true,
                Polarity.Positive, 1, new MzRange(0, 1000), "f", MZAnalyzerType.Orbitrap, 1, 1, null, "scan=1"), 500, 1, "f", cp);
            PeptideWithSetModifications Form(string sequence) => new Protein(sequence, "P_" + sequence)
                .Digest(new DigestionParams(protease: "top-down"), new List<Modification>(), new List<Modification>()).First();
            SpectralMatch Match(params string[] sequences)
            {
                var psm = new PeptideSpectralMatch(Form(sequences[0]), 0, 10, 0, scan, cp, new List<Omics.Fragmentation.MatchedFragmentIon>());
                foreach (string sequence in sequences.Skip(1))
                {
                    psm.AddOrReplace(Form(sequence), 10, 0, true, new List<Omics.Fragmentation.MatchedFragmentIon>());
                }
                psm.ResolveAllAmbiguities();
                return psm;
            }

            SpectralMatch ab = Match("PEPTIDE", "PEPTIDEK");
            SpectralMatch ba = Match("PEPTIDEK", "PEPTIDE");
            SpectralMatch cd = Match("ACDEFGH", "ACDEFGHK");
            SpectralMatch plain = Match("MNPQRST");

            Assert.That(ab.FullSequence, Is.Null);
            Assert.That(cd.FullSequence, Is.Null);
            Assert.That(SearchTask.ProteoformKey(ab), Is.Not.EqualTo(SearchTask.ProteoformKey(cd)));
            Assert.That(SearchTask.ProteoformKey(ab), Is.EqualTo(SearchTask.ProteoformKey(ba)));
            Assert.That(SearchTask.ProteoformKey(plain), Is.EqualTo(plain.FullSequence));
        }

        /// <summary>
        /// The parent filter (#3) reads proteoform-level FDR when FdrAnalysisEngine assigned it, and the PSM-level
        /// values when the proteoform-level block is still the unassigned all-2 sentinel (any PSM that was not its
        /// proteoform's representative, including every pipe-ambiguous match after the first).
        /// </summary>
        [Test]
        public void ParentFilter_FallsBackToPsmLevelOnlyWhenProteoformLevelIsUnassigned()
        {
            var cp = new CommonParameters();
            var scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new[] { 100.0 }, new[] { 1.0 }, false), 1, 2, true,
                Polarity.Positive, 1, new MzRange(0, 1000), "f", MZAnalyzerType.Orbitrap, 1, 1, null, "scan=1"), 500, 1, "f", cp);
            var form = new Protein("PEPTIDE", "U").Digest(new DigestionParams(protease: "top-down"),
                new List<Modification>(), new List<Modification>()).First();
            var psm = new PeptideSpectralMatch(form, 0, 10, 0, scan, cp, new List<Omics.Fragmentation.MatchedFragmentIon>());
            psm.PsmFdrInfo = new FdrInfo { QValue = 0.01, QValueNotch = 0.01 };

            psm.PeptideFdrInfo = new FdrInfo();
            Assert.That(SearchTask.ParentFdrInfo(psm), Is.SameAs(psm.PsmFdrInfo));
            Assert.That(SearchTask.PassesTruncationParentFilter(psm, 0.10), Is.True);

            psm.PeptideFdrInfo = new FdrInfo { QValue = 0.5, QValueNotch = 0.5 };
            Assert.That(SearchTask.ParentFdrInfo(psm), Is.SameAs(psm.PeptideFdrInfo));
            Assert.That(SearchTask.PassesTruncationParentFilter(psm, 0.10), Is.False);
        }

        /// <summary>The search's WriteDecoys / WriteContaminants decide which rows reach the truncation files (#17).</summary>
        [TestCase(true, true, 3)]
        [TestCase(false, true, 2)]
        [TestCase(true, false, 2)]
        [TestCase(false, false, 1)]
        public void RowsToWrite_HonorsWriteDecoysAndWriteContaminants(bool writeDecoys, bool writeContaminants, int expectedRows)
        {
            var cp = new CommonParameters();
            var scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new[] { 100.0 }, new[] { 1.0 }, false), 1, 2, true,
                Polarity.Positive, 1, new MzRange(0, 1000), "f", MZAnalyzerType.Orbitrap, 1, 1, null, "scan=1"), 500, 1, "f", cp);
            SpectralMatch Psm(string accession, bool isDecoy, bool isContaminant)
            {
                var protein = new Protein("PEPTIDE", accession, isDecoy: isDecoy, isContaminant: isContaminant);
                var form = protein.Digest(new DigestionParams(protease: "top-down"), new List<Modification>(), new List<Modification>()).First();
                var psm = new PeptideSpectralMatch(form, 0, 10, 0, scan, cp, new List<Omics.Fragmentation.MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                return psm;
            }

            var parameters = new SearchParameters { WriteDecoys = writeDecoys, WriteContaminants = writeContaminants };
            List<SpectralMatch> rows = SearchTask.TruncationRowsToWrite(
                new[] { Psm("T", false, false), Psm("DECOY_T", true, false), Psm("C", false, true), null }, parameters);

            Assert.That(rows, Has.Count.EqualTo(expectedRows));
            Assert.That(rows.Any(p => p.IsDecoy), Is.EqualTo(writeDecoys));
            Assert.That(rows.Any(p => p.IsContaminant), Is.EqualTo(writeContaminants));
        }

        /// <summary>The perf log is benchmarking instrumentation, not a search setting: it is not in the search TOML.</summary>
        [Test]
        public void PerfLogPath_IsNotPartOfTheSearchToml()
        {
            var task = new SearchTask { TruncationPerfLogPath = @"C:\somewhere\perf_log.tsv" };

            string tomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TruncationNoPerfLog.toml");
            Toml.WriteFile(task, tomlPath, MetaMorpheusTask.tomlConfig);
            try
            {
                Assert.That(File.ReadAllText(tomlPath), Does.Not.Contain("PerfLog"));
                Assert.That(File.ReadAllText(tomlPath), Does.Contain("AddTruncations = false"));
            }
            finally
            {
                File.Delete(tomlPath);
            }
        }

        [Test]
        public void PerfLogger_ParsesRunFolderConvention()
        {
            var (phase, dataset, label) = PerfLogger.ParseRunFolderName("2026-05-22_Phase4_KaulichSPE_snippetFromFile");
            Assert.That(phase, Is.EqualTo("Phase4"));
            Assert.That(dataset, Is.EqualTo("KaulichSPE"));
            Assert.That(label, Is.EqualTo("snippetFromFile"));

            // Non-conforming names degrade gracefully.
            var (phase2, dataset2, label2) = PerfLogger.ParseRunFolderName("TruncationE2E");
            Assert.That(phase2, Is.EqualTo("Manual"));
            Assert.That(dataset2, Is.EqualTo(""));
            Assert.That(label2, Is.EqualTo(""));
        }

        [Test]
        public void PerfLogger_WritesHeaderOnceThenAppends()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "perf_log_unit.tsv");
            if (File.Exists(path)) File.Delete(path);
            try
            {
                PerfLogger.Append(path, new TruncationPerfMetrics { NPsmsEmitted = 5, NTruncationsNterm = 3 });
                PerfLogger.Append(path, new TruncationPerfMetrics { NPsmsEmitted = 9, NTruncationsCterm = 2 });

                string[] lines = File.ReadAllLines(path);
                Assert.That(lines.Length, Is.EqualTo(3), "Expected one header + two data rows.");
                Assert.That(lines[0].Split('\t'), Does.Contain("n_psms_emitted"));
                int col = System.Array.IndexOf(lines[0].Split('\t'), "n_psms_emitted");
                Assert.That(lines[1].Split('\t')[col], Is.EqualTo("5"));
                Assert.That(lines[2].Split('\t')[col], Is.EqualTo("9"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void PerfLogger_EmptyOrNullFolder_ReturnsManualDefaults()
        {
            foreach (string input in new[] { "", null })
            {
                var (phase, dataset, label) = PerfLogger.ParseRunFolderName(input);
                Assert.That(phase, Is.EqualTo("Manual"));
                Assert.That(dataset, Is.EqualTo(""));
                Assert.That(label, Is.EqualTo(""));
            }
        }

        [Test]
        public void PerfLogger_SanitizesTabsAndNewlinesInStringCells()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "perf_log_sanitize.tsv");
            if (File.Exists(path)) File.Delete(path);
            try
            {
                // A raw tab/newline in a string field would otherwise split into extra columns / extra rows.
                PerfLogger.Append(path, new TruncationPerfMetrics { RunLabel = "a\tb\nc", DatasetTag = "x\ty" });

                string[] lines = File.ReadAllLines(path);
                int headerCols = lines[0].Split('\t').Length;
                Assert.That(lines.Length, Is.EqualTo(2));                          // header + one row (no extra line from \n)
                Assert.That(lines[1].Split('\t').Length, Is.EqualTo(headerCols));  // no extra columns from \t
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// A search TOML with AddTruncations = true runs the truncation search inside the search: the run list holds one
        /// task, and both truncation files, with the standard psmtsv header, are written beside the search's AllPSMs. Every
        /// confident search PSM is inherited exactly once as full-length, on its own precursor, with the same Essential
        /// Sequence as AllPSMs, and none of those precursors is also reported as a truncation (#4a). Uses the tiny top-down
        /// fixture, and the perf-log path.
        /// </summary>
        [Test]
        public void SearchWithAddTruncations_WritesTruncationOutputsInTheSearchFolder()
        {
            var (data, db, toml) = TopDownFixture();
            string outDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TruncationInSearch");
            if (Directory.Exists(outDirectory)) Directory.Delete(outDirectory, true);
            Directory.CreateDirectory(outDirectory);
            try
            {
                string addTruncationsToml = Path.Combine(outDirectory, "AddTruncations.toml");
                string tomlText = File.ReadAllText(toml);
                Assert.That(tomlText, Does.Contain("AddTruncations = false"));
                File.WriteAllText(addTruncationsToml, tomlText.Replace("AddTruncations = false", "AddTruncations = true"));
                var searchTask = Toml.ReadFile<SearchTask>(addTruncationsToml, MetaMorpheusTask.tomlConfig);
                Assert.That(searchTask.CommonParameters.AddTruncations, Is.True);
                string perfLogPath = Path.Combine(outDirectory, "perf_log.tsv");
                searchTask.TruncationPerfLogPath = perfLogPath;

                var engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task1-SearchTask", searchTask) },
                    new List<string> { data }, new List<DbForTask> { new DbForTask(db, false) }, outDirectory);
                Assert.DoesNotThrow(() => engine.Run());

                string searchFolder = Path.Combine(outDirectory, "Task1-SearchTask");
                Assert.That(Directory.GetDirectories(outDirectory).Select(Path.GetFileName),
                    Is.EquivalentTo(new[] { "Task1-SearchTask", "Task Settings" }), "the truncation search is not a second task");
                AssertTruncationOutputsWellFormed(searchFolder);

                string psmsPath = Path.Combine(searchFolder, SearchTask.TruncatedPsmsFileName);
                string[] psmLines = File.ReadAllLines(psmsPath);
                Assert.That(psmLines.Length, Is.GreaterThan(1), "Expected at least one pooled PSM row.");
                Assert.That(File.ReadAllLines(Path.Combine(searchFolder, SearchTask.TruncatedProteoformsFileName)).Length, Is.GreaterThan(1),
                    "Expected at least one pooled proteoform row.");

                var pass1 = ReadTsv(Path.Combine(searchFolder, "AllPSMs.psmtsv"));
                var truncation = ReadTsv(psmsPath);
                var confident = pass1.Where(r => r["Decoy/Contaminant/Target"] == "T" && double.Parse(r["QValue"], CultureInfo.InvariantCulture) <= 0.01).ToList();
                Assert.That(confident, Is.Not.Empty);
                foreach (var row in confident)
                {
                    var samePrecursor = truncation.Where(t => t["Scan Number"] == row["Scan Number"]
                        && Math.Abs(double.Parse(t["Precursor Mass"], CultureInfo.InvariantCulture) - double.Parse(row["Precursor Mass"], CultureInfo.InvariantCulture)) < 0.01).ToList();
                    Assert.That(samePrecursor.Count(t => t["Description"] == TruncationPass3.FullLength), Is.EqualTo(1),
                        $"scan {row["Scan Number"]} should be inherited once as full-length");
                    Assert.That(samePrecursor, Has.Count.EqualTo(1), $"scan {row["Scan Number"]} also carries a truncation");
                    Assert.That(samePrecursor[0]["Essential Sequence"], Is.EqualTo(row["Essential Sequence"]));
                }

                // The perf log got one row (header + data row) counting the PSM file's rows, and the CandidateRanks side file.
                string[] perfLines = File.ReadAllLines(perfLogPath);
                Assert.That(perfLines.Length, Is.EqualTo(2));
                int nPsmsCol = System.Array.IndexOf(perfLines[0].Split('\t'), "n_psms_emitted");
                Assert.That(int.Parse(perfLines[1].Split('\t')[nPsmsCol]), Is.EqualTo(psmLines.Length - 1));
                Assert.That(File.Exists(Path.Combine(searchFolder, "CandidateRanks.tsv")), Is.True);
            }
            finally
            {
                if (Directory.Exists(outDirectory)) Directory.Delete(outDirectory, true);
            }
        }

        /// <summary>Without AddTruncations the search writes no truncation files.</summary>
        [Test]
        public void SearchWithoutAddTruncations_WritesNoTruncationFiles()
        {
            var (data, db, toml) = TopDownFixture();
            string outDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TruncationSearchOnly");
            if (Directory.Exists(outDirectory)) Directory.Delete(outDirectory, true);
            try
            {
                var searchTask = Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig);
                Assert.That(searchTask.CommonParameters.AddTruncations, Is.False);
                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task1-SearchTask", searchTask) },
                    new List<string> { data }, new List<DbForTask> { new DbForTask(db, false) }, outDirectory).Run();

                Assert.That(Directory.GetFiles(outDirectory, "AllPSMs.psmtsv", SearchOption.AllDirectories), Is.Not.Empty);
                Assert.That(Directory.GetFiles(outDirectory, "AllTruncated*", SearchOption.AllDirectories), Is.Empty);
            }
            finally
            {
                if (Directory.Exists(outDirectory)) Directory.Delete(outDirectory, true);
            }
        }

        /// <summary>
        /// A search match with no FDR assigned cannot pass the parent filter (#3), so it seeds no parent; with no parents
        /// the truncation search says so and still writes both result files, header-only (#16).
        /// </summary>
        [Test]
        public void UnscoredMatchSeedsNoParent_AndTheFilesAreHeaderOnly()
        {
            var (data, _, _) = TopDownFixture();
            string outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TruncationNoParents");
            if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            Directory.CreateDirectory(outDir);
            var warnings = new List<string>();
            EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += onWarn;
            try
            {
                var cp = new CommonParameters();
                var scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new[] { 100.0 }, new[] { 1.0 }, false), 1, 2, true,
                    Polarity.Positive, 1, new MzRange(0, 1000), "f", MZAnalyzerType.Orbitrap, 1, 1, null, "scan=1"), 500, 1, data, cp);
                var form = new Protein("PEPTIDE", "U").Digest(new DigestionParams(protease: "top-down"),
                    new List<Modification>(), new List<Modification>()).First();
                var unscored = new PeptideSpectralMatch(form, 0, 10, 0, scan, cp, new List<Omics.Fragmentation.MatchedFragmentIon>());

                new SearchTask().RunTruncationSearch(new List<SpectralMatch> { unscored, null }, outDir, new List<string> { data },
                    "Task1-SearchTask", new FileSpecificParameters[1], new List<Modification>());

                Assert.That(warnings, Has.Some.Contains("found no proteoforms to seed it"));
                string header = SpectralMatch.GetTabSeparatedHeader().TrimEnd();
                foreach (string file in new[] { SearchTask.TruncatedPsmsFileName, SearchTask.TruncatedProteoformsFileName })
                {
                    Assert.That(File.ReadAllLines(Path.Combine(outDir, file)).Select(l => l.TrimEnd()), Is.EqualTo(new[] { header }), file);
                }
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= onWarn;
                if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            }
        }

        /// <summary>AddTruncations no longer expands the database with truncation products; the truncation search replaces that.</summary>
        [Test]
        public void DatabaseLoading_WithAddTruncations_DoesNotExpandTheDatabase()
        {
            var (_, db, _) = TopDownFixture();
            List<Protein> proteins = DatabaseLoadingEngine.LoadProteinDb(db, true, UsefulProteomicsDatabases.DecoyType.None, new List<string>(),
                false, out _, out _, new CommonParameters(addTruncations: true)).ToList();

            Assert.That(proteins, Is.Not.Empty);
            Assert.That(proteins.SelectMany(p => p.TruncationProducts), Is.Empty);
        }

        // ---------- helpers ----------

        private static (string data, string db, string toml) TopDownFixture()
        {
            string dir = TestContext.CurrentContext.TestDirectory;
            return (Path.Combine(dir, "TopDownTestData", "slicedTDYeast.mzML"),
                    Path.Combine(dir, "TestData", "smalldb.fasta"),
                    Path.Combine(dir, "TopDownTestData", "TopDownSearchToml.toml"));
        }

        /// <summary>Rows of a psmtsv as column-name -> cell dictionaries.</summary>
        private static List<Dictionary<string, string>> ReadTsv(string path)
        {
            string[] lines = File.ReadAllLines(path);
            string[] header = lines[0].Split('\t');
            return lines.Skip(1)
                .Select(line => line.Split('\t'))
                .Select(cells => header.Select((name, i) => (name, cell: i < cells.Length ? cells[i] : ""))
                    .GroupBy(c => c.name).ToDictionary(g => g.Key, g => g.First().cell))
                .ToList();
        }

        private static void AssertTruncationOutputsWellFormed(string searchFolder)
        {
            string psms = Path.Combine(searchFolder, SearchTask.TruncatedPsmsFileName);
            string proteoforms = Path.Combine(searchFolder, SearchTask.TruncatedProteoformsFileName);
            Assert.That(File.Exists(psms), Is.True, "AllTruncatedPSMs.psmtsv not written");
            Assert.That(File.Exists(proteoforms), Is.True, "AllTruncatedProteoforms.psmtsv not written");
            string header = SpectralMatch.GetTabSeparatedHeader().TrimEnd();
            Assert.That(File.ReadAllLines(psms)[0].TrimEnd(), Is.EqualTo(header));
            Assert.That(File.ReadAllLines(proteoforms)[0].TrimEnd(), Is.EqualTo(header));
        }
    }
}
