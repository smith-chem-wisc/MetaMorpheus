using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Readers;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Tests for the opt-in SDRF-Proteomics output: <see cref="SearchParameters.WriteSdrf"/>, the
    /// validation that refuses a run without sample metadata, and the file the search writes.
    ///
    /// The division of labour these tests assume: everything that reasons ABOUT SDRF -- the format,
    /// the vocabulary, structural validation -- lives in mzLib and is tested there against a curated
    /// corpus of 1,236 real documents. What is testable only from here is the ADAPTER: whether
    /// MetaMorpheus hands mzLib the right facts, and whether opting in means what it claims to mean.
    /// So these assert on values in named columns, not on SDRF's rules.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class SdrfOutputTests
    {
        private const string SdrfFileName = "experiment.sdrf.tsv";

        #region Opting in and out

        /// <summary>
        /// The default is off, and off writes nothing. This is what makes the requirement that
        /// opting in imposes (see the validation tests below) tolerable at all.
        /// </summary>
        [Test]
        public static void NoSdrfIsWritten_WhenTheSearchDidNotAskForOne()
        {
            string folder = SetUpIsolatedRun(nameof(NoSdrfIsWritten_WhenTheSearchDidNotAskForOne),
                out string spectraPath, out DbForTask database);

            // Deliberately NO ExperimentalDesign.tsv: a search that did not ask for an SDRF must not
            // acquire a new prerequisite because this feature exists.
            var task = BuildSearchTask(writeSdrf: false);
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(output);

            task.RunTask(output, new List<DbForTask> { database }, new List<string> { spectraPath }, "no-sdrf");

            Assert.That(File.Exists(Path.Combine(output, SdrfFileName)), Is.False,
                SdrfFileName + " must not be written when WriteSdrf is false.");
            Assert.That(File.Exists(Path.Combine(output, "AllPSMs.psmtsv")), Is.True,
                "The search itself should still have produced results.");

            Directory.Delete(folder, true);
        }

        [Test]
        public static void WriteSdrfDefaultsToOff()
        {
            Assert.That(new SearchParameters().WriteSdrf, Is.False,
                "SDRF output is opt-in; defaulting it on would make every existing search acquire a " +
                "new prerequisite.");
        }

        #endregion

        #region Validation happens before the search, not after it (D17)

        /// <summary>
        /// Opting in does NOT block the search when sample metadata is absent (D21).
        ///
        /// The earlier design refused the run. It was stricter than the specification -- which marks
        /// organism part required but explicitly permits the reserved words -- and stricter than the
        /// community, a fifth of whose curated cells are one. It was also self-defeating: a user who
        /// is blocked unticks the box, and then there is no file at all.
        ///
        /// So the run proceeds, the file is written, and the gap is REPORTED rather than hidden.
        /// That last part is what keeps this honest rather than merely permissive.
        /// </summary>
        [Test]
        public static void OptingInWithoutAnExperimentalDesign_StillRunsAndStillWrites()
        {
            string folder = SetUpIsolatedRun(nameof(OptingInWithoutAnExperimentalDesign_StillRunsAndStillWrites),
                out string spectraPath, out DbForTask database);

            var task = BuildSearchTask(writeSdrf: true);
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(output);

            Assert.DoesNotThrow(() =>
                task.RunTask(output, new List<DbForTask> { database }, new List<string> { spectraPath }, "sdrf-no-design"));

            Assert.That(File.Exists(Path.Combine(output, "AllPSMs.psmtsv")), Is.True,
                "The search itself must not be blocked by missing sample metadata.");

            string sdrfPath = Path.Combine(output, SdrfFileName);
            Assert.That(File.Exists(sdrfPath), Is.True, "The SDRF is still written.");

            var document = new SdrfDocument(sdrfPath);
            document.LoadResults();

            var missing = SdrfValidator.RequiredColumns.Where(c => !document.Header.Contains(c)).ToList();
            Assert.That(missing, Is.Empty,
                "Every required column must be present even with no design file: " + string.Join(", ", missing));
            Assert.That(document.Results.Single()["characteristics[organism part]"], Is.EqualTo("not available"),
                "Absent metadata is stated with the reserved word, not invented and not omitted.");

            // With no design nobody established the replicate or fraction numbers, so the file must not
            // claim 1 (mzLib #1378 made them nullable for exactly this).
            foreach (string column in new[]
                     {
                         "characteristics[biological replicate]", "comment[technical replicate]",
                         "comment[fraction identifier]"
                     })
                Assert.That(document.Results.Single()[column], Is.EqualTo("not available"),
                    column + " is unknown without a design, and must say so.");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// A design file that exists but cannot be parsed does not block either -- it degrades to
        /// the same position as no design at all, and says so.
        ///
        /// Note this is not the only consequence of a malformed design: quantification independently
        /// refuses to run on one. SDRF is not the component that should be enforcing that.
        /// </summary>
        [Test]
        public static void OptingInWithAnUnusableExperimentalDesign_StillRunsAndStillWrites()
        {
            string folder = SetUpIsolatedRun(nameof(OptingInWithAnUnusableExperimentalDesign_StillRunsAndStillWrites),
                out string spectraPath, out DbForTask database);

            // Four cells where five are required -- the same malformed shape the calibration tests
            // pin, written against this run's own spectra file name.
            File.WriteAllLines(
                Path.Combine(Path.GetDirectoryName(spectraPath)!, GlobalVariables.ExperimentalDesignFileName),
                new[]
                {
                    "FileName\tCondition\tBiorep\tFraction\tTechrep",
                    Path.GetFileName(spectraPath) + "\tcondition\t1\t1"
                });

            var task = BuildSearchTask(writeSdrf: true);
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(output);

            Assert.DoesNotThrow(() =>
                task.RunTask(output, new List<DbForTask> { database }, new List<string> { spectraPath }, "sdrf-bad-design"));

            Assert.That(File.Exists(Path.Combine(output, SdrfFileName)), Is.True);

            Directory.Delete(folder, true);
        }

        #endregion

        #region Values that must not be asserted when they are not known

        /// <summary>
        /// The acquisition method is READ, not assumed.
        ///
        /// It was hardcoded to DDA. MetaMorpheus does not support DIA searching today, so nothing
        /// was visibly wrong -- but a hardcoded term is wrong the day support arrives, and wrong in
        /// the worst way: the column comes out 100% filled with a false CV term, so SdrfCoverage
        /// cannot see it and a reader has no reason to doubt it.
        /// </summary>
        [Test]
        public static void AcquisitionMethodIsReadFromTheSearch_NotHardcoded()
        {
            Assert.That(InvokeAcquisitionMethod(null)?.Name, Is.EqualTo("Data-dependent acquisition"),
                "No DIA parameters means DDA, which is every MetaMorpheus search today.");

            var dia = new CommonParameters(diaParameters: DiaParams(DIAanalysisType.DIA));
            Assert.That(InvokeAcquisitionMethod(dia)?.Accession, Is.EqualTo("PRIDE:0000450"),
                "A DIA search must say Data-independent acquisition.");

            var isd = new CommonParameters(diaParameters: DiaParams(DIAanalysisType.ISD));
            Assert.That(InvokeAcquisitionMethod(isd), Is.Null,
                "In-source decay has no acquisition-method term. Leave it unresolved rather than " +
                "borrowing the nearest-looking one.");
        }

        /// <summary>
        /// A labelled search does not get to claim it was label free.
        ///
        /// SDRF wants one row per sample per channel, and the writer emits one row per file. Saying
        /// "label free sample" for a SILAC or TMT run is a confident falsehood; until rows are
        /// expanded per channel, the reserved word is the truth.
        /// </summary>
        [Test]
        public static void LabelledSearchesDoNotClaimToBeLabelFree()
        {
            Assert.That(InvokeLabel(new SearchParameters())?.Name, Is.EqualTo("label free sample"));

            Assert.That(InvokeLabel(new SearchParameters { DoMultiplexQuantification = true }), Is.Null,
                "Isobaric labelling: rows are not yet expanded per channel.");

            var silac = new SearchParameters
            {
                SilacLabels = new List<SilacLabel> { new('K', 'a', "C{6}H{12}N{2}O{1}", 6.020129) }
            };
            Assert.That(InvokeLabel(silac), Is.Null,
                "SILAC: which sample carries which label is not something the search knows.");
        }

        #endregion

        #region What the search actually writes

        /// <summary>
        /// The whole feature, end to end: opt in with a usable design and an SDRF lands with the
        /// results, readable by the library that will later pool it.
        /// </summary>
        [Test]
        public static void AnSdrfIsWrittenWithTheResults_AndCanBeReadBack()
        {
            string output = RunSearchWritingSdrf(nameof(AnSdrfIsWrittenWithTheResults_AndCanBeReadBack),
                out string folder, out string spectraPath);

            string sdrfPath = Path.Combine(output, SdrfFileName);
            Assert.That(File.Exists(sdrfPath), Is.True,
                SdrfFileName + " should be written into the search's own output folder.");

            var document = new SdrfDocument(sdrfPath);
            document.LoadResults();

            Assert.That(document.Results.Count, Is.EqualTo(1),
                "One row per spectra file, and this run had one file.");
            Assert.That(document.Header.Count, Is.GreaterThan(0));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// The adapter's actual job: put the facts THIS search used into the right columns. A row
        /// that reported the defaults rather than the run would be worse than no row.
        /// </summary>
        [Test]
        public static void TheWrittenSdrfCarriesThisSearchsOwnAssayParameters()
        {
            string output = RunSearchWritingSdrf(nameof(TheWrittenSdrfCarriesThisSearchsOwnAssayParameters),
                out string folder, out string spectraPath);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();
            SdrfRow row = document.Results.Single();

            Assert.That(row["comment[cleavage agent details]"], Does.Contain("Trypsin").IgnoreCase,
                "The search digested with trypsin, so the row has to say so.");
            Assert.That(row["comment[software]"], Does.Contain("MetaMorpheus"),
                "The file must record what produced it.");
            Assert.That(row["technology type"], Is.EqualTo("proteomic profiling by mass spectrometry"));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// A per-file override reaches the row for THAT file, and only that file.
        ///
        /// MetaMorpheus lets a spectra file carry its own .toml overriding protease, tolerances and
        /// dissociation type, and the adapter looks those up per file. When that lookup misses, the
        /// failure is silent: every row reports the task-level values, the document still validates,
        /// and the fill rate is unchanged. Only two files with DIFFERENT parameters can see it --
        /// which is why every other end-to-end case here, all single-file, missed a lookup that
        /// compared a bare file name against a stored full path and therefore never matched.
        /// </summary>
        [Test]
        public static void APerFileParameterOverrideReachesThatFilesRow()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "SdrfOutput_" + nameof(APerFileParameterOverrideReachesThatFilesRow));
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);

            string source = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string taskLevelFile = Path.Combine(folder, "tryptic.mzml");
            string overriddenFile = Path.Combine(folder, "aspN.mzml");
            File.Copy(source, taskLevelFile, true);
            File.Copy(source, overriddenFile, true);

            // A file-specific toml is named after its spectra file and sits beside it. Only the
            // second file gets one, so the two rows must not come out the same.
            File.WriteAllLines(Path.Combine(folder, "aspN.toml"), new[] { "Protease = \"Asp-N\"" });

            var database = new DbForTask(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"), false);

            ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo>
            {
                new(taskLevelFile, "condition", 0, 0, 0),
                new(overriddenFile, "condition", 1, 0, 0)
            });

            var task = BuildSearchTask(writeSdrf: true);
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(output);

            task.RunTask(output, new List<DbForTask> { database },
                new List<string> { taskLevelFile, overriddenFile }, "sdrf-file-specific");

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            Assert.That(document.Results.Count, Is.EqualTo(2), "One row per spectra file.");

            SdrfRow RowFor(string path) => document.Results
                .Single(r => r["comment[data file]"] == Path.GetFileName(path));

            Assert.That(RowFor(taskLevelFile)["comment[cleavage agent details]"],
                Does.Contain("Trypsin").IgnoreCase,
                "The file with no toml of its own is digested with the task-level protease.");
            Assert.That(RowFor(overriddenFile)["comment[cleavage agent details]"],
                Does.Contain("Asp-N").IgnoreCase,
                "This file was digested with Asp-N by its own toml, so its row has to say Asp-N -- " +
                "reporting the task-level protease here is exactly the flattening this column exists " +
                "to prevent.");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// comment[data file] names the acquired file (sdrf D46). In a search-only run the search reads that
        /// same file, so there is no comment[searched data file] column.
        /// </summary>
        [Test]
        public static void ASearchOnlySdrfNamesTheAcquiredFileAndNoSearchedFile()
        {
            string output = RunSearchWritingSdrf(nameof(ASearchOnlySdrfNamesTheAcquiredFileAndNoSearchedFile),
                out string folder, out string spectraPath);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();
            SdrfRow row = document.Results.Single();

            Assert.That(row["comment[data file]"], Is.EqualTo(Path.GetFileName(spectraPath)));
            Assert.That(document.Header, Does.Not.Contain("comment[searched data file]"),
                "Nothing was calibrated, so the searched file is the acquired one.");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Calibrate -> Search: the search reads the -calib derivative. comment[data file] still names the
        /// acquired file, extension and all, and comment[searched data file] names the derivative (sdrf D46,
        /// pcruzparri's review of #2816). The instrument comes from the search's own load of the file.
        /// </summary>
        [Test]
        public static void AfterCalibrationTheSdrfNamesTheAcquiredFileAndTheCalibratedFileItSearched()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SdrfOutput_Calib");
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            string spectraPath = Path.Combine(folder, "sample1.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML"), spectraPath, true);
            string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta");
            ExperimentalDesign.WriteExperimentalDesignToFile(
                new List<SpectraFileInfo> { new(spectraPath, "condition", 0, 0, 0) });

            var engine = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Task1-Calibrate", new CalibrationTask()), ("Task2-Search", BuildSearchTask(writeSdrf: true)) },
                new List<string> { spectraPath }, new List<DbForTask> { new DbForTask(database, false) }, folder);
            engine.Run();

            var document = new SdrfDocument(Path.Combine(folder, "Task2-Search", SdrfFileName));
            document.LoadResults();
            SdrfRow row = document.Results.Single();

            Assert.That(row["comment[data file]"], Is.EqualTo("sample1.mzML"), "the acquired file, as deposited");
            Assert.That(row["comment[searched data file]"], Is.EqualTo("sample1-calib.mzML"), "the file the search read");
            Assert.That(row["comment[instrument]"], Does.Contain("LTQ Orbitrap Velos"),
                "read from the search's own load of the file, not a second read");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Every column SDRF-Proteomics requires is present, and mzLib's own validator passes the
        /// document the search wrote.
        ///
        /// A column that is MISSING is invisible to SdrfCoverage, which measures the fill rate of
        /// columns that exist, so an absent column has no fill rate to report. SdrfValidator does
        /// check, and asking it directly -- rather than a copy of its column list -- means this
        /// test cannot drift from the rules that actually govern.
        /// </summary>
        [Test]
        public static void TheWrittenSdrfCarriesEveryRequiredColumn()
        {
            string output = RunSearchWritingSdrf(nameof(TheWrittenSdrfCarriesEveryRequiredColumn),
                out string folder, out string spectraPath);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            var missing = SdrfValidator.RequiredColumns.Where(c => !document.Header.Contains(c)).ToList();

            Assert.That(missing, Is.Empty,
                "SDRF-Proteomics requires these columns of every document, and they are absent from " +
                "the one this search wrote: " + string.Join(", ", missing) + ". A document missing a " +
                "required column cannot be pooled across experiments, which is the only reason to " +
                "write one.");

            SdrfValidationResult validation = SdrfValidator.Validate(document);
            Assert.That(validation.Errors, Is.Empty,
                "mzLib's validator rejects the SDRF this search wrote: " + validation);

            // mzLib demotes these to warnings because some curated files lack them, but the specification
            // requires them: disease and cell type among them, written `not available` when unknown.
            Assert.That(SdrfValidator.RecommendedColumns.Where(c => !document.Header.Contains(c)), Is.Empty);
            Assert.That(document.Results[0]["characteristics[disease]"], Is.EqualTo("not available"));
            Assert.That(document.Results[0]["comment[label]"], Is.EqualTo("label free sample"),
                "written bare (sdrf D25, SdrfLabelForm.Bare): what the community writes, and quantms crashes on the accessioned form");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// The organism is the first TARGET protein's NCBI taxonomy id, as an NCBITaxon term.
        ///
        /// Driven through the resolver rather than a search because no search in this suite can
        /// reach it: MetaMorpheus hands LoadProteinFasta explicit UniProt regexes, which skips mzLib's
        /// header detection and with it the OX= regex, so a FASTA search never populates
        /// NcbiTaxonomyId today. #2782 records that and deliberately leaves it unchanged, because
        /// parsing OX= moves other output. A decoy listed first must not decide the organism.
        /// </summary>
        [Test]
        public static void TheOrganismIsTheTargetsOneTaxonNeverADecoysOrAContaminants()
        {
            var taxon = new List<DatabaseReference>
            {
                new(Protein.NcbiTaxonomyDatabaseReferenceType, "559292", new List<Tuple<string, string>>())
            };
            var decoyTaxon = new List<DatabaseReference>
            {
                new(Protein.NcbiTaxonomyDatabaseReferenceType, "9606", new List<Tuple<string, string>>())
            };

            var task = new PostSearchAnalysisTask
            {
                Parameters = new PostSearchAnalysisParameters
                {
                    BioPolymerList = new List<IBioPolymer>
                    {
                        new Protein("PEPTIDE", "DECOY_P1", organism: "Homo sapiens", isDecoy: true, databaseReferences: decoyTaxon),
                        // MetaMorpheusContaminants.xml carries NCBI Taxonomy 9913 (Alexander-Sol's review of #2816).
                        new Protein("PEPTIDER", "CONTAM_P02769", organism: "Bos taurus", isContaminant: true,
                            databaseReferences: new List<DatabaseReference> { new(Protein.NcbiTaxonomyDatabaseReferenceType, "9913", new List<Tuple<string, string>>()) }),
                        new Protein("PEPTIDEK", "P38266", organism: "Saccharomyces cerevisiae", databaseReferences: taxon)
                    }
                }
            };

            var organism = (CvParam)typeof(PostSearchAnalysisTask)
                .GetMethod("ResolveOrganismFromSearchDatabase", BindingFlags.NonPublic | BindingFlags.Instance)!
                .Invoke(task, null);

            Assert.That(organism?.Accession, Is.EqualTo("NCBITaxon:559292"));
            Assert.That(organism?.Name, Is.EqualTo("Saccharomyces cerevisiae"));

            // Two target organisms: not one of them at random, but none.
            task.Parameters.BioPolymerList.Add(new Protein("PEPTIDEKK", "P99999", organism: "Homo sapiens", databaseReferences: decoyTaxon));
            Assert.That(typeof(PostSearchAnalysisTask)
                .GetMethod("ResolveOrganismFromSearchDatabase", BindingFlags.NonPublic | BindingFlags.Instance)!
                .Invoke(task, null), Is.Null);
        }

        #endregion

        #region A failed SDRF never costs the search

        /// <summary>
        /// The SDRF writer is the one writer in PostSearchAnalysisTask that is wrapped: it is metadata
        /// about results, so a failure writing it must leave the finished search standing and leave a
        /// crash report behind, not throw. Here the file cannot be created because a folder already
        /// holds its name.
        /// </summary>
        [Test]
        public static void ASearchSurvivesAnSdrfThatCannotBeWritten()
        {
            string folder = SetUpIsolatedRun(nameof(ASearchSurvivesAnSdrfThatCannotBeWritten),
                out string spectraPath, out DbForTask database);

            var task = BuildSearchTask(writeSdrf: true);
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(Path.Combine(output, SdrfFileName));

            Assert.DoesNotThrow(() =>
                task.RunTask(output, new List<DbForTask> { database }, new List<string> { spectraPath }, "sdrf-blocked"));

            Assert.That(File.Exists(Path.Combine(output, "AllPSMs.psmtsv")), Is.True,
                "The search's own results are written regardless.");
            Assert.That(File.Exists(Path.Combine(output, "SdrfWriter_crash.txt")), Is.True,
                "The failure is reported, not swallowed silently.");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// The instrument comes from the search's own load of each file; the SDRF never opens a data file again
        /// (Alexander-Sol's review of #2816: a RAW's SourceFile hashes the whole file). A file the search did
        /// not record, even a readable one, resolves to no instrument.
        /// </summary>
        [Test]
        public static void TheInstrumentIsTheOneTheSearchReadAndNoFileIsOpenedAgain()
        {
            string readable = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML");
            var velos = new CvParam("MS", "MS:1001742", "LTQ Orbitrap Velos", "");
            var task = new PostSearchAnalysisTask
            {
                Parameters = new PostSearchAnalysisParameters
                {
                    InstrumentModelsByFile = new Dictionary<string, CvParam>(StringComparer.OrdinalIgnoreCase) { ["a.raw"] = velos }
                }
            };
            var resolve = typeof(PostSearchAnalysisTask).GetMethod("ResolveInstrument", BindingFlags.NonPublic | BindingFlags.Instance)!;

            Assert.That(resolve.Invoke(task, new object[] { "a.raw" }), Is.EqualTo(velos));
            Assert.That(resolve.Invoke(task, new object[] { readable }), Is.Null, "not read again, though it could be");
        }

        /// <summary>
        /// Not every data file names its instrument (nbollis's review of #2817). An MGF carries no instrument model,
        /// so the search records none, and the row says `not available` rather than failing or guessing.
        /// </summary>
        [Test]
        public static void AFileThatNamesNoInstrumentWritesNotAvailable()
        {
            string folder = SetUpIsolatedRun(nameof(AFileThatNamesNoInstrumentWritesNotAvailable), out _, out DbForTask database);
            string mgfPath = Path.Combine(folder, "ok.mgf");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "ok.mgf"), mgfPath, true);
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(output);

            BuildSearchTask(writeSdrf: true).RunTask(output, new List<DbForTask> { database }, new List<string> { mgfPath }, "sdrf");

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();
            Assert.That(document.Results.Single()["comment[instrument]"], Is.EqualTo("not available"));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// An ExperimentalDesign.tsv another program holds open (Excel) costs the SDRF its design, not the SDRF:
        /// the read warns and describes the search only, as the pre-run warning promised.
        /// </summary>
        [Test]
        public static void ALockedDesignFileLeavesTheSdrfDescribingTheSearchOnly()
        {
            string folder = SetUpIsolatedRun(nameof(ALockedDesignFileLeavesTheSdrfDescribingTheSearchOnly), out string spectraPath, out _);
            ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { new(spectraPath, "condition", 0, 0, 0) });
            string designPath = Path.Combine(folder, GlobalVariables.ExperimentalDesignFileName);
            var task = new PostSearchAnalysisTask { Parameters = new PostSearchAnalysisParameters { CurrentRawFileList = new List<string> { spectraPath } } };
            var read = typeof(PostSearchAnalysisTask).GetMethod("ReadExperimentalDesignIfPresent", BindingFlags.NonPublic | BindingFlags.Instance)!;

            using (new FileStream(designPath, FileMode.Open, FileAccess.ReadWrite, FileShare.None))
            {
                object byStem = null;
                Assert.DoesNotThrow(() => byStem = read.Invoke(task, null));
                Assert.That((System.Collections.IDictionary)byStem!, Is.Empty);
            }
            Assert.That(((System.Collections.IDictionary)read.Invoke(task, null)!).Count, Is.EqualTo(1), "read once released");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// A search that describes no spectra file writes no SDRF and says so, rather than handing the
        /// builder an empty table.
        /// </summary>
        [Test]
        public static void ASearchWithNoSpectraFilesWritesNoSdrfAndSaysSo()
        {
            string output = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(ASearchWithNoSpectraFilesWritesNoSdrfAndSaysSo));
            if (Directory.Exists(output)) Directory.Delete(output, true);
            Directory.CreateDirectory(output);
            var task = new PostSearchAnalysisTask
            {
                Parameters = new PostSearchAnalysisParameters
                {
                    OutputFolder = output,
                    SearchTaskId = "no-files",
                    SearchParameters = new SearchParameters(),
                    CurrentRawFileList = new List<string>()
                }
            };

            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                typeof(PostSearchAnalysisTask).GetMethod("WriteSdrf", BindingFlags.NonPublic | BindingFlags.Instance)!.Invoke(task, null);
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }

            Assert.That(warnings.Any(w => w.Contains("No spectra files to describe")), Is.True, string.Join(" | ", warnings));
            Assert.That(File.Exists(Path.Combine(output, SdrfFileName)), Is.False);
            Assert.That(File.Exists(Path.Combine(output, "SdrfWriter_crash.txt")), Is.False, "an empty search is not a crash");

            Directory.Delete(output, true);
        }

        #endregion

        #region Warnings before the run

        /// <summary>
        /// A labelled search is told before it starts that comment[label] will not be filled in, so
        /// the gap is not first discovered in the written file.
        ///
        /// A turnover search that sets only StartTurnoverLabel / EndTurnoverLabel is warned too: those are
        /// copied into SilacLabels only after this check runs.
        /// </summary>
        [TestCase(false)]
        [TestCase(true)]
        public static void ALabelledSearchIsWarnedThatItsLabelWillNotBeFilledIn(bool turnoverOnly)
        {
            string folder = SetUpIsolatedRun(nameof(ALabelledSearchIsWarnedThatItsLabelWillNotBeFilledIn) + turnoverOnly,
                out string spectraPath, out _);

            var label = new SilacLabel('K', 'a', "C{6}H{12}N{2}O{1}", 6.020129);
            var task = new SearchTask
            {
                SearchParameters = turnoverOnly
                    ? new SearchParameters { WriteSdrf = true, EndTurnoverLabel = label }
                    : new SearchParameters { WriteSdrf = true, SilacLabels = new List<SilacLabel> { label } }
            };

            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                typeof(SearchTask)
                    .GetMethod("WarnAboutSdrfGaps", BindingFlags.NonPublic | BindingFlags.Instance)!
                    .Invoke(task, new object[] { new List<string> { spectraPath } });
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }

            Assert.That(warnings.Any(w => w.Contains("comment[label]")), Is.True, string.Join(" | ", warnings));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// The pre-run SDRF check only warns, so it must never be what stops a search. An
        /// ExperimentalDesign.tsv held open by another program (Excel locks what it opens) cannot be
        /// read; the search should start anyway, told why the design could not be checked.
        /// </summary>
        [Test]
        public static void AnUnreadableDesignIsAWarningNotAFailedSearch()
        {
            string folder = SetUpIsolatedRun(nameof(AnUnreadableDesignIsAWarningNotAFailedSearch),
                out string spectraPath, out _);
            string designPath = Path.Combine(Path.GetDirectoryName(spectraPath)!,
                GlobalVariables.ExperimentalDesignFileName);
            File.WriteAllText(designPath, "FileName\tCondition\tBiorep\tFraction\tTechrep\n");

            var task = new SearchTask { SearchParameters = new SearchParameters { WriteSdrf = true } };

            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                using (new FileStream(designPath, FileMode.Open, FileAccess.ReadWrite, FileShare.None))
                {
                    Assert.DoesNotThrow(() => typeof(SearchTask)
                        .GetMethod("WarnAboutSdrfGaps", BindingFlags.NonPublic | BindingFlags.Instance)!
                        .Invoke(task, new object[] { new List<string> { spectraPath } }));
                }
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }

            Assert.That(warnings.Any(w => w.Contains(GlobalVariables.ExperimentalDesignFileName)
                                          && w.Contains("could not be read")),
                Is.True, string.Join(" | ", warnings));

            Directory.Delete(folder, true);
        }

        #endregion

        #region Isobaric searches

        private static readonly string[] Tmt11Channels =
            { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" };

        /// <summary>
        /// A channel is a PRIDE term, looked up rather than spelled. The accessions are the check:
        /// PRIDE numbers 127C before 127N, so a resolver that walked the vocabulary in accession
        /// order instead of by name would transpose them and still produce plausible-looking cells.
        /// </summary>
        [TestCase(IsobaricMassTagType.TMT11, "126", "PRIDE:0000516")]
        [TestCase(IsobaricMassTagType.TMT11, "127N", "PRIDE:0000519")]
        [TestCase(IsobaricMassTagType.TMT11, "127C", "PRIDE:0000518")]
        [TestCase(IsobaricMassTagType.TMT6, "127", "PRIDE:0000517")]
        [TestCase(IsobaricMassTagType.TMT18, "135N", "PRIDE:0000670")]
        [TestCase(IsobaricMassTagType.iTRAQ8, "121", "PRIDE:0000538")]
        // 131 depends on the kit: TMT131 (PRIDE:0000529) for TMT6/TMT10, TMT131N (PRIDE:0000580, a TMTpro
        // term) for TMT11 and TMTpro, as sdrf-pipelines' channel_map.yaml spells them.
        [TestCase(IsobaricMassTagType.TMT10, " 131n ", "PRIDE:0000529")]
        [TestCase(IsobaricMassTagType.TMT6, "131", "PRIDE:0000529")]
        [TestCase(IsobaricMassTagType.TMT11, "131N", "PRIDE:0000580")]
        [TestCase(IsobaricMassTagType.TMT11, "131C", "PRIDE:0000581")]
        [TestCase(IsobaricMassTagType.TMT18, "131N", "PRIDE:0000580")]
        public static void AnIsobaricChannelResolvesToItsPrideTerm(IsobaricMassTagType tagType, string channel, string accession)
        {
            Assert.That(InvokeChannelLabel(tagType, channel)?.Accession, Is.EqualTo(accession));
        }

        /// <summary>
        /// A TMT10 file's label set is TMT10PLEX's, ending in TMT131. With TMT131N, a TMTpro term, the
        /// smallest plex holding every label is tmt11plex, and tools that infer the plex from the labels
        /// (sdrf-pipelines' _infer_plex) configure the wrong one (pcruzparri, Alexander-Sol, #2817 review).
        /// </summary>
        [Test]
        public static void ATmt10FileWritesTheTmt10LabelSet()
        {
            var tmt10 = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.TMT10)!;
            var file = new TmtFileInfo(@"C:\data\tmt10.raw", "Plex1", 1, 1, tmt10
                .Select((tag, i) => new TmtPlexAnnotation { Tag = tag, SampleName = "S" + (i + 1), BiologicalReplicate = 1 })
                .ToList());

            var labels = InvokeChannelRows(file, IsobaricMassTagType.TMT10).Select(r => r.Sample.Label?.Name).ToList();

            Assert.That(labels, Is.EqualTo(new[]
            {
                "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C",
                "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131"
            }));
        }

        /// <summary>
        /// PRIDE defines no DiLeu channels, and a channel name that is not a real channel resolves to
        /// nothing rather than to a term someone assembled from a prefix.
        /// </summary>
        [TestCase(IsobaricMassTagType.diLeu12, "115a")]
        [TestCase(IsobaricMassTagType.TMT11, "999N")]
        [TestCase(IsobaricMassTagType.TMT11, "")]
        public static void AChannelPrideDoesNotDefineResolvesToNothing(IsobaricMassTagType tagType, string channel)
        {
            Assert.That(InvokeChannelLabel(tagType, channel), Is.Null);
        }

        /// <summary>
        /// The whole isobaric feature, end to end, on the TMT11 fixture MultiplexQuantificationTests
        /// uses. A TMT search with a TmtDesign.txt writes one row per sample per channel: the source
        /// name is the design's sample, the label is the channel's PRIDE term, the rows run in
        /// reporter m/z order, and EVERY annotated channel gets a row -- including the Empty one,
        /// whose absence would shorten the label set infer_tmtplex reads and make the round trip
        /// impossible by construction (QuantProject 002 §5).
        ///
        /// The document is then handed to mzLib's two independent readers of SDRF. The validator
        /// checks the row key (source name + assay name + label) is unique; the quantification auditor
        /// has to recognise a channel-level design from the file alone, which is what a later reader
        /// projecting it back into a design will depend on.
        /// </summary>
        [Test]
        public static void ATmtSearchWritesOneRowPerSamplePerChannel()
        {
            string root = RunTmtSearchWritingSdrf("SdrfOutput_TmtChannels", writeDesign: true,
                out string output, out _);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            var expectedChannels = Tmt11Channels.ToList();
            Assert.That(document.Results.Count, Is.EqualTo(expectedChannels.Count),
                "One row per annotated channel of the one file this run searched -- all eleven, so " +
                "the channel count in the file is the plex size.");

            for (int i = 0; i < expectedChannels.Count; i++)
            {
                SdrfRow row = document.Results[i];
                Assert.That(row["comment[label]"], Is.EqualTo("TMT" + expectedChannels[i]),
                    "Rows follow the plex's reporter m/z order, not the design file's row order, and " +
                    "the label is written bare (D25): it is what the community writes, and quantms " +
                    "crashes on the accessioned form.");
                Assert.That(row["source name"], Is.EqualTo("Sample" + (i + 1)),
                    "The source name is the sample the design put in this channel.");
                Assert.That(row["comment[data file]"], Is.EqualTo("VA084TQ_6.mzML"));
                Assert.That(row["assay name"], Is.EqualTo("run VA084TQ_6"),
                    "Every channel of one file is one assay.");
            }

            Assert.That(document.Results.Any(r => r["comment[label]"].Contains("TMT131C")), Is.True,
                "The Empty channel is a real channel of a real plex: it gets a row, so the set of " +
                "labels in this file is the whole plex and not ten elevenths of it.");
            Assert.That(document.Header.Contains("characteristics[sample type]"), Is.False,
                "Deferred, not unknown: the design states each channel's sample type (131C is empty " +
                "here), but characteristics[sample type] waits for QuantProject's M4 to put the concept " +
                "in mzLib (MAP-08). When it lands a column appears; no row or source name changes.");

            SdrfValidationResult validation = SdrfValidator.Validate(document);
            Assert.That(validation.Errors, Is.Empty, "mzLib's validator rejects the TMT SDRF: " + validation);

            SdrfQuantAudit audit = SdrfQuantAuditor.Audit(document);
            Assert.That(audit.Kind, Is.EqualTo(SdrfQuantKind.ChannelLevel), audit.ToReport());
            Assert.That(audit.Channels, Has.Count.EqualTo(expectedChannels.Count), audit.ToReport());
            Assert.That(audit.DuplicateFileLabelPairs, Is.Empty, audit.ToReport());

            Directory.Delete(root, true);
        }

        /// <summary>
        /// An empty channel whose design row names no sample still gets a row, with a source name
        /// that is present and distinct.
        ///
        /// This is the case that writing every channel creates. A design is allowed to leave an empty
        /// channel's sample name blank, and `source name` is the one column REQ-2 keys incoming sample
        /// blocks on (D27) -- so a blank there would either collide two empty channels into one sample
        /// or be written as a reserved word, and both are worse than naming the channel after the file
        /// and tag it actually is.
        ///
        /// What it deliberately does NOT do is spell "empty" in the source name. That fact belongs in
        /// characteristics[sample type] when M4 lands, and putting it here would mean rewriting source
        /// names -- and breaking every join made against them -- on the day it does.
        /// </summary>
        [Test]
        public static void AnUnnamedEmptyChannelIsStillNamedInTheDocument()
        {
            string root = RunTmtSearchWritingSdrf("SdrfOutput_TmtUnnamedEmpty", writeDesign: true,
                out string output, out _, emptyChannelIsUnnamed: true);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            Assert.That(document.Results.Count, Is.EqualTo(Tmt11Channels.Length),
                "Every annotated channel is written whether or not its sample is named.");

            SdrfRow empty = document.Results.Single(r => r["comment[label]"].Contains("TMT131C"));
            Assert.That(empty["source name"], Is.EqualTo("Plex1 131C"),
                "The plex and the tag are facts already in the design, and together they are unique.");
            Assert.That(document.Results.Select(r => r["source name"]).Distinct().Count(),
                Is.EqualTo(Tmt11Channels.Length),
                "No two channels of one file may share a source name.");

            SdrfValidationResult validation = SdrfValidator.Validate(document);
            Assert.That(validation.Errors, Is.Empty, "mzLib's validator rejects the document: " + validation);

            Directory.Delete(root, true);
        }

        /// <summary>
        /// Without a TmtDesign.txt there is no channel-to-sample map, so the SDRF falls back to one row
        /// per file with the label unresolved -- and the user is told so BEFORE the search, naming
        /// the TMT design file rather than ExperimentalDesign.tsv, which an isobaric search never reads.
        /// </summary>
        [Test]
        public static void ATmtSearchWithoutItsDesignDescribesEachFileOnce_AndSaysWhy()
        {
            string root = RunTmtSearchWritingSdrf("SdrfOutput_TmtNoDesign", writeDesign: false,
                out string output, out List<string> warnings);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            Assert.That(document.Results.Count, Is.EqualTo(1), "One row for the one file, as for any search.");
            Assert.That(document.Results.Single()["comment[label]"], Is.EqualTo("not available"),
                "No design means no channel to name, and a TMT run is never called label free.");

            Assert.That(warnings.Any(w => w.Contains(GlobalVariables.TmtExperimentalDesignFileName)
                                          && w.Contains("without its channels or samples")), Is.True,
                "The fallback is announced, naming the file that would fix it: " + string.Join(" | ", warnings));
            // Matched on the file name, not on "no ...": the post-run warning starts "No ExperimentalDesign.tsv".
            Assert.That(warnings.Where(w => w.Contains(GlobalVariables.ExperimentalDesignFileName)), Is.Empty,
                "An isobaric search does not read ExperimentalDesign.tsv, before or after the run, so its " +
                "absence is not reported.");

            Directory.Delete(root, true);
        }

        /// <summary>
        /// A TmtDesign.txt the reader refuses is no better than none: the SDRF falls back to one row
        /// per file, and the warning before the search says the design is the reason, not its absence.
        /// </summary>
        [Test]
        public static void ATmtSearchWithAnUnusableDesignDescribesEachFileOnce_AndSaysWhy()
        {
            string root = RunTmtSearchWritingSdrf("SdrfOutput_TmtUnusableDesign", writeDesign: true,
                out string output, out List<string> warnings, designIsUnusable: true);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            Assert.That(document.Results.Count, Is.EqualTo(1), "A refused design gives no channels to expand.");
            Assert.That(document.Results.Single()["comment[label]"], Is.EqualTo("not available"));

            Assert.That(warnings.Any(w => w.Contains("cannot be used as it stands")
                                          && w.Contains("not a sample type")), Is.True,
                "The warning names the design's own error: " + string.Join(" | ", warnings));
            Assert.That(warnings.Where(w => w.Contains(GlobalVariables.ExperimentalDesignFileName)), Is.Empty,
                "The fallback does not turn to ExperimentalDesign.tsv either.");

            Directory.Delete(root, true);
        }

        /// <summary>
        /// PRIDE has no DiLeu channel terms, so a DiLeu channel row could not say which channel it is:
        /// N rows on one file with no label is what mzLib's auditor reads as a label-free file holding
        /// N samples. A DiLeu file is therefore described once, like a TMT file with no design, and the
        /// warning before the search says so.
        /// </summary>
        [Test]
        public static void ADiLeuSearchIsWarnedThatEachFileIsDescribedOnce()
        {
            var warnings = WarningsBeforeAnIsobaricSearch(nameof(ADiLeuSearchIsWarnedThatEachFileIsDescribedOnce),
                "DiLeu-12plex on K", spectra => new[] { $"{spectra}	Plex1	Sample1	115a	CondA	1	1	1	study sample" });

            Assert.That(warnings, Has.Count.EqualTo(1), string.Join(" | ", warnings));
            Assert.That(warnings.Single(), Does.Contain("diLeu12").And.Contain("comment[label]")
                .And.Contain("without their channels or samples"));
        }

        /// <summary>
        /// The gate that decides channel rows, for the three cases that must fall back to one row per
        /// file: a tag type PRIDE has no channel terms for, a plex nobody annotated (the placeholder row
        /// TmtExperimentalDesign.Write emits, which Read accepts without error), and a channel the
        /// search's plex does not have. Without the gate, the first gives unlabelled channel rows, the
        /// second drops the file from the SDRF, and the third writes a TMT6 label on a TMT11 search.
        /// </summary>
        [TestCase(IsobaricMassTagType.diLeu12, "115a", "not diLeu12")]
        [TestCase(IsobaricMassTagType.TMT11, null, "no annotated channels")]
        [TestCase(IsobaricMassTagType.TMT11, "127", "not a channel of TMT11")]
        public static void ChannelRowsAreRefusedWhenTheyCouldNotSayWhichChannelIsWhich(
            IsobaricMassTagType tagType, string tag, string reason)
        {
            var annotations = tag is null
                ? new List<TmtPlexAnnotation>()
                : new List<TmtPlexAnnotation> { new() { Tag = tag, SampleName = "S1", BiologicalReplicate = 1 } };
            var file = new TmtFileInfo(@"C:\data\run1.raw", "Plex1", 1, 1, annotations);

            Assert.That(InvokeChannelRowsUnusableReason(file, tagType), Does.Contain(reason));
        }

        [Test]
        public static void ChannelRowsAreAllowedForAnAnnotatedPlexOfTheSearchedTag()
        {
            var file = new TmtFileInfo(@"C:\data\run1.raw", "Plex1", 1, 1, new List<TmtPlexAnnotation>
            {
                new() { Tag = "126", SampleName = "S1", BiologicalReplicate = 1 },
                new() { Tag = "127N", SampleName = "S2", BiologicalReplicate = 1 }
            });

            Assert.That(InvokeChannelRowsUnusableReason(file, IsobaricMassTagType.TMT11), Is.Null);
        }

        /// <summary>
        /// A plex with no annotated channels, a channel off the search's plex, and a multiplex label
        /// MetaMorpheus does not recognise all mean one row per file. Each is named before the search,
        /// when it is still cheap to fix.
        /// </summary>
        [Test]
        public static void AnUnannotatedPlexIsWarnedBeforeTheSearch()
        {
            var warnings = WarningsBeforeAnIsobaricSearch(nameof(AnUnannotatedPlexIsWarnedBeforeTheSearch),
                "TMT11", spectra => new[] { spectra + "\tPlex1" + new string('\t', 7) });

            Assert.That(warnings, Has.Count.EqualTo(1), string.Join(" | ", warnings));
            Assert.That(warnings.Single(), Does.Contain("no annotated channels")
                .And.Contain("without their channels or samples"));
        }

        [Test]
        public static void AChannelOffTheSearchedPlexIsWarnedBeforeTheSearch()
        {
            var warnings = WarningsBeforeAnIsobaricSearch(nameof(AChannelOffTheSearchedPlexIsWarnedBeforeTheSearch),
                "TMT11", spectra => new[] { $"{spectra}	Plex1	Sample1	127	CondA	1	1	1	study sample" });

            Assert.That(warnings, Has.Count.EqualTo(1), string.Join(" | ", warnings));
            Assert.That(warnings.Single(), Does.Contain("'127'").And.Contain("not a channel of TMT11"));
        }

        [Test]
        public static void AnUnrecognisedMultiplexLabelIsWarnedBeforeTheSearch()
        {
            var warnings = WarningsBeforeAnIsobaricSearch(nameof(AnUnrecognisedMultiplexLabelIsWarnedBeforeTheSearch),
                "Nonsense on K", spectra => new[] { $"{spectra}	Plex1	Sample1	126	CondA	1	1	1	study sample" });

            Assert.That(warnings, Has.Count.EqualTo(1), string.Join(" | ", warnings));
            Assert.That(warnings.Single(), Does.Contain("Nonsense on K")
                .And.Contain("without its channels or samples"));
        }

        /// <summary>
        /// After the search, the same unrecognised multiplex label leaves no isobaric design to write
        /// channel rows from, even with a usable TmtDesign.txt beside the spectra: with no tag type
        /// there is no plex to check the channels against, so the file is described once.
        /// </summary>
        [Test]
        public static void AnUnrecognisedMultiplexLabelReadsNoIsobaricDesign()
        {
            string folder = SetUpIsolatedRun(nameof(AnUnrecognisedMultiplexLabelReadsNoIsobaricDesign), out string spectraPath, out _);
            File.WriteAllLines(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName),
                new[] { TmtExperimentalDesign.Header, $"{spectraPath}	Plex1	Sample1	126	CondA	1	1	1	study sample" });
            var task = new PostSearchAnalysisTask
            {
                Parameters = new PostSearchAnalysisParameters
                {
                    CurrentRawFileList = new List<string> { spectraPath },
                    SearchParameters = new SearchParameters { DoMultiplexQuantification = true, MultiplexModId = "Nonsense on K" }
                }
            };

            var args = new object[] { null };
            var design = typeof(PostSearchAnalysisTask)
                .GetMethod("ReadIsobaricDesignIfPresent", BindingFlags.NonPublic | BindingFlags.Instance)!
                .Invoke(task, args);

            Assert.That(design, Is.Null);
            Assert.That(args[0], Is.Null, "no tag type");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// An unnamed channel belongs to the plex, and every fraction of a plex carries the same
        /// channel-to-sample map. So its source name is built from the plex, not the file: with the
        /// file stem, one unnamed channel of a two-fraction plex became two samples downstream.
        /// </summary>
        [Test]
        public static void AnUnnamedChannelHasOneSourceNameAcrossTheFractionsOfItsPlex()
        {
            var annotations = new List<TmtPlexAnnotation>
            {
                new() { Tag = "126", SampleName = "S1", BiologicalReplicate = 1 },
                new() { Tag = "131C", SampleName = "", BiologicalReplicate = 1, SampleType = TmtSampleType.Empty }
            };
            var fraction1 = new TmtFileInfo(@"C:\data\plexA_F01.raw", "PlexA", 1, 1, annotations);
            var fraction2 = new TmtFileInfo(@"C:\data\plexA_F02.raw", "PlexA", 2, 1, annotations);

            string UnnamedSourceName(TmtFileInfo file) => InvokeChannelRows(file, IsobaricMassTagType.TMT11)
                .Single(r => r.Sample.Label?.Name == "TMT131C").Sample.SourceName;

            Assert.That(UnnamedSourceName(fraction1), Is.EqualTo("PlexA 131C"));
            Assert.That(UnnamedSourceName(fraction2), Is.EqualTo(UnnamedSourceName(fraction1)),
                "Both fractions measure the same channel of the same plex, so it is one sample.");

            var noPlex = new TmtFileInfo(@"C:\data\plexA_F01.raw", "", 1, 1, annotations);
            Assert.That(UnnamedSourceName(noPlex), Is.EqualTo("plexA_F01 131C"),
                "With no plex name the file stem is the only scope left, and it is still unique.");
        }

        /// <summary>
        /// TmtExperimentalDesign.Read checks sample names only within a plex, so a two-plex design can
        /// name each plex's channels S1, S2, ... for different material (pcruzparri, #2817 review). Source
        /// name is the sample key, so such a name is scoped to its plex; a bridge, named the same in both
        /// plexes with the same condition and replicate, is one sample and keeps its one name.
        /// </summary>
        [Test]
        public static void ASampleNameReusedInAnotherPlexForADifferentSampleIsScopedToItsPlex()
        {
            var plex1 = new TmtFileInfo(@"C:\data\p1.raw", "Plex1", 1, 1, new List<TmtPlexAnnotation>
            {
                new() { Tag = "126", SampleName = "S1", Condition = "Ctrl", BiologicalReplicate = 1 },
                new() { Tag = "127N", SampleName = "Bridge", Condition = "Pool", BiologicalReplicate = 1 }
            });
            var plex2 = new TmtFileInfo(@"C:\data\p2.raw", "Plex2", 1, 1, new List<TmtPlexAnnotation>
            {
                new() { Tag = "126", SampleName = "S1", Condition = "Drug", BiologicalReplicate = 2 },
                new() { Tag = "127N", SampleName = "Bridge", Condition = "Pool", BiologicalReplicate = 1 }
            });

            var reused = InvokeSampleNamesReusedForDifferentSamples(plex1, plex2);
            Assert.That(reused, Is.EquivalentTo(new[] { "S1" }), "Only the name that names two samples.");

            string SourceName(TmtFileInfo file, string label) => InvokeChannelRows(file, IsobaricMassTagType.TMT11, plex1, plex2)
                .Single(r => r.Sample.Label?.Name == label).Sample.SourceName;

            Assert.That(SourceName(plex1, "TMT126"), Is.EqualTo("Plex1 S1"));
            Assert.That(SourceName(plex2, "TMT126"), Is.EqualTo("Plex2 S1"),
                "Two samples, two source names.");
            Assert.That(SourceName(plex1, "TMT127N"), Is.EqualTo("Bridge"));
            Assert.That(SourceName(plex2, "TMT127N"), Is.EqualTo("Bridge"),
                "A bridge is one sample measured in both plexes, so it keeps one source name.");

            Assert.That(InvokeSampleNamesReusedForDifferentSamples(plex1), Is.Empty,
                "Plex1 names each channel differently, so no name in it names two samples.");
        }

        /// <summary>
        /// Read checks uniqueness within a plex on (sample, biological replicate, fraction, technical
        /// replicate), not on the name, so one plex can name every channel after its cell line with a
        /// different biological replicate on each (what Auto-fill writes). Adding the plex alone left all
        /// those channels one source name, so they take the channel as well, and the warning says the name
        /// repeats inside a plex rather than across plexes (pcruzparri, #2817 review).
        /// </summary>
        [Test]
        public static void ASampleNameRepeatedOnChannelsOfOnePlexGetsOneSourceNamePerChannel()
        {
            var annotations = new List<TmtPlexAnnotation>
            {
                new() { Tag = "126", SampleName = "HeLa", Condition = "Ctrl", BiologicalReplicate = 1 },
                new() { Tag = "127N", SampleName = "HeLa", Condition = "Ctrl", BiologicalReplicate = 2 },
                new() { Tag = "127C", SampleName = "HeLa", Condition = "Drug", BiologicalReplicate = 3 },
                new() { Tag = "128N", SampleName = "Pool", Condition = "Pool", BiologicalReplicate = 1 }
            };
            var fraction1 = new TmtFileInfo(@"C:\data\p1_F01.raw", "Plex1", 1, 1, annotations);
            var fraction2 = new TmtFileInfo(@"C:\data\p1_F02.raw", "Plex1", 2, 1, annotations);

            var reused = InvokeSampleNamesReusedForDifferentSamples(fraction1, fraction2);
            Assert.That(reused, Is.EquivalentTo(new[] { "HeLa" }));

            List<string> SourceNames(TmtFileInfo file) => InvokeChannelRows(file, IsobaricMassTagType.TMT11, fraction1, fraction2)
                .Select(r => r.Sample.SourceName).ToList();

            Assert.That(SourceNames(fraction1),
                Is.EqualTo(new[] { "Plex1 HeLa 126", "Plex1 HeLa 127N", "Plex1 HeLa 127C", "Pool" }),
                "Three samples, three source names; Pool names one sample and keeps its name.");
            Assert.That(SourceNames(fraction2), Is.EqualTo(SourceNames(fraction1)),
                "Every fraction of the plex measures the same channels, so the names are the same.");

            string warning = (string)typeof(PostSearchAnalysisTask)
                .GetMethod("ReusedSampleNamesWarning", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { reused, InvokeSampleNamesRepeatedWithinAPlex(fraction1, fraction2) });
            Assert.That(warning, Does.Contain("more than one channel of the same plex")
                .And.Contain("'HeLa'").And.Contain("<plex> <sample name> <channel>")
                .And.Not.Contain("in more than one plex"));
        }

        /// <summary>The within-plex repeat is named before the search too, in its own words.</summary>
        [Test]
        public static void ASampleNameRepeatedOnChannelsOfOnePlexIsWarnedBeforeTheSearch()
        {
            var warnings = WarningsBeforeAnIsobaricSearch(nameof(ASampleNameRepeatedOnChannelsOfOnePlexIsWarnedBeforeTheSearch),
                "TMT11", spectra => new[]
                {
                    $"{spectra}	Plex1	HeLa	126	Ctrl	1	1	1	study sample",
                    $"{spectra}	Plex1	HeLa	127N	Drug	2	1	1	study sample"
                });

            Assert.That(warnings, Has.Count.EqualTo(1), string.Join(" | ", warnings));
            Assert.That(warnings.Single(), Does.Contain("'HeLa'").And.Contain("more than one channel of the same plex")
                .And.Not.Contain("in more than one plex"));
        }

        /// <summary>The same reuse is named before the search, when it is still cheap to rename.</summary>
        [Test]
        public static void ASampleNameReusedInAnotherPlexIsWarnedBeforeTheSearch()
        {
            var warnings = WarningsBeforeAnIsobaricSearch(nameof(ASampleNameReusedInAnotherPlexIsWarnedBeforeTheSearch),
                "TMT11", spectra =>
                {
                    string second = Path.Combine(Path.GetDirectoryName(spectra)!, "plex2.mzML");
                    return new[]
                    {
                        $"{spectra}	Plex1	S1	126	Ctrl	1	1	1	study sample",
                        $"{spectra}	Plex1	Bridge	127N	Pool	1	1	1	reference",
                        $"{second}	Plex2	S1	126	Drug	2	1	1	study sample",
                        $"{second}	Plex2	Bridge	127N	Pool	1	1	1	reference"
                    };
                }, "plex2.mzML");

            Assert.That(warnings, Has.Count.EqualTo(1), string.Join(" | ", warnings));
            Assert.That(warnings.Single(), Does.Contain("'S1'").And.Not.Contain("'Bridge'")
                .And.Contain("<plex> <sample name>"));
        }

        /// <summary>
        /// mzLib's SdrfSampleBlock.BySourceName groups rows by source name ignoring case, so "Liver" in one plex
        /// and "liver" in another are one sample downstream. For different samples they are therefore a reused
        /// name, and each is scoped to its plex.
        /// </summary>
        [Test]
        public static void SampleNamesThatDifferOnlyInCaseAreOneName()
        {
            var rows = IsobaricSdrfRows(nameof(SampleNamesThatDifferOnlyInCaseAreOneName),
                new[] { "p1.mzML", "p2.mzML" }, folder => new[]
                {
                    TmtDesignRow(folder, "p1.mzML", "Plex1", "Liver", "126", "Ctrl", 1, 1),
                    TmtDesignRow(folder, "p2.mzML", "Plex2", "liver", "126", "Drug", 2, 1)
                }, out _);

            Assert.That(rows.Select(r => r.Sample.SourceName), Is.EqualTo(new[] { "Plex1 Liver", "Plex2 liver" }));
            AssertEachSourceNameIsOneSample(rows);
        }

        /// <summary>
        /// A source name the SDRF makes up is checked against every name in the document (Alexander-Sol,
        /// #2817 review). The unnamed channel 126 of Plex1 is called "Plex1 126", and nothing stops the design
        /// giving another sample exactly that name. Both would then be one source name, and anything keyed on
        /// it (REQ-2, mzLib's SdrfSampleBlock.BySourceName) would merge a Ctrl sample with a Drug one. The name
        /// the user wrote is kept; the made-up one is numbered, the same way in every fraction of the plex.
        /// </summary>
        [Test]
        public static void AMadeUpSourceNameIsNeverANameTheDesignGivesAnotherSample()
        {
            var rows = IsobaricSdrfRows(nameof(AMadeUpSourceNameIsNeverANameTheDesignGivesAnotherSample),
                new[] { "plex1_F1.mzML", "plex1_F2.mzML" }, folder => new[] { 1, 2 }.SelectMany(fraction => new[]
                {
                    TmtDesignRow(folder, $"plex1_F{fraction}.mzML", "Plex1", "", "126", "Ctrl", 1, fraction),
                    TmtDesignRow(folder, $"plex1_F{fraction}.mzML", "Plex1", "Plex1 126", "127N", "Drug", 2, fraction)
                }), out var warnings);

            string SourceName(int fraction, string label) => rows
                .Single(r => r.Assay.DataFileName == $"plex1_F{fraction}.mzML" && r.Sample.Label?.Name == label)
                .Sample.SourceName;

            Assert.That(SourceName(1, "TMT127N"), Is.EqualTo("Plex1 126"), "The name the user wrote is written as given.");
            Assert.That(SourceName(1, "TMT126"), Is.EqualTo("Plex1 126 (2)"),
                "The unnamed channel is another sample, so it takes the next free name.");
            Assert.That(SourceName(2, "TMT126"), Is.EqualTo(SourceName(1, "TMT126")),
                "Both fractions measure the same channel of the same plex, so it is one sample.");
            Assert.That(SourceName(2, "TMT127N"), Is.EqualTo(SourceName(1, "TMT127N")));
            AssertEachSourceNameIsOneSample(rows);

            Assert.That(warnings.Any(w => w.Contains("'Plex1 126' is written 'Plex1 126 (2)'")), Is.True,
                "The numbering is announced: " + string.Join(" | ", warnings));
        }

        /// <summary>
        /// Every kind of made-up source name is checked, not only an unnamed channel's: a name scoped to its
        /// plex because another plex uses it for a different sample ("Plex2 S1"), a name scoped to its channel
        /// because one plex gives it to several samples ("Plex1 HeLa 128N"), and the stem a file described
        /// once is named after ("Ref"). Names are compared ignoring case, as mzLib's BySourceName compares them.
        /// </summary>
        [Test]
        public static void EveryKindOfMadeUpSourceNameIsCheckedAgainstTheDesign()
        {
            var rows = IsobaricSdrfRows(nameof(EveryKindOfMadeUpSourceNameIsCheckedAgainstTheDesign),
                new[] { "p1.mzML", "p2.mzML", "Ref.mzML" }, folder => new[]
                {
                    TmtDesignRow(folder, "p1.mzML", "Plex1", "S1", "126", "Ctrl", 1, 1),
                    TmtDesignRow(folder, "p1.mzML", "Plex1", "plex2 s1", "127N", "Ctrl", 2, 1),
                    TmtDesignRow(folder, "p1.mzML", "Plex1", "HeLa", "128N", "Ctrl", 3, 1),
                    TmtDesignRow(folder, "p1.mzML", "Plex1", "HeLa", "128C", "Drug", 4, 1),
                    TmtDesignRow(folder, "p2.mzML", "Plex2", "S1", "126", "Drug", 5, 1),
                    TmtDesignRow(folder, "p2.mzML", "Plex2", "Plex1 HeLa 128N", "127N", "Drug", 6, 1),
                    TmtDesignRow(folder, "p2.mzML", "Plex2", "Ref", "127C", "Pool", 1, 1),
                    // A plex nobody annotated, so its file is described once, under its stem.
                    Path.Combine(folder, "Ref.mzML") + "\tPlex3" + new string('\t', 7)
                }, out _);

            string SourceName(string file, string label) => rows
                .Single(r => r.Assay.DataFileName == file && r.Sample.Label?.Name == label)
                .Sample.SourceName;

            Assert.That(SourceName("p1.mzML", "TMT126"), Is.EqualTo("Plex1 S1"));
            Assert.That(SourceName("p1.mzML", "TMT127N"), Is.EqualTo("plex2 s1"));
            Assert.That(SourceName("p2.mzML", "TMT126"), Is.EqualTo("Plex2 S1 (2)"), "Not the user's 'plex2 s1'.");
            Assert.That(SourceName("p2.mzML", "TMT127N"), Is.EqualTo("Plex1 HeLa 128N"));
            Assert.That(SourceName("p1.mzML", "TMT128N"), Is.EqualTo("Plex1 HeLa 128N (2)"));
            Assert.That(SourceName("p1.mzML", "TMT128C"), Is.EqualTo("Plex1 HeLa 128C"));
            Assert.That(SourceName("p2.mzML", "TMT127C"), Is.EqualTo("Ref"));
            Assert.That(rows.Single(r => r.Assay.DataFileName == "Ref.mzML").Sample.SourceName, Is.EqualTo("Ref (2)"),
                "A file described once is not the sample the design calls Ref.");
            AssertEachSourceNameIsOneSample(rows);
        }

        /// <summary>The numbering is named before the search too, when the user can still choose the names.</summary>
        [Test]
        public static void AMadeUpSourceNameTheDesignAlreadyUsesIsWarnedBeforeTheSearch()
        {
            var warnings = WarningsBeforeAnIsobaricSearch(nameof(AMadeUpSourceNameTheDesignAlreadyUsesIsWarnedBeforeTheSearch),
                "TMT11", spectra => new[]
                {
                    $"{spectra}\tPlex1\t\t126\tCtrl\t1\t1\t1\tstudy sample",
                    $"{spectra}\tPlex1\tPlex1 126\t127N\tDrug\t2\t1\t1\tstudy sample"
                });

            Assert.That(warnings, Has.Count.EqualTo(1), string.Join(" | ", warnings));
            Assert.That(warnings.Single(), Does.Contain("'Plex1 126' is written 'Plex1 126 (2)'")
                .And.Contain(GlobalVariables.TmtExperimentalDesignFileName));
        }

        /// <summary>
        /// The rows BuildSdrfRows gives for a TMT11 search of <paramref name="fileNames"/>, with a TmtDesign.txt
        /// of <paramref name="designRows"/> (given the folder) beside them, and the warnings it raised. No
        /// search is run and the spectra files need not exist: only the design is read. With
        /// <paramref name="inputSdrf"/>, its rows are the search's input SDRF. <paramref name="multiplexModId"/> is the kit.
        /// </summary>
        private static List<SdrfRowInput> IsobaricSdrfRows(string testName, string[] fileNames,
            Func<string, IEnumerable<string>> designRows, out List<string> warnings, IEnumerable<string[]> inputSdrf = null,
            string multiplexModId = "TMT11")
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SdrfOutput_" + testName);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            File.WriteAllLines(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName),
                new[] { TmtExperimentalDesign.Header }.Concat(designRows(folder)));

            var task = new PostSearchAnalysisTask
            {
                CommonParameters = new CommonParameters(),
                Parameters = new PostSearchAnalysisParameters
                {
                    CurrentRawFileList = fileNames.Select(n => Path.Combine(folder, n)).ToList(),
                    SearchParameters = new SearchParameters { DoMultiplexQuantification = true, MultiplexModId = multiplexModId }
                }
            };

            // BuildSdrfRows's one argument is the input SDRF, scoped to the searched files as WriteSdrf scopes it.
            SdrfSearchScoping scopedInput = null;
            if (inputSdrf is not null)
            {
                string inputPath = Path.Combine(folder, "input.sdrf.tsv");
                File.WriteAllLines(inputPath, inputSdrf.Select(cells => string.Join("\t", cells)));
                task.Parameters.SearchParameters.InputSdrfPath = inputPath;
                scopedInput = SdrfSearchScope.Restrict(new SdrfDocument(inputPath), task.Parameters.CurrentRawFileList);
            }

            var captured = new List<string>();
            warnings = captured;
            EventHandler<StringEventArgs> handler = (o, e) => captured.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                return ((IEnumerable<SdrfRowInput>)typeof(PostSearchAnalysisTask)
                        .GetMethod("BuildSdrfRows", BindingFlags.NonPublic | BindingFlags.Instance)!
                        .Invoke(task, new object[] { scopedInput }))
                    .ToList();
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
                Directory.Delete(folder, true);
            }
        }

        /// <summary>One TmtDesign.txt row, technical replicate 1, for a study sample.</summary>
        private static string TmtDesignRow(string folder, string file, string plex, string sampleName, string channel,
            string condition, int biologicalReplicate, int fraction) =>
            string.Join("\t", Path.Combine(folder, file), plex, sampleName, channel, condition, biologicalReplicate,
                fraction, 1, "study sample");

        /// <summary>
        /// No source name, compared ignoring case as mzLib's SdrfSampleBlock.BySourceName compares it, is
        /// given to two samples, where a sample is a (condition, biological replicate).
        /// </summary>
        private static void AssertEachSourceNameIsOneSample(IEnumerable<SdrfRowInput> rows)
        {
            var merged = rows
                .GroupBy(r => r.Sample.SourceName, StringComparer.OrdinalIgnoreCase)
                .Where(g => g.Select(r => (r.Sample.FactorValue, r.Sample.BiologicalReplicate)).Distinct().Count() > 1)
                .Select(g => g.Key)
                .ToList();
            Assert.That(merged, Is.Empty, "These source names each name more than one sample.");
        }

        /// <summary>
        /// End to end: a TmtDesign.txt whose only row for the file is the channel-less placeholder
        /// (the GUI's normal state while a design is half-filled). The file used to vanish from the
        /// SDRF, and with every file in that state no SDRF was written at all.
        /// </summary>
        [Test]
        public static void ATmtSearchWhosePlexHasNoAnnotatedChannelsStillDescribesTheFile()
        {
            string root = RunTmtSearchWritingSdrf("SdrfOutput_TmtPlaceholderOnly", writeDesign: true,
                out string output, out List<string> warnings, placeholderOnly: true);

            string path = Path.Combine(output, SdrfFileName);
            Assert.That(File.Exists(path), Is.True, "An SDRF is written: " + string.Join(" | ", warnings));
            var document = new SdrfDocument(path);
            document.LoadResults();

            Assert.That(document.Results.Count, Is.EqualTo(1), "The file is described once.");
            Assert.That(document.Results.Single()["comment[data file]"], Is.EqualTo("VA084TQ_6.mzML"));
            Assert.That(document.Results.Single()["comment[label]"], Is.EqualTo("not available"));
            Assert.That(warnings.Any(w => w.Contains("no annotated channels")), Is.True, string.Join(" | ", warnings));

            Directory.Delete(root, true);
        }

        /// <summary>
        /// The warnings SearchTask gives before an isobaric search with SDRF output on, for a
        /// TmtDesign.txt of <paramref name="designRows"/> (given the spectra path) beside the spectra.
        /// </summary>
        private static List<string> WarningsBeforeAnIsobaricSearch(string testName, string multiplexModId,
            Func<string, string[]> designRows, params string[] otherSpectraFileNames) =>
            WarningsBeforeAnIsobaricSearch(testName, multiplexModId, designRows, null, otherSpectraFileNames);

        /// <summary>
        /// As above, and with <paramref name="inputSdrf"/>, its rows are the search's input SDRF. The spectra file is
        /// PrunedDbSpectra.mzml, so an input row for it names PrunedDbSpectra.raw, say.
        /// </summary>
        private static List<string> WarningsBeforeAnIsobaricSearch(string testName, string multiplexModId,
            Func<string, string[]> designRows, IEnumerable<string[]> inputSdrf, params string[] otherSpectraFileNames)
        {
            string folder = SetUpIsolatedRun(testName, out string spectraPath, out _);
            var spectraPaths = new List<string> { spectraPath }
                .Concat(otherSpectraFileNames.Select(n => Path.Combine(folder, n))).ToList();
            File.WriteAllLines(Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName),
                new[] { TmtExperimentalDesign.Header }.Concat(designRows(spectraPath)));

            var task = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WriteSdrf = true,
                    DoMultiplexQuantification = true,
                    MultiplexModId = multiplexModId
                }
            };
            if (inputSdrf is not null)
            {
                string inputPath = Path.Combine(folder, "input.sdrf.tsv");
                File.WriteAllLines(inputPath, inputSdrf.Select(cells => string.Join("\t", cells)));
                task.SearchParameters.InputSdrfPath = inputPath;
            }

            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                typeof(SearchTask)
                    .GetMethod("WarnAboutSdrfGaps", BindingFlags.NonPublic | BindingFlags.Instance)!
                    .Invoke(task, new object[] { spectraPaths });
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
                Directory.Delete(folder, true);
            }

            return warnings;
        }

        /// <summary>
        /// Runs the TMT11 fixture MultiplexQuantificationTests uses, with SDRF output on, in a folder of
        /// its own. With <paramref name="writeDesign"/>, a TmtDesign.txt goes beside the spectra: the
        /// last channel is Empty, and the rows are written in reverse so document order cannot pass
        /// for reporter m/z order. With <paramref name="designIsUnusable"/>, every row's Sample Type is
        /// one the design reader refuses. With <paramref name="inputSdrf"/>, its rows are the search's input SDRF.
        /// Returns the root folder to delete.
        /// </summary>
        private static string RunTmtSearchWritingSdrf(string folderName, bool writeDesign,
            out string output, out List<string> warnings, bool designIsUnusable = false,
            bool emptyChannelIsUnnamed = false, bool placeholderOnly = false, IEnumerable<string[]> inputSdrf = null)
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, folderName);
            if (Directory.Exists(root)) Directory.Delete(root, true);
            string dataFolder = Path.Combine(root, "data");
            Directory.CreateDirectory(dataFolder);

            string mzml = Path.Combine(dataFolder, "VA084TQ_6.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\VA084TQ_6.mzML"), mzml);

            if (writeDesign)
            {
                bool IsEmptyChannel(int i) => i == Tmt11Channels.Length - 1;
                IEnumerable<string> designRows = Tmt11Channels
                    .Select((tag, i) =>
                        $"{mzml}	Plex1	" +
                        // A design may leave an empty channel's sample name blank -- there is no
                        // sample to name -- and a real one often does.
                        (emptyChannelIsUnnamed && IsEmptyChannel(i) ? "" : $"Sample{i + 1}") +
                        $"	{tag}	Cond{(i % 2 == 0 ? "A" : "B")}	{i / 2 + 1}	1	1	" +
                        (designIsUnusable ? "not a sample type"
                            : IsEmptyChannel(i) ? "empty" : "study sample"))
                    .Reverse();
                if (placeholderOnly)
                    // The channel-less row TmtExperimentalDesign.Write emits for a plex nobody has annotated.
                    designRows = new[] { mzml + "\tPlex1" + new string('\t', 7) };
                File.WriteAllLines(Path.Combine(dataFolder, GlobalVariables.TmtExperimentalDesignFileName),
                    new[] { TmtExperimentalDesign.Header }.Concat(designRows));
            }

            var searchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\TMT-Task1-SearchTaskconfig.toml"),
                MetaMorpheusTask.tomlConfig);
            searchTask.SearchParameters.WriteSdrf = true;
            searchTask.SearchParameters.DoParsimony = true;
            if (inputSdrf is not null)
            {
                string inputPath = Path.Combine(dataFolder, "input.sdrf.tsv");
                File.WriteAllLines(inputPath, inputSdrf.Select(cells => string.Join("\t", cells)));
                searchTask.SearchParameters.InputSdrfPath = inputPath;
            }

            output = Path.Combine(root, "out");
            Directory.CreateDirectory(output);
            string fasta = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\mouseTmt.fasta");

            var captured = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => captured.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                searchTask.RunTask(output, new List<DbForTask> { new(fasta, false) }, new List<string> { mzml }, "tmt-sdrf");
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }

            warnings = captured;
            return root;
        }

        #endregion

        #region The dataset this search re-analyses

        /// <summary>
        /// A reanalysis names the dataset it searched. The accession is the join key for pooling, and
        /// no spectra file carries it, so it has to come from the task.
        /// </summary>
        [Test]
        public static void TheProteomeXchangeAccessionIsWrittenWhenSupplied()
        {
            string folder = SetUpIsolatedRun(nameof(TheProteomeXchangeAccessionIsWrittenWhenSupplied),
                out string spectraPath, out DbForTask database);

            var task = BuildSearchTask(writeSdrf: true);
            task.SearchParameters.ProteomeXchangeAccession = " PXD012345 ";
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(output);

            task.RunTask(output, new List<DbForTask> { database }, new List<string> { spectraPath }, "sdrf-pxd");

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            Assert.That(document.Results.Single()["comment[proteomexchange accession number]"],
                Is.EqualTo("PXD012345"), "Written trimmed, into the spec's own column.");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Data that was never deposited has no accession, and the column is left out rather than
        /// filled with a reserved word that would read as "deposited, accession unknown".
        /// </summary>
        [Test]
        public static void NoAccessionColumnIsWrittenWhenNoneWasSupplied()
        {
            string output = RunSearchWritingSdrf(nameof(NoAccessionColumnIsWrittenWhenNoneWasSupplied),
                out string folder, out _);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            Assert.That(document.Header, Does.Not.Contain("comment[proteomexchange accession number]"));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// A pipeline sets the accession in the task file, so it has to survive a TOML round trip.
        /// </summary>
        [Test]
        public static void TheAccessionSurvivesTheTaskFile()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SdrfOutput_AccessionToml");
            Directory.CreateDirectory(folder);
            string tomlPath = Path.Combine(folder, "SearchTask.toml");

            var task = BuildSearchTask(writeSdrf: true);
            task.SearchParameters.ProteomeXchangeAccession = "PXD012345";
            Toml.WriteFile(task, tomlPath, MetaMorpheusTask.tomlConfig);

            var reread = Toml.ReadFile<SearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
            Assert.That(reread.SearchParameters.ProteomeXchangeAccession, Is.EqualTo("PXD012345"));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Any accession that is set is named before the run, well-formed or not: it lives in task
        /// settings, which are reused, so a valid accession left over from another dataset would
        /// otherwise reach this run's SDRF unseen. A malformed one is also called malformed. Neither
        /// is refused: the user's value is still written, because a warning costs nothing and a
        /// refusal would cost the whole file.
        /// </summary>
        [TestCase("PXD12345", true, true)]
        [TestCase("MSV000012345", true, true)]
        [TestCase("PXD012345", true, false)]
        [TestCase(" PXD012345 ", true, false)]
        // The spec's value type is PXD\d+: seven digits will come. RPXD is PRIDE's reprocessed-dataset
        // accession, legitimate for a re-analysis. Case is not the user's mistake worth a warning.
        [TestCase("PXD1234567", true, false)]
        [TestCase("RPXD012345", true, false)]
        [TestCase("pxd012345", true, false)]
        [TestCase(null, false, false)]
        [TestCase("", false, false)]
        // An isobaric search too: merged with #2817, an early return in its multiplex branch skipped
        // this warning (Alexander-Sol, #2817 and #2870 review).
        [TestCase("PXD012345", true, false, true)]
        [TestCase("PXD12345", true, true, true)]
        public static void ASetAccessionIsNamedBeforeTheRun(string accession, bool named, bool malformed,
            bool multiplex = false)
        {
            string folder = SetUpIsolatedRun(nameof(ASetAccessionIsNamedBeforeTheRun) + named + malformed + multiplex +
                accession?.Trim(), out string spectraPath, out _);
            ExperimentalDesign.WriteExperimentalDesignToFile(
                new List<SpectraFileInfo> { new(spectraPath, "condition", 0, 0, 0) });

            var task = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WriteSdrf = true,
                    ProteomeXchangeAccession = accession,
                    DoMultiplexQuantification = multiplex
                }
            };

            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                typeof(SearchTask)
                    .GetMethod("WarnAboutSdrfGaps", BindingFlags.NonPublic | BindingFlags.Instance)!
                    .Invoke(task, new object[] { new List<string> { spectraPath } });
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }

            var accessionWarnings = warnings.Where(w => w.Contains("ProteomeXchange accession")).ToList();
            Assert.That(accessionWarnings.Count, Is.EqualTo(named ? 1 : 0), string.Join(" | ", warnings));
            if (named)
            {
                Assert.That(accessionWarnings[0], Does.Contain("'" + accession.Trim() + "'"));
                Assert.That(accessionWarnings[0], Does.Contain("clear ProteomeXchangeAccession"));
                Assert.That(accessionWarnings[0].Contains("not the form of a ProteomeXchange accession"), Is.EqualTo(malformed));
            }

            Directory.Delete(folder, true);
        }

        #endregion

        #region An input SDRF's sample columns are copied (sdrf D50, QuantProject 037 PR-B)

        /// <summary>The sample columns <see cref="LabelFreeInputSdrf"/> gives, in the order the copy must keep.</summary>
        private static readonly string[] InputSampleColumns =
        {
            "source name", "characteristics[organism]", "characteristics[organism part]",
            "characteristics[organism part]", "characteristics[individual]", "characteristics[age]",
            "characteristics[sex]", "characteristics[disease]", "characteristics[cell type]",
            "characteristics[biological replicate]", "material type", "comment[sample preparation batch]",
            "comment[file uri]", "factor value[disease]", "factor value[treatment]"
        };

        /// <summary>
        /// A curated label-free SDRF for a.raw, b.raw and c.raw: a term with a key mzLib does not know (TA=x), a
        /// repeated column, `not available`, a `pooled` replicate, a blank comment, an accession, and one misspelled
        /// column (characteristic[mouse line]) that no copy carries.
        /// </summary>
        private static List<string[]> LabelFreeInputSdrf(params string[] dataFiles)
        {
            var rows = new List<string[]>
            {
                new[]
                {
                    "source name", "characteristics[organism]", "characteristics[organism part]",
                    "characteristics[organism part]", "characteristics[individual]", "characteristics[age]",
                    "characteristics[sex]", "characteristics[disease]", "characteristics[cell type]",
                    "characteristics[biological replicate]", "material type", "assay name", "technology type",
                    "comment[data file]", "comment[label]", "comment[proteomexchange accession number]",
                    "comment[sample preparation batch]", "comment[file uri]", "characteristic[mouse line]",
                    "factor value[disease]", "factor value[treatment]"
                }
            };
            var samples = new Dictionary<string, string[]>
            {
                ["a.raw"] = new[]
                {
                    "donor A liver", "NT=Mus musculus;AC=NCBITaxon:10090", "NT=liver;AC=UBERON:0002107;TA=x",
                    "NT=hepatic lobule;AC=UBERON:0004647", "donor A", "12M", "female", "NT=normal;AC=PATO:0000461",
                    "not available", "1", "tissue", "run a", "proteomic profiling by mass spectrometry", "a.raw",
                    "label free sample", "PXD012345", "batch 1", "", "C57BL/6", "normal", "vehicle"
                },
                ["b.raw"] = new[]
                {
                    "pool QC", "NT=Mus musculus;AC=NCBITaxon:10090", "NT=liver;AC=UBERON:0002107",
                    "NT=hepatic lobule;AC=UBERON:0004647", "pool", "not available", "not available",
                    "NT=normal;AC=PATO:0000461", "not available", "pooled", "tissue", "run b",
                    "proteomic profiling by mass spectrometry", "b.raw", "label free sample", "PXD012345", "batch 2",
                    "https://example.org/b.raw", "C57BL/6", "normal", "none"
                },
                ["c.raw"] = new[]
                {
                    "donor C liver", "NT=Mus musculus;AC=NCBITaxon:10090", "NT=liver;AC=UBERON:0002107",
                    "NT=hepatic lobule;AC=UBERON:0004647", "donor C", "14M", "male", "NT=normal;AC=PATO:0000461",
                    "not available", "2", "tissue", "run c", "proteomic profiling by mass spectrometry", "c.raw",
                    "label free sample", "PXD012345", "batch 1", "", "C57BL/6", "normal", "vehicle"
                }
            };
            rows.AddRange((dataFiles.Length == 0 ? samples.Keys.ToArray() : dataFiles).Select(f => samples[f]));
            return rows;
        }

        /// <summary>
        /// The whole copy (B6.1, B6.3, B6.5, B6.6), end to end. Two files are searched against an input SDRF that also
        /// describes a third: the two rows written carry the input's sample columns byte for byte, in the input's order,
        /// with the repeat repeated and the reserved words as written; c.raw, not searched, is left out. The assay half
        /// is still this search's own, and the input's accession is written because the task set none (B8).
        /// </summary>
        [Test]
        public static void AnInputSdrfsSampleColumnsAreCopiedCellForCell()
        {
            string output = RunSearchWithInputSdrf(nameof(AnInputSdrfsSampleColumnsAreCopiedCellForCell),
                LabelFreeInputSdrf(), out string folder, out var warnings);

            var input = new SdrfDocument(Path.Combine(folder, "input.sdrf.tsv"));
            input.LoadResults();
            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            Assert.That(document.Results.Select(r => r["source name"]), Is.EquivalentTo(new[] { "donor A liver", "pool QC" }),
                "One row per searched file, named as the input names it; c.raw was not searched, so its row is left out.");

            var copied = new HashSet<string>(InputSampleColumns, StringComparer.Ordinal);
            Assert.That(document.Header.Where(copied.Contains), Is.EqualTo(InputSampleColumns),
                "The input's sample columns, in the input's order, the repeated one twice.");

            foreach (SdrfRow row in document.Results)
            {
                SdrfRow source = input.Results.Single(r => r["source name"] == row["source name"]);
                foreach (string column in copied)
                    Assert.That(row.All(column), Is.EqualTo(source.All(column)),
                        $"'{column}' of '{row["source name"]}' is copied byte for byte.");
            }

            SdrfRow donor = document.Results.Single(r => r["source name"] == "donor A liver");
            SdrfRow pool = document.Results.Single(r => r["source name"] == "pool QC");
            Assert.That(donor.All("characteristics[organism part]"),
                Is.EqualTo(new[] { "NT=liver;AC=UBERON:0002107;TA=x", "NT=hepatic lobule;AC=UBERON:0004647" }),
                "A repeated column is written twice, in order, and a key mzLib does not know (TA=x) survives.");
            Assert.That(pool["characteristics[biological replicate]"], Is.EqualTo("pooled"));
            Assert.That(pool["characteristics[sex]"], Is.EqualTo("not available"));
            Assert.That(donor["comment[file uri]"], Is.EqualTo(""), "A blank comment is written blank, never invented.");

            Assert.That(donor["comment[data file]"], Is.EqualTo("a.mzml"), "The assay half is this search's: the file it read.");
            Assert.That(donor["comment[cleavage agent details]"], Does.Contain("Trypsin").IgnoreCase);
            Assert.That(donor["comment[label]"], Is.EqualTo("label free sample"));
            Assert.That(donor["comment[proteomexchange accession number]"], Is.EqualTo("PXD012345"),
                "The task set no accession, so the input's is written (B8).");

            Assert.That(SdrfValidator.Validate(document).Errors, Is.Empty, "mzLib's validator rejects the copied SDRF.");

            string notCarried = warnings.SingleOrDefault(w => w.Contains("not carried"));
            Assert.That(notCarried, Does.Contain("characteristic[mouse line]"),
                "An input column the written SDRF lacks is named (B9). " + string.Join(" | ", warnings));
            Assert.That(notCarried, Does.Not.Contain("comment[searched data file]").And.Not.Contain("comment[label]"));
            Assert.That(warnings.Any(w => w.Contains("since MetaMorpheus reads no input SDRF")), Is.False,
                "The coverage report must not say MetaMorpheus read no input SDRF when it did (B5).");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// The task's own accession wins over the input's (B8): the input's is only the fallback for an empty one.
        /// Checked on the builder options alone, without a search.
        /// </summary>
        [TestCase(null, "PXD012345")]
        [TestCase(" ", "PXD012345")]
        [TestCase(" PXD000001 ", "PXD000001")]
        public static void TheInputsAccessionIsWrittenOnlyWhenTheTaskSetsNone(string taskAccession, string written)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SdrfOutput_InputAccession");
            Directory.CreateDirectory(folder);
            string inputPath = Path.Combine(folder, "input.sdrf.tsv");
            File.WriteAllLines(inputPath, LabelFreeInputSdrf("a.raw").Select(cells => string.Join("\t", cells)));
            var scope = SdrfSearchScope.Restrict(new SdrfDocument(inputPath), new[] { "a.mzml" });

            string accession = (string)typeof(PostSearchAnalysisTask)
                .GetMethod("AccessionToWrite", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { taskAccession, scope });

            Assert.That(accession, Is.EqualTo(written));
            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Calibrate -> Search (B6.2): the input names the acquired sample1.mzML and the search reads
        /// sample1-calib.mzML. The row is still found, its sample columns copied, and the derivative named.
        /// </summary>
        [Test]
        public static void AfterCalibrationTheInputRowOfTheAcquiredFileIsCopied()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SdrfOutput_InputCalib");
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            string spectraPath = Path.Combine(folder, "sample1.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML"), spectraPath, true);
            string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta");
            ExperimentalDesign.WriteExperimentalDesignToFile(
                new List<SpectraFileInfo> { new(spectraPath, "condition", 0, 0, 0) });

            var input = LabelFreeInputSdrf("a.raw");
            input[1][Array.IndexOf(input[0], "comment[data file]")] = "sample1.mzML";
            string inputPath = Path.Combine(folder, "input.sdrf.tsv");
            File.WriteAllLines(inputPath, input.Select(cells => string.Join("\t", cells)));

            var search = BuildSearchTask(writeSdrf: true);
            search.SearchParameters.InputSdrfPath = inputPath;
            new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Task1-Calibrate", new CalibrationTask()), ("Task2-Search", search) },
                new List<string> { spectraPath }, new List<DbForTask> { new DbForTask(database, false) }, folder).Run();

            var document = new SdrfDocument(Path.Combine(folder, "Task2-Search", SdrfFileName));
            document.LoadResults();
            SdrfRow row = document.Results.Single();

            Assert.That(row["source name"], Is.EqualTo("donor A liver"));
            Assert.That(row["characteristics[age]"], Is.EqualTo("12M"));
            Assert.That(row["comment[data file]"], Is.EqualTo("sample1.mzML"));
            Assert.That(row["comment[searched data file]"], Is.EqualTo("sample1-calib.mzML"));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// A searched file the input has no row for, or (label-free) more than one row for, cannot be described from
        /// the input, so no SDRF is written (B2.3, B6.4). The crash record names the file and the search's own
        /// results stand. A missing row is also named before the run.
        /// </summary>
        [TestCase(0)]
        [TestCase(2)]
        public static void ASearchedFileWithoutExactlyOneInputRowWritesNoSdrfButKeepsTheSearch(int rowsForB)
        {
            var input = LabelFreeInputSdrf("a.raw");
            for (int i = 0; i < rowsForB; i++)
            {
                var row = LabelFreeInputSdrf("b.raw")[1];
                row[0] = "pool QC " + (i + 1);
                input.Add(row);
            }

            string output = RunSearchWithInputSdrf(nameof(ASearchedFileWithoutExactlyOneInputRowWritesNoSdrfButKeepsTheSearch) + rowsForB,
                input, out string folder, out var warnings);

            Assert.That(File.Exists(Path.Combine(output, SdrfFileName)), Is.False);
            Assert.That(File.Exists(Path.Combine(output, "AllPSMs.psmtsv")), Is.True, "The search's results stand.");
            string crash = File.ReadAllText(Path.Combine(output, "SdrfWriter_crash.txt"));
            Assert.That(crash, Does.Contain("b.mzml"));
            Assert.That(crash, Does.Not.Contain("a.mzml"), "a.mzml has its one row; only the problem is listed.");
            if (rowsForB == 0)
                Assert.That(warnings.Any(w => w.Contains("b.mzml") && w.Contains("no " + SdrfFileName)), Is.True,
                    "A file with no input row is named before the search. " + string.Join(" | ", warnings));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// TMT (B6.7): each channel row takes the input row of its file whose comment[label] is the channel's label,
        /// whatever order the input lists them in and whatever the case. A channel the input has no row for refuses
        /// the SDRF, naming the channel.
        /// </summary>
        [TestCase(false)]
        [TestCase(true)]
        public static void ATmtChannelRowCopiesTheInputRowWithItsLabel(bool inputLacksAChannel)
        {
            var input = new List<string[]>
            {
                new[]
                {
                    "source name", "characteristics[organism]", "characteristics[individual]",
                    "characteristics[biological replicate]", "assay name", "technology type", "comment[data file]",
                    "comment[label]", "factor value[treatment]"
                }
            };
            input.AddRange(Tmt11Channels
                .Select((tag, i) => new[]
                {
                    "donor " + (i + 1), "NT=Mus musculus;AC=NCBITaxon:10090", "D" + (i + 1), (i + 1).ToString(),
                    "run VA084TQ_6", "proteomic profiling by mass spectrometry", "VA084TQ_6.raw",
                    tag == "127N" ? "tmt127n" : "TMT" + tag, i % 2 == 0 ? "drug" : "vehicle"
                })
                .Where(r => !(inputLacksAChannel && r[7] == "TMT131C"))
                .Reverse());

            string root = RunTmtSearchWritingSdrf("SdrfOutput_TmtInput" + inputLacksAChannel, writeDesign: true,
                out string output, out _, inputSdrf: input);

            if (inputLacksAChannel)
            {
                Assert.That(File.Exists(Path.Combine(output, SdrfFileName)), Is.False);
                Assert.That(File.ReadAllText(Path.Combine(output, "SdrfWriter_crash.txt")), Does.Contain("TMT131C"));
                Directory.Delete(root, true);
                return;
            }

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();
            Assert.That(document.Results.Count, Is.EqualTo(Tmt11Channels.Length));
            for (int i = 0; i < Tmt11Channels.Length; i++)
            {
                SdrfRow row = document.Results.Single(r => r["comment[label]"] == "TMT" + Tmt11Channels[i]);
                Assert.That(row["source name"], Is.EqualTo("donor " + (i + 1)),
                    "The input's sample for this channel, not TmtDesign.txt's Sample" + (i + 1) + ".");
                Assert.That(row["characteristics[individual]"], Is.EqualTo("D" + (i + 1)));
                Assert.That(row["factor value[treatment]"], Is.EqualTo(i % 2 == 0 ? "drug" : "vehicle"));
            }

            Directory.Delete(root, true);
        }

        /// <summary>
        /// An isobaric input SDRF with one row per (data file, label, source name), each a different individual.
        /// </summary>
        private static List<string[]> TmtInputSdrf(params (string DataFile, string Label, string SourceName)[] channels)
        {
            var rows = new List<string[]>
            {
                new[]
                {
                    "source name", "characteristics[organism]", "characteristics[individual]",
                    "characteristics[biological replicate]", "assay name", "technology type", "comment[data file]",
                    "comment[label]", "factor value[treatment]"
                }
            };
            rows.AddRange(channels.Select((c, i) => new[]
            {
                c.SourceName, "NT=Mus musculus;AC=NCBITaxon:10090", "D" + (i + 1), (i + 1).ToString(),
                "run " + Path.GetFileNameWithoutExtension(c.DataFile), "proteomic profiling by mass spectrometry",
                c.DataFile, c.Label, i % 2 == 0 ? "drug" : "vehicle"
            }));
            return rows;
        }

        /// <summary>
        /// One plex whose TmtDesign.txt makes #2817 make up source names: channel 126 is unnamed, so it would be
        /// "Plex1 126", which the design gives 127N; and "HeLa" (Ctrl) and "hela" (Drug), one name when case is
        /// ignored, name two samples, so each would be scoped to its channel.
        /// </summary>
        private static readonly (string Name, string Channel, string Condition)[] DesignThatMakesUpSourceNames =
        {
            ("", "126", "Ctrl"), ("Plex1 126", "127N", "Drug"), ("HeLa", "128N", "Ctrl"), ("hela", "128C", "Drug")
        };

        /// <summary>
        /// With an input SDRF every source name is the input's (sdrf D50), so the SDRF makes none up from TmtDesign.txt,
        /// and neither of #2817's warnings about made-up names is true: no name is scoped to its plex or channel, and
        /// none is numbered because another sample has it. Without an input SDRF the same design gives both warnings,
        /// which is what makes their absence mean something.
        /// </summary>
        [TestCase(false)]
        [TestCase(true)]
        public static void WithAnInputSdrfNoSourceNameIsMadeUpSoNoneIsScopedOrNumbered(bool withInputSdrf)
        {
            var rows = IsobaricSdrfRows(nameof(WithAnInputSdrfNoSourceNameIsMadeUpSoNoneIsScopedOrNumbered) + withInputSdrf,
                new[] { "p1.mzML" },
                folder => DesignThatMakesUpSourceNames.Select((d, i) =>
                    TmtDesignRow(folder, "p1.mzML", "Plex1", d.Name, d.Channel, d.Condition, i + 1, 1)),
                out var warnings,
                withInputSdrf
                    ? TmtInputSdrf(DesignThatMakesUpSourceNames.Select((d, i) => ("p1.raw", "TMT" + d.Channel, "donor " + (i + 1))).ToArray())
                    : null);

            string SourceName(string channel) => rows.Single(r => r.Sample.Label?.Name == "TMT" + channel).Sample.SourceName;
            string allWarnings = string.Join(" | ", warnings);
            bool scoped = warnings.Any(w => w.Contains("more than one channel of the same plex"));
            bool numbered = warnings.Any(w => w.Contains("'Plex1 126' is written 'Plex1 126 (2)'"));

            if (!withInputSdrf)
            {
                Assert.That(SourceName("126"), Is.EqualTo("Plex1 126 (2)"));
                Assert.That(SourceName("128C"), Is.EqualTo("Plex1 hela 128C"));
                Assert.That(scoped && numbered, Is.True, "The design alone makes names up: " + allWarnings);
                return;
            }

            Assert.That(DesignThatMakesUpSourceNames.Select(d => SourceName(d.Channel)),
                Is.EqualTo(new[] { "donor 1", "donor 2", "donor 3", "donor 4" }), "Every source name is the input's.");
            Assert.That(scoped, Is.False, "No name is scoped to its channel: " + allWarnings);
            Assert.That(numbered, Is.False, "No made-up name is written, so none is numbered: " + allWarnings);
        }

        /// <summary>The same two warnings are not given before the search either, when an input SDRF gives the names.</summary>
        [TestCase(false)]
        [TestCase(true)]
        public static void WithAnInputSdrfNoSourceNameIsWarnedAsScopedOrNumberedBeforeTheSearch(bool withInputSdrf)
        {
            var warnings = WarningsBeforeAnIsobaricSearch(
                nameof(WithAnInputSdrfNoSourceNameIsWarnedAsScopedOrNumberedBeforeTheSearch) + withInputSdrf, "TMT11",
                spectra => DesignThatMakesUpSourceNames
                    .Select((d, i) => $"{spectra}\tPlex1\t{d.Name}\t{d.Channel}\t{d.Condition}\t{i + 1}\t1\t1\tstudy sample")
                    .ToArray(),
                withInputSdrf
                    ? TmtInputSdrf(DesignThatMakesUpSourceNames
                        .Select((d, i) => ("PrunedDbSpectra.raw", "TMT" + d.Channel, "donor " + (i + 1))).ToArray())
                    : null);

            string allWarnings = string.Join(" | ", warnings);
            bool scoped = warnings.Any(w => w.Contains("more than one channel of the same plex"));
            bool numbered = warnings.Any(w => w.Contains("'Plex1 126' is written 'Plex1 126 (2)'"));

            if (!withInputSdrf)
            {
                Assert.That(scoped && numbered, Is.True, "The design alone makes names up: " + allWarnings);
                return;
            }

            Assert.That(scoped, Is.False, "No name is scoped to its channel: " + allWarnings);
            Assert.That(numbered, Is.False, "No made-up name is written, so none is numbered: " + allWarnings);
            Assert.That(warnings.Any(w => w.Contains("input SDRF")), Is.False,
                "The input describes every channel, so it gives no warning of its own: " + allWarnings);
        }

        /// <summary>A TmtDesign.txt row per channel for the spectra file, one study sample each, all usable.</summary>
        private static string[] TmtDesignOf(string spectra, string[] channels, string sampleType = "study sample") =>
            channels.Select((c, i) => $"{spectra}\tPlex1\tS{i + 1}\t{c}\tCtrl\t{i + 1}\t1\t1\t{sampleType}").ToArray();

        /// <summary>
        /// A channel's input row is found by mzLib's one channel/label rule (SdrfIsobaricLabels.TryGetChannel, sdrf D51),
        /// not by comparing the label as written bare (QuantProject 067 R2). An accessioned cell, which 6 of the corpus's
        /// 86 isobaric SDRFs write, and TMT10's last channel written TMT131N, as 2 of its 25 TMT10 SDRFs write it, each
        /// find their row. The written label is still the search's own.
        /// </summary>
        [TestCase("accessioned126", "NT=TMT126;AC=PRIDE:0000516", "TMT131N")]
        [TestCase("accessioned131", "tmt126", "NT=TMT131;AC=PRIDE:0000529")]
        public static void AChannelsInputRowIsFoundWhicheverWayItsLabelIsWritten(string caseName, string label126,
            string label131)
        {
            var rows = IsobaricSdrfRows(nameof(AChannelsInputRowIsFoundWhicheverWayItsLabelIsWritten) + caseName,
                new[] { "p1.mzML" },
                folder => new[]
                {
                    TmtDesignRow(folder, "p1.mzML", "Plex1", "S1", "126", "Ctrl", 1, 1),
                    TmtDesignRow(folder, "p1.mzML", "Plex1", "S2", "131N", "Drug", 2, 1)
                },
                out _, TmtInputSdrf(("p1.raw", label126, "donor 1"), ("p1.raw", label131, "donor 2")), "TMT10");

            Assert.That(rows.Select(r => (r.Sample.Label?.Name, r.Sample.SourceName)),
                Is.EquivalentTo(new[] { ("TMT126", "donor 1"), ("TMT131", "donor 2") }));
        }

        /// <summary>
        /// A row whose label names no channel of the searched kit is named in the refusal, since it, not a missing row,
        /// is why the channel has none (QuantProject 067 R2): here a TMTpro label in a TMT10 search.
        /// </summary>
        [Test]
        public static void AChannelWhoseRowIsLabelledOutsideTheKitIsRefusedNamingThatLabel()
        {
            var refusal = Assert.Throws<InvalidDataException>(() => IsobaricSdrfRows(
                nameof(AChannelWhoseRowIsLabelledOutsideTheKitIsRefusedNamingThatLabel), new[] { "p1.mzML" },
                folder => new[]
                {
                    TmtDesignRow(folder, "p1.mzML", "Plex1", "S1", "126", "Ctrl", 1, 1),
                    TmtDesignRow(folder, "p1.mzML", "Plex1", "S2", "131N", "Drug", 2, 1)
                },
                out _, TmtInputSdrf(("p1.raw", "TMT126", "donor 1"), ("p1.raw", "TMT134N", "donor 2")), "TMT10"));

            Assert.That(refusal!.Message, Does.Contain("p1.mzML channel TMT131 has no row")
                .And.Contain("'TMT134N'").And.Contain("TMT10"));
        }

        /// <summary>
        /// Before a TMT search whose TmtDesign.txt cannot give channel rows (QuantProject 067 R1): the SDRF can then
        /// describe the file only once, so an input written a row per channel is refused at write time. With such an
        /// input, a warning names the file and its rows and says no SDRF will be written, and none promises that the
        /// SDRF will describe the file once. With one input row for the file the SDRF is written, and nothing says not.
        /// </summary>
        [TestCase("unusable design", true)]
        [TestCase("unusable design", false)]
        [TestCase("unknown tag", true)]
        [TestCase("placeholder", true)]
        public static void WithoutChannelRowsAnInputWithARowPerChannelIsNamedBeforeTheRun(string why, bool rowPerChannel)
        {
            var channels = new[] { "126", "127N" };
            var warnings = WarningsBeforeAnIsobaricSearch(
                nameof(WithoutChannelRowsAnInputWithARowPerChannelIsNamedBeforeTheRun) + why.Replace(" ", "") + rowPerChannel,
                why == "unknown tag" ? "NotAnIsobaricTag" : "TMT11",
                spectra => why == "placeholder"
                    ? new[] { spectra + "\tPlex1" + new string('\t', 7) }
                    : TmtDesignOf(spectra, channels, why == "unusable design" ? "not a sample type" : "study sample"),
                TmtInputSdrf((rowPerChannel ? channels : channels.Take(1))
                    .Select((c, i) => ("PrunedDbSpectra.raw", "TMT" + c, "donor " + (i + 1))).ToArray()));

            string allWarnings = string.Join(" | ", warnings);
            bool refusalNamed = warnings.Any(w => w.Contains("no " + SdrfFileName)
                && w.Contains("PrunedDbSpectra.mzml has 2 rows in it", StringComparison.OrdinalIgnoreCase));
            Assert.That(refusalNamed, Is.EqualTo(rowPerChannel), allWarnings);
            Assert.That(warnings.Any(w => w.Contains("SDRF will describe")), Is.False,
                "With an input SDRF the file is described by its input row, which may be refused: " + allWarnings);
        }

        /// <summary>
        /// Before the run only the input rows the SDRF copies are checked for blank sample facts (QuantProject 067 R3). A
        /// row for a channel TmtDesign.txt does not annotate is never written (QP-S6), so its blank is not named; a blank
        /// in a row that is copied still is.
        /// </summary>
        [TestCase(false)]
        [TestCase(true)]
        public static void OnlyTheInputRowsTheSdrfCopiesAreCheckedForBlanksBeforeTheRun(bool blankInACopiedRow)
        {
            var input = TmtInputSdrf(("PrunedDbSpectra.raw", "TMT126", "donor 1"),
                ("PrunedDbSpectra.raw", "TMT127N", "donor 2"), ("PrunedDbSpectra.raw", "TMT128N", "donor 3"));
            // Row 2 is 127N's, which the design annotates; row 3 is 128N's, which it does not.
            input[blankInACopiedRow ? 2 : 3][Array.IndexOf(input[0], "characteristics[individual]")] = "";

            var warnings = WarningsBeforeAnIsobaricSearch(
                nameof(OnlyTheInputRowsTheSdrfCopiesAreCheckedForBlanksBeforeTheRun) + blankInACopiedRow, "TMT11",
                spectra => TmtDesignOf(spectra, new[] { "126", "127N" }), input);

            string allWarnings = string.Join(" | ", warnings);
            Assert.That(warnings.Any(w => w.Contains("blank")), Is.EqualTo(blankInACopiedRow), allWarnings);
            if (blankInACopiedRow)
                Assert.That(warnings.Single(w => w.Contains("blank")), Does.Contain("donor 2"));
        }

        /// <summary>
        /// Before the run, the input rows the writer would refuse are named (G51), not only the files with no row: a
        /// channel the input has no row for, and a label-free file with two.
        /// </summary>
        [Test]
        public static void AChannelWithoutAnInputRowAndAFileWithTwoAreNamedBeforeTheRun()
        {
            var tmt = WarningsBeforeAnIsobaricSearch(nameof(AChannelWithoutAnInputRowAndAFileWithTwoAreNamedBeforeTheRun),
                "TMT11", spectra => TmtDesignOf(spectra, new[] { "126", "127N" }),
                TmtInputSdrf(("PrunedDbSpectra.raw", "TMT126", "donor 1")));
            Assert.That(tmt.Any(w => w.Contains("no " + SdrfFileName) && w.Contains("channel TMT127N has no row")), Is.True,
                string.Join(" | ", tmt));

            var input = LabelFreeInputSdrf();
            var second = (string[])input[1].Clone();
            second[0] = "donor A again";
            input.Add(second);
            var labelFree = WarningsBeforeASearchWithAnInputSdrf(
                nameof(AChannelWithoutAnInputRowAndAFileWithTwoAreNamedBeforeTheRun) + "LabelFree", input);
            Assert.That(labelFree.Any(w => w.Contains("no " + SdrfFileName) && w.Contains("a.mzml has 2 rows in it")),
                Is.True, string.Join(" | ", labelFree));
        }

        /// <summary>A pipeline sets the input SDRF in the task file (B1, B6.8), so it survives a TOML round trip.</summary>
        [Test]
        public static void TheInputSdrfPathSurvivesTheTaskFile()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SdrfOutput_InputSdrfToml");
            Directory.CreateDirectory(folder);
            string tomlPath = Path.Combine(folder, "SearchTask.toml");

            Assert.That(new SearchParameters().InputSdrfPath, Is.Null, "Off unless a pipeline sets it.");
            var task = BuildSearchTask(writeSdrf: true);
            task.SearchParameters.InputSdrfPath = @"C:\data\PXD012345.sdrf.tsv";
            Toml.WriteFile(task, tomlPath, MetaMorpheusTask.tomlConfig);

            var reread = Toml.ReadFile<SearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
            Assert.That(reread.SearchParameters.InputSdrfPath, Is.EqualTo(@"C:\data\PXD012345.sdrf.tsv"));

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// Before the run (B7): every blank sample fact in the rows the search will use is named, with its file, source
        /// name and column, because the writer refuses one and the SDRF would be lost. A blank comment is not named: it
        /// is written blank. Rows for files the search does not read are not checked.
        /// </summary>
        [Test]
        public static void ABlankSampleFactInTheInputIsNamedBeforeTheRun()
        {
            var input = LabelFreeInputSdrf();
            input[1][Array.IndexOf(input[0], "characteristics[age]")] = "";
            input[2][Array.IndexOf(input[0], "factor value[treatment]")] = " ";
            input[3][Array.IndexOf(input[0], "characteristics[sex]")] = "";

            var warnings = WarningsBeforeASearchWithAnInputSdrf(nameof(ABlankSampleFactInTheInputIsNamedBeforeTheRun), input);

            string blanks = warnings.Single(w => w.Contains("blank"));
            Assert.That(blanks, Does.Contain("a.mzml").And.Contain("donor A liver").And.Contain("characteristics[age]"));
            Assert.That(blanks, Does.Contain("b.mzml").And.Contain("pool QC").And.Contain("factor value[treatment]"));
            Assert.That(blanks, Does.Not.Contain("comment[file uri]"), "a.raw's blank comment is written blank.");
            Assert.That(blanks, Does.Not.Contain("characteristics[sex]"), "c.raw is not searched.");
            Assert.That(blanks, Does.Contain("no " + SdrfFileName));
        }

        /// <summary>
        /// Before the run (B4): an input SDRF that cannot be read is named, with the consequence. And with an input SDRF,
        /// a missing ExperimentalDesign.tsv no longer warns that the sample columns will say `not available`: the input
        /// gives them, and only the technical replicate and fraction are lost.
        /// </summary>
        [Test]
        public static void AnUnreadableInputSdrfIsNamedBeforeTheRun_AndNoDesignIsNoLongerASampleGap()
        {
            var warnings = WarningsBeforeASearchWithAnInputSdrf(
                nameof(AnUnreadableInputSdrfIsNamedBeforeTheRun_AndNoDesignIsNoLongerASampleGap), inputSdrf: null);

            Assert.That(warnings.Any(w => w.Contains("missing.sdrf.tsv") && w.Contains("no " + SdrfFileName)), Is.True,
                string.Join(" | ", warnings));
            Assert.That(warnings.Any(w => w.Contains("sample columns will say 'not available'")), Is.False,
                string.Join(" | ", warnings));
            Assert.That(warnings.Any(w => w.Contains(GlobalVariables.ExperimentalDesignFileName)
                                          && w.Contains("technical replicate and fraction")), Is.True,
                string.Join(" | ", warnings));

            var readable = WarningsBeforeASearchWithAnInputSdrf(
                nameof(AnUnreadableInputSdrfIsNamedBeforeTheRun_AndNoDesignIsNoLongerASampleGap) + "Readable",
                LabelFreeInputSdrf());
            Assert.That(readable.Any(w => w.Contains("no " + SdrfFileName)), Is.False,
                "A readable input that covers every searched file with no blank sample fact gives no input warning. " +
                string.Join(" | ", readable));
        }

        /// <summary>
        /// Searches a.mzml and b.mzml (copies of the pruned spectra) with SDRF output on, <paramref name="inputSdrf"/>
        /// as the input SDRF, and no ExperimentalDesign.tsv. Returns the output folder.
        /// </summary>
        private static string RunSearchWithInputSdrf(string testName, IEnumerable<string[]> inputSdrf,
            out string folder, out List<string> warnings)
        {
            folder = SetUpIsolatedRun(testName, out string spectraPath, out DbForTask database);
            var searched = new List<string> { Path.Combine(folder, "a.mzml"), Path.Combine(folder, "b.mzml") };
            foreach (string path in searched)
                File.Copy(spectraPath, path, true);

            string inputPath = Path.Combine(folder, "input.sdrf.tsv");
            File.WriteAllLines(inputPath, inputSdrf.Select(cells => string.Join("\t", cells)));

            var task = BuildSearchTask(writeSdrf: true);
            task.SearchParameters.InputSdrfPath = inputPath;
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(output);

            var captured = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => captured.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                task.RunTask(output, new List<DbForTask> { database }, searched, "sdrf-input");
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }

            warnings = captured;
            return output;
        }

        /// <summary>
        /// The warnings SearchTask gives before a label-free search of a.mzml and b.mzml with SDRF output on and an input
        /// SDRF of <paramref name="inputSdrf"/>, or a path that does not exist when it is null. No ExperimentalDesign.tsv.
        /// </summary>
        private static List<string> WarningsBeforeASearchWithAnInputSdrf(string testName, IEnumerable<string[]> inputSdrf)
        {
            string folder = SetUpIsolatedRun(testName, out _, out _);
            string inputPath = Path.Combine(folder, inputSdrf is null ? "missing.sdrf.tsv" : "input.sdrf.tsv");
            if (inputSdrf is not null)
                File.WriteAllLines(inputPath, inputSdrf.Select(cells => string.Join("\t", cells)));

            var task = new SearchTask { SearchParameters = new SearchParameters { WriteSdrf = true, InputSdrfPath = inputPath } };
            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                typeof(SearchTask)
                    .GetMethod("WarnAboutSdrfGaps", BindingFlags.NonPublic | BindingFlags.Instance)!
                    .Invoke(task, new object[] { new List<string> { Path.Combine(folder, "a.mzml"), Path.Combine(folder, "b.mzml") } });
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
                Directory.Delete(folder, true);
            }

            return warnings;
        }

        #endregion

        #region Helpers

        /// <summary>
        /// Both resolvers are private statics on the adapter. Reaching them by reflection follows
        /// the idiom PostSearchAnalysisTaskTests already uses for this class's private writers, and
        /// is far cheaper than driving a whole DIA or SILAC search to observe one cell.
        /// </summary>
        /// <summary>
        /// Only the analysis type matters to the SDRF; the XIC/grouping machinery does not, so it is
        /// left null rather than stood up.
        /// </summary>
        private static DIAparameters DiaParams(DIAanalysisType type) =>
            new(type, null, null, null, default);

        private static CvParam InvokeAcquisitionMethod(CommonParameters common) =>
            (CvParam)typeof(PostSearchAnalysisTask)
                .GetMethod("ResolveAcquisitionMethod", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { common });

        /// <summary>
        /// A modification is identified by (ModificationType, IdWithMotif), and the lookup must key
        /// on both. Keyed on the id alone, a search that fixes one mod gets every known mod sharing
        /// its id -- a wrong mod, or the same column written twice.
        /// </summary>
        [Test]
        public static void AModificationIsResolvedByItsTypeAsWellAsItsId()
        {
            var shared = GlobalVariables.AllModsKnown
                .GroupBy(m => m.IdWithMotif)
                .FirstOrDefault(g => g.Select(m => m.ModificationType).Distinct().Count() > 1);
            Assert.That(shared, Is.Not.Null,
                "The known mods no longer contain two types sharing an id; this test needs one.");

            Modification wanted = shared!.First();
            var resolved = (IReadOnlyList<Modification>)typeof(PostSearchAnalysisTask)
                .GetMethod("ResolveModifications", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { new List<(string, string)> { (wanted.ModificationType, wanted.IdWithMotif) } });

            Assert.That(resolved, Is.Not.Empty);
            Assert.That(resolved.Select(m => m.ModificationType).Distinct().ToList(), Is.EqualTo(new[] { wanted.ModificationType }),
                $"Asked for '{wanted.IdWithMotif}' of type '{wanted.ModificationType}' only, and " +
                $"'{wanted.IdWithMotif}' is also known under another type.");
        }

        private static CvParam InvokeLabel(SearchParameters searchParameters) =>
            (CvParam)typeof(PostSearchAnalysisTask)
                .GetMethod("ResolveLabel", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { searchParameters });

        private static CvParam InvokeChannelLabel(IsobaricMassTagType tagType, string channel) =>
            (CvParam)typeof(PostSearchAnalysisTask)
                .GetMethod("ResolveChannelLabel", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { tagType, channel });

        private static string InvokeChannelRowsUnusableReason(TmtFileInfo file, IsobaricMassTagType tagType) =>
            (string)typeof(PostSearchAnalysisTask)
                .GetMethod("ChannelRowsUnusableReason", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { file, tagType });

        /// <summary>
        /// The channel rows of <paramref name="file"/>, with source names resolved over
        /// <paramref name="document"/>: every file of the search's channel rows, or just this one.
        /// </summary>
        private static List<SdrfRowInput> InvokeChannelRows(TmtFileInfo file, IsobaricMassTagType tagType,
            params TmtFileInfo[] document)
        {
            var args = new object[] { document.Length > 0 ? document : new[] { file }, new List<string>(), null };
            var sourceNames = typeof(PostSearchAnalysisTask)
                .GetMethod("IsobaricSourceNames", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, args);
            return ((IEnumerable<SdrfRowInput>)typeof(PostSearchAnalysisTask)
                    .GetMethod("BuildChannelRows", BindingFlags.NonPublic | BindingFlags.Static)!
                    .Invoke(null, new[] { file, tagType, null, new SdrfAssay { DataFileName = "run1.raw", AssayName = "run run1" },
                        sourceNames }))
                .ToList();
        }

        private static HashSet<string> InvokeSampleNamesRepeatedWithinAPlex(params TmtFileInfo[] files) =>
            (HashSet<string>)typeof(PostSearchAnalysisTask)
                .GetMethod("SampleNamesRepeatedWithinAPlex", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { files });

        private static HashSet<string> InvokeSampleNamesReusedForDifferentSamples(params TmtFileInfo[] files) =>
            (HashSet<string>)typeof(PostSearchAnalysisTask)
                .GetMethod("SampleNamesReusedForDifferentSamples", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { files });

        /// <summary>
        /// A search whose parameters are cheap but not degenerate. Notch/parsimony settings match
        /// the other end-to-end search tests so this runs in the same time they do.
        /// </summary>
        private static SearchTask BuildSearchTask(bool writeSdrf) => new SearchTask
        {
            SearchParameters = new SearchParameters
            {
                DoParsimony = true,
                SearchType = SearchType.Classic,
                SearchTarget = true,
                DecoyType = DecoyType.None,
                WriteSdrf = writeSdrf
            }
        };

        /// <summary>
        /// Copies the smallest committed spectra/database pair into a folder of this test's own.
        ///
        /// The copy matters. ExperimentalDesign.tsv has to sit BESIDE the spectra file, so a test
        /// that used the shared TestData folder in place would be writing a design file other tests
        /// can see -- and the feature under test changes behaviour based on whether that file
        /// exists, so a stray one turns a real failure green.
        /// </summary>
        private static string SetUpIsolatedRun(string testName, out string spectraPath, out DbForTask database)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SdrfOutput_" + testName);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);

            string spectraSource = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            spectraPath = Path.Combine(folder, "PrunedDbSpectra.mzml");
            File.Copy(spectraSource, spectraPath, true);

            database = new DbForTask(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"), false);

            return folder;
        }

        /// <summary>
        /// The happy path shared by the write tests: isolated folder, a valid single-file
        /// experimental design beside the spectra, SDRF requested. Returns the output folder.
        /// </summary>
        private static string RunSearchWritingSdrf(string testName, out string folder, out string spectraPath)
        {
            folder = SetUpIsolatedRun(testName, out spectraPath, out DbForTask database);

            // Written through the API rather than by hand so the file's shape is whatever
            // MetaMorpheus itself considers correct, not whatever this test believes it to be.
            ExperimentalDesign.WriteExperimentalDesignToFile(
                new List<SpectraFileInfo> { new(spectraPath, "condition", 0, 0, 0) });

            var task = BuildSearchTask(writeSdrf: true);
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(output);

            task.RunTask(output, new List<DbForTask> { database }, new List<string> { spectraPath }, "sdrf");
            return output;
        }

        #endregion
    }
}
