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
            Directory.CreateDirectory(output);
            var task = new PostSearchAnalysisTask
            {
                Parameters = new PostSearchAnalysisParameters
                {
                    OutputFolder = output,
                    SearchTaskId = "no-files",
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
        /// </summary>
        [Test]
        public static void ALabelledSearchIsWarnedThatItsLabelWillNotBeFilledIn()
        {
            string folder = SetUpIsolatedRun(nameof(ALabelledSearchIsWarnedThatItsLabelWillNotBeFilledIn),
                out string spectraPath, out _);

            var task = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WriteSdrf = true,
                    SilacLabels = new List<SilacLabel> { new('K', 'a', "C{6}H{12}N{2}O{1}", 6.020129) }
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
        [TestCase(IsobaricMassTagType.TMT10, " 131n ", "PRIDE:0000580")]
        public static void AnIsobaricChannelResolvesToItsPrideTerm(IsobaricMassTagType tagType, string channel, string accession)
        {
            Assert.That(InvokeChannelLabel(tagType, channel)?.Accession, Is.EqualTo(accession));
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
            Assert.That(warnings.Any(w => w.Contains("no " + GlobalVariables.ExperimentalDesignFileName)), Is.False,
                "An isobaric search does not read ExperimentalDesign.tsv, so its absence is not reported.");

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
                "TMT11-plex on K", spectra => new[] { spectra + "\tPlex1" + new string('\t', 7) });

            Assert.That(warnings, Has.Count.EqualTo(1), string.Join(" | ", warnings));
            Assert.That(warnings.Single(), Does.Contain("no annotated channels")
                .And.Contain("without their channels or samples"));
        }

        [Test]
        public static void AChannelOffTheSearchedPlexIsWarnedBeforeTheSearch()
        {
            var warnings = WarningsBeforeAnIsobaricSearch(nameof(AChannelOffTheSearchedPlexIsWarnedBeforeTheSearch),
                "TMT11-plex on K", spectra => new[] { $"{spectra}	Plex1	Sample1	127	CondA	1	1	1	study sample" });

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
            Func<string, string[]> designRows)
        {
            string folder = SetUpIsolatedRun(testName, out string spectraPath, out _);
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
                Directory.Delete(folder, true);
            }

            return warnings;
        }

        /// <summary>
        /// Runs the TMT11 fixture MultiplexQuantificationTests uses, with SDRF output on, in a folder of
        /// its own. With <paramref name="writeDesign"/>, a TmtDesign.txt goes beside the spectra: the
        /// last channel is Empty, and the rows are written in reverse so document order cannot pass
        /// for reporter m/z order. With <paramref name="designIsUnusable"/>, every row's Sample Type is
        /// one the design reader refuses. Returns the root folder to delete.
        /// </summary>
        private static string RunTmtSearchWritingSdrf(string folderName, bool writeDesign,
            out string output, out List<string> warnings, bool designIsUnusable = false,
            bool emptyChannelIsUnnamed = false, bool placeholderOnly = false)
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

        private static List<SdrfRowInput> InvokeChannelRows(TmtFileInfo file, IsobaricMassTagType tagType) =>
            ((IEnumerable<SdrfRowInput>)typeof(PostSearchAnalysisTask)
                .GetMethod("BuildChannelRows", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { file, tagType, null, new SdrfAssay { DataFileName = "run1.raw", AssayName = "run run1" } }))
            .ToList();

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
