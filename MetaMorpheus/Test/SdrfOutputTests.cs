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
using Omics.Modifications;
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
        /// comment[data file] names the ORIGINAL acquisition, never a calibrated or averaged
        /// derivative this run happened to produce. An SDRF describes data as acquired; pointing it
        /// at an intermediate makes the row unjoinable to the deposited dataset.
        /// </summary>
        [Test]
        public static void TheWrittenSdrfNamesTheOriginalDataFile()
        {
            string output = RunSearchWritingSdrf(nameof(TheWrittenSdrfNamesTheOriginalDataFile),
                out string folder, out string spectraPath);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();
            SdrfRow row = document.Results.Single();

            Assert.That(row["comment[data file]"], Is.EqualTo(Path.GetFileName(spectraPath)));
            Assert.That(row["comment[data file]"], Does.Not.Contain("-calib"),
                "A calibrated intermediate is not the acquisition the SDRF describes.");

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
        /// reporter m/z order, and an Empty channel gets no row.
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

            var expectedChannels = Tmt11Channels.Take(Tmt11Channels.Length - 1).ToList();
            Assert.That(document.Results.Count, Is.EqualTo(expectedChannels.Count),
                "One row per annotated, non-empty channel of the one file this run searched.");

            for (int i = 0; i < expectedChannels.Count; i++)
            {
                SdrfRow row = document.Results[i];
                Assert.That(row["comment[label]"], Does.Contain("NT=TMT" + expectedChannels[i] + ";"),
                    "Rows follow the plex's reporter m/z order, not the design file's row order.");
                Assert.That(row["source name"], Is.EqualTo("Sample" + (i + 1)),
                    "The source name is the sample the design put in this channel.");
                Assert.That(row["comment[data file]"], Is.EqualTo("VA084TQ_6.mzML"));
                Assert.That(row["assay name"], Is.EqualTo("run VA084TQ_6"),
                    "Every channel of one file is one assay.");
            }

            Assert.That(document.Results.Any(r => r["comment[label]"].Contains("TMT131C")), Is.False,
                "The Empty channel holds no sample and gets no row.");

            SdrfValidationResult validation = SdrfValidator.Validate(document);
            Assert.That(validation.Errors, Is.Empty, "mzLib's validator rejects the TMT SDRF: " + validation);

            SdrfQuantAudit audit = SdrfQuantAuditor.Audit(document);
            Assert.That(audit.Kind, Is.EqualTo(SdrfQuantKind.ChannelLevel), audit.ToReport());
            Assert.That(audit.Channels, Has.Count.EqualTo(expectedChannels.Count), audit.ToReport());
            Assert.That(audit.DuplicateFileLabelPairs, Is.Empty, audit.ToReport());

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
        /// Runs the TMT11 fixture MultiplexQuantificationTests uses, with SDRF output on, in a folder of
        /// its own. With <paramref name="writeDesign"/>, a TmtDesign.txt goes beside the spectra: the
        /// last channel is Empty, and the rows are written in reverse so document order cannot pass
        /// for reporter m/z order. Returns the root folder to delete.
        /// </summary>
        private static string RunTmtSearchWritingSdrf(string folderName, bool writeDesign,
            out string output, out List<string> warnings)
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, folderName);
            if (Directory.Exists(root)) Directory.Delete(root, true);
            string dataFolder = Path.Combine(root, "data");
            Directory.CreateDirectory(dataFolder);

            string mzml = Path.Combine(dataFolder, "VA084TQ_6.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\VA084TQ_6.mzML"), mzml);

            if (writeDesign)
            {
                var designRows = Tmt11Channels
                    .Select((tag, i) => $"{mzml}\tPlex1\tSample{i + 1}\t{tag}\tCond{(i % 2 == 0 ? "A" : "B")}\t{i / 2 + 1}\t1\t1\t" +
                                        (i == Tmt11Channels.Length - 1 ? "empty" : "study sample"))
                    .Reverse();
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

        private static CvParam InvokeLabel(SearchParameters searchParameters) =>
            (CvParam)typeof(PostSearchAnalysisTask)
                .GetMethod("ResolveLabel", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { searchParameters });

        private static CvParam InvokeChannelLabel(IsobaricMassTagType tagType, string channel) =>
            (CvParam)typeof(PostSearchAnalysisTask)
                .GetMethod("ResolveChannelLabel", BindingFlags.NonPublic | BindingFlags.Static)!
                .Invoke(null, new object[] { tagType, channel });

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
