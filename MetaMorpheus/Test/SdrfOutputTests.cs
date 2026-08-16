using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using NUnit.Framework;
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

        /// <summary>
        /// The columns SDRF-Proteomics requires of every document.
        ///
        /// DUPLICATED FROM mzLib ON PURPOSE, and the duplication is the point rather than an
        /// oversight. mzLib's SdrfValidator -- which owns this list -- is internal: it ships as
        /// NuGet, and a type goes public only when something outside its assembly already calls it.
        /// Nothing here does, so a MetaMorpheus test cannot ask mzLib whether the file it just
        /// wrote is valid; it can only read the header back and check.
        ///
        /// The accepted cost of that decision is that this list can drift from the one that
        /// actually governs. If mzLib changes it, this fails and should be updated to match --
        /// never the other way round.
        /// </summary>
        private static readonly string[] RequiredSdrfColumns =
        {
            "source name",
            "assay name",
            "technology type",
            "characteristics[organism]",
            "characteristics[organism part]",
            "characteristics[biological replicate]",
            "comment[data file]",
            "comment[instrument]",
            "comment[label]",
            "comment[cleavage agent details]",
            "comment[technical replicate]",
            "comment[fraction identifier]",
            "comment[proteomics data acquisition method]"
        };

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
        /// The decision this feature rests on: asking for an SDRF is an assertion that the sample
        /// metadata exists, and if it does not the run is refused BEFORE a spectrum is read.
        ///
        /// Failing late would be the same defect wearing a different coat -- a three-hour search
        /// thrown away, or worse, a document padded with "not available" that passes every validator
        /// and answers no question.
        /// </summary>
        [Test]
        public static void OptingInWithoutAnExperimentalDesign_RefusesTheRun()
        {
            string folder = SetUpIsolatedRun(nameof(OptingInWithoutAnExperimentalDesign_RefusesTheRun),
                out string spectraPath, out DbForTask database);

            var task = BuildSearchTask(writeSdrf: true);
            string output = Path.Combine(folder, "TaskOutput");
            Directory.CreateDirectory(output);

            var exception = Assert.Throws<MetaMorpheusException>(() =>
                task.RunTask(output, new List<DbForTask> { database }, new List<string> { spectraPath }, "sdrf-no-design"));

            Assert.That(exception.Message, Does.Contain(GlobalVariables.ExperimentalDesignFileName),
                "The message has to name the file the user is expected to produce.");

            // The point of validating at task level rather than at write time: nothing was searched.
            Assert.That(File.Exists(Path.Combine(output, "AllPSMs.psmtsv")), Is.False,
                "The run must be refused BEFORE the search, not after it.");

            Directory.Delete(folder, true);
        }

        /// <summary>
        /// A design file that exists but cannot be parsed is refused too. Present-but-broken is the
        /// more dangerous case: it looks satisfied.
        /// </summary>
        [Test]
        public static void OptingInWithAnUnusableExperimentalDesign_RefusesTheRun()
        {
            string folder = SetUpIsolatedRun(nameof(OptingInWithAnUnusableExperimentalDesign_RefusesTheRun),
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

            var exception = Assert.Throws<MetaMorpheusException>(() =>
                task.RunTask(output, new List<DbForTask> { database }, new List<string> { spectraPath }, "sdrf-bad-design"));

            Assert.That(exception.Message, Does.Contain("cannot be used as it stands"),
                "A malformed design must be reported as malformed, not as missing.");
            Assert.That(File.Exists(Path.Combine(output, "AllPSMs.psmtsv")), Is.False,
                "The run must be refused BEFORE the search, not after it.");

            Directory.Delete(folder, true);
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
        /// Every column SDRF-Proteomics requires is present.
        ///
        /// This is the assertion the rest of the stack cannot make on its own, and that is exactly
        /// why it is here. A column that is MISSING is invisible to both instruments built to catch
        /// thin metadata: SdrfCoverage measures the fill rate of columns that exist, so an absent
        /// column has no fill rate to report, and mzLib's validator -- which does check -- is
        /// internal and unreachable from MetaMorpheus. So nothing between the search and a mined
        /// corpus notices.
        /// </summary>
        [Test]
        public static void TheWrittenSdrfCarriesEveryRequiredColumn()
        {
            string output = RunSearchWritingSdrf(nameof(TheWrittenSdrfCarriesEveryRequiredColumn),
                out string folder, out string spectraPath);

            var document = new SdrfDocument(Path.Combine(output, SdrfFileName));
            document.LoadResults();

            var missing = RequiredSdrfColumns.Where(c => !document.Header.Contains(c)).ToList();

            Assert.That(missing, Is.Empty,
                "SDRF-Proteomics requires these columns of every document, and they are absent from " +
                "the one this search wrote: " + string.Join(", ", missing) + ". A document missing a " +
                "required column cannot be pooled across experiments, which is the only reason to " +
                "write one.");

            Directory.Delete(folder, true);
        }

        #endregion

        #region Helpers

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
