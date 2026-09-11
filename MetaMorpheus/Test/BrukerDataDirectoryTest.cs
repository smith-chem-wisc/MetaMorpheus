using System.IO;
using EngineLayer;
using EngineLayer.Util;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Covers the single predicate every entry point (CMD, main GUI, MetaDraw) now uses to decide whether a ".d" folder
    /// holds Bruker data MetaMorpheus can read. Folders are synthesized rather than committed - the classifier only looks
    /// at which files are present, so real vendor data would add megabytes and prove nothing extra.
    /// </summary>
    [TestFixture]
    public static class BrukerDataDirectoryTest
    {
        private static string _scratch;

        [OneTimeSetUp]
        public static void SetUp()
        {
            _scratch = Path.Combine(TestContext.CurrentContext.TestDirectory, "BrukerDataDirectoryTest");
            if (Directory.Exists(_scratch))
            {
                Directory.Delete(_scratch, true);
            }
            Directory.CreateDirectory(_scratch);
        }

        [OneTimeTearDown]
        public static void TearDown()
        {
            if (Directory.Exists(_scratch))
            {
                Directory.Delete(_scratch, true);
            }
        }

        /// <summary>
        /// Creates a ".d" folder containing the named files, each with placeholder content. The classifier never opens
        /// them, so their contents are irrelevant.
        /// </summary>
        private static string MakeDotDFolder(string folderName, params string[] innerFileNames)
        {
            string folder = Path.Combine(_scratch, folderName);
            Directory.CreateDirectory(folder);
            foreach (string innerFileName in innerFileNames)
            {
                File.WriteAllText(Path.Combine(folder, innerFileName), "placeholder");
            }
            return folder;
        }

        [Test]
        public static void QTofFolderWithBafIsValid()
        {
            Assert.That(BrukerDataDirectory.IsValid(MakeDotDFolder("qtof.d", "analysis.baf")), Is.True);
        }

        [Test]
        public static void TimsTofFolderWithTdfAndSidecarIsValid()
        {
            Assert.That(BrukerDataDirectory.IsValid(MakeDotDFolder("tims.d", "analysis.tdf", "analysis.tdf_bin")), Is.True);
        }

        /// <summary>
        /// mzLib's TimsTofFileReader needs the .tsf_bin sidecar just as it needs .tdf_bin, so a .tsf folder counts as
        /// readable Bruker data. MetaMorpheus rejected these outright before, even though mzLib could open them.
        /// </summary>
        [Test]
        public static void TimsTofFolderWithTsfAndSidecarIsValid()
        {
            Assert.That(BrukerDataDirectory.IsValid(MakeDotDFolder("tsf.d", "analysis.tsf", "analysis.tsf_bin")), Is.True);
        }

        [Test]
        [TestCase("no_bin.d", "analysis.tdf")]
        [TestCase("no_tsf_bin.d", "analysis.tsf")]
        public static void TimsTofFolderMissingItsBinarySidecarIsNotValid(string folderName, string innerFileName)
        {
            // mzLib classifies these as timsTOF from the SQLite file alone, but the reader throws without the sidecar.
            Assert.That(BrukerDataDirectory.IsValid(MakeDotDFolder(folderName, innerFileName)), Is.False);
        }

        [Test]
        public static void EmptyDotDFolderIsNotValid()
        {
            Assert.That(BrukerDataDirectory.IsValid(MakeDotDFolder("empty.d")), Is.False);
        }

        [Test]
        public static void FolderWithUnrelatedContentsIsNotValid()
        {
            Assert.That(BrukerDataDirectory.IsValid(MakeDotDFolder("junk.d", "notes.txt", "spectra.mzML")), Is.False);
        }

        [Test]
        public static void FolderNotEndingInDotDIsNotValid()
        {
            Assert.That(BrukerDataDirectory.IsValid(MakeDotDFolder("plain_folder", "analysis.baf")), Is.False);
        }

        [Test]
        public static void NonexistentAndNullPathsAreNotValid()
        {
            Assert.That(BrukerDataDirectory.IsValid(Path.Combine(_scratch, "does_not_exist.d")), Is.False);
            Assert.That(BrukerDataDirectory.IsValid(null), Is.False);
        }

        [Test]
        public static void IsValidIsCaseInsensitiveOnTheDotDSuffix()
        {
            string folder = MakeDotDFolder("upper.D", "analysis.baf");
            Assert.That(BrukerDataDirectory.IsValid(folder), Is.True);
        }

        [Test]
        [TestCase("analysis.baf")]
        [TestCase("analysis.tdf")]
        [TestCase("analysis.tdf_bin")]
        [TestCase("analysis.tsf")]
        [TestCase("analysis.tsf_bin")]
        public static void InnerFileRedirectsToItsParentDotDFolder(string innerFileName)
        {
            // one folder that satisfies every flavour, so a single fixture covers all five inner files
            string folder = MakeDotDFolder("redirect_" + innerFileName.Replace('.', '_') + ".d",
                "analysis.baf", "analysis.tdf", "analysis.tdf_bin", "analysis.tsf", "analysis.tsf_bin");

            Assert.That(BrukerDataDirectory.TryGetParentDotDFolder(Path.Combine(folder, innerFileName), out string dotDFolder), Is.True);
            Assert.That(dotDFolder, Is.EqualTo(folder));
        }

        [Test]
        public static void InnerFileInAnUnreadableFolderDoesNotRedirect()
        {
            // analysis.tdf with no sidecar - the folder is not readable, so the caller must keep the original path
            string folder = MakeDotDFolder("unreadable.d", "analysis.tdf");

            Assert.That(BrukerDataDirectory.TryGetParentDotDFolder(Path.Combine(folder, "analysis.tdf"), out string dotDFolder), Is.False);
            Assert.That(dotDFolder, Is.Null);
        }

        [Test]
        public static void NonBrukerFileDoesNotRedirect()
        {
            string folder = MakeDotDFolder("mixed.d", "analysis.baf", "spectra.mzML");

            Assert.That(BrukerDataDirectory.TryGetParentDotDFolder(Path.Combine(folder, "spectra.mzML"), out string dotDFolder), Is.False);
            Assert.That(dotDFolder, Is.Null);
        }

        [Test]
        public static void NullAndEmptyPathsDoNotRedirect()
        {
            Assert.That(BrukerDataDirectory.TryGetParentDotDFolder(null, out _), Is.False);
            Assert.That(BrukerDataDirectory.TryGetParentDotDFolder(string.Empty, out _), Is.False);
        }

        [Test]
        [TestCase(".baf", true)]
        [TestCase(".TDF", true)]
        [TestCase(".tdf_bin", true)]
        [TestCase(".tsf", true)]
        [TestCase(".tsf_bin", true)]
        [TestCase(".d", false)]
        [TestCase(".mzml", false)]
        [TestCase(null, false)]
        public static void IsInnerFileExtensionRecognizesTheBrukerInnerFiles(string extension, bool expected)
        {
            Assert.That(BrukerDataDirectory.IsInnerFileExtension(extension), Is.EqualTo(expected));
        }

        /// <summary>
        /// The regression that motivated this fixture: a duplicated assignment in SetUpGlobalVariables silently dropped
        /// ".baf" from the accepted list, which left MetaDraw unable to see or open Bruker qTOF data even though every
        /// other piece of the feature was in place. Every Bruker inner file must survive into the accepted list, because
        /// CMD and MetaDraw both gate on it before the redirect ever runs.
        /// </summary>
        [Test]
        public static void EveryBrukerExtensionIsAnAcceptedSpectraFormat()
        {
            GlobalVariables.SetUpGlobalVariables();

            Assert.That(GlobalVariables.AcceptedSpectraFormats, Contains.Item(BrukerDataDirectory.DotD));
            foreach (string extension in BrukerDataDirectory.InnerFileExtensions)
            {
                Assert.That(GlobalVariables.AcceptedSpectraFormats, Contains.Item(extension),
                    $"{extension} must be an accepted spectra format or MetaDraw and the CMD directory scan will reject it");
            }
        }

        /// <summary>
        /// Accepted formats are compared against ToLowerInvariant()'d extensions at every call site, so an upper-case
        /// entry in the list would never match.
        /// </summary>
        [Test]
        public static void AcceptedSpectraFormatsAreLowerCase()
        {
            GlobalVariables.SetUpGlobalVariables();

            foreach (string extension in GlobalVariables.AcceptedSpectraFormats)
            {
                Assert.That(extension, Is.EqualTo(extension.ToLowerInvariant()));
            }
        }
    }
}
