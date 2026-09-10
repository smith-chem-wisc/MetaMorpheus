using EngineLayer;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    /// <summary>
    /// The custom-file contract, enforced. Every file a user is expected to edit must be seeded with a
    /// template when it is absent and left strictly alone when it is present -- that second half is the
    /// one that matters, and the one #2752 broke.
    /// </summary>
    [TestFixture]
    public static class CustomDataFileTests
    {
        private static string _dir;

        [SetUp]
        public static void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "CustomDataFileTests", Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public static void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        /// <summary>
        /// The rule the whole class exists for: an existing file is never touched, whatever is in it --
        /// including content the template would never have produced.
        /// </summary>
        [Test]
        public static void EnsureExists_NeverRewritesAnExistingFile()
        {
            string path = Path.Combine(_dir, "already_here.tsv");
            const string userContent = "the user typed this\tand meant it\n";
            File.WriteAllText(path, userContent);

            CustomDataFile.EnsureExists(path, () => "TEMPLATE THAT MUST NOT APPEAR", "test");

            Assert.That(File.ReadAllText(path), Is.EqualTo(userContent));
        }

        /// <summary>
        /// An empty file is still the user's file. Seeding it would silently replace a file someone had
        /// truncated on purpose, and "it was empty so I overwrote it" is exactly the reasoning that loses
        /// data.
        /// </summary>
        [Test]
        public static void EnsureExists_TreatsAnEmptyExistingFileAsTheUsers()
        {
            string path = Path.Combine(_dir, "empty.tsv");
            File.WriteAllText(path, string.Empty);

            CustomDataFile.EnsureExists(path, () => "TEMPLATE", "test");

            Assert.That(File.ReadAllText(path), Is.Empty);
        }

        [Test]
        public static void EnsureExists_SeedsWhenAbsent_AndCreatesTheDirectory()
        {
            string path = Path.Combine(_dir, "nested", "seeded.tsv");
            Assert.That(File.Exists(path), Is.False, "precondition");

            CustomDataFile.EnsureExists(path, () => "Name\tValue", "test");

            Assert.That(File.ReadAllText(path), Is.EqualTo("Name\tValue"));
        }

        /// <summary>
        /// The template must not be built when the file already exists -- reading an embedded resource on
        /// every startup for a file that is already there is waste, and the deferred Func is the only
        /// thing preventing it.
        /// </summary>
        [Test]
        public static void EnsureExists_DoesNotBuildTheTemplateWhenTheFileExists()
        {
            string path = Path.Combine(_dir, "present.tsv");
            File.WriteAllText(path, "x");
            bool built = false;

            CustomDataFile.EnsureExists(path, () => { built = true; return "y"; }, "test");

            Assert.That(built, Is.False);
        }

        [Test]
        public static void EnsureExists_ReportsAFailureRatherThanLeavingTheFileMissing()
        {
            string path = Path.Combine(_dir, "boom.tsv");

            var ex = Assert.Throws<MetaMorpheusException>(() =>
                CustomDataFile.EnsureExists(path, () => throw new InvalidOperationException("no template"), "test thing"));

            Assert.Multiple(() =>
            {
                Assert.That(ex.Message, Does.Contain("test thing"), "the message has to say which file");
                Assert.That(ex.Message, Does.Contain(path));
                Assert.That(File.Exists(path), Is.False);
            });
        }

        /// <summary>
        /// The template is the shipped file's banner and header with the data dropped. Keeping the banner
        /// is the point: the user opens a file that documents its own format, which is what makes the
        /// custom-protease template worth copying everywhere else.
        /// </summary>
        [Test]
        public static void BannerAndHeader_KeepsCommentsAndHeader_AndDropsEveryDataRow()
        {
            string shipped = Path.Combine(_dir, "shipped.tsv");
            File.WriteAllLines(shipped, new[]
            {
                "# what this file is",
                "# how to edit it",
                "Name\tSites\tMass",
                "Trypsin\tKR\t0",
                "LysC\tK\t0",
            });

            string template = CustomDataFile.BannerAndHeaderFromFile(shipped, "Name\t");

            var lines = template.Split(new[] { "\r\n", "\n" }, StringSplitOptions.RemoveEmptyEntries);
            Assert.Multiple(() =>
            {
                Assert.That(lines, Is.EqualTo(new[] { "# what this file is", "# how to edit it", "Name\tSites\tMass" }));
                Assert.That(template, Does.Not.Contain("Trypsin"), "a data row would become a fake custom entry");
                Assert.That(template, Does.Not.Contain("LysC"));
            });
        }

        /// <summary>
        /// Every shipped file whose custom counterpart we seed must actually contain the header we look
        /// for. If a header is renamed upstream this silently produces a data-row-free but header-free
        /// template, and the user's first custom entry lands in a file the parser rejects.
        /// </summary>
        [Test]
        public static void BannerAndHeader_FindsTheHeaderInTheRealShippedCrosslinkerFile()
        {
            string shipped = Path.Combine(GlobalVariables.DataDir, "Data", "Crosslinkers.tsv");
            Assume.That(File.Exists(shipped), "shipped crosslinker file is required for this test");

            string template = CustomDataFile.BannerAndHeaderFromFile(shipped, "Name\t");

            Assert.Multiple(() =>
            {
                Assert.That(template, Does.StartWith("Name\t"));
                Assert.That(template.Split(new[] { "\r\n", "\n" }, StringSplitOptions.RemoveEmptyEntries).Length,
                    Is.EqualTo(1), "LoadCrosslinkers skips only line 1, so this template must be exactly the header");
                Assert.That(template, Does.Not.Contain("DSSO"), "no shipped crosslinker may leak into the custom file");
            });
        }

        /// <summary>
        /// Every custom file is seeded by startup, and none of them is a file the installer or the build
        /// also writes -- so listing them here is the inventory that a new custom file has to join.
        /// </summary>
        /// <remarks>
        /// Runs the startup itself rather than reading whatever the process happens to have lying around.
        /// Test order is not fixed, and CustomAminoAcidsTest deletes CustomAminoAcids.txt when it finishes
        /// ("Delete so it doesn't crash the next time") without seeding it again -- so asserting on ambient
        /// state passes or fails depending on what ran first. Re-running SetUpGlobalVariables is what the
        /// application does on launch, and is the pattern LoadDigestionAgentTest already uses to put the
        /// globals back.
        /// </remarks>
        [Test]
        public static void EveryKnownCustomFileIsSeededByStartup()
        {
            GlobalVariables.SetUpGlobalVariables();
            string d = GlobalVariables.DataDir;
            var expected = new Dictionary<string, string>
            {
                ["custom proteases"] = GlobalVariables.CustomProteasePath,
                ["custom rnases"] = GlobalVariables.CustomRnasePath,
                ["custom monosaccharides"] = GlobalVariables.CustomMonosaccharidePath,
                ["custom O-glycan database"] = GlobalVariables.CustomOGlycanDatabasePath,
                ["custom crosslinkers"] = Path.Combine(d, "Data", "CustomCrosslinkers.tsv"),
                ["custom modifications"] = Path.Combine(d, "Mods", "CustomModifications.txt"),
                ["custom RNA modifications"] = Path.Combine(d, "Mods", "RnaCustomModifications.txt"),
                ["custom amino acids"] = Path.Combine(d, "CustomAminoAcids", "CustomAminoAcids.txt"),
            };

            var missing = expected.Where(p => !File.Exists(p.Value)).Select(p => p.Key).ToList();
            Assert.That(missing, Is.Empty, "not seeded by SetUpGlobalVariables: " + string.Join(", ", missing));
        }

        /// <summary>
        /// A hand-edited crosslinker file picks up blank lines and notes. Both used to reach
        /// ParseCrosslinkerFromString and surface as an unhandled IndexOutOfRangeException during startup,
        /// before any window opened, with nothing naming the file.
        /// </summary>
        [Test]
        public static void LoadCrosslinkers_SkipsBlankAndCommentLines()
        {
            string path = Path.Combine(_dir, "CustomCrosslinkers.tsv");
            File.WriteAllLines(path, new[]
            {
                "Name\tCrosslinkAminoAcid\tCrosslinkerAminoAcid2\tCleavable\tDissociationType\tCrosslinkerTotalMass\tCrosslinkerShortMass\tCrosslinkerLongMass\tQuenchMassH2O\tQuenchMassNH2\tQuenchMassTris",
                "# a note the user left themselves",
                "",
                "MyLinker\tK\tK\tT\tCID|HCD\t158.0038\t54.01056\t85.982635\t176.0143\t175.0303\t279.0777",
                "   ",
            });

            var loaded = Crosslinker.LoadCrosslinkers(path).ToList();

            Assert.That(loaded.Select(p => p.CrosslinkerName), Is.EqualTo(new[] { "MyLinker" }));
        }

        /// <summary>
        /// A genuinely malformed row still has to fail -- but as a message naming the file and the line,
        /// not as IndexOutOfRangeException from inside a LINQ iterator.
        /// </summary>
        [Test]
        public static void LoadCrosslinkers_NamesTheFileAndLineForAShortRow()
        {
            string path = Path.Combine(_dir, "CustomCrosslinkers.tsv");
            File.WriteAllLines(path, new[]
            {
                "Name\tCrosslinkAminoAcid\tCrosslinkerAminoAcid2\tCleavable\tDissociationType\tCrosslinkerTotalMass\tCrosslinkerShortMass\tCrosslinkerLongMass\tQuenchMassH2O\tQuenchMassNH2\tQuenchMassTris",
                "OopsOnlyThree\tK\tK",
            });

            var ex = Assert.Throws<MetaMorpheusException>(() => Crosslinker.LoadCrosslinkers(path).ToList());

            Assert.Multiple(() =>
            {
                Assert.That(ex.Message, Does.Contain("CustomCrosslinkers.tsv"));
                Assert.That(ex.Message, Does.Contain("Line 2"));
                Assert.That(ex.Message, Does.Contain("OopsOnlyThree"), "the offending line is what the user has to fix");
            });
        }
    }
}
