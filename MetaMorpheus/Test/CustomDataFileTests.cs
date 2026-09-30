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
    [NonParallelizable]
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
                    Is.EqualTo(1), "the shipped file has no banner, so the template is exactly the header");
                Assert.That(template, Does.Not.Contain("DSSO"), "no shipped crosslinker may leak into the custom file");
            });
        }

        /// <summary>
        /// Every custom file is seeded by startup, and none of them is a file the installer or the build
        /// also writes -- so listing them here is the inventory that a new custom file has to join.
        /// </summary>
        /// <remarks>
        /// Deletes the files first, so the assertion is about what startup creates and not about what
        /// earlier runs left in the test's DataDir. The originals are put back afterwards, since other
        /// tests edit some of them.
        /// </remarks>
        [Test]
        public static void EveryKnownCustomFileIsSeededByStartup()
        {
            string d = GlobalVariables.DataDir;
            var expected = new Dictionary<string, string>
            {
                ["custom proteases"] = GlobalVariables.CustomProteasePath,
                ["custom rnases"] = GlobalVariables.CustomRnasePath,
                ["custom monosaccharides"] = GlobalVariables.CustomMonosaccharidePath,
                ["custom O-glycan database"] = GlobalVariables.CustomOGlycanDatabasePath,
                ["custom N-glycan database"] = GlobalVariables.CustomNGlycanDatabasePath,
                ["custom crosslinkers"] = Path.Combine(d, "Data", "CustomCrosslinkers.tsv"),
                ["custom modifications"] = Path.Combine(d, "Mods", "CustomModifications.txt"),
                ["custom RNA modifications"] = Path.Combine(d, "Mods", "RnaCustomModifications.txt"),
                ["custom amino acids"] = Path.Combine(d, "CustomAminoAcids", "CustomAminoAcids.txt"),
            };

            var originals = expected.Values.Where(File.Exists).ToDictionary(p => p, File.ReadAllBytes);
            try
            {
                foreach (var path in expected.Values)
                {
                    File.Delete(path);
                }

                GlobalVariables.SetUpGlobalVariables();

                var missing = expected.Where(p => !File.Exists(p.Value)).Select(p => p.Key).ToList();
                Assert.That(missing, Is.Empty, "not seeded by SetUpGlobalVariables: " + string.Join(", ", missing));
            }
            finally
            {
                foreach (var original in originals)
                {
                    File.WriteAllBytes(original.Key, original.Value);
                }
                GlobalVariables.SetUpGlobalVariables();
            }
        }

        /// <summary>
        /// A note or a blank line above the header must not push the header into the data rows, where it
        /// passes the column count and fails on double.Parse("CrosslinkerTotalMass").
        /// </summary>
        [Test]
        public static void LoadCrosslinkers_FindsTheHeaderBelowANoteOrBlankLine()
        {
            string path = Path.Combine(_dir, "CustomCrosslinkers.tsv");
            File.WriteAllLines(path, new[]
            {
                "# a note the user left above the header",
                "",
                CrosslinkerHeader,
                "MyLinker\tK\tK\tT\tCID|HCD\t158.0038\t54.01056\t85.982635\t176.0143\t175.0303\t279.0777",
            });

            var loaded = Crosslinker.LoadCrosslinkers(path).ToList();

            Assert.That(loaded.Select(p => p.CrosslinkerName), Is.EqualTo(new[] { "MyLinker" }));
        }

        /// <summary>
        /// The same file through the real startup, which is where the crash surfaced: before any window,
        /// with nothing naming the file.
        /// </summary>
        [Test]
        public static void Startup_SurvivesACustomCrosslinkerFileWithANoteAboveTheHeader()
        {
            string path = Path.Combine(GlobalVariables.DataDir, "Data", "CustomCrosslinkers.tsv");
            byte[] original = File.Exists(path) ? File.ReadAllBytes(path) : null;
            try
            {
                File.WriteAllLines(path, new[]
                {
                    "# my linkers",
                    CrosslinkerHeader,
                    "MyStartupLinker\tK\tK\tT\tCID|HCD\t158.0038\t54.01056\t85.982635\t176.0143\t175.0303\t279.0777",
                });

                Assert.DoesNotThrow(GlobalVariables.SetUpGlobalVariables);
                Assert.That(GlobalVariables.Crosslinkers.Select(c => c.CrosslinkerName), Does.Contain("MyStartupLinker"));
            }
            finally
            {
                if (original == null)
                {
                    File.Delete(path);
                }
                else
                {
                    File.WriteAllBytes(path, original);
                }
                GlobalVariables.SetUpGlobalVariables();
            }
        }

        /// <summary>
        /// A blank line inside the shipped banner, or an indented comment, used to end the scan before the
        /// header and hand back a headerless or empty template.
        /// </summary>
        [TestCase("# first", "", "# after a blank line")]
        [TestCase("  # an indented comment", "# a plain one")]
        public static void BannerAndHeader_KeepsTheHeaderPastBlankAndIndentedCommentLines(params string[] banner)
        {
            string shipped = Path.Combine(_dir, "shipped.tsv");
            File.WriteAllLines(shipped, banner.Concat(new[] { "Name\tSites", "Trypsin\tKR" }));

            string template = CustomDataFile.BannerAndHeaderFromFile(shipped, "Name\t");

            var lines = template.Split(new[] { "\r\n", "\n" }, StringSplitOptions.None);
            Assert.Multiple(() =>
            {
                Assert.That(lines.Take(banner.Length), Is.EqualTo(banner), "the banner is kept as it was");
                Assert.That(lines[banner.Length], Is.EqualTo("Name\tSites"));
                Assert.That(template, Does.Not.Contain("Trypsin"));
            });
        }

        /// <summary>
        /// A shipped file whose header has gone, or moved below the data, is reported rather than turned
        /// into a template that has no header.
        /// </summary>
        [Test]
        public static void BannerAndHeader_ReportsAShippedFileWithNoHeaderBeforeTheData()
        {
            string shipped = Path.Combine(_dir, "shipped.tsv");
            File.WriteAllLines(shipped, new[] { "# banner", "Trypsin\tKR", "Name\tSites" });

            var ex = Assert.Throws<MetaMorpheusException>(() => CustomDataFile.BannerAndHeaderFromFile(shipped, "Name\t"));

            Assert.That(ex.Message, Does.Contain("Name"));
        }

        /// <summary>
        /// An empty template would be written as an empty file, which rule 1 then protects forever.
        /// </summary>
        [TestCase("")]
        [TestCase("  \r\n")]
        public static void EnsureExists_RefusesAnEmptyTemplate(string template)
        {
            string path = Path.Combine(_dir, "blank.tsv");

            var ex = Assert.Throws<MetaMorpheusException>(() => CustomDataFile.EnsureExists(path, () => template, "test thing"));

            Assert.Multiple(() =>
            {
                Assert.That(ex.Message, Does.Contain("test thing"));
                Assert.That(File.Exists(path), Is.False, "nothing is written, so the next startup tries again");
            });
        }

        private const string CrosslinkerHeader =
            "Name\tCrosslinkAminoAcid\tCrosslinkerAminoAcid2\tCleavable\tDissociationType\tCrosslinkerTotalMass\tCrosslinkerShortMass\tCrosslinkerLongMass\tQuenchMassH2O\tQuenchMassNH2\tQuenchMassTris";

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

        /// <summary>
        /// A hand-written file with no header row at all must not lose a crosslinker. Detecting the header
        /// by position (the original) or by "first line that is not a note" (its first fix) both consumed a
        /// real row: with a '#' note on top, the note was skipped as a note and MyLinker1 was eaten as the
        /// header. Silently, which is worse than the crash the note used to cause. The header is now the
        /// first line that actually LOOKS like the header.
        /// </summary>
        [Test]
        public static void LoadCrosslinkers_KeepsEveryRowOfAHeaderlessFile()
        {
            string path = Path.Combine(_dir, "CustomCrosslinkers.tsv");
            File.WriteAllLines(path, new[]
            {
                "# my own crosslinkers",
                "MyLinker1\tK\tK\tT\tCID\t158.0\t0\t0\t1\t1\t1",
                "MyLinker2\tK\tK\tT\tCID\t159.0\t0\t0\t1\t1\t1",
            });

            var loaded = Crosslinker.LoadCrosslinkers(path).ToList();

            Assert.That(loaded.Select(p => p.CrosslinkerName), Is.EqualTo(new[] { "MyLinker1", "MyLinker2" }));
        }

        /// <summary>
        /// A row with the right column count but an unreadable mass is the last route to a raw .NET crash
        /// out of SetUpGlobalVariables: the column-count guard passes it and double.Parse throws
        /// FormatException naming neither the file nor the line. It now follows the same failure contract as
        /// every other custom-file reader.
        /// </summary>
        [Test]
        public static void LoadCrosslinkers_NamesTheFileAndLineForARowItCannotParse()
        {
            string path = Path.Combine(_dir, "CustomCrosslinkers.tsv");
            File.WriteAllLines(path, new[]
            {
                CrosslinkerHeader,
                "BadMass\tK\tK\tT\tCID\tnot-a-number\t0\t0\t1\t1\t1",
            });

            var ex = Assert.Throws<MetaMorpheusException>(() => Crosslinker.LoadCrosslinkers(path).ToList());

            Assert.Multiple(() =>
            {
                Assert.That(ex.Message, Does.Contain("CustomCrosslinkers.tsv"));
                Assert.That(ex.Message, Does.Contain("Line 2"));
                Assert.That(ex.Message, Does.Contain("BadMass"), "the offending line is what the user has to fix");
            });
        }

        /// <summary>
        /// Seeding writes beside the destination and moves the finished file into place, so a write that
        /// fails leaves nothing behind. It matters more here than it looks: rule 1 protects whatever sits at
        /// the destination on every later startup, so a partial file written there would be permanent.
        /// </summary>
        [Test]
        public static void EnsureExists_LeavesNoPartialFileWhenItCannotFinish()
        {
            // a directory where the file should go: the template builds and writes, and the move fails
            string path = Path.Combine(_dir, "blocked.tsv");
            Directory.CreateDirectory(path);

            Assert.Throws<MetaMorpheusException>(() =>
                CustomDataFile.EnsureExists(path, () => "a perfectly good template", "test thing"));

            Assert.That(File.Exists(path + ".tmp"), Is.False, "a half-written file must not be left behind");
        }
    }
}
