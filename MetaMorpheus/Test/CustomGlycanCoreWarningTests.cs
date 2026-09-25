using EngineLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
using System;
using System.IO;

namespace Test
{
    /// <summary>
    /// A custom glycan that does not start from HexNAc is warned about, not refused: nearly every O-glycan
    /// starts from GalNAc and every N-glycan from GlcNAc, but O-mannose, O-fucose and O-glucose glycans are
    /// real, and one of the shipped O-glycan databases holds four of them.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // registers custom monosaccharides, which is process-wide state
    public static class CustomGlycanCoreWarningTests
    {
        private static string _dir;

        [SetUp]
        public static void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "CustomGlycanCoreWarningTests", Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public static void TearDown()
        {
            Glycan.ResetCustomMonosaccharides();
            if (Directory.Exists(_dir))
            {
                Directory.Delete(_dir, true);
            }
        }

        private static string Path_(string name) => Path.Combine(_dir, name);

        /// <summary>
        /// A structure is judged by its root, the monosaccharide on the amino acid. A composition has no
        /// order, so it is judged by whether it holds a HexNAc at all -- "Hex(1)HexNAc(1)" starts from
        /// HexNAc as far as anyone can tell from the counts.
        /// </summary>
        [TestCase("(H)", true)]
        [TestCase("(H(N))", true)]
        [TestCase("(F(H))", true)]
        [TestCase("Hex(1)", true)]
        [TestCase("Hex(2)Fuc(1)", true)]
        [TestCase("(N)", false)]
        [TestCase("(N(H(A)))", false)]
        [TestCase("HexNAc(1)", false)]
        [TestCase("Hex(1)HexNAc(1)", false)]
        public static void AGlycanThatDoesNotStartWithHexNAcIsWarnedAbout(string glycan, bool warned)
        {
            foreach (bool isOGlycan in new[] { true, false })
            {
                string warning = GlycanDatabase.CoreWarning(glycan, isOGlycan);

                if (warned)
                {
                    Assert.That(warning, Does.Contain(glycan).And.Contain("HexNAc"), $"isOGlycan={isOGlycan}");
                }
                else
                {
                    Assert.That(warning, Is.Null, $"isOGlycan={isOGlycan}");
                }
            }
        }

        /// <summary>The O- and N-glycan messages name the core that is missing, which differs.</summary>
        [Test]
        public static void TheWarningNamesTheCoreForTheGlycanClass()
        {
            Assert.That(GlycanDatabase.CoreWarning("(H)", true), Does.Contain("GalNAc").And.Contain("O-mannose"));
            Assert.That(GlycanDatabase.CoreWarning("(H)", false), Does.Contain("GlcNAc"));
        }

        /// <summary>
        /// A warning, not a refusal: the glycan is written, and the warning is there for the window to show.
        /// </summary>
        [TestCase("(H(H))", true)]
        [TestCase("Hex(2)", false)]
        public static void AGlycanThatDoesNotStartWithHexNAcIsStillAdded(string glycan, bool isOGlycan)
        {
            string path = Path_("custom.gdb");

            Assert.DoesNotThrow(() => GlycanDatabase.PersistCustomGlycan(glycan, path, isOGlycan));

            Assert.That(File.ReadAllLines(path), Is.EqualTo(new[] { glycan }));
            Assert.That(GlycanDatabase.CoreWarning(glycan, isOGlycan), Is.Not.Null);
        }

        /// <summary>
        /// What is not a glycan is the validator's and the loader's business; this has nothing to add.
        /// </summary>
        [TestCase("")]
        [TestCase("# a comment")]
        [TestCase("not a glycan")]
        [TestCase("Nope(1)")]
        public static void SomethingThatIsNotAGlycanIsNotWarnedAbout(string text)
        {
            Assert.That(GlycanDatabase.CoreWarning(text, true), Is.Null);
        }

        /// <summary>
        /// The same check over a whole database, for startup: every offending line named with its line
        /// number, in one message; comments, annotations and tab columns handled as the loader handles them.
        /// </summary>
        [Test]
        public static void EveryGlycanInADatabaseThatDoesNotStartWithHexNAcIsNamedWithItsLine()
        {
            string path = Path_("OGlycan_Custom.gdb");
            File.WriteAllLines(path, new[]
            {
                "# Hex(1) in a comment is not a glycan",
                "HexNAc(1)Hex(1)",
                "Hex(1) # O-mannose",
                "",
                "Hex(1)Fuc(1)\t308.11",
            });

            string warning = GlycanDatabase.CoreWarningsFor(path, true);

            Assert.That(warning, Does.Contain("'OGlycan_Custom.gdb'"));
            Assert.That(warning, Does.Contain("line 3: Hex(1)").And.Contain("line 5: Hex(1)Fuc(1)"));
            Assert.That(warning, Does.Not.Contain("line 1").And.Not.Contain("line 2"));
        }

        [Test]
        public static void ADatabaseThatAllStartsWithHexNAcOrIsMissingGivesNoWarning()
        {
            string path = Path_("NGlycan_Custom.gdb");
            File.WriteAllLines(path, new[] { "(N(N(H(H)(H))))", "(N(N(H)))" });

            Assert.That(GlycanDatabase.CoreWarningsFor(path, false), Is.Null);
            Assert.That(GlycanDatabase.CoreWarningsFor(Path_("missing.gdb"), false), Is.Null);
        }

        /// <summary>The seeded templates are all banner, and say nothing.</summary>
        [TestCase("EngineLayer.Glycan_Mods.OGlycan_Custom.gdb", true)]
        [TestCase("EngineLayer.Glycan_Mods.NGlycan_Custom.gdb", false)]
        public static void TheSeededTemplatesGiveNoWarning(string resource, bool isOGlycan)
        {
            string path = Path_("template.gdb");
            File.WriteAllText(path, CustomDataFile.EmbeddedText(typeof(GlobalVariables).Assembly, resource));

            Assert.That(GlycanDatabase.CoreWarningsFor(path, isOGlycan), Is.Null);
        }

        /// <summary>
        /// A code declared in MonosaccharidesCustom.tsv is judged the same way: it is not HexNAc, so a
        /// structure rooted on it is warned about.
        /// </summary>
        [Test]
        public static void AStructureRootedOnACustomMonosaccharideIsWarnedAbout()
        {
            Glycan.RegisterCustomMonosaccharide("HexA", 'U', (int)Math.Round(176.03209 * 1E5), null);

            Assert.That(GlycanDatabase.CoreWarning("(U(N))", true), Is.Not.Null);
        }
    }
}
