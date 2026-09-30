using EngineLayer;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;

namespace Test
{
    /// <summary>
    /// A custom glycan that does not start from its core is warned about when it is added, not refused:
    /// nearly every O-glycan starts from GalNAc and every N-glycan from the GlcNAc-GlcNAc chitobiose core,
    /// but O-mannose, O-fucose, O-glucose and O-xylose glycans are real, one of the shipped O-glycan
    /// databases holds four of them, and an endoglycosidase leaves a single GlcNAc on the asparagine.
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
        /// order, so it is judged by its HexNAc count. An O-glycan needs one HexNAc to start from; an
        /// N-glycan needs the two of the chitobiose core -- in a structure, a HexNAc root with a HexNAc
        /// branch on it.
        /// </summary>
        [TestCase("(H)", true, true)]
        [TestCase("(H(N))", true, true)]
        [TestCase("(F(H))", true, true)]
        [TestCase("Hex(1)", true, true)]
        [TestCase("Hex(2)Fuc(1)", true, true)]
        [TestCase("(X(H(H)))", true, true)]
        [TestCase("Xylose(1)Hex(2)", true, true)]
        [TestCase("(H(N(N)))", true, true)]
        [TestCase("(N)", false, true)]
        [TestCase("(N(H(A)))", false, true)]
        [TestCase("(N(F))", false, true)]
        [TestCase("(N(H(N)))", false, true)]
        [TestCase("HexNAc(1)", false, true)]
        [TestCase("Hex(1)HexNAc(1)", false, true)]
        [TestCase("HexNAc(1)Fuc(1)", false, true)]
        [TestCase("(N(N))", false, false)]
        [TestCase("(N(N(H(H)(H))))", false, false)]
        [TestCase("(N(F)(N(H)))", false, false)]
        [TestCase("HexNAc(2)", false, false)]
        [TestCase("HexNAc(2)Hex(5)", false, false)]
        public static void AGlycanThatDoesNotStartFromItsCoreIsWarnedAbout(string glycan, bool warnedAsOGlycan, bool warnedAsNGlycan)
        {
            foreach (var (isOGlycan, warned) in new[] { (true, warnedAsOGlycan), (false, warnedAsNGlycan) })
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
            Assert.That(GlycanDatabase.CoreWarning("(H)", false), Does.Contain("GlcNAc").And.Contain("chitobiose"));
        }

        /// <summary>
        /// A single-HexNAc N-glycan is warned about for the second HexNAc it lacks, and the message says why
        /// the user might still mean it: an endoglycosidase leaves just one GlcNAc on the asparagine.
        /// </summary>
        [TestCase("(N(H(H)))", "does not begin with two HexNAc")]
        [TestCase("HexNAc(1)Hex(3)", "fewer than two HexNAc")]
        public static void AnNGlycanWithOneHexNAcIsWarnedAboutTheSecond(string glycan, string says)
        {
            string warning = GlycanDatabase.CoreWarning(glycan, false);

            Assert.That(warning, Does.Contain(says).And.Contain("endoglycosidase"));
        }

        /// <summary>
        /// O-xylose is one of the exceptions too: every proteoglycan attaches through it (Ser-Xyl-Gal-Gal-GlcA).
        /// A GAG linker still warns, because it does not start from HexNAc, but the message must not tell the
        /// user it is outside the known exceptions.
        /// </summary>
        [Test]
        public static void TheOGlycanWarningListsOXyloseAmongTheExceptions()
        {
            Assert.That(GlycanDatabase.CoreWarning("(X(H(H)))", true), Does.Contain("O-xylose"));
        }

        /// <summary>
        /// A warning, not a refusal: with no one to ask, the glycan is written.
        /// </summary>
        [TestCase("(H(H))", true)]
        [TestCase("Hex(2)", false)]
        [TestCase("HexNAc(1)Hex(3)", false)]
        public static void AGlycanThatDoesNotStartFromItsCoreIsStillAdded(string glycan, bool isOGlycan)
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
            Assert.That(GlycanDatabase.CoreWarning(text, false), Is.Null);
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

        // ---------------------------------------------------------------------------------------------
        // Asking before the write: PersistCustomGlycan's confirm
        // ---------------------------------------------------------------------------------------------

        /// <summary>
        /// An entry that is going to be refused anyway is refused without being asked about first: the
        /// question comes after every check that can refuse it, so the user is never asked to confirm a
        /// glycan and then told it cannot be added.
        /// </summary>
        [TestCase("Hex(1)", "Hex(1)", false, "already contains")]           // duplicate
        [TestCase("(H)", "Hex(1)", true, "format")]                         // structure into a composition file
        [TestCase("Hex(1)", "(N(H))", true, "format")]                      // composition into a structure file
        [TestCase("Hex(0)", null, false, "no monosaccharides")]
        [TestCase("()", null, true, "no monosaccharides")]
        [TestCase("Hex(41)", null, false, "at most")]
        public static void AnEntryThatIsRefusedIsNotAskedAbout(string entry, string alreadyInFile, bool isOGlycan, string refusal)
        {
            string path = Path_("custom.gdb");
            if (alreadyInFile != null)
            {
                File.WriteAllLines(path, new[] { alreadyInFile });
            }
            var asked = new List<string>();

            var ex = Assert.Throws<MetaMorpheusException>(() => GlycanDatabase.PersistCustomGlycan(entry, path, isOGlycan,
                question => { asked.Add(question); return true; }));

            Assert.That(ex.Message, Does.Contain(refusal));
            Assert.That(asked, Is.Empty);
        }

        /// <summary>
        /// A glycan that will be accepted and does not start from its core is asked about, with the warning
        /// text; yes writes it.
        /// </summary>
        [Test]
        public static void AGlycanWithoutItsCoreIsAskedAboutAndWrittenOnYes()
        {
            string path = Path_("custom.gdb");
            var asked = new List<string>();

            GlycanDatabase.PersistCustomGlycan("Hex(2)", path, true, question => { asked.Add(question); return true; });

            Assert.That(asked, Is.EqualTo(new[] { GlycanDatabase.CoreWarning("Hex(2)", true) }));
            Assert.That(File.ReadAllLines(path), Is.EqualTo(new[] { "Hex(2)" }));
        }

        /// <summary>No writes nothing: not into a new file, and not onto an existing one.</summary>
        [Test]
        public static void AGlycanWithoutItsCoreIsNotWrittenOnNo()
        {
            string newPath = Path_("new.gdb");
            string existingPath = Path_("existing.gdb");
            File.WriteAllLines(existingPath, new[] { "# mine", "HexNAc(1)Hex(1)" });
            string before = File.ReadAllText(existingPath);

            Assert.That(GlycanDatabase.PersistCustomGlycan("Hex(2)", newPath, true, _ => false), Is.Null);
            Assert.That(GlycanDatabase.PersistCustomGlycan("Hex(2)", existingPath, true, _ => false), Is.Null);

            Assert.That(File.Exists(newPath), Is.False);
            Assert.That(File.ReadAllText(existingPath), Is.EqualTo(before));
        }

        /// <summary>A glycan that starts from its core has nothing to ask about.</summary>
        [TestCase("(N(H(A)))", true)]
        [TestCase("HexNAc(2)Hex(5)", false)]
        public static void AGlycanThatStartsFromItsCoreIsNotAskedAbout(string glycan, bool isOGlycan)
        {
            string path = Path_("custom.gdb");
            bool asked = false;

            GlycanDatabase.PersistCustomGlycan(glycan, path, isOGlycan, _ => asked = true);

            Assert.That(asked, Is.False);
            Assert.That(File.ReadAllLines(path), Is.EqualTo(new[] { glycan }));
        }
    }
}
