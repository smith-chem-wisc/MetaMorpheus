using EngineLayer;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;

namespace Test
{
    /// <summary>
    /// The custom O-glycan database: the template a user is handed, the loader tolerance that lets a
    /// documented database be read at all, and the validation that stands between what they type and
    /// what the search will later have to parse.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // registers custom monosaccharides, which is process-wide state
    public static class CustomGlycanDatabaseTests
    {
        private static string _dir;

        [SetUp]
        public static void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "CustomGlycanDatabaseTests", Guid.NewGuid().ToString("N"));
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

        private static string EmbeddedTemplate() => CustomDataFile.EmbeddedText(
            typeof(GlobalVariables).Assembly, "EngineLayer.Glycan_Mods.OGlycan_Custom.gdb");

        // ---------------------------------------------------------------------------------------------
        // The template
        // ---------------------------------------------------------------------------------------------

        /// <summary>
        /// The template is embedded under the name the seeding code asks for. If the EmbeddedResource entry
        /// in EngineLayer.csproj is dropped or the file is renamed, this is what says so.
        /// </summary>
        [Test]
        public static void TheOGlycanTemplateIsEmbeddedUnderTheExpectedName()
        {
            var names = typeof(GlobalVariables).Assembly.GetManifestResourceNames();
            Assert.That(names, Does.Contain("EngineLayer.Glycan_Mods.OGlycan_Custom.gdb"),
                "the <EmbeddedResource Include=\"Glycan_Mods\\OGlycan_Custom.gdb\" /> entry in EngineLayer.csproj is missing");
        }

        /// <summary>
        /// Every line of the shipped template is a comment or blank -- it documents the format and
        /// contributes no glycans. A data row slipped into it would be silently searched by every user who
        /// never opened the file.
        /// </summary>
        [Test]
        public static void TheOGlycanTemplateHasNoDataRows()
        {
            string[] lines = EmbeddedTemplate().Split(new[] { "\r\n", "\n" }, StringSplitOptions.None);

            for (int i = 0; i < lines.Length; i++)
            {
                string line = lines[i];
                Assert.That(string.IsNullOrWhiteSpace(line) || line.TrimStart().StartsWith("#"), Is.True,
                    $"line {i + 1} of the embedded O-glycan template is a data row: \"{line}\"");
            }
        }

        /// <summary>
        /// A freshly seeded database loads, and loads to nothing. This is the case that runs on every user's
        /// first launch: GlobalVariables.LoadGlycans reads every database in the list eagerly, so a template
        /// the loaders could not read would throw inside SetUpGlobalVariables, before any window opened.
        /// </summary>
        [Test]
        public static void TheSeededTemplateLoadsToZeroGlycansWithoutThrowing()
        {
            string path = Path_("OGlycan_Custom.gdb");
            CustomDataFile.EnsureExists(path, EmbeddedTemplate, "custom O-glycan database");

            Assert.That(File.Exists(path), Is.True);
            Assert.That(GlycanDatabase.LoadGlycan(path, false, true).ToList(), Is.Empty);
        }

        // ---------------------------------------------------------------------------------------------
        // Loader tolerance
        // ---------------------------------------------------------------------------------------------

        /// <summary>
        /// A composition database that opens with a banner is still read as a composition database. The
        /// format is sniffed from the first line, so before this the banner decided the format and the whole
        /// file went to the structure parser.
        /// </summary>
        [Test]
        public static void ACommentBannerDoesNotChangeWhichParserReadsTheFile()
        {
            string path = Path_("commented_composition.gdb");
            File.WriteAllLines(path, new[]
            {
                "# my N-glycans",
                "#",
                "",
                "HexNAc(2)Hex(5)",
                "HexNAc(2)Hex(3)Fuc(1)",
            });

            var glycans = GlycanDatabase.LoadGlycan(path, false, false).ToList();

            // two motifs (Nxs, Nxt) per composition
            Assert.That(glycans.Count, Is.EqualTo(4));
            Assert.That(glycans.Select(g => Glycan.GetKindString(g.Kind)).Distinct().Count(), Is.EqualTo(2));
        }

        /// <summary>
        /// A comment that quotes a composition -- which the template does, to show the format -- is a comment.
        /// LoadKindGlycan used to skip only lines that lacked "Hex", so this line reached String2Kind and
        /// threw a KeyNotFoundException naming neither the file nor the line.
        /// </summary>
        [Test]
        public static void ACommentThatQuotesACompositionIsNotReadAsAGlycan()
        {
            string path = Path_("comment_quotes_composition.gdb");
            File.WriteAllLines(path, new[]
            {
                "HexNAc(2)Hex(5)",
                "# for example HexNAc(2)Hex(3)Fuc(1)",
            });

            var glycans = GlycanDatabase.LoadGlycan(path, false, false).ToList();

            Assert.That(glycans.Count, Is.EqualTo(2), "the commented composition was read as a glycan");
        }

        /// <summary>
        /// The structure parser skips comments and blanks too. Before this it ran ValidateStructureLine over
        /// the '#' and threw.
        /// </summary>
        [Test]
        public static void TheStructureParserSkipsCommentsAndBlankLines()
        {
            string path = Path_("commented_structure.gdb");
            File.WriteAllLines(path, new[]
            {
                "# my O-glycans",
                "",
                "(N)",
                "   # indented comment",
                "(N(H))",
            });

            var glycans = GlycanDatabase.LoadGlycan(path, false, true).ToList();

            Assert.That(glycans.Count, Is.EqualTo(4)); // two motifs (S, T) per structure
        }

        /// <summary>
        /// A structure with a character no monosaccharide claims is refused by name, naming the file and the
        /// line -- it used to be a bare FormatException that said neither.
        /// </summary>
        [Test]
        public static void AnUnknownMonosaccharideCodeNamesTheFileAndTheLine()
        {
            string path = Path_("bad_structure.gdb");
            File.WriteAllLines(path, new[] { "# banner", "(N(H))", "(N(Q))" });

            var ex = Assert.Throws<MetaMorpheusException>(() => GlycanDatabase.LoadGlycan(path, false, true).ToList());

            Assert.Multiple(() =>
            {
                Assert.That(ex.Message, Does.Contain("bad_structure.gdb"));
                Assert.That(ex.Message, Does.Contain("line 3"));
                Assert.That(ex.Message, Does.Contain("'Q'"));
            });
        }

        // ---------------------------------------------------------------------------------------------
        // Adding a glycan
        // ---------------------------------------------------------------------------------------------

        [Test]
        public static void AStructureIsAppendedAndReadsBack()
        {
            string path = Path_("OGlycan_Custom.gdb");
            CustomDataFile.EnsureExists(path, EmbeddedTemplate, "custom O-glycan database");

            GlycanDatabase.PersistCustomGlycan("(N(H(A)))", path, true);

            var glycans = GlycanDatabase.LoadGlycan(path, false, true).ToList();
            Assert.That(glycans.Count, Is.EqualTo(2)); // S and T
            Assert.That(File.ReadAllText(path), Does.Contain("(N(H(A)))"));
        }

        [Test]
        public static void ACompositionIsAppendedAndReadsBack()
        {
            string path = Path_("composition.gdb");
            GlycanDatabase.PersistCustomGlycan("HexNAc(2)Hex(5)", path, false);

            var glycans = GlycanDatabase.LoadGlycan(path, false, false).ToList();
            Assert.That(glycans.Count, Is.EqualTo(2)); // Nxs and Nxt
        }

        /// <summary>
        /// A database is read entirely as one format, so an entry in the other one would change how every
        /// line already in the file is interpreted. It is refused rather than mixed in.
        /// </summary>
        [Test]
        public static void AnEntryInTheOtherFormatIsRefused()
        {
            string path = Path_("structure_db.gdb");
            File.WriteAllLines(path, new[] { "# banner", "(N(H))" });

            var ex = Assert.Throws<MetaMorpheusException>(
                () => GlycanDatabase.PersistCustomGlycan("HexNAc(2)Hex(5)", path, true));

            Assert.That(ex.Message, Does.Contain("structure"));
            Assert.That(ex.Message, Does.Contain("composition"));
        }

        /// <summary>
        /// The format of a database that holds only a banner is undecided, so the first glycan added to a
        /// freshly seeded file may be in either format.
        /// </summary>
        [Test]
        public static void AFreshlySeededDatabaseAcceptsEitherFormat()
        {
            string structurePath = Path_("seeded_structure.gdb");
            string compositionPath = Path_("seeded_composition.gdb");
            CustomDataFile.EnsureExists(structurePath, EmbeddedTemplate, "custom O-glycan database");
            CustomDataFile.EnsureExists(compositionPath, EmbeddedTemplate, "custom O-glycan database");

            Assert.DoesNotThrow(() => GlycanDatabase.PersistCustomGlycan("(N(H))", structurePath, true));
            Assert.DoesNotThrow(() => GlycanDatabase.PersistCustomGlycan("HexNAc(2)Hex(5)", compositionPath, true));
        }

        [Test]
        public static void TheSameGlycanIsNotAddedTwice()
        {
            string path = Path_("dupe.gdb");
            GlycanDatabase.PersistCustomGlycan("(N(H))", path, true);

            var ex = Assert.Throws<MetaMorpheusException>(() => GlycanDatabase.PersistCustomGlycan("(N(H))", path, true));
            Assert.That(ex.Message, Does.Contain("already contains"));
        }

        /// <summary>
        /// A hand-edited .gdb very often has no final newline -- the shipped NGlycan_ForNoSearch.gdb test
        /// fixture does not. Without the guard the new glycan is glued onto the end of the last one and both
        /// are lost.
        /// </summary>
        [Test]
        public static void AppendingToAFileWithNoTrailingNewlineDoesNotCorruptTheLastEntry()
        {
            string path = Path_("no_trailing_newline.gdb");
            File.WriteAllText(path, "(N)\r\n(N(H))"); // deliberately unterminated

            GlycanDatabase.PersistCustomGlycan("(N(A))", path, true);

            var entries = File.ReadAllLines(path).Where(l => l.Trim().Length > 0).ToList();
            Assert.That(entries, Is.EqualTo(new List<string> { "(N)", "(N(H))", "(N(A))" }));
        }

        [Test]
        public static void AnUnknownMonosaccharideNameIsRefusedBeforeItReachesTheFile()
        {
            string path = Path_("unknown_name.gdb");

            var ex = Assert.Throws<MetaMorpheusException>(
                () => GlycanDatabase.PersistCustomGlycan("HexA(1)Hex(1)", path, false));

            Assert.Multiple(() =>
            {
                Assert.That(ex.Message, Does.Contain("HexA"));
                Assert.That(ex.Message, Does.Contain("MonosaccharidesCustom.tsv"));
                Assert.That(File.Exists(path), Is.False, "the file was written despite the entry being refused");
            });
        }

        /// <summary>
        /// Once a monosaccharide is declared, it is legal in both formats -- which is the promise the
        /// MonosaccharidesCustom.tsv banner already makes about this file.
        /// </summary>
        [Test]
        public static void ADeclaredCustomMonosaccharideIsAcceptedInBothFormats()
        {
            Glycan.ResetCustomMonosaccharides();
            Glycan.RegisterCustomMonosaccharide("HexA", 'U', (int)Math.Round(176.03209 * 1E5), null);

            string compositionPath = Path_("custom_sugar_composition.gdb");
            string structurePath = Path_("custom_sugar_structure.gdb");

            Assert.DoesNotThrow(() => GlycanDatabase.PersistCustomGlycan("HexNAc(1)HexA(1)", compositionPath, true));
            Assert.DoesNotThrow(() => GlycanDatabase.PersistCustomGlycan("(N(U))", structurePath, true));

            Assert.That(GlycanDatabase.LoadGlycan(structurePath, false, true).ToList(), Is.Not.Empty);
        }

        [Test]
        public static void UnbalancedParenthesesAreRefused()
        {
            string path = Path_("unbalanced.gdb");

            Assert.Multiple(() =>
            {
                Assert.That(Assert.Throws<MetaMorpheusException>(
                    () => GlycanDatabase.PersistCustomGlycan("(N(H)", path, true)).Message, Does.Contain("balance"));
                Assert.That(Assert.Throws<MetaMorpheusException>(
                    () => GlycanDatabase.PersistCustomGlycan("(N(H)))", path, true)).Message, Does.Contain("balance"));
            });
        }

        [Test]
        public static void TextThatIsNeitherFormatIsRefused()
        {
            string path = Path_("nonsense.gdb");

            var ex = Assert.Throws<MetaMorpheusException>(
                () => GlycanDatabase.PersistCustomGlycan("a nice glycan please", path, true));
            Assert.That(ex.Message, Does.Contain("neither a structure"));
        }

        [Test]
        public static void ACommentIsRefusedRatherThanSilentlyIgnored()
        {
            string path = Path_("comment_entry.gdb");

            var ex = Assert.Throws<MetaMorpheusException>(
                () => GlycanDatabase.PersistCustomGlycan("# (N(H))", path, true));
            Assert.That(ex.Message, Does.Contain("comment"));
        }

        [Test]
        public static void ARepeatedMonosaccharideInOneCompositionIsRefused()
        {
            string path = Path_("repeated.gdb");

            var ex = Assert.Throws<MetaMorpheusException>(
                () => GlycanDatabase.PersistCustomGlycan("Hex(1)Hex(2)", path, false));
            Assert.That(ex.Message, Does.Contain("more than once"));
        }

        // ---------------------------------------------------------------------------------------------
        // Startup wiring, and the guard on an empty database
        // ---------------------------------------------------------------------------------------------

        /// <summary>
        /// The custom database sits at the DataDir root, so unlike the shipped ones it is not found by the
        /// directory sweep. If it is not added to the list by name it exists but is never offered.
        /// </summary>
        [Test]
        public static void StartupOffersTheCustomOGlycanDatabase()
        {
            GlobalVariables.SetUpGlobalVariables();

            Assert.Multiple(() =>
            {
                Assert.That(File.Exists(GlobalVariables.CustomOGlycanDatabasePath), Is.True);
                Assert.That(GlobalVariables.OGlycanDatabasePaths, Does.Contain(GlobalVariables.CustomOGlycanDatabasePath));
            });
        }

        /// <summary>
        /// Searching a database with no glycans in it used to build an empty box array without complaint and
        /// then throw "Sequence contains no elements" from inside the parallel search loop -- or, if no scan
        /// reached that branch, quietly return nothing. It matters more now that a user can be handed an
        /// empty database: a freshly seeded one is all banner until they add a glycan.
        /// </summary>
        [Test]
        public static void SearchingAnEmptyGlycanDatabaseIsRefusedByName()
        {
            string path = Path_("OGlycan_Empty.gdb");
            CustomDataFile.EnsureExists(path, EmbeddedTemplate, "custom O-glycan database");
            GlobalVariables.OGlycanDatabasePaths.Add(path);

            try
            {
                var ex = Assert.Throws<MetaMorpheusException>(() => new GlycoSearchEngine(
                    new List<GlycoSpectralMatch>[0], new Ms2ScanWithSpecificMass[0],
                    new List<PeptideWithSetModifications>(), null, null, 0,
                    new CommonParameters(), null, "OGlycan_Empty.gdb", "NGlycan.gdb",
                    GlycoSearchType.OGlycanSearch, 30, 3, false, null));

                Assert.Multiple(() =>
                {
                    Assert.That(ex.Message, Does.Contain("OGlycan_Empty.gdb"));
                    Assert.That(ex.Message, Does.Contain("no glycans"));
                });
            }
            finally
            {
                GlobalVariables.OGlycanDatabasePaths.Remove(path);
            }
        }

        /// <summary>
        /// A database that is selected but no longer in the folder threw the same unhelpful
        /// "Sequence contains no elements" from the path lookup. It now says which one is missing.
        /// </summary>
        [Test]
        public static void SelectingAGlycanDatabaseThatIsNotThereIsRefusedByName()
        {
            var ex = Assert.Throws<MetaMorpheusException>(() => new GlycoSearchEngine(
                new List<GlycoSpectralMatch>[0], new Ms2ScanWithSpecificMass[0],
                new List<PeptideWithSetModifications>(), null, null, 0,
                new CommonParameters(), null, "NoSuchDatabase.gdb", "NGlycan.gdb",
                GlycoSearchType.OGlycanSearch, 30, 3, false, null));

            Assert.That(ex.Message, Does.Contain("NoSuchDatabase.gdb"));
        }
    }
}
