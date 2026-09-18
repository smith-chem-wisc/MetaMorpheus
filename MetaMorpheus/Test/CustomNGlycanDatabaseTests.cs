using EngineLayer;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    /// <summary>
    /// The custom N-glycan database. The loader tolerance and the validate-and-append it rests on are
    /// covered by <see cref="CustomGlycanDatabaseTests"/> and shared; what is N-specific is the template,
    /// the seeding, and that a glycan added to it is read back on the N-glycan motifs.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // registers custom monosaccharides, which is process-wide state
    public static class CustomNGlycanDatabaseTests
    {
        private static string _dir;

        [SetUp]
        public static void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "CustomNGlycanDatabaseTests", Guid.NewGuid().ToString("N"));
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
            typeof(GlobalVariables).Assembly, "EngineLayer.Glycan_Mods.NGlycan_Custom.gdb");

        /// <summary>
        /// The template is embedded under the name the seeding code asks for. If the EmbeddedResource entry
        /// in EngineLayer.csproj is dropped or the file is renamed, this is what says so.
        /// </summary>
        [Test]
        public static void TheNGlycanTemplateIsEmbeddedUnderTheExpectedName()
        {
            var names = typeof(GlobalVariables).Assembly.GetManifestResourceNames();
            Assert.That(names, Does.Contain("EngineLayer.Glycan_Mods.NGlycan_Custom.gdb"),
                "the <EmbeddedResource Include=\"Glycan_Mods\\NGlycan_Custom.gdb\" /> entry in EngineLayer.csproj is missing");
        }

        /// <summary>
        /// Every line of the shipped template is a comment or blank -- it documents the format and
        /// contributes no glycans. A data row slipped into it would be silently searched by every user who
        /// never opened the file.
        /// </summary>
        [Test]
        public static void TheNGlycanTemplateHasNoDataRows()
        {
            string[] lines = EmbeddedTemplate().Split(new[] { "\r\n", "\n" }, StringSplitOptions.None);

            for (int i = 0; i < lines.Length; i++)
            {
                string line = lines[i];
                Assert.That(string.IsNullOrWhiteSpace(line) || line.TrimStart().StartsWith("#"), Is.True,
                    $"line {i + 1} of the embedded N-glycan template is a data row: \"{line}\"");
            }
        }

        /// <summary>
        /// A freshly seeded database loads, and loads to nothing. GlobalVariables.LoadGlycans reads every
        /// database in the list eagerly, so a template the loaders could not read would throw inside
        /// SetUpGlobalVariables, before any window opened.
        /// </summary>
        [Test]
        public static void TheSeededTemplateLoadsToZeroGlycansWithoutThrowing()
        {
            string path = Path_("NGlycan_Custom.gdb");
            CustomDataFile.EnsureExists(path, EmbeddedTemplate, "custom N-glycan database");

            Assert.That(File.Exists(path), Is.True);
            Assert.That(GlycanDatabase.LoadGlycan(path, false, false).ToList(), Is.Empty);
        }

        /// <summary>
        /// A composition added to the seeded template reads back on both N-glycan sequons, with the
        /// composition it was written with.
        /// </summary>
        [Test]
        public static void ACompositionAddedToTheSeededTemplateReadsBackOnBothSequons()
        {
            string path = Path_("NGlycan_Custom.gdb");
            CustomDataFile.EnsureExists(path, EmbeddedTemplate, "custom N-glycan database");

            GlycanDatabase.PersistCustomGlycan("HexNAc(4)Hex(5)Fuc(1)NeuAc(1)", path, false);

            var glycans = GlycanDatabase.LoadGlycan(path, false, false).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(glycans.Count, Is.EqualTo(2));
                Assert.That(glycans.Select(g => g.Target.ToString()), Is.EquivalentTo(new[] { "Nxs", "Nxt" }));
                Assert.That(glycans.Select(g => Glycan.GetKindString(g.Kind)).Distinct().Single(),
                    Is.EqualTo(Glycan.GetKindString(GlycanDatabase.String2Kind("HexNAc(4)Hex(5)Fuc(1)NeuAc(1)"))));
            });
        }

        /// <summary>
        /// The N-glycan template invites composition format, but structure format is legal in it too -- the
        /// loader picks the parser from the first data line, not from which database it is.
        /// </summary>
        [Test]
        public static void AStructureIsAlsoLegalInAnNGlycanDatabase()
        {
            string path = Path_("NGlycan_structure.gdb");
            CustomDataFile.EnsureExists(path, EmbeddedTemplate, "custom N-glycan database");

            GlycanDatabase.PersistCustomGlycan("(N(N(H(H)(H))))", path, false);

            var glycans = GlycanDatabase.LoadGlycan(path, false, false).ToList();
            Assert.That(glycans, Is.Not.Empty);
        }

        /// <summary>
        /// The custom database sits at the DataDir root, so unlike the shipped ones it is not found by the
        /// directory sweep. If it is not added to the list by name it exists but is never offered.
        /// </summary>
        [Test]
        public static void StartupOffersTheCustomNGlycanDatabase()
        {
            GlobalVariables.SetUpGlobalVariables();

            Assert.Multiple(() =>
            {
                Assert.That(File.Exists(GlobalVariables.CustomNGlycanDatabasePath), Is.True);
                Assert.That(GlobalVariables.NGlycanDatabasePaths, Does.Contain(GlobalVariables.CustomNGlycanDatabasePath));
            });
        }

        /// <summary>
        /// The empty-database guard covers the N-glycan search too -- that branch loads its database through
        /// a separate call, and a freshly seeded N-glycan database is exactly the empty case.
        /// </summary>
        [Test]
        public static void SearchingAnEmptyNGlycanDatabaseIsRefusedByName()
        {
            string path = Path_("NGlycan_Empty.gdb");
            CustomDataFile.EnsureExists(path, EmbeddedTemplate, "custom N-glycan database");
            GlobalVariables.NGlycanDatabasePaths.Add(path);

            try
            {
                var ex = Assert.Throws<MetaMorpheusException>(() => new GlycoSearchEngine(
                    new List<GlycoSpectralMatch>[0], new Ms2ScanWithSpecificMass[0],
                    new List<PeptideWithSetModifications>(), null, null, 0,
                    new CommonParameters(), null, "OGlycan.gdb", "NGlycan_Empty.gdb",
                    GlycoSearchType.NGlycanSearch, 30, 3, false, null));

                Assert.Multiple(() =>
                {
                    Assert.That(ex.Message, Does.Contain("NGlycan_Empty.gdb"));
                    Assert.That(ex.Message, Does.Contain("no glycans"));
                });
            }
            finally
            {
                GlobalVariables.NGlycanDatabasePaths.Remove(path);
            }
        }
    }
}
