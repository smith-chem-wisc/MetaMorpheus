using Nett;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The individual-glycan selection has to survive a task being saved and reopened, and a task
    /// saved before the selection existed has to keep behaving as it did.
    /// </summary>
    [TestFixture]
    public class GlycoSearchParametersTomlTests
    {
        private string _tomlPath;

        [SetUp]
        public void SetUp()
        {
            _tomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoSelectedGlycans.toml");
        }

        [TearDown]
        public void TearDown()
        {
            if (File.Exists(_tomlPath))
            {
                File.Delete(_tomlPath);
            }
        }

        private GlycoSearchTask WriteAndReadBack(GlycoSearchTask task)
        {
            Toml.WriteFile(task, _tomlPath, MetaMorpheusTask.tomlConfig);
            return Toml.ReadFile<GlycoSearchTask>(_tomlPath, MetaMorpheusTask.tomlConfig);
        }

        [Test]
        public void SelectedGlycans_DefaultsToEmptyMeaningTheWholeDatabase()
        {
            var task = new GlycoSearchTask();

            Assert.That(task._glycoSearchParameters.SelectedGlycans, Is.Not.Null);
            Assert.That(task._glycoSearchParameters.SelectedGlycans, Is.Empty);
        }

        [Test]
        public void SelectedGlycans_RoundTripsThroughToml()
        {
            var task = new GlycoSearchTask();
            task._glycoSearchParameters.SelectedGlycans = new List<(string, string)>
            {
                ("OGlycan.gdb", "N1 on S"),
                ("OGlycan.gdb", "H1N1 on T"),
                ("NGlycan.gdb", "H5N2 on Nxs"),
            };

            var loaded = WriteAndReadBack(task);

            Assert.That(loaded._glycoSearchParameters.SelectedGlycans,
                Is.EqualTo(task._glycoSearchParameters.SelectedGlycans));
        }

        [Test]
        public void SelectedGlycans_RoundTripsASelectionSpanningTwoDatabases()
        {
            var task = new GlycoSearchTask();
            task._glycoSearchParameters.SelectedGlycans = new List<(string, string)>
            {
                ("OGlycan.gdb", "N1 on S"),
                ("NGlycan.gdb", "H5N2 on Nxs"),
            };

            var loaded = WriteAndReadBack(task);

            // The database name is half the key, and it is what the engine matches on to decide which
            // database a chosen glycan belongs to -- so it has to survive, not just the glycan id.
            Assert.That(loaded._glycoSearchParameters.SelectedGlycans.Select(s => s.Item1),
                Is.EquivalentTo(new[] { "OGlycan.gdb", "NGlycan.gdb" }));
            Assert.That(loaded._glycoSearchParameters.SelectedGlycans.Select(s => s.Item2),
                Is.EquivalentTo(new[] { "N1 on S", "H5N2 on Nxs" }));
        }

        [Test]
        public void SelectedGlycans_RoundTripsAnEmptySelection()
        {
            var task = new GlycoSearchTask();
            task._glycoSearchParameters.SelectedGlycans = new List<(string, string)>();

            var loaded = WriteAndReadBack(task);

            Assert.That(loaded._glycoSearchParameters.SelectedGlycans, Is.Empty);
        }

        [Test]
        public void SelectedGlycans_KeepsTheDatabaseFileNamesAlongsideIt()
        {
            // The tree narrows a database; it does not replace choosing one. Both have to persist.
            var task = new GlycoSearchTask();
            task._glycoSearchParameters.OGlycanDatabasefile = "OGlycan_withIsobaric.gdb";
            task._glycoSearchParameters.SelectedGlycans = new List<(string, string)>
            {
                ("OGlycan_withIsobaric.gdb", "N1 on S"),
            };

            var loaded = WriteAndReadBack(task);

            Assert.That(loaded._glycoSearchParameters.OGlycanDatabasefile, Is.EqualTo("OGlycan_withIsobaric.gdb"));
            Assert.That(loaded._glycoSearchParameters.SelectedGlycans.Single().Item1, Is.EqualTo("OGlycan_withIsobaric.gdb"));
        }

        [Test]
        public void ATaskFileWrittenBeforeSelectedGlycansExistedStillLoads()
        {
            // The compatibility guarantee, tested the way it will actually be met: every GlycoSearch
            // toml already on a user's disk has no SelectedGlycans line at all. It must load, and it
            // must load as "nothing selected", which the engine reads as the whole database.
            var task = new GlycoSearchTask();
            task._glycoSearchParameters.SelectedGlycans = new List<(string, string)>
            {
                ("OGlycan.gdb", "N1 on S"),
            };
            Toml.WriteFile(task, _tomlPath, MetaMorpheusTask.tomlConfig);

            var withoutTheLine = File.ReadAllLines(_tomlPath)
                .Where(line => !line.TrimStart().StartsWith("SelectedGlycans"))
                .ToArray();
            Assert.That(withoutTheLine.Length, Is.LessThan(File.ReadAllLines(_tomlPath).Length),
                "the property should have been written, otherwise this test proves nothing");
            File.WriteAllLines(_tomlPath, withoutTheLine);

            var loaded = Toml.ReadFile<GlycoSearchTask>(_tomlPath, MetaMorpheusTask.tomlConfig);

            Assert.That(loaded._glycoSearchParameters.SelectedGlycans, Is.Not.Null);
            Assert.That(loaded._glycoSearchParameters.SelectedGlycans, Is.Empty);
        }
    }
}
