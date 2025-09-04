using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Media;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using NUnit.Framework;
using Omics;
using Proteomics;
using Readers;

namespace Test.MetaDraw
{
    public class DummyMetaDrawLogic : MetaDrawLogic
    {
        public DummyMetaDrawLogic()
        {
            AllSpectralMatches = new List<Readers.SpectrumMatchFromTsv>();
        }
    }

    [TestFixture]
    public class BioPolymerTabViewModelTests
    {
        private class DummyBioPolymer(string accession, string name, string baseSeq)
            : Protein(baseSeq, accession, name: name);

        [Test]
        public void Constructor_InitializesProperties()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic, "C:\\Export");
            Assert.That(vm.AllGroups, Is.Not.Null);
            Assert.That(vm.FilteredGroups, Is.Not.Null);
            Assert.That(vm.CoverageMapViewModel, Is.Not.Null);
            Assert.That(vm.ExportDirectory, Is.EqualTo("C:\\Export"));
            Assert.That(vm.LoadDatabaseCommand, Is.Not.Null);
            Assert.That(vm.ResetDatabaseCommand, Is.Not.Null);
            Assert.That(vm.ExportImageCommand, Is.Not.Null);
        }

        [Test]
        public void IsDatabaseLoaded_PropertyChanged()
        {
            var vm = new BioPolymerTabViewModel(new DummyMetaDrawLogic());
            bool called = false;
            vm.PropertyChanged += (s, e) => { if (e.PropertyName == "IsDatabaseLoaded") called = true; };
            vm.IsDatabaseLoaded = true;
            Assert.That(vm.IsDatabaseLoaded, Is.True);
            Assert.That(called, Is.True);
        }

        [Test]
        public void DatabasePath_And_DatabaseName_PropertyChanged()
        {
            var vm = new BioPolymerTabViewModel(new DummyMetaDrawLogic());
            bool pathChanged = false, nameChanged = false;
            vm.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName == "DatabasePath") pathChanged = true;
                if (e.PropertyName == "DatabaseName") nameChanged = true;
            };
            vm.DatabasePath = "C:\\db.fasta";
            Assert.That(vm.DatabasePath, Is.EqualTo("C:\\db.fasta"));
            Assert.That(vm.DatabaseName, Is.Not.Null.And.Not.Empty);
            Assert.That(pathChanged, Is.True);
            Assert.That(nameChanged, Is.True);
        }

        [Test]
        public void SearchText_PropertyChanged_And_FilteredGroups()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            var group1 = new BioPolymerGroupViewModel("A1", "Alpha", "ABC", new BioPolymerCoverageResultModel[0]);
            var group2 = new BioPolymerGroupViewModel("B2", "Beta", "DEF", new BioPolymerCoverageResultModel[0]);
            vm.AllGroups.Add(group1);
            vm.AllGroups.Add(group2);

            bool changed = false;
            vm.PropertyChanged += (s, e) => { if (e.PropertyName == "SearchText") changed = true; };

            vm.SearchText = "Alpha";
            Assert.That(changed, Is.True);
            Assert.That(vm.FilteredGroups.Count, Is.EqualTo(1));
            Assert.That(vm.FilteredGroups[0], Is.EqualTo(group1));

            vm.SearchText = "B2";
            Assert.That(vm.FilteredGroups.Count, Is.EqualTo(1));
            Assert.That(vm.FilteredGroups[0], Is.EqualTo(group2));

            vm.SearchText = "";
            Assert.That(vm.FilteredGroups.Count, Is.EqualTo(2));
        }

        [Test]
        public void SelectedGroup_PropertyChanged_And_CoverageMapViewModel()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            var group = new BioPolymerGroupViewModel("A1", "Alpha", "ABC", new BioPolymerCoverageResultModel[0]);
            bool changed = false;
            vm.PropertyChanged += (s, e) => { if (e.PropertyName == "SelectedGroup") changed = true; };
            vm.SelectedGroup = group;
            Assert.That(vm.SelectedGroup, Is.EqualTo(group));
            Assert.That(vm.CoverageMapViewModel.Group, Is.EqualTo(group));
            Assert.That(changed, Is.True);
        }

        [Test]
        public void ResetDatabase_ClearsState()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            vm.DatabasePath = "C:\\db.fasta";
            vm.DatabaseName = "db";
            vm.AllGroups.Add(new BioPolymerGroupViewModel("A1", "Alpha", "ABC", new BioPolymerCoverageResultModel[0]));
            vm.FilteredGroups.Add(new BioPolymerGroupViewModel("B2", "Beta", "DEF", new BioPolymerCoverageResultModel[0]));

            vm.GetType().GetMethod("ResetDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)
                .Invoke(vm, null);

            Assert.That(vm.DatabasePath, Is.Empty);
            Assert.That(vm.DatabaseName, Is.Empty);
            Assert.That(vm.AllGroups.Count, Is.EqualTo(0));
            Assert.That(vm.FilteredGroups.Count, Is.EqualTo(0));
        }

        [Test]
        public void ExportDirectory_CreatesDirectoryIfNotExists()
        {
            var logic = new DummyMetaDrawLogic();
            var tempDir = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "BioPolymerTabViewModelTestDir");
            if (System.IO.Directory.Exists(tempDir))
                System.IO.Directory.Delete(tempDir, true);

            var vm = new BioPolymerTabViewModel(logic, tempDir);
            bool changed = false;
            vm.PropertyChanged += (s, e) => { if (e.PropertyName == "ExportDirectory") changed = true; };
            var dir = vm.ExportDirectory;
            Assert.That(System.IO.Directory.Exists(dir), Is.True);
            Assert.That(changed, Is.False); // Only fires on set, not get

            vm.ExportDirectory = tempDir + "2";
            Assert.That(vm.ExportDirectory, Is.EqualTo(tempDir + "2"));
            Assert.That(changed, Is.True);

            System.IO.Directory.Delete(tempDir, true);
            if (System.IO.Directory.Exists(tempDir + "2"))
                System.IO.Directory.Delete(tempDir + "2", true);
        }

        [Test]
        public void ProcessSpectralMatches_AddsGroups()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);

            // Setup _allBioPolymers with a dummy
            var dummyBioPolymer = new DummyBioPolymer("ACC", "Protein", "ABC");
            var allBioPolymersField = typeof(BioPolymerTabViewModel).GetField("_allBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            allBioPolymersField.SetValue(vm, new Dictionary<string, IBioPolymer> { { "ACC", dummyBioPolymer } });

            var match = new DummySpectralmatch();
            var matches = new List<Readers.SpectrumMatchFromTsv> { match };
            vm.ProcessSpectralMatches(matches);

            Assert.That(vm.AllGroups.Count, Is.EqualTo(1));
            Assert.That(vm.AllGroups[0].Accession, Is.EqualTo("ACC"));
        }

        [Test]
        public void LoadDatabase_DoesNothing_WhenPathIsNullOrEmpty()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            
            var allBioPolymersField = typeof(BioPolymerTabViewModel)
                .GetField("_allBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            allBioPolymersField.SetValue(vm, new Dictionary<string, IBioPolymer>());
            var method = typeof(BioPolymerTabViewModel)
                .GetMethod("LoadDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            method.Invoke(vm, null);
            Assert.That(vm.IsDatabaseLoaded, Is.False);
        }

        [Test]
        public void LoadDatabase_CatchesException_AndDoesNotThrow()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            // Set an invalid path to force exception
            vm.DatabasePath = "Z:\\this\\does\\not\\exist.fasta";
            var method = typeof(BioPolymerTabViewModel)
                .GetMethod("LoadDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            Assert.DoesNotThrow(() => method.Invoke(vm, null));
            Assert.That(vm.IsDatabaseLoaded, Is.False);
        }

        [Test]
        public void LoadDatabase_NonExistentDatabase_DoesNotSetLoadedToTrue()
        {
            var logic = new DummyMetaDrawLogic();
            logic.AllSpectralMatches = new List<SpectrumMatchFromTsv>
            {
                new DummySpectralmatch(accession: "ACC", baseSeq: "ABC"),
                new DummySpectralmatch(accession: "ACC2", baseSeq: "DEF")
            };
            var vm = new BioPolymerTabViewModel(logic);

            vm.DatabasePath = "C:\\db.fasta";
            var method = typeof(BioPolymerTabViewModel)
                .GetMethod("LoadDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            method.Invoke(vm, null);
            Assert.That(vm.IsDatabaseLoaded, Is.False);
            Assert.That(vm.AllGroups.Count, Is.EqualTo(0));
        }

        [Test]
        public void LoadDatabase_EmptyDatabase_DoesNotSetLoadedToTrue()
        {
            var logic = new DummyMetaDrawLogic();
            logic.AllSpectralMatches = new List<SpectrumMatchFromTsv>
            {
                new DummySpectralmatch(accession: "ACC", baseSeq: "ABC"),
                new DummySpectralmatch(accession: "ACC2", baseSeq: "DEF")
            };
            var vm = new BioPolymerTabViewModel(logic);

            var dbPath = Path.Combine(Path.GetTempPath(), "empty.fasta");
            File.WriteAllText(dbPath, string.Empty);

            vm.DatabasePath = dbPath;
            var method = typeof(BioPolymerTabViewModel)
                .GetMethod("LoadDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            method.Invoke(vm, null);
            Assert.That(vm.IsDatabaseLoaded, Is.False);
            Assert.That(vm.AllGroups.Count, Is.EqualTo(0));
        }

        [Test]
        public void LoadDatabase_RealDatabase_LoadsWithoutIssue()
        {
            // Prepare Files
            if (!EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpQValue, out var searchTestCase))
                Assert.Fail();

            var dbPath = searchTestCase.DatabaseList.First().FilePath;
            var smPath = Directory.GetFiles(searchTestCase.OutputDirectory, "*.psmtsv", SearchOption.AllDirectories).First();
            var matches = Readers.SpectrumMatchTsvReader.ReadTsv(smPath, out _);

            // Build View Model
            var logic = new DummyMetaDrawLogic();
            logic.AllSpectralMatches = matches;
            var vm = new BioPolymerTabViewModel(logic);

            // Test
            vm.DatabasePath = dbPath;
            var method = typeof(BioPolymerTabViewModel)
                .GetMethod("LoadDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            method.Invoke(vm, null);
            Assert.That(vm.IsDatabaseLoaded, Is.True);
        }

        [Test]
        public void ProcessSpectralMatches_HandlesSingleAccessionSingleLocalization()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            var allBioPolymersField = typeof(BioPolymerTabViewModel)
                .GetField("_allBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            allBioPolymersField.SetValue(vm, new Dictionary<string, IBioPolymer>
            {
                { "ACC", new DummyBioPolymer("ACC", "Protein", "ABC") }
            });
            var match = new DummySpectralmatch(accession: "ACC", baseSeq: "ABC", missedCleavages: "0", startAndEnd: "[1 to 3]");
            vm.ProcessSpectralMatches(new List<SpectrumMatchFromTsv> { match });
            Assert.That(vm.AllGroups.Count, Is.EqualTo(1));
            var group = vm.AllGroups[0];
            Assert.That(group.CoverageResults.Count, Is.EqualTo(1));
            Assert.That(group.CoverageResults[0].CoverageType, Is.EqualTo(BioPolymerCoverageType.Unique));
        }

        [Test]
        public void ProcessSpectralMatches_HandlesSingleAccessionSingleLocalization_MissedCleavage()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            var allBioPolymersField = typeof(BioPolymerTabViewModel)
                .GetField("_allBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            allBioPolymersField.SetValue(vm, new Dictionary<string, IBioPolymer>
            {
                { "ACC", new DummyBioPolymer("ACC", "Protein", "ABC") }
            });
            var match = new DummySpectralmatch(accession: "ACC", baseSeq: "ABC", missedCleavages: "1", startAndEnd: "[1 to 3]");
            vm.ProcessSpectralMatches(new List<SpectrumMatchFromTsv> { match });
            Assert.That(vm.AllGroups.Count, Is.EqualTo(1));
            var group = vm.AllGroups[0];
            Assert.That(group.CoverageResults[0].CoverageType, Is.EqualTo(BioPolymerCoverageType.UniqueMissedCleavage));
        }

        [Test]
        public void ProcessSpectralMatches_HandlesMultipleAccessions()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            var allBioPolymersField = typeof(BioPolymerTabViewModel)
                .GetField("_allBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            allBioPolymersField.SetValue(vm, new Dictionary<string, IBioPolymer>
            {
                { "ACC", new DummyBioPolymer("ACC", "Protein", "ABC") },
                { "ACC2", new DummyBioPolymer("ACC2", "Protein2", "ABC") }
            });
            var match = new DummySpectralmatch(accession: "ACC|ACC2", baseSeq: "ABC", missedCleavages: "0|0", startAndEnd: "[1 to 3]|[1 to 3]");
            vm.ProcessSpectralMatches(new List<SpectrumMatchFromTsv> { match });
            Assert.That(vm.AllGroups.Count, Is.EqualTo(2));
            Assert.That(vm.AllGroups.All(g => g.CoverageResults[0].CoverageType == BioPolymerCoverageType.Shared));
        }

        [Test]
        public void ProcessSpectralMatches_HandlesSingleAccessionMultipleLocalizations()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            var allBioPolymersField = typeof(BioPolymerTabViewModel)
                .GetField("_allBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            allBioPolymersField.SetValue(vm, new Dictionary<string, IBioPolymer>
            {
                { "ACC", new DummyBioPolymer("ACC", "Protein", "ABC") }
            });
            var match = new DummySpectralmatch(accession: "ACC", baseSeq: "ABC", missedCleavages: "0", startAndEnd: "[1 to 2]|[2 to 3]");
            vm.ProcessSpectralMatches(new List<SpectrumMatchFromTsv> { match });
            Assert.That(vm.AllGroups.Count, Is.EqualTo(1));
            var group = vm.AllGroups[0];
            Assert.That(group.CoverageResults.Count, Is.EqualTo(2));
            Assert.That(group.CoverageResults.All(r => r.CoverageType == BioPolymerCoverageType.TandemRepeat));
        }

        [Test]
        public void ProcessSpectralMatches_HandlesMultipleAccessionsSingleLocalization()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            var allBioPolymersField = typeof(BioPolymerTabViewModel)
                .GetField("_allBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            allBioPolymersField.SetValue(vm, new Dictionary<string, IBioPolymer>
            {
                { "ACC", new DummyBioPolymer("ACC", "Protein", "ABC") },
                { "ACC2", new DummyBioPolymer("ACC2", "Protein2", "ABC") }
            });
            var match = new DummySpectralmatch(accession: "ACC|ACC2", baseSeq: "ABC", missedCleavages: "0|0", startAndEnd: "[1 to 3]");
            vm.ProcessSpectralMatches(new List<SpectrumMatchFromTsv> { match });
            Assert.That(vm.AllGroups.Count, Is.EqualTo(2));
            Assert.That(vm.AllGroups.All(g => g.CoverageResults[0].CoverageType == BioPolymerCoverageType.Shared));
        }

        [Test]
        public void ExportImage_DoesNothing_IfCoverageDrawingIsNull()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);
            vm.CoverageMapViewModel.CoverageDrawing = null;
            var method = typeof(BioPolymerTabViewModel)
                .GetMethod("ExportImage", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            Assert.DoesNotThrow(() => method.Invoke(vm, null));
        }

        [Test]
        public void ExportImage_ExportsImage_WhenCoverageDrawingIsPresent()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic, Path.GetTempPath());
            // Setup a dummy DrawingImage
            var drawingVisual = new DrawingVisual();
            using (var dc = drawingVisual.RenderOpen())
            {
                dc.DrawRectangle(Brushes.Red, null, new System.Windows.Rect(0, 0, 100, 100));
            }
            var drawingImage = new DrawingImage(drawingVisual.Drawing);
            vm.CoverageMapViewModel.CoverageDrawing = drawingImage;
            // Setup a selected group for export path
            vm.SelectedGroup = new BioPolymerGroupViewModel("ACC", "Protein", "ABC", new List<BioPolymerCoverageResultModel>());
            var method = typeof(BioPolymerTabViewModel)
                .GetMethod("ExportImage", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            MessageBoxHelper.SuppressMessageBoxes = true;
            Assert.DoesNotThrow(() => method.Invoke(vm, null));
            // Check that file was created
            var expectedPath = Path.Combine(vm.ExportDirectory, "ACC_SequenceCoverage.Pdf");
            Assert.That(File.Exists(expectedPath), Is.True);
            File.Delete(expectedPath);
        }

        [Test]
        public void ResetDatabase_ClearsAllBioPolymers()
        {
            // Prepare Files
            if (!EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpQValue, out var searchTestCase))
                Assert.Fail();

            var dbPath = searchTestCase.DatabaseList.First().FilePath;
            var smPath = Directory.GetFiles(searchTestCase.OutputDirectory, "*.psmtsv", SearchOption.AllDirectories).First();
            var matches = Readers.SpectrumMatchTsvReader.ReadTsv(smPath, out _);

            // Build View Model
            var logic = new DummyMetaDrawLogic();
            logic.AllSpectralMatches = matches;
            var vm = new BioPolymerTabViewModel(logic);

            vm.DatabasePath = dbPath;
            var method = typeof(BioPolymerTabViewModel)
                .GetMethod("LoadDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)!
                .Invoke(vm, null);

            var allBioPolymersField = typeof(BioPolymerTabViewModel)
                .GetField("_allBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

            // Call ResetDatabase
            vm.GetType().GetMethod("ResetDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)!
                .Invoke(vm, null);

            // Assert all cleared
            Assert.That(vm.DatabasePath, Is.Empty);
            Assert.That(vm.DatabaseName, Is.Empty);
            Assert.That(vm.AllGroups.Count, Is.EqualTo(0));
            Assert.That(vm.FilteredGroups.Count, Is.EqualTo(0));
            var clearedDict = allBioPolymersField.GetValue(vm) as Dictionary<string, IBioPolymer>;
            Assert.That(clearedDict, Is.Not.Null);
            Assert.That(clearedDict.Count, Is.EqualTo(0));
        }

        [Test]
        public void SelectedGroup_Setter_UpdatesGroupAndCoverageDrawing()
        {
            var logic = new DummyMetaDrawLogic();
            var vm = new BioPolymerTabViewModel(logic);

            // Create two groups with different sequences
            var group1 = new BioPolymerGroupViewModel("A1", "Alpha", "ABC", new[] {
                new BioPolymerCoverageResultModel(new DummySpectralmatch(), "ABC", 1, 3, BioPolymerCoverageType.Unique)
            });
            var group2 = new BioPolymerGroupViewModel("B2", "Beta", "DEF", new[] {
                new BioPolymerCoverageResultModel(new DummySpectralmatch(), "DEF", 1, 3, BioPolymerCoverageType.Shared)
            });

            // Set first group and capture the plot
            vm.SelectedGroup = group1;
            var firstDrawing = vm.CoverageMapViewModel.CoverageDrawing;

            // Set second group and capture the plot
            vm.SelectedGroup = group2;
            var secondDrawing = vm.CoverageMapViewModel.CoverageDrawing;

            // Assert that the selected group and CoverageMapViewModel.Group are updated
            Assert.That(vm.SelectedGroup, Is.EqualTo(group2));
            Assert.That(vm.CoverageMapViewModel.Group, Is.EqualTo(group2));

            // Assert that the plot changes when the group changes
            Assert.That(firstDrawing, Is.Not.Null);
            Assert.That(secondDrawing, Is.Not.Null);
            Assert.That(secondDrawing, Is.Not.EqualTo(firstDrawing));
        }

    }
}
