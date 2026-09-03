using EngineLayer;
using EngineLayer.GlycoSearch;
using GuiFunctions;
using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;

namespace Test.GuiTests
{
    [TestFixture]
    public class GlycanSelectionViewModelTests
    {
        /// <summary>
        /// A composition-format glycan, built the way GlycanDatabase builds one from a .gdb line.
        /// </summary>
        private static Glycan MakeGlycan(string composition, string motif, GlycanType type = GlycanType.O_glycan)
        {
            return new Glycan(GlycanDatabase.String2Kind(composition), motif, type);
        }

        private static Dictionary<string, List<Glycan>> TwoDatabases()
        {
            return new Dictionary<string, List<Glycan>>
            {
                ["OGlycan.gdb"] = new List<Glycan>
                {
                    MakeGlycan("HexNAc(1)", "S"),
                    MakeGlycan("HexNAc(1)Hex(1)", "S"),
                    MakeGlycan("HexNAc(1)Hex(1)NeuAc(1)", "T"),
                },
                ["NGlycan.gdb"] = new List<Glycan>
                {
                    MakeGlycan("HexNAc(2)Hex(5)", "Nxs", GlycanType.N_glycan),
                    MakeGlycan("HexNAc(2)Hex(6)", "Nxt", GlycanType.N_glycan),
                },
            };
        }

        [Test]
        public void Constructor_GroupsGlycansByDatabase()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());

            Assert.That(viewModel.Databases.Count, Is.EqualTo(2));
            Assert.That(viewModel.Databases.Select(d => d.DisplayName), Is.EqualTo(new[] { "OGlycan.gdb", "NGlycan.gdb" }));
            Assert.That(viewModel.Databases[0].Children.Count, Is.EqualTo(3));
            Assert.That(viewModel.Databases[1].Children.Count, Is.EqualTo(2));
        }

        [Test]
        public void Constructor_NothingIsCheckedByDefault()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());

            Assert.That(viewModel.Databases.SelectMany(d => d.Children).Any(c => c.Use), Is.False);
            Assert.That(viewModel.ToSelectedGlycans(), Is.Empty);
        }

        [Test]
        public void Constructor_DatabaseWithNoGlycansIsUncheckedNotIndeterminate()
        {
            // ModTypeForTreeViewModel takes `bad` as its second argument, not `use`, so Use is never
            // initialised and a bool? defaults to null -- which renders as the indeterminate box. An
            // empty database would otherwise be the only group that looked partly selected.
            var databases = new Dictionary<string, List<Glycan>>
            {
                ["OGlycan_Custom.gdb"] = new List<Glycan>(),
            };

            var viewModel = new GlycanSelectionViewModel(databases);

            Assert.That(viewModel.Databases.Single().Use, Is.False);
        }

        [Test]
        public void CheckingAGlycan_MakesItsDatabaseIndeterminate()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());
            var database = viewModel.Databases.First(d => d.DisplayName == "OGlycan.gdb");

            database.Children[0].Use = true;

            Assert.That(database.Use, Is.Null);
        }

        [Test]
        public void CheckingEveryGlycan_MakesItsDatabaseChecked()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());
            var database = viewModel.Databases.First(d => d.DisplayName == "OGlycan.gdb");

            foreach (var child in database.Children)
            {
                child.Use = true;
            }

            // Also the re-entrancy case: settling the parent to true cascades back down to every
            // child, each raising PropertyChanged again. Without the guard this never returns.
            Assert.That(database.Use, Is.True);
        }

        [Test]
        public void UncheckingTheLastGlycan_ReturnsItsDatabaseToUnchecked()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());
            var database = viewModel.Databases.First(d => d.DisplayName == "OGlycan.gdb");

            database.Children[0].Use = true;
            database.Children[0].Use = false;

            Assert.That(database.Use, Is.False);
        }

        [Test]
        public void Summary_ReportsNothingCheckedWithTheTotalAvailable()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());

            Assert.That(viewModel.Summary, Does.Contain("none checked"));
            Assert.That(viewModel.Summary, Does.Contain("5 entries available"));
        }

        [Test]
        public void Summary_CountsCheckedEntriesAndDatabases()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());

            viewModel.Databases[0].Children[0].Use = true;
            Assert.That(viewModel.Summary, Does.Contain("1 of 5 entries checked, across 1 database"));

            viewModel.Databases[1].Children[0].Use = true;
            Assert.That(viewModel.Summary, Does.Contain("2 of 5 entries checked, across 2 databases"));
        }

        [Test]
        public void Summary_RaisesPropertyChanged()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());
            bool propertyChangedFired = false;
            viewModel.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(viewModel.Summary)) propertyChangedFired = true; };

            viewModel.Databases[0].Children[0].Use = true;

            Assert.That(propertyChangedFired, Is.True);
        }

        [Test]
        public void ToSelectedGlycans_ProjectsDatabaseNameAndIdWithMotif()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());
            var glycan = viewModel.Databases[0].Children[1];
            glycan.Use = true;

            var selected = viewModel.ToSelectedGlycans();

            Assert.That(selected.Count, Is.EqualTo(1));
            Assert.That(selected[0].Item1, Is.EqualTo("OGlycan.gdb"));
            Assert.That(selected[0].Item2, Is.EqualTo(glycan.ModName));
        }

        [Test]
        public void ToSelectedGlycans_SpansMultipleDatabases()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases());
            viewModel.Databases[0].Children[0].Use = true;
            viewModel.Databases[1].Children[1].Use = true;

            var selected = viewModel.ToSelectedGlycans();

            Assert.That(selected.Select(s => s.Item1), Is.EquivalentTo(new[] { "OGlycan.gdb", "NGlycan.gdb" }));
        }

        [Test]
        public void Constructor_RestoresASavedSelection()
        {
            var databases = TwoDatabases();
            var wanted = databases["NGlycan.gdb"][1].IdWithMotif;

            var viewModel = new GlycanSelectionViewModel(databases, new[] { ("NGlycan.gdb", wanted) });

            var restored = viewModel.Databases.First(d => d.DisplayName == "NGlycan.gdb");
            Assert.That(restored.Children.Single(c => c.Use).ModName, Is.EqualTo(wanted));
            Assert.That(restored.Use, Is.Null, "one of two checked should read as indeterminate");
        }

        [Test]
        public void Constructor_RestoredSelectionRoundTripsThroughToSelectedGlycans()
        {
            var databases = TwoDatabases();
            var saved = new List<(string, string)>
            {
                ("OGlycan.gdb", databases["OGlycan.gdb"][0].IdWithMotif),
                ("NGlycan.gdb", databases["NGlycan.gdb"][0].IdWithMotif),
            };

            var viewModel = new GlycanSelectionViewModel(databases, saved);

            Assert.That(viewModel.ToSelectedGlycans(), Is.EquivalentTo(saved));
        }

        [Test]
        public void Constructor_MarksAGlycanTheDatabaseNoLongerHasAsBad()
        {
            var databases = TwoDatabases();

            var viewModel = new GlycanSelectionViewModel(databases, new[] { ("OGlycan.gdb", "H9N9 on S") });

            var database = viewModel.Databases.First(d => d.DisplayName == "OGlycan.gdb");
            var unknown = database.Children.Single(c => c.ModName == "H9N9 on S");
            Assert.That(unknown.Use, Is.True);
            Assert.That(unknown.ToolTipStuff, Is.EqualTo("UNKNOWN GLYCAN!"));
        }

        [Test]
        public void Constructor_InventsTheGroupWhenTheWholeDatabaseIsGone()
        {
            var viewModel = new GlycanSelectionViewModel(TwoDatabases(), new[] { ("Retired.gdb", "H9N9 on S") });

            var invented = viewModel.Databases.Single(d => d.DisplayName == "Retired.gdb");
            Assert.That(invented.Children.Single().ModName, Is.EqualTo("H9N9 on S"));
            Assert.That(invented.Use, Is.True);
        }

        [Test]
        public void Constructor_AcceptsNoGlycansAndNoSelection()
        {
            var viewModel = new GlycanSelectionViewModel(null, null);

            Assert.That(viewModel.Databases, Is.Empty);
            Assert.That(viewModel.ToSelectedGlycans(), Is.Empty);
            Assert.That(viewModel.Summary, Does.Contain("0 entries available"));
        }

        [Test]
        [TestCase("H5N4A2", "Hex(5)HexNAc(4)NeuAc(2)")]
        [TestCase("N1", "HexNAc(1)")]
        [TestCase("N1F1", "HexNAc(1)Fuc(1)")]
        public void ExpandComposition_SpellsOutTheMonosaccharides(string composition, string expected)
        {
            Assert.That(GlycanSelectionViewModel.ExpandComposition(composition), Is.EqualTo(expected));
        }

        [Test]
        public void ExpandComposition_TreatsAMissingCountAsOne()
        {
            Assert.That(GlycanSelectionViewModel.ExpandComposition("N"), Is.EqualTo("HexNAc(1)"));
        }

        [Test]
        [TestCase(null)]
        [TestCase("")]
        [TestCase("   ")]
        public void ExpandComposition_ReturnsEmptyForNothing(string composition)
        {
            Assert.That(GlycanSelectionViewModel.ExpandComposition(composition), Is.Empty);
        }

        [Test]
        public void ExpandComposition_PassesAnUnknownCodeThrough()
        {
            Assert.That(GlycanSelectionViewModel.ExpandComposition("Z2"), Is.EqualTo("Z(2)"));
        }

        [Test]
        public void GlycanLabel_CarriesIdentifierExpandedCompositionAndMass()
        {
            var glycan = MakeGlycan("HexNAc(1)", "S");

            var label = GlycanSelectionViewModel.GlycanLabel(glycan);

            // The whole point of the label: the identifier alone is a composition code, so searching
            // "HexNAc" would find nothing without the expansion.
            Assert.That(label, Does.Contain(glycan.IdWithMotif));
            Assert.That(label, Does.Contain("HexNAc(1)"));
            Assert.That(label, Does.Contain("Da"));
        }

        [Test]
        public void GlycanLabel_ShowsMassDividedOutOfItsScaledInteger()
        {
            var glycan = MakeGlycan("HexNAc(1)", "S");

            Assert.That(GlycanSelectionViewModel.GlycanLabel(glycan),
                Does.Contain((glycan.Mass / 1E5).ToString("F2", System.Globalization.CultureInfo.InvariantCulture)));
        }

        [Test]
        public void GlycanToolTip_ReportsCompositionExpansionMassAndType()
        {
            var glycan = MakeGlycan("HexNAc(2)Hex(5)", "Nxs", GlycanType.N_glycan);

            var toolTip = GlycanSelectionViewModel.GlycanToolTip(glycan);

            Assert.That(toolTip, Does.Contain("Composition: " + glycan.Composition));
            // Expansion follows the Kind[] slot order (Hex, HexNAc, ...), not the order the .gdb line
            // happens to write the monosaccharides in.
            Assert.That(toolTip, Does.Contain("Hex(5)HexNAc(2)"));
            Assert.That(toolTip, Does.Contain("Mass:"));
            Assert.That(toolTip, Does.Contain(GlycanType.N_glycan.ToString()));
        }

        [Test]
        public void GlycanToolTip_OmitsTheStructureLineForACompositionGlycan()
        {
            // Struc is only populated for structure-format databases.
            var glycan = MakeGlycan("HexNAc(1)", "S");

            Assert.That(GlycanSelectionViewModel.GlycanToolTip(glycan), Does.Not.Contain("Structure:"));
        }
    }
}
