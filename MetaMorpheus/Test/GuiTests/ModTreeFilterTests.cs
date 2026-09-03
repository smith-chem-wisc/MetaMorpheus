using GuiFunctions;
using NUnit.Framework;
using System.Collections.ObjectModel;
using System.Linq;

namespace Test.GuiTests
{
    [TestFixture]
    public class ModTreeFilterTests
    {
        private static ObservableCollection<ModTypeForTreeViewModel> TwoGroups()
        {
            var collection = new ObservableCollection<ModTypeForTreeViewModel>();

            var biological = new ModTypeForTreeViewModel("Common Biological", false);
            biological.Children.Add(new ModForTreeViewModel("tip", false, "Acetylation on K", false, biological));
            biological.Children.Add(new ModForTreeViewModel("tip", false, "ADP-ribosylation on S", false, biological));
            collection.Add(biological);

            var glycans = new ModTypeForTreeViewModel("OGlycan.gdb", false);
            glycans.Children.Add(new ModForTreeViewModel("tip", false, "N1 on S", false, glycans, "N1 on S — HexNAc(1) — 203.08 Da"));
            glycans.Children.Add(new ModForTreeViewModel("tip", false, "H1N1 on T", false, glycans, "H1N1 on T — Hex(1)HexNAc(1) — 365.13 Da"));
            collection.Add(glycans);

            return collection;
        }

        [Test]
        public void Filter_KeepsOnlyGroupsHoldingAMatch()
        {
            var filtered = ModTreeFilter.Filter(TwoGroups(), "ribosyl");

            Assert.That(filtered.Count, Is.EqualTo(1));
            Assert.That(filtered.Single().DisplayName, Is.EqualTo("Common Biological"));
        }

        [Test]
        public void Filter_KeepsOnlyMatchingChildren()
        {
            var filtered = ModTreeFilter.Filter(TwoGroups(), "ribosyl");

            Assert.That(filtered.Single().Children.Count, Is.EqualTo(1));
            Assert.That(filtered.Single().Children.Single().ModName, Is.EqualTo("ADP-ribosylation on S"));
        }

        [Test]
        public void Filter_IsCaseInsensitive()
        {
            Assert.That(ModTreeFilter.Filter(TwoGroups(), "ACETYL").Single().Children.Single().ModName,
                Is.EqualTo("Acetylation on K"));
        }

        [Test]
        public void Filter_MatchesModNameByDefault()
        {
            // The default search text is the identifier, which is what the modification trees rely on.
            var filtered = ModTreeFilter.Filter(TwoGroups(), "H1N1");

            Assert.That(filtered.Single().DisplayName, Is.EqualTo("OGlycan.gdb"));
        }

        [Test]
        public void Filter_ByDefaultCannotFindTextThatIsOnlyInTheDisplayName()
        {
            // This is exactly why the overload exists: a glycan identifier is a composition code, so
            // "HexNAc" is absent from ModName even though the row displays it.
            Assert.That(ModTreeFilter.Filter(TwoGroups(), "HexNAc"), Is.Empty);
        }

        [Test]
        public void Filter_WithADisplayNameSelectorFindsTheExpandedComposition()
        {
            var filtered = ModTreeFilter.Filter(TwoGroups(), "HexNAc", m => m.DisplayName);

            Assert.That(filtered.Single().DisplayName, Is.EqualTo("OGlycan.gdb"));
            Assert.That(filtered.Single().Children.Count, Is.EqualTo(2));
        }

        [Test]
        public void Filter_WithADisplayNameSelectorFindsAMass()
        {
            var filtered = ModTreeFilter.Filter(TwoGroups(), "203.08", m => m.DisplayName);

            Assert.That(filtered.Single().Children.Single().ModName, Is.EqualTo("N1 on S"));
        }

        [Test]
        public void Filter_ExpandsTheGroupsItReturns()
        {
            Assert.That(ModTreeFilter.Filter(TwoGroups(), "ribosyl").Single().Expanded, Is.True);
        }

        [Test]
        public void Filter_SharesChildInstancesWithTheMasterCollection()
        {
            // What makes filtering non-destructive: ticking a row while filtered has to be visible in
            // the master collection, because that is what the read-back reads.
            var master = TwoGroups();

            var filtered = ModTreeFilter.Filter(master, "ribosyl");
            filtered.Single().Children.Single().Use = true;

            Assert.That(master[0].Children.Single(c => c.ModName == "ADP-ribosylation on S").Use, Is.True);
        }

        [Test]
        public void Filter_DoesNotMutateTheMasterCollection()
        {
            var master = TwoGroups();

            ModTreeFilter.Filter(master, "ribosyl");

            Assert.That(master.Count, Is.EqualTo(2));
            Assert.That(master[0].Children.Count, Is.EqualTo(2));
        }

        [Test]
        public void Filter_CopiesGroupCheckStateWithoutClobberingTheChildren()
        {
            // The clone's Use is written while its Children is still empty, because the Use setter
            // cascades. Writing it after the children were added would set every visible row to true.
            var master = TwoGroups();
            master[0].Use = true;

            var filtered = ModTreeFilter.Filter(master, "acetyl");

            Assert.That(filtered.Single().Use, Is.True);
            Assert.That(master[0].Children.Single(c => c.ModName == "ADP-ribosylation on S").Use, Is.True,
                "cascade from the master group, not from the filtered clone");
        }

        [Test]
        public void Filter_ReturnsNothingWhenNothingMatches()
        {
            Assert.That(ModTreeFilter.Filter(TwoGroups(), "no such modification"), Is.Empty);
        }

        [Test]
        public void Filter_WithAnEmptyKeyKeepsEverything()
        {
            var filtered = ModTreeFilter.Filter(TwoGroups(), "");

            Assert.That(filtered.Count, Is.EqualTo(2));
            Assert.That(filtered.Sum(g => g.Children.Count), Is.EqualTo(4));
        }

        [Test]
        public void DefaultSearchText_IsTheModName()
        {
            var parent = new ModTypeForTreeViewModel("group", false);
            var mod = new ModForTreeViewModel("tip", false, "Oxidation on M", false, parent, "a different label");

            Assert.That(ModTreeFilter.DefaultSearchText(mod), Is.EqualTo("Oxidation on M"));
        }
    }
}
