using EngineLayer;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    /// <summary>
    /// GlycoSearchEngine narrowing a glycan database to the individual glycans the user checked.
    /// </summary>
    /// <remarks>
    /// Marked NonParallelizable because the subset lands on GlycanBox's static GlobalOGlycans, which
    /// GlobalVariables never resets between fixtures. TearDown puts the whole database back so a
    /// later fixture asserting a baked-in glycan count does not fail on this one's leftovers.
    /// </remarks>
    [TestFixture]
    [NonParallelizable]
    public class GlycoSearchEngineSelectedGlycansTests
    {
        private const string OGlycanDatabase = "OGlycan.gdb";
        private const string NGlycanDatabase = "NGlycan.gdb";

        private static string PathOf(string databaseFileName)
        {
            return GlobalVariables.OGlycanDatabasePaths.Concat(GlobalVariables.NGlycanDatabasePaths)
                .First(p => Path.GetFileName(p) == databaseFileName);
        }

        private static Glycan[] WholeODatabase()
        {
            return GlycanDatabase.LoadGlycan(PathOf(OGlycanDatabase), true, true).ToArray();
        }

        /// <summary>
        /// Builds the engine and returns nothing: everything under test is the effect its constructor
        /// has on GlycanBox's statics, which is where the subset has to be applied. It is applied
        /// before the array is assigned and before any box is built, because a box identifies a
        /// glycan by its POSITION in that array.
        /// </summary>
        private static void RunEngineConstructor(GlycoSearchType searchType, List<(string, string)> selectedGlycans, int maxOGlycanNum = 2)
        {
            var commonParameters = new CommonParameters(dissociationType: DissociationType.HCD, trimMsMsPeaks: false);

            _ = new GlycoSearchEngine(
                new List<GlycoSpectralMatch>[0],
                new Ms2ScanWithSpecificMass[0],
                new List<PeptideWithSetModifications>(),
                null,
                null,
                0,
                commonParameters,
                null,
                OGlycanDatabase,
                NGlycanDatabase,
                searchType,
                glycoSearchTopNum: 30,
                maxOGlycanNum: maxOGlycanNum,
                oxoniumIonFilter: false,
                nestedIds: null,
                selectedGlycans: selectedGlycans);
        }

        [TearDown]
        public void TearDown()
        {
            // Leave the statics holding the whole databases, the way start-up and every other fixture
            // expect to find them. N_O_GlycanSearch settles both sides plus the box arrays.
            RunEngineConstructor(GlycoSearchType.N_O_GlycanSearch, null);
            RunEngineConstructor(GlycoSearchType.OGlycanSearch, null);
        }

        [Test]
        public void NoSelection_SearchesTheWholeDatabase()
        {
            RunEngineConstructor(GlycoSearchType.OGlycanSearch, null);

            Assert.That(GlycanBox.GlobalOGlycans.Length, Is.EqualTo(WholeODatabase().Length));
        }

        [Test]
        public void EmptySelection_SearchesTheWholeDatabase()
        {
            // This is the compatibility case: it is what every task written before the feature says.
            RunEngineConstructor(GlycoSearchType.OGlycanSearch, new List<(string, string)>());

            Assert.That(GlycanBox.GlobalOGlycans.Length, Is.EqualTo(WholeODatabase().Length));
        }

        [Test]
        public void ASelection_NarrowsTheDatabaseToTheChosenGlycans()
        {
            var chosen = WholeODatabase().Select(g => g.IdWithMotif).Distinct().Take(3).ToList();
            var selection = chosen.Select(id => (OGlycanDatabase, id)).ToList();

            RunEngineConstructor(GlycoSearchType.OGlycanSearch, selection);

            Assert.That(GlycanBox.GlobalOGlycans.Length, Is.LessThan(WholeODatabase().Length));
            Assert.That(GlycanBox.GlobalOGlycans.Select(g => g.IdWithMotif).Distinct(), Is.EquivalentTo(chosen));
        }

        [Test]
        public void ASelectionNamingNothingThatExists_FallsBackToTheWholeDatabase()
        {
            // Better than handing BuildOGlycanBoxes an empty array, whose failure surfaces much later
            // and elsewhere -- GlycanBoxes.First().Mass inside the parallel search loop.
            var selection = new List<(string, string)> { (OGlycanDatabase, "H99N99 on S") };

            RunEngineConstructor(GlycoSearchType.OGlycanSearch, selection);

            Assert.That(GlycanBox.GlobalOGlycans.Length, Is.EqualTo(WholeODatabase().Length));
        }

        [Test]
        public void ASelectionNamingOnlyTheOtherDatabase_LeavesThisOneWhole()
        {
            // Choosing individual N-glycans must not silently narrow the O-glycan side too.
            var selection = new List<(string, string)> { (NGlycanDatabase, "H5N2 on Nxs") };

            RunEngineConstructor(GlycoSearchType.OGlycanSearch, selection);

            Assert.That(GlycanBox.GlobalOGlycans.Length, Is.EqualTo(WholeODatabase().Length));
        }

        [Test]
        public void ANarrowedDatabase_BuildsFewerBoxes()
        {
            // The point of the feature. BuildOGlycanBoxes is combinations-with-repetition over the
            // whole array, so the box count is what a selection actually collapses.
            RunEngineConstructor(GlycoSearchType.OGlycanSearch, null);
            int boxesForTheWholeDatabase = GlycanBox.OGlycanBoxes.Length;

            var chosen = WholeODatabase().Select(g => g.IdWithMotif).Distinct().Take(2).ToList();
            RunEngineConstructor(GlycoSearchType.OGlycanSearch, chosen.Select(id => (OGlycanDatabase, id)).ToList());

            Assert.That(GlycanBox.OGlycanBoxes.Length, Is.LessThan(boxesForTheWholeDatabase));
        }

        [Test]
        public void ANarrowedDatabase_BuildsBoxesThatPointOnlyAtChosenGlycans()
        {
            // The reason the filter is applied before GlobalOGlycans is assigned. A box identifies a
            // glycan by its index into that array, so filtering after the assignment would leave every
            // box pointing at a different glycan than the one the user picked -- silently.
            var chosen = WholeODatabase().Select(g => g.IdWithMotif).Distinct().Take(3).ToList();

            RunEngineConstructor(GlycoSearchType.OGlycanSearch, chosen.Select(id => (OGlycanDatabase, id)).ToList());

            Assert.That(GlycanBox.OGlycanBoxes, Is.Not.Empty);
            foreach (var box in GlycanBox.OGlycanBoxes)
            {
                foreach (var id in box.ModIds)
                {
                    Assert.That(id, Is.InRange(0, GlycanBox.GlobalOGlycans.Length - 1));
                    Assert.That(chosen, Does.Contain(GlycanBox.GlobalOGlycans[id].IdWithMotif));
                }
            }
        }

        [Test]
        public void NOGlycanSearch_NarrowsTheNGlycanSideToo()
        {
            var wholeN = GlycanDatabase.LoadGlycan(PathOf(NGlycanDatabase), true, false).ToArray();
            var chosen = wholeN.Select(g => g.IdWithMotif).Distinct().Take(4).ToList();

            RunEngineConstructor(GlycoSearchType.N_O_GlycanSearch, chosen.Select(id => (NGlycanDatabase, id)).ToList());

            Assert.That(GlycanBox.GlobalNGlycans.Values.Select(g => g.IdWithMotif).Distinct(), Is.EquivalentTo(chosen));
            Assert.That(GlycanBox.GlobalNGlycans.Count, Is.LessThan(wholeN.Length));
        }

        [Test]
        public void NOGlycanSearch_ReIndexesTheNarrowedNGlycansDenselyAndNegatively()
        {
            // N-glycans are keyed by negative index to tell them apart from O-glycans, and the keys
            // are handed out in load order. Narrowing has to produce -1..-m with no gaps, or a box
            // referring to a key that was skipped resolves to nothing.
            var wholeN = GlycanDatabase.LoadGlycan(PathOf(NGlycanDatabase), true, false).ToArray();
            var chosen = wholeN.Select(g => g.IdWithMotif).Distinct().Take(4).ToList();

            RunEngineConstructor(GlycoSearchType.N_O_GlycanSearch, chosen.Select(id => (NGlycanDatabase, id)).ToList());

            var keys = GlycanBox.GlobalNGlycans.Keys.OrderByDescending(k => k).ToList();
            Assert.That(keys, Is.EqualTo(Enumerable.Range(1, keys.Count).Select(i => -i).ToList()));
        }

        [Test]
        public void NOGlycanSearch_NarrowingOnlyTheNSideLeavesTheOSideWhole()
        {
            var wholeN = GlycanDatabase.LoadGlycan(PathOf(NGlycanDatabase), true, false).ToArray();
            var chosen = wholeN.Select(g => g.IdWithMotif).Distinct().Take(2).ToList();

            RunEngineConstructor(GlycoSearchType.N_O_GlycanSearch, chosen.Select(id => (NGlycanDatabase, id)).ToList());

            Assert.That(GlycanBox.GlobalOGlycans.Length, Is.EqualTo(WholeODatabase().Length));
        }
    }
}
