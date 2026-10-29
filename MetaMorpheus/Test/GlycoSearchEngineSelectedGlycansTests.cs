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
        private static GlycoSearchEngine RunEngineConstructor(GlycoSearchType searchType, List<(string, string)> selectedGlycans, int maxOGlycanNum = 2)
        {
            var commonParameters = new CommonParameters(dissociationType: DissociationType.HCD, trimMsMsPeaks: false);

            return new GlycoSearchEngine(
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
            Assert.That(GlycanBox.GlobalOGlycans.Select(g => g.IdWithMotif), Is.EquivalentTo(chosen));
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
            Assert.That(keys, Is.EqualTo(Enumerable.Range(1, chosen.Count).Select(i => -i).ToList()));
        }

        [Test]
        public void NGlycanSearch_ASelection_NarrowsTheNGlycansTheEngineScores()
        {
            // An N-only search keeps its glycans on the engine, not in a GlycanBox static, so the N+O tests
            // above cannot see whether this path applies the selection.
            var wholeN = GlycanDatabase.LoadGlycan(PathOf(NGlycanDatabase), true, false).ToArray();
            var chosen = wholeN.Select(g => g.IdWithMotif).Distinct().Take(2).ToList();

            var engine = RunEngineConstructor(GlycoSearchType.NGlycanSearch, chosen.Select(id => (NGlycanDatabase, id)).ToList());

            Assert.That(engine.NGlycans.Select(g => g.IdWithMotif), Is.EquivalentTo(chosen));
            Assert.That(engine.NGlycans.Length, Is.LessThan(wholeN.Length));
        }

        /// <summary>
        /// Runs an O-glycan search task over one small spectra file with the given selection, and
        /// returns the manuscript prose it wrote and every warning it raised.
        /// </summary>
        private static (string Prose, List<string> Warnings) RunOGlycoTask(List<(string, string)> selectedGlycans, GlycoSearchType searchType = GlycoSearchType.OGlycanSearch)
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TESTGlycoSelectionProse");
            string dataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "N_O_glycoWithFileSpecific");
            string proteinDatabase = Path.Combine(dataFolder, "FourMucins_NoSigPeps_FASTA.fasta");
            string spectraFile = Path.Combine(dataFolder, "2019_09_16_StcEmix_35trig_EThcD25_rep1_4999_5500.mzML");

            var warnings = new List<string>();
            System.EventHandler<StringEventArgs> onWarn = (_, e) => { lock (warnings) { warnings.Add(e.S); } };
            TaskLayer.MetaMorpheusTask.WarnHandler += onWarn;
            MetaMorpheusEngine.WarnHandler += onWarn;
            try
            {
                Directory.CreateDirectory(outputFolder);
                var task = new TaskLayer.GlycoSearchTask
                {
                    CommonParameters = new CommonParameters(dissociationType: DissociationType.HCD, ms2childScanDissociationType: DissociationType.EThcD),
                    _glycoSearchParameters = new TaskLayer.GlycoSearchParameters
                    {
                        OGlycanDatabasefile = OGlycanDatabase,
                        NGlycanDatabasefile = NGlycanDatabase,
                        GlycoSearchType = searchType,
                        MaximumOGlycanAllowed = 2,
                        SelectedGlycans = selectedGlycans,
                    }
                };
                task.RunTask(outputFolder, new List<EngineLayer.DatabaseLoading.DbForTask> { new EngineLayer.DatabaseLoading.DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

                return (File.ReadAllText(Path.Combine(outputFolder, "AutoGeneratedManuscriptProse.txt")), warnings);
            }
            finally
            {
                TaskLayer.MetaMorpheusTask.WarnHandler -= onWarn;
                MetaMorpheusEngine.WarnHandler -= onWarn;
                if (Directory.Exists(outputFolder))
                {
                    Directory.Delete(outputFolder, true);
                }
            }
        }

        [Test]
        public void TaskRun_ASelectionNamingAnEntryNoLongerInTheDatabase_WarnsAndSaysSoInTheProse()
        {
            // The checked glycan the file has lost is skipped; the search must not stay silent about it.
            var kept = WholeODatabase().Select(g => g.IdWithMotif).First();
            var selection = new List<(string, string)> { (OGlycanDatabase, kept), (OGlycanDatabase, "H99N99 on S") };

            var (prose, warnings) = RunOGlycoTask(selection);

            Assert.That(warnings.Where(w => w.Contains("H99N99 on S")).ToList(), Has.Count.EqualTo(1),
                "one warning naming the missing entry, raised once rather than per partition or per file");
            Assert.That(warnings.Single(w => w.Contains("H99N99 on S")), Does.Contain(OGlycanDatabase));
            Assert.That(prose, Does.Contain("H99N99 on S"));
        }

        [Test]
        public void TaskRun_ASelectionNamingNothingThatExists_WarnsThatTheWholeDatabaseWasSearched()
        {
            var selection = new List<(string, string)> { (OGlycanDatabase, "H99N99 on S") };

            var (prose, warnings) = RunOGlycoTask(selection);

            var warning = warnings.Single(w => w.Contains("H99N99 on S"));
            Assert.That(warning, Does.Contain(OGlycanDatabase));
            Assert.That(warning, Does.Contain("whole database"));
            Assert.That(prose, Does.Contain("whole database"));
        }

        [Test]
        public void TaskRun_ASelection_IsNamedInTheMethodsProse()
        {
            var whole = WholeODatabase().Select(g => g.IdWithMotif).Distinct().ToList();
            var chosen = whole.Take(2).ToList();

            var (prose, _) = RunOGlycoTask(chosen.Select(id => (OGlycanDatabase, id)).ToList());

            Assert.That(prose, Does.Contain($"The O-glycan database: {OGlycanDatabase} (2 of {whole.Count} entries selected: {chosen[0]}, {chosen[1]})"));
        }

        [Test]
        public void TaskRun_NGlycanSearch_NamesTheNGlycanDatabaseInTheMethodsProse()
        {
            // The N-glycan line used to print the O-glycan database's file name.
            var wholeN = GlycanDatabase.LoadGlycan(PathOf(NGlycanDatabase), false, false).Select(g => g.IdWithMotif).Distinct().ToList();
            var chosen = wholeN.Take(2).ToList();

            var (prose, _) = RunOGlycoTask(chosen.Select(id => (NGlycanDatabase, id)).ToList(), GlycoSearchType.NGlycanSearch);

            Assert.That(prose, Does.Contain($"The N-glycan database: {NGlycanDatabase} (2 of {wholeN.Count} entries selected: {chosen[0]}, {chosen[1]})"));
            Assert.That(prose, Does.Not.Contain("The N-glycan database: " + OGlycanDatabase));
        }

        [Test]
        public void TaskRun_NoSelection_LeavesTheMethodsProseAsItWas()
        {
            var (prose, _) = RunOGlycoTask(new List<(string, string)>());

            Assert.That(prose, Does.Contain($"The O-glycan database: {OGlycanDatabase}\n"));
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
