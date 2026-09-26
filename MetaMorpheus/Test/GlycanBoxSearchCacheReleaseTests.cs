using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Locks in that a glyco search lets go of what it built on its glycan boxes before post-search analysis.
    /// </summary>
    /// <remarks>
    /// During a search each glycan box builds, on first read, its child boxes, its motif count and its localization cache. The boxes
    /// stay reachable after the search through static fields (GlycanBox.OGlycanBoxes, NOGlycanBoxes and GlycoSpectralMatch.GlycanBoxes,
    /// which the result writer reads), so on an N+O search those caches kept about ten million small objects alive. The PEP model's
    /// FastTree trainer forces three full garbage collections per fit, 24 per task, and each had to walk them: on a six-file N+O
    /// search PEP spent 16 s of its 22.5 s in those pauses. Every cache rebuilds identically if it is read again, so dropping them
    /// cannot change results.
    /// </remarks>
    [TestFixture]
    public class GlycanBoxSearchCacheReleaseTests
    {
        private const int MaxOGlycans = 2;

        private static GlycanBox[] BuildNOBoxes()
        {
            string oglycanPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "OGlycan.gdb");
            string nglycanPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "NGlycan_NOBoxesTesting.gdb");

            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(oglycanPath, true, true).ToArray();
            GlycanBox.GlobalNGlycans = new Dictionary<int, Glycan>();
            int nGlycanId = -1;
            foreach (var nGlycan in GlycanDatabase.LoadGlycan(nglycanPath, true, false))
            {
                GlycanBox.GlobalNGlycans.Add(nGlycanId--, nGlycan);
            }
            return GlycanBox.BuildNOGlycanBoxes(MaxOGlycans, false, GlycanBox.DefaultMaximumGlycanBoxMass).OrderBy(p => p.Mass).ToArray();
        }

        [Test]
        [NonParallelizable] // writes the static glycan lists
        public static void ReleasingDropsWhatTheSearchBuiltAndItRebuildsTheSame()
        {
            GlycanBox[] boxes = BuildNOBoxes();
            GlycanBox box = boxes.Last(b => b.NumberOfMods == MaxOGlycans + 1);
            GlycanBox untouched = boxes.First();

            GlycanBox[] childrenBefore = box.ChildGlycanBoxes;
            Assert.That(box.GetMotifCount(), Is.Not.Null);
            Assert.That(box.GetLocalizationCache(), Is.Not.Null);
            Assert.That((box.HasBuiltChildGlycanBoxes, box.HasBuiltMotifCount, box.HasBuiltLocalizationCache), Is.EqualTo((true, true, true)),
                "the search's reads build all three");

            GlycanBox.ReleaseSearchCaches(boxes);

            Assert.That((box.HasBuiltChildGlycanBoxes, box.HasBuiltMotifCount, box.HasBuiltLocalizationCache), Is.EqualTo((false, false, false)),
                "releasing drops all three");
            Assert.That((untouched.HasBuiltChildGlycanBoxes, untouched.HasBuiltMotifCount, untouched.HasBuiltLocalizationCache), Is.EqualTo((false, false, false)));

            GlycanBox[] childrenAfter = box.ChildGlycanBoxes;
            Assert.That(childrenAfter, Is.Not.SameAs(childrenBefore), "read again, the child boxes are rebuilt");
            Assert.That(childrenAfter.Select(c => (string.Join(",", c.ModIds), c.Mass, c.TargetDecoy)),
                Is.EqualTo(childrenBefore.Select(c => (string.Join(",", c.ModIds), c.Mass, c.TargetDecoy))), "and are the same child boxes");
            Assert.That(box.GetMotifCount(), Is.Not.Null);
            Assert.That(box.GetLocalizationCache(), Is.Not.Null);
        }

        [Test]
        public static void ReleasingKeepsAnAssignedChildArrayAndIgnoresNoBoxes()
        {
            // no glycan ids, so the box needs no glycan database loaded
            var box = new GlycanBox(new int[0], true);
            var assigned = new GlycanBox[] { new GlycanBox(new int[0], true) };
            box.ChildGlycanBoxes = assigned;

            GlycanBox.ReleaseSearchCaches(new[] { box });
            Assert.That(box.ChildGlycanBoxes, Is.SameAs(assigned), "a box made with a constructor cannot rebuild an assigned array, so it keeps it");

            Assert.DoesNotThrow(() => GlycanBox.ReleaseSearchCaches(null));
        }

        /// <summary>
        /// After a real N+O glyco search, no glycan box still holds anything the search built on it.
        /// </summary>
        [Test]
        [NonParallelizable] // the glyco search writes process-wide glycan state
        public static void AGlycoSearchReleasesItsBoxCachesBeforePostSearchAnalysis()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestGlycanBoxSearchCacheRelease");
            Directory.CreateDirectory(outputFolder);
            var db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\Q9C0Y4.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\yeast_glycan_25170.mgf");
            var task = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\NGlycanSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            task._glycoSearchParameters.GlycoSearchType = GlycoSearchType.N_O_GlycanSearch;
            task._glycoSearchParameters.NGlycanDatabasefile = "NGlycan_ForNoSearch.gdb";
            task._glycoSearchParameters.OGlycanDatabasefile = "OGlycan.gdb";
            task._glycoSearchParameters.MaximumOGlycanAllowed = 1;
            string nGlycanDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "NGlycan_ForNoSearch.gdb");
            if (!GlobalVariables.NGlycanDatabasePaths.Contains(nGlycanDatabase))
            {
                GlobalVariables.NGlycanDatabasePaths.Add(nGlycanDatabase);
            }

            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, outputFolder).Run();

            Assert.That(File.Exists(Path.Combine(outputFolder, "Task", "AllPSMs.psmtsv")), Is.True, "the search must have run to the end");
            Assert.That(GlycanBox.NOGlycanBoxes, Is.Not.Empty);
            var holding = GlycanBox.NOGlycanBoxes.Where(b => b.HasBuiltChildGlycanBoxes || b.HasBuiltMotifCount || b.HasBuiltLocalizationCache).ToList();
            Assert.That(holding, Is.Empty, $"{holding.Count} of {GlycanBox.NOGlycanBoxes.Length} glycan boxes still hold caches the search built");

            Directory.Delete(outputFolder, true);
        }
    }
}
