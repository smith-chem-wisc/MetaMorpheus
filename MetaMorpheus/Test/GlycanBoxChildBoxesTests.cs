using EngineLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace Test
{
    /// <summary>
    /// Locks in how a glycan box gets its child boxes, the sub-combinations the localization graph walks.
    /// </summary>
    /// <remarks>
    /// An N+O search with 43 O-glycans, 42 N-glycans and up to 3 O-glycans per peptide makes about 650,000 boxes, each with
    /// about 16 child boxes. Building every child box up front took most of the 14 s spent in
    /// <see cref="GlycanBox.BuildNOGlycanBoxes(int, bool, double)"/> on a semi-tryptic search and kept about 10 million
    /// objects alive for the whole search, but the localization graph only reads the child boxes of boxes that match a
    /// precursor mass. These tests pin that the builders leave child boxes unbuilt, that reading them gives exactly what the
    /// explicit child builders give, and that every read returns the same array (the localization cache and graph compare
    /// the array by reference).
    /// </remarks>
    [TestFixture]
    public class GlycanBoxChildBoxesTests
    {
        private const int MaxOGlycans = 2;

        [SetUp]
        public void LoadGlycanDatabases()
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
        }

        /// <summary>
        /// The N+O builder must not build any child boxes. Reading one box's child boxes builds that box's only.
        /// </summary>
        [Test]
        public static void BuildNOGlycanBoxes_LeavesChildBoxesUnbuiltUntilRead()
        {
            GlycanBox[] boxes = GlycanBox.BuildNOGlycanBoxes(MaxOGlycans, true, GlycanBox.DefaultMaximumGlycanBoxMass).ToArray();
            Assert.That(boxes, Is.Not.Empty);
            Assert.That(boxes.Count(b => b.HasBuiltChildGlycanBoxes), Is.EqualTo(0), "no box should build its child boxes before they are read");

            _ = boxes[boxes.Length / 2].ChildGlycanBoxes;

            Assert.That(boxes.Count(b => b.HasBuiltChildGlycanBoxes), Is.EqualTo(1), "reading one box's child boxes must build only that box's");
        }

        /// <summary>
        /// Same rule for the O-glycan-only builder.
        /// </summary>
        [Test]
        public static void BuildOGlycanBoxes_LeavesChildBoxesUnbuiltUntilRead()
        {
            GlycanBox[] boxes = GlycanBox.BuildOGlycanBoxes(MaxOGlycans, true, GlycanBox.DefaultMaximumGlycanBoxMass).ToArray();
            Assert.That(boxes, Is.Not.Empty);
            Assert.That(boxes.Count(b => b.HasBuiltChildGlycanBoxes), Is.EqualTo(0));

            _ = boxes[0].ChildGlycanBoxes;

            Assert.That(boxes.Count(b => b.HasBuiltChildGlycanBoxes), Is.EqualTo(1));
        }

        /// <summary>
        /// Every target N+O box's child boxes, built on read, are exactly what <see cref="GlycanBox.BulidChildNOBoxes"/> gives
        /// for that box: same order, ids, masses, compositions and target flag. Holds before and after child boxes became
        /// lazy, so the search sees the same child boxes either way.
        /// </summary>
        [Test]
        public static void NOBoxChildBoxes_AreExactlyTheExplicitChildBuildersOutput()
        {
            GlycanBox[] boxes = GlycanBox.BuildNOGlycanBoxes(MaxOGlycans, false, GlycanBox.DefaultMaximumGlycanBoxMass).ToArray();
            Assert.That(boxes.Any(b => b.NGlycanId != 0 && b.OGlycanIds != null && b.OGlycanIds.Length == MaxOGlycans), Is.True,
                "the fixture must include boxes with both an N-glycan and the maximum number of O-glycans");

            foreach (GlycanBox box in boxes)
            {
                GlycanBox[] expected = GlycanBox.BulidChildNOBoxes(box.NumberOfMods, box.ModIds, box.TargetDecoy).ToArray();
                AssertSameChildBoxes(box.ChildGlycanBoxes, expected, box.GlycanIdString);
            }
        }

        /// <summary>
        /// Same check for O-glycan-only boxes against <see cref="GlycanBox.BuildChildOGlycanBoxes"/>.
        /// </summary>
        [Test]
        public static void OBoxChildBoxes_AreExactlyTheExplicitChildBuildersOutput()
        {
            GlycanBox[] boxes = GlycanBox.BuildOGlycanBoxes(MaxOGlycans, false, GlycanBox.DefaultMaximumGlycanBoxMass).ToArray();
            Assert.That(boxes, Is.Not.Empty);

            foreach (GlycanBox box in boxes)
            {
                GlycanBox[] expected = GlycanBox.BuildChildOGlycanBoxes(box.NumberOfMods, box.ModIds, box.TargetDecoy).ToArray();
                AssertSameChildBoxes(box.ChildGlycanBoxes, expected, box.GlycanIdString);
            }
        }

        /// <summary>
        /// Every read, including reads racing on many threads (a task shares one set of boxes across all files and
        /// threads), returns one and the same array. <see cref="LocalizationGraph"/> only uses the box's localization cache
        /// when the child box array it was given is that same reference.
        /// </summary>
        [Test]
        public static void ChildGlycanBoxes_EveryReadReturnsTheSameArray()
        {
            GlycanBox[] boxes = GlycanBox.BuildNOGlycanBoxes(MaxOGlycans, false, GlycanBox.DefaultMaximumGlycanBoxMass).ToArray();
            GlycanBox box = boxes.Last(b => b.NumberOfMods == MaxOGlycans + 1);

            var seen = new GlycanBox[64][];
            Parallel.For(0, seen.Length, new ParallelOptions { MaxDegreeOfParallelism = 16 }, i => seen[i] = box.ChildGlycanBoxes);

            Assert.That(seen.All(children => children != null && ReferenceEquals(children, seen[0])), Is.True);
            Assert.That(box.ChildGlycanBoxes, Is.SameAs(seen[0]));
        }

        /// <summary>
        /// A box made directly with a constructor (child boxes themselves, and boxes in tests) still has no child boxes, and
        /// an array assigned to <see cref="GlycanBox.ChildGlycanBoxes"/> is returned as assigned.
        /// </summary>
        [Test]
        public static void ChildGlycanBoxes_ConstructedBoxHasNoneAndAnAssignedArrayIsKept()
        {
            var box = new GlycanBox(new[] { 0, 1 }, true);
            Assert.That(box.ChildGlycanBoxes, Is.Null);
            Assert.That(box.HasBuiltChildGlycanBoxes, Is.False);

            GlycanBox[] assigned = GlycanBox.BuildChildOGlycanBoxes(box.NumberOfMods, box.ModIds).ToArray();
            box.ChildGlycanBoxes = assigned;
            Assert.That(box.ChildGlycanBoxes, Is.SameAs(assigned));

            GlycanBox built = GlycanBox.BuildOGlycanBoxes(1, false, GlycanBox.DefaultMaximumGlycanBoxMass).First();
            built.ChildGlycanBoxes = assigned;
            Assert.That(built.ChildGlycanBoxes, Is.SameAs(assigned), "an assigned array must win over building on read");
        }

        private static void AssertSameChildBoxes(GlycanBox[] actual, GlycanBox[] expected, string parent)
        {
            Assert.That(actual, Is.Not.Null, $"box [{parent}] has no child boxes");
            Assert.That(actual.Length, Is.EqualTo(expected.Length), $"child box count for box [{parent}]");
            for (int i = 0; i < expected.Length; i++)
            {
                Assert.That(actual[i].ModIds, Is.EqualTo(expected[i].ModIds), $"child {i} ids of box [{parent}]");
                Assert.That(actual[i].OGlycanIds, Is.EqualTo(expected[i].OGlycanIds), $"child {i} O-glycan ids of box [{parent}]");
                Assert.That(actual[i].NGlycanId, Is.EqualTo(expected[i].NGlycanId), $"child {i} N-glycan id of box [{parent}]");
                Assert.That(actual[i].Mass, Is.EqualTo(expected[i].Mass), $"child {i} mass of box [{parent}]");
                Assert.That(actual[i].Kind, Is.EqualTo(expected[i].Kind), $"child {i} composition of box [{parent}]");
                Assert.That(actual[i].TargetDecoy, Is.EqualTo(expected[i].TargetDecoy), $"child {i} target flag of box [{parent}]");
            }
        }
    }
}
