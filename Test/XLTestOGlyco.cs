using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using MzLibUtil;
using Nett;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Annotations;

namespace Test
{
    [TestFixture]
    public static class XLTestOGlyco
    {
        [Test]
        public static void OGlycoTest_LoadGlycanBox()
        {
            var GlycanBoxes = Glycan.BuildGlycanBoxes(GlobalVariables.OGlycans.ToList(), 3);
            Assert.AreEqual(GlycanBoxes.Count(), 454);
        }

        [Test]
        public static void OGlycoTest_GetK()
        {
            List<int> input = new List<int> {1, 2, 3, 4, 5 };

            //Combination test
            var kcombs = GlycoPeptides.GetKCombs(input, 3);

            Assert.AreEqual(kcombs.Count(), 10);

            var allcombs = GlycoPeptides.GetKCombs(input, 5);

            Assert.AreEqual(allcombs.Count(), 1);

            //Combination test with repetition
            var kcombs_rep = GlycoPeptides.GetKCombsWithRept(input, 3);

            Assert.AreEqual(kcombs_rep.Count(), 35);

            var allcombs_rep = GlycoPeptides.GetKCombsWithRept(input, 5);

            Assert.AreEqual(allcombs_rep.Count(), 126);

            //Permutation
            var kperm = GlycoPeptides.GetPermutations(input, 3);

            Assert.AreEqual(kperm.Count(), 60);

            var allperm = GlycoPeptides.GetPermutations(input, 5).ToList();

            Assert.AreEqual(allperm.Count(), 120);

            //


        }
    }
}
