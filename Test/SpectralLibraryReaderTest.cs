using NUnit.Framework;
using System.IO;
using System;
using System.Linq;
using EngineLayer.spectralLibrarySearch;
using MzLibUtil;
using TaskLayer;
using System.Collections.Generic;
using Proteomics.Fragmentation;

namespace Test
{
    [TestFixture]
    public static class SpectralLibraryTest
    {
        [Test]
        public static void SpectralReaderWithoutDecpyTest()
        {
            //
            var testLibraryWithoutDecoy = new SpectralLibraryReader(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\myPrositLib.msp")).SpectralLibraryDictionary;
            var test1 = testLibraryWithoutDecoy["ALAVDGAGKPGAEE/2"];
            Assert.AreEqual(test1.Charge_state, 2);
            PpmTolerance newTolerance = new PpmTolerance(50);
            Assert.That(newTolerance.Within(test1.PrecursorMz, 642.825));
            Assert.AreEqual(test1.MatchedFragmentIons.Count, 31);
            Assert.That(newTolerance.Within(test1.MatchedFragmentIons[0].Mz, 148.06044));
            Assert.That(newTolerance.Within(test1.MatchedFragmentIons[0].Intensity, 0.03711248));
            Assert.AreEqual(test1.MatchedFragmentIons[0].Charge, 1);
        }
    }
}
