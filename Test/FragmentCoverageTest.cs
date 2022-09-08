using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
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

namespace Test
{
    [TestFixture]
    public static class FragmentCoverageTest
    {

        [Test]
        public static void TestFragmentCoverage()
        {
            var testPSMs = PsmTsvReader.ReadTsv(@"TestData\SequenceCoverageTestPSM.psmtsv", out var warnings);

            foreach (var ion in testPSMs[0].MatchedIons)
            {
                Console.WriteLine(ion.NeutralTheoreticalProduct.ProductType);
            }
        }
    }
}