using EngineLayer;
using EngineLayer.Indexing;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class IndexEngineTest
    {
        [Test]
        public static void TestIndexEngine()
        {
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();

            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            Protease p = new Protease("Custom Protease2", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(p.Name, p);
            CommonParameters CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: new DigestionParams(protease: p.Name, minPeptideLength: 1));

            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, CommonParameters, 30000, false, new List<string>());

            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(5, results.PeptideIndex.Count);

            var digestedList = proteinList[0].Digest(CommonParameters.DigestionParams, new List<Modification>(), variableModifications).ToList();

            Assert.AreEqual(5, digestedList.Count);
            foreach (var fdfd in digestedList)
            {
                Assert.Contains(fdfd, results.PeptideIndex);
            }
        }

        [Test]
        public static void TestIndexEngineWithWeirdSeq()
        {
            var proteinList = new List<Protein> { new Protein("MQXQ", null) };
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();

            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
            foreach (var mod in fixedModifications)
            {
                modsDictionary.Add(mod, 0);
            }
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            Protease protease = new Protease("Custom Protease", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    minPeptideLength: 1,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
                scoreCutoff: 1);

            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, CommonParameters, 30000, false, new List<string>());

            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(1, results.PeptideIndex.Count);

            Assert.IsNaN(results.PeptideIndex[0].MonoisotopicMass);
            Assert.AreEqual(30000000 + 1, results.FragmentIndex.Length);
        }
    }
}