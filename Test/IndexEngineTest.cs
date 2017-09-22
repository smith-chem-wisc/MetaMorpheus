using EngineLayer;
using EngineLayer.Indexing;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class IndexEngineTest
    {
        #region Public Methods

        [Test]
        public static void TestIndexEngine()
        {
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
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

            CommonParameters CommonParameters = new CommonParameters
            {
                DigestionParams = new DigestionParams
                {
                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                    MinPeptideLength = null,
                },
                ConserveMemory = false,
                ScoreCutoff = 1,
            };
            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.BnoB1ions, ProductType.Y }, 1, true, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.TotalPartitions, new List<string>());

            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(5, results.PeptideIndex.Count);

            var digestedList = proteinList[0].Digest(CommonParameters.DigestionParams, new List<ModificationWithMass>()).ToList();

            Assert.AreEqual(5, digestedList.Count);
            foreach (var fdfd in digestedList)
            {
                var dfdfse = fdfd.GetPeptidesWithSetModifications(CommonParameters.DigestionParams, variableModifications).ToList();
                Assert.AreEqual(1, dfdfse.Count);
                foreach (var kjdfk in dfdfse)
                {
                    Assert.Contains(kjdfk.CompactPeptide(TerminusType.None), results.PeptideIndex);
                }
            }
        }

        [Test]
        public static void TestIndexEngineWithWeirdSeq()
        {
            var proteinList = new List<Protein> { new Protein("MQXQ", null) };
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
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

            CommonParameters CommonParameters = new CommonParameters
            {
                DigestionParams = new DigestionParams
                {
                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                    MinPeptideLength = null,
                    InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain
                },
                ConserveMemory = false,
                ScoreCutoff = 1,
            };
            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.BnoB1ions, ProductType.Y }, 1, true, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.TotalPartitions, new List<string>());

            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(1, results.PeptideIndex.Count);

            Assert.IsNaN(results.PeptideIndex[0].MonoisotopicMassIncludingFixedMods);
            Assert.AreEqual(0, results.FragmentIndex.Length);
        }

        #endregion Public Methods
    }
}