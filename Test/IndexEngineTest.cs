using EngineLayer;
using EngineLayer.Indexing;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class IndexEngineTest
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

            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, protease, InitiatorMethionineBehavior.Variable, 2, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, null, false);
            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(5, results.PeptideIndex.Count);

            var digestedList = proteinList[0].Digest(protease, 2, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>(), false).ToList();

            Assert.AreEqual(5, digestedList.Count);
            foreach (var fdfd in digestedList)
            {
                var dfdfse = fdfd.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();
                Assert.AreEqual(1, dfdfse.Count);
                foreach (var kjdfk in dfdfse)
                {
                    Assert.Contains(kjdfk.CompactPeptide, results.PeptideIndex);
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

            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, protease, InitiatorMethionineBehavior.Retain, 2, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, null, false);
            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(1, results.PeptideIndex.Count);

            Assert.IsNaN(results.PeptideIndex[0].MonoisotopicMassIncludingFixedMods);
            Assert.AreEqual(2, results.FragmentIndexDict.Count);
        }

        #endregion Public Methods

    }
}