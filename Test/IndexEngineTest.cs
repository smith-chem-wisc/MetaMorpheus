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
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null, null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, false, false, null) };
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, protease, InitiatorMethionineBehavior.Variable, 2, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, null);
            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(5, results.PeptideIndex.Count);

            var listOfPeptides = results.PeptideIndex.Select(b => string.Join("", b.BaseSequence.Select(c => char.ConvertFromUtf32(c)))).ToList();

            Assert.Contains("MNNNK", listOfPeptides);
            Assert.Contains("NNNK", listOfPeptides);
            Assert.Contains("QQQ", listOfPeptides);
            Assert.Contains("MNNNKQQQ", listOfPeptides);
            Assert.Contains("NNNKQQQ", listOfPeptides);
        }

        [Test]
        public static void TestIndexEngineWithWeirdSeq()
        {
            var proteinList = new List<Protein> { new Protein("MQXQ", null, null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, false, false, null) };
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var engine = new IndexingEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, protease, InitiatorMethionineBehavior.Retain, 2, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, null);
            var results = (IndexingResults)engine.Run();

            Assert.AreEqual(1, results.PeptideIndex.Count);

            Assert.IsNaN(results.PeptideIndex[0].MonoisotopicMassIncludingFixedMods);
            Assert.AreEqual(2, results.FragmentIndexDict.Count);
        }

        #endregion Public Methods

    }
}