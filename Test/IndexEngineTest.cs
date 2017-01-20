using InternalLogicEngineLayer;
using NUnit.Framework;
using OldInternalLogic;
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
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false) };
            var variableModifications = new List<MorpheusModification>();
            var fixedModifications = new List<MorpheusModification>();
            var localizeableModifications = new List<MorpheusModification>();
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var engine = new IndexEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, protease, InitiatorMethionineBehavior.Variable, 2, 4096);
            var results = (IndexResults)engine.Run();

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
            var proteinList = new List<Protein> { new Protein("MQXQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false) };
            var variableModifications = new List<MorpheusModification>();
            var fixedModifications = new List<MorpheusModification>();
            var localizeableModifications = new List<MorpheusModification>();
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var engine = new IndexEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, protease, InitiatorMethionineBehavior.Retain, 2, 4096);
            var results = (IndexResults)engine.Run();

            Assert.AreEqual(1, results.PeptideIndex.Count);

            Assert.IsNaN(results.PeptideIndex[0].MonoisotopicMass);
            Assert.AreEqual(2, results.FragmentIndexDict.Count);
        }

        #endregion Public Methods

    }
}