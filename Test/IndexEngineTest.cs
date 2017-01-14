using InternalLogicEngineLayer;
using NUnit.Framework;
using System.Collections.Generic;
using OldInternalLogic;
using System.Linq;

namespace Test
{
	[TestFixture]
	public class IndexEngineTest
	{

		[Test]
		public void TestIndexEngine()
		{
			var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null, null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false) };
			var variableModifications = new List<MorpheusModification>();
			var fixedModifications = new List<MorpheusModification>();
			var localizeableModifications = new List<MorpheusModification>();
			var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), Terminus.C, CleavageSpecificity.Full, null, null, null);

			var engine = new IndexEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, protease);
			var results = (IndexResults)engine.Run();

			Assert.AreEqual(5, results.peptideIndex.Count);

			var listOfPeptides = results.peptideIndex.Select(b => string.Join("", b.BaseSequence.Select(c => char.ConvertFromUtf32(c)))).ToList();

			Assert.Contains("MNNNK", listOfPeptides);
			Assert.Contains("NNNK", listOfPeptides);
			Assert.Contains("QQQ", listOfPeptides);
			Assert.Contains("MNNNKQQQ", listOfPeptides);
			Assert.Contains("NNNKQQQ", listOfPeptides);

		}


	}
}