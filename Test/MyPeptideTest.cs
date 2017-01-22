using NUnit.Framework;
using OldInternalLogic;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class MyPeptideTest
    {

        #region Public Methods

        [Test]
        public static void TestGoodPeptide()
        {
            var prot = new Protein("MNNNKQQQQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).ToList();

            Assert.AreEqual(2, ye.Count);

			List<MorpheusModification> variableModifications = new List<MorpheusModification>();
			var pep1 = ye[0].GetPeptideWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

			var pep2 = ye[1].GetPeptideWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).ToList();

            Assert.AreEqual(2, ye.Count);

			List<MorpheusModification> variableModifications = new List<MorpheusModification>();
			var pep1 = ye[0].GetPeptideWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

			var pep2 = ye[1].GetPeptideWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.FastSortedProductMasses(new List<ProductType> { ProductType.Y });
            Assert.IsTrue(cool[0] > 0);
            Assert.IsNaN(cool[1]);
            Assert.IsNaN(cool[2]);
        }

        #endregion Public Methods

    }
}