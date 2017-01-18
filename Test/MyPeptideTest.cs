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
            var prot = new Protein("MNNNKQQQQ", null, null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).ToList();

            Assert.AreEqual(2, ye.Count);

            var pep1 = new PeptideWithSetModifications(ye[0], new Dictionary<int, MorpheusModification>());
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = new PeptideWithSetModifications(ye[1], new Dictionary<int, MorpheusModification>());
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null, null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).ToList();

            Assert.AreEqual(2, ye.Count);

            var pep1 = new PeptideWithSetModifications(ye[0], new Dictionary<int, MorpheusModification>());
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = new PeptideWithSetModifications(ye[1], new Dictionary<int, MorpheusModification>());
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.FastSortedProductMasses(new List<ProductType> { ProductType.Y });
            Assert.IsTrue(cool[0] > 0);
            Assert.IsNaN(cool[1]);
            Assert.IsNaN(cool[2]);
        }

        #endregion Public Methods

    }
}