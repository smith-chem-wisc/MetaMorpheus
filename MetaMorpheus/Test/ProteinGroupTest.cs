using EngineLayer;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    [TestFixture]
    public static class ProteinGroupTest
    {
        [Test]
        public static void TestProteinGroupEquals()
        {
            Protein prot1 = new Protein("MEDEEK", "prot1");
            List<Protein> proteinList1 = new List<Protein> { prot1 };
            ProteinGroup proteinGroup1 = new ProteinGroup(new HashSet<Protein>(proteinList1),
                new HashSet<PeptideWithSetModifications>(), new HashSet<PeptideWithSetModifications>());

            List<Protein> proteinList2 = new List<Protein> { prot1 };
            ProteinGroup proteinGroup2 = new ProteinGroup(new HashSet<Protein>(proteinList2),
                new HashSet<PeptideWithSetModifications>(), new HashSet<PeptideWithSetModifications>());

            Assert.IsTrue(proteinGroup1.Equals(proteinGroup2));

            Protein prot3 = new Protein("EDEEK", "prot3");
            List<Protein> proteinList3 = new List<Protein> { prot3 };
            ProteinGroup proteinGroup3 = new ProteinGroup(new HashSet<Protein>(proteinList3),
                new HashSet<PeptideWithSetModifications>(), new HashSet<PeptideWithSetModifications>());

            Assert.IsFalse(proteinGroup1.Equals(proteinGroup3));

            object myObject = new object();
            Assert.IsFalse(proteinGroup1.Equals(myObject));

        }
    }
}
