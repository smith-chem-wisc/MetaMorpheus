using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
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
using MzLibUtil;
using Nett;
using EngineLayer.GlycoSearch;
using MathNet.Numerics.LinearRegression;

namespace Test
{
    [TestFixture]
    public class XLTestGlycan
    {
        #region Basic Glycan Test
        [Test]
        public static void GlycanTest_Node()
        {
            Node node_left = new Node('H');
            Node node_right = new Node('F');
            Node node_middle = new Node('A');
            Node node_root = new Node('N', node_left, node_right, node_middle);
            node_root.Level = 1;
            string node2string = node_root.ToString();
            Assert.That(node_root.Level == 1);
            Assert.That(node2string == "N");              
        }

        [Test]
        public static void GlyTest_GetKindString()
        {
            byte[] kind = new byte[Glycan.SugarLength];
            kind[0] = 3; kind[1] = 4; kind[4] = 1; 
            string kindString = Glycan.GetKindString(kind);
            Assert.AreEqual("H3N4F1", kindString);
        }

        [Test]
        public static void GlyTest_DistinguishGlycans()
        {
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(H(H(H)))(H(H(H(H)))(H(H(H)))))))", 0);
            Glycan glycan2 = Glycan.Struct2Glycan("(N(N(H(H(H))(H(H(H))(H(H(H(H(H)))))))))", 0);

            var test = Glycan.Equals(glycan, glycan2);
            Assert.AreEqual(test, true);

            //TO DO: Test the glycan ions. 
            Glycan glycan3 = Glycan.Struct2Glycan("(N(F)(N(H(H(N(H(N(H(N(H))))))(N(H(N(H(N(F)(H(G))))))))(H(N(H(N(H(N(H(A)))))))(N(F)(H(N(F)(H(N(H)(F))))))))))", 8086);
            Glycan glycan4 = Glycan.Struct2Glycan("(N(F)(N(H(H(N(H(N(H(N(H))))))(N(H(N(H(N(F)(H(A))))))))(H(N(H(N(H(N(H(G)))))))(N(F)(H(N(F)(H(N(H)(F))))))))))", 8087);
        }

        [Test]
        public static void GlyTest_BisectHexNAc()
        {
            //The node here is for check the structure of the glycan. 
            Node node = Glycan.Struct2Node("(N(N(H(N)(H(N)(N))(H(N(H))))))"); //This glycan has a bisect hexnac 
            Assert.That(node.LeftChild.LeftChild.MiddleChild != null);

            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(N)(H(N)(N))(H(N(H))))))", 0);
            Assert.AreEqual(glycan.Ions.Count, 18);
        }

        [Test]
        public static void GlyTest_GlycanDecoy()
        {
            Glycan glycan = Glycan.Struct2Glycan("(N(N(H(N)(H(N)(N))(H(N(H))))))", 0);
            var test = Glycan.BuildTargetDecoyGlycans(new Glycan[] { glycan });
            Assert.AreEqual(test.Last().Decoy, true);
            foreach (var ion in test.Last().Ions)
            {
                Assert.AreEqual(ion.IonMass + ion.LossIonMass, test.Last().Mass);
            }
        }
        #endregion

        #region NGlycan Test

        [Test]
        public static void GlyTest_NGlycanCompositionFragments()
        {
            var kind = GlycanDatabase.String2Kind("HexNAc(3)Hex(4)Fuc(2)NeuAc(1)");

            var ions = GlycanDatabase.NGlycanCompositionFragments(kind);

            Glycan glycan = Glycan.Struct2Glycan("(N(F)(N(H(H)(H(N(F)(H(A)))))))", 0);

            var ionMass = ions.Select(p => p.IonMass).ToList();

            var glycanIonmass = glycan.Ions.Select(p => p.IonMass).ToList();

            var overlap = glycanIonmass.Intersect(ionMass).Count();

            Assert.That(overlap == 13);

        }

        [Test]
        public static void GlyTest_NGlycanCompositionFragments2()
        {
            //The purpose of this test is to test the code will run if glycans contain 'Xylose'. 
            //However, there is no strong evidence that '103' is the right answer. Since we don't know the structure.
            var kind = GlycanDatabase.String2Kind("HexNAc(5)Hex(8)Fuc(4)NeuAc(2)Xylose(2)");
            var x = GlycanDatabase.NGlycanCompositionFragments(kind);
            Assert.That(x.Count() == 103);
        }

        #endregion

        #region OGlycan Test
        //OGlycan test
        [Test]
        public static void OGlycoTest_LoadGlycanBox()
        {
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => p.Contains("OGlycan.gdb")).First(), true, true).ToArray();
            var OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(3).OrderBy(p => p.Mass).ToArray();
            Assert.AreEqual(OGlycanBoxes.Count(), 454);

            //To build the DecoyOGlycanBox, the number of glycan box will double.
            //However, how to build and applied the decoy glycan box is still unclear and require more investigation.
            var OGlycanBoxesWithDeocy = GlycanBox.BuildOGlycanBoxes(3, true).OrderBy(p => p.Mass).ToArray();
            Assert.That(OGlycanBoxesWithDeocy.Count() == 908);
        }

        [Test]
        public static void OGlycoTest_GetK()
        {
            List<int> input = new List<int> { 1, 2, 3, 4, 5 };

            //Combination test
            var kcombs = Glycan.GetKCombs(input, 3);

            Assert.AreEqual(kcombs.Count(), 10);

            var allcombs = Glycan.GetKCombs(input, 5);

            Assert.AreEqual(allcombs.Count(), 1);

            //Combination test with repetition
            var kcombs_rep = Glycan.GetKCombsWithRept(input, 3);

            Assert.AreEqual(kcombs_rep.Count(), 35);

            var allcombs_rep = Glycan.GetKCombsWithRept(input, 5);

            Assert.AreEqual(allcombs_rep.Count(), 126);

            //Permutation test
            var kperm = Glycan.GetPermutations(input, 3);

            Assert.AreEqual(kperm.Count(), 60);

            var allperm = Glycan.GetPermutations(input, 5).ToList();

            Assert.AreEqual(allperm.Count(), 120);

            //Permutation test with repetition
            var kperm_rep = Glycan.GetPermutationsWithRept(input, 3);

            Assert.AreEqual(kperm_rep.Count(), 125);
        }

        [Test]
        public static void OGlycoTest_OGlycanChildIons()
        {
            var glycan = GlycanBox.GlobalOGlycans[5];

            Assert.That(glycan.Ions.Count == 5);

            var kind = glycan.Kind;

            var glycanIons = GlycanDatabase.OGlycanCompositionCombinationChildIons(kind);

            Assert.That(glycanIons.Count() == 6);

            var coreIons = GlycanDatabase.OGlycanCompositionFragments(kind);
            Assert.That(coreIons.Count() == 6);
        }

        [Test]
        public static void OGlycoTest_GetKPerWithDuplicate()
        {
            List<int> input = new List<int> { 3, 5, 2, 7 };
            int[] ids = new int[3] { 2, 2, 3 };
            var perWithDuplicate = GlycoPeptides.GetPermutations(input, ids);
            var allPermutation = Glycan.GetPermutations(input, ids.Length);
            Assert.That(perWithDuplicate.Count() == allPermutation.Count() / 2);
        }

        #endregion
    }
}
