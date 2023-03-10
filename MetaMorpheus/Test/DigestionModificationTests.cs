using EngineLayer;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    class DigestionModificationTests
    {
        [Test]
        public static void ProteaseModTest()
        {
            Loaders.LoadElements();
            string subFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DigestionTest");
            Directory.CreateDirectory(subFolder);

            string databasePath1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "ProteaseModTest.fasta");           
            var protDic = ProteaseDictionary.LoadProteaseDictionary(Path.Combine(GlobalVariables.DataDir, @"ProteolyticDigestion", @"proteases.tsv"), GlobalVariables.ProteaseMods);

            DigestionParams param = new DigestionParams(protease: "CNBr", maxMissedCleavages: 0);
            var proteinList = ProteinDbLoader.LoadProteinFasta(databasePath1, true, DecoyType.None, false, out List<string> errors);
            var protein = proteinList[0];
            var peptides = protein.Digest(param, new List<Modification>(), new List<Modification>()).ToList();
            Assert.AreEqual(2, peptides.Count());
            Assert.AreNotEqual(peptides[0].FullSequence, peptides[1].FullSequence);
            Assert.AreEqual(882.39707781799996, peptides[0].MonoisotopicMass);
            Assert.AreEqual(930.400449121, peptides[1].MonoisotopicMass);
           
        }

       
    }
}
