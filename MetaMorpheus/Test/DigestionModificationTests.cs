using EngineLayer;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Modifications;
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
            Assert.That(peptides.Count(), Is.EqualTo(2));
            Assert.That(peptides[0].FullSequence, !Is.EqualTo(peptides[1].FullSequence));
            Assert.That(peptides[0].MonoisotopicMass, Is.EqualTo(882.39707781799996));
            Assert.That(peptides[1].MonoisotopicMass, Is.EqualTo(930.400449121));

        }

       
    }
}
