using EngineLayer;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
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
            Assert.AreEqual(2, peptides.Count());
            Assert.AreNotEqual(peptides[0].FullSequence, peptides[1].FullSequence);
            Assert.AreEqual(882.39707781799996, peptides[0].MonoisotopicMass);
            Assert.AreEqual(930.400449121, peptides[1].MonoisotopicMass);
           
        }

        [Test]
        public static void EcoliDbTest()
        {
            Loaders.LoadElements();
            string subFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DigestionTest");
            Directory.CreateDirectory(subFolder);

            string databasePath1 = @"C:\Users\Alex\Documents\Proteomes\Escheria_coli_K12.fasta";
            var protDic = ProteaseDictionary.LoadProteaseDictionary(Path.Combine(GlobalVariables.DataDir, @"ProteolyticDigestion", @"proteases.tsv"), GlobalVariables.ProteaseMods);

            DigestionParams param = new DigestionParams();
            var proteinList = ProteinDbLoader.LoadProteinFasta(databasePath1, true, DecoyType.None, false, out List<string> errors);
            //var protein = proteinList[0];

            GlobalVariables.SetUpGlobalVariables();
            var methOx = GlobalVariables.AllModsKnown.First(mod => mod.IdWithMotif == "Oxidation on M");

            List<PeptideWithSetModifications> allPeptides = new();
            foreach(var protein in proteinList)
            {
                allPeptides.AddRange(protein.Digest(param, new List<Modification>(), new List<Modification> { methOx }));
            }

            List<PeptideWithSetModifications> uniquePeptides = allPeptides.DistinctBy(p => p.FullSequence).ToList();

            

            string humanDbPath = @"C:\Users\Alex\Documents\Proteomes\HumanProteome.fasta";
            proteinList = ProteinDbLoader.LoadProteinFasta(humanDbPath, true, DecoyType.None, false, out errors);

            List<PeptideWithSetModifications> humanPeptides = new();
            foreach (var protein in proteinList)
            {
                humanPeptides.AddRange(protein.Digest(param, new List<Modification>(), new List<Modification> { methOx }));
                
            }

            List<PeptideWithSetModifications> uniqueHumanPeptides = humanPeptides.DistinctBy(p => p.FullSequence).ToList();
            int placeholder = 0;
        }

        [Test]
        public static void UpsPeptideCount()
        {
            Loaders.LoadElements();
            string subFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DigestionTest");
            Directory.CreateDirectory(subFolder);

            string databasePath1 = @"C:\Users\Alex\Documents\Proteomes\ups1-ups2-sequences_corrected.fasta";
            var protDic = ProteaseDictionary.LoadProteaseDictionary(Path.Combine(GlobalVariables.DataDir, @"ProteolyticDigestion", @"proteases.tsv"), GlobalVariables.ProteaseMods);

            DigestionParams param = new DigestionParams();
            param.MaxLength = 25;

            var proteinList = ProteinDbLoader.LoadProteinFasta(databasePath1, true, DecoyType.None, false, out List<string> errors);
            //var protein = proteinList[0];

            GlobalVariables.SetUpGlobalVariables();
            var methOx = GlobalVariables.AllModsKnown.First(mod => mod.IdWithMotif == "Oxidation on M");

            List<PeptideWithSetModifications> allPeptides = new();
            foreach (var protein in proteinList)
            {
                allPeptides.AddRange(protein.Digest(param, new List<Modification>(), new List<Modification> {  }));
            }

            List<PeptideWithSetModifications> uniquePeptides = allPeptides.DistinctBy(p => p.FullSequence).ToList();

            var longestPep = uniquePeptides.MaxBy(pep => pep.BaseSequence.Length);
            var shortestPep = uniquePeptides.MinBy(pep => pep.BaseSequence.Length);

            string humanDbPath = @"C:\Users\Alex\Documents\Proteomes\HumanProteome.fasta";
            proteinList = ProteinDbLoader.LoadProteinFasta(humanDbPath, true, DecoyType.None, false, out errors);

            List<PeptideWithSetModifications> humanPeptides = new();
            foreach (var protein in proteinList)
            {
                humanPeptides.AddRange(protein.Digest(param, new List<Modification>(), new List<Modification> { methOx }));

            }

            List<PeptideWithSetModifications> uniqueHumanPeptides = humanPeptides.DistinctBy(p => p.FullSequence).ToList();
            int placeholder = 0;
        }

        [Test]
        public static void MultiDbTest()
        {
            Loaders.LoadElements();
            string parentFolder = @"C:\Users\Alex\Documents\Proteomes\MBR Proteomes March 2024";
            DirectoryInfo d = new DirectoryInfo(parentFolder);
            var proteomeFiles = d.GetFiles();
            Dictionary<string, int> dbFileToUniquePeptideCountDict = new();

            var protDic = ProteaseDictionary.LoadProteaseDictionary(Path.Combine(GlobalVariables.DataDir, @"ProteolyticDigestion", @"proteases.tsv"), GlobalVariables.ProteaseMods);
            DigestionParams param = new DigestionParams();
            GlobalVariables.SetUpGlobalVariables();
            var methOx = GlobalVariables.AllModsKnown.First(mod => mod.IdWithMotif == "Oxidation on M");

            foreach (FileInfo file in proteomeFiles)
            {
                var proteinList = ProteinDbLoader.LoadProteinFasta(file.FullName, true, DecoyType.None, false, out List<string> errors);

                List<PeptideWithSetModifications> allPeptides = new();
                foreach (var protein in proteinList)
                {
                    allPeptides.AddRange(protein.Digest(param, new List<Modification>(), new List<Modification> { methOx }));
                }

                List<PeptideWithSetModifications> uniquePeptides = allPeptides.DistinctBy(p => p.FullSequence).ToList();
                dbFileToUniquePeptideCountDict.Add(file.Name, uniquePeptides.Count);
            }
            
            using (StreamWriter sw = new StreamWriter(Path.Combine(parentFolder, "UniquePeptideCount.txt")))
            {
                foreach(var kvp in dbFileToUniquePeptideCountDict)
                {
                    sw.WriteLine(kvp.Key + ": " + kvp.Value.ToString());
                }
            }

        }

    }
}
