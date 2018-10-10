using EngineLayer;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using TaskLayer;
using Proteomics.ProteolyticDigestion;
using Proteomics.Fragmentation;

namespace Test
{
    [TestFixture]
    public static class SearchTaskTest
    {
        /// <summary>
        /// Tests each type of mass difference acceptor type to make sure values are assigned properly
        /// </summary>
        [Test]
        public static void MassDiffAceptorTest()
        {
            SearchTask searchTask = new SearchTask();
            var result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, searchTask.SearchParameters.MassDiffAcceptorType, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("1mm"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.TwoMM, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("2mm"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.ThreeMM, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("3mm"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.ModOpen, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("-187andUp"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Open, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("OpenSearch"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "custom ppmAroundZero 4");
            Assert.That(result.FileNameAddition.Equals("4ppmAroundZero"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Exact, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("5ppmAroundZero"));
        }

        /// <summary>
        /// Tests to make sure custom mass difference acceptor inputs are parsed properly
        /// </summary>
        [Test]
        public static void ParseSearchModeTest()
        {
            SearchTask searchTask = new SearchTask();
            var result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom dot 5 ppm 0,1.0029,2.0052");
            Assert.That(result.NumNotches == 3);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom dot 5 da 0,1.0029,2.0052");
            Assert.That(result.NumNotches == 3);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom interval [0;5],[0;5]");
            Assert.That(result.NumNotches == 1);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom OpenSearch 5");
            Assert.That(result.FileNameAddition.Equals("OpenSearch"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom daltonsAroundZero 5");
            Assert.That(result.FileNameAddition.Equals("5daltonsAroundZero"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom ppmAroundZero 5");
            Assert.That(result.FileNameAddition.Equals("5ppmAroundZero"));

            Assert.That(() => SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom Test 5"),
                Throws.TypeOf<MetaMorpheusException>());
        }

        /// <summary>
        /// Ensures semispecific search runs and outputs properly
        /// </summary>
        [Test]
        public static void SemiSpecificTest()
        {
            List<FragmentationTerminus> terminiToTest = new List<FragmentationTerminus>
            {
                FragmentationTerminus.N,
                FragmentationTerminus.C
            };
            foreach (FragmentationTerminus fragTerm in terminiToTest)
            {
                SearchTask searchTask = new SearchTask()
                {
                    SearchParameters = new SearchParameters
                    {
                        SearchType = SearchType.NonSpecific,
                        LocalFdrCategories = new List<FdrCategory>
                        {
                            FdrCategory.FullySpecific,
                            FdrCategory.SemiSpecific
                        }
                    },
                    CommonParameters = new CommonParameters(scoreCutoff: 4, addCompIons: true,
                    digestionParams: new DigestionParams(searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: fragTerm))
                };

                string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSameSettingsOutput");
                string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
                string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
                DbForTask db = new DbForTask(myDatabase, false);

                List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("TestSemiSpecific", searchTask) };

                var engine = new EverythingRunnerEngine(taskList, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, Environment.CurrentDirectory);
                engine.Run();

                string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSemiSpecific\AllPSMs.psmtsv");
                var output = File.ReadAllLines(outputPath);
                Assert.That(output.Length == 11); //if N is only producing 10 lines, then the c is not being searched with it.
            }
        }
               
        /// <summary>
        /// Tests that normalization in a search task works properly with an Experimental Design file read in,
        /// and crashes when that file is absent
        /// </summary>
        [Test]
        public static void PostSearchNormalizeTest()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    Normalize = true
                },
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestNormalization");
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\ExperimentalDesign.tsv");
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine("FileName\tCondition\tBiorep\tFraction\tTechrep");
                output.WriteLine("PrunedDbSpectra" + "\t" + "condition" + "\t" + "1" + "\t" + "1" + "\t" + "1");
            }
            DbForTask db = new DbForTask(myDatabase, false);

            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");

            File.Delete(filePath);

            Assert.That(() => searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal"),
               Throws.TypeOf<MetaMorpheusException>());
        }

        /// <summary>
        /// Test that we don't get a crash if protein groups are not constructed
        /// </summary>
        [Test]
        public static void ProteinGroupsNoParsimonyTest()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = false
                },
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestProteinGroupsNoParsimony");

            DbForTask db = new DbForTask(myDatabase, false);

            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");
        }

        /// <summary>
        /// Test ensures pruned databases are written when contaminant DB is searched
        /// </summary>
        [Test]
        public static void PrunedDbWithContaminantsTest()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true
                },
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestNormalization");
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\ExperimentalDesign.tsv");

            // contaminant DB
            DbForTask db = new DbForTask(myDatabase, true);

            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");

            Assert.That(File.ReadAllLines(Path.Combine(folderPath, @"DbForPrunedDbproteinPruned.xml")).Length > 0);
            Assert.That(File.ReadAllLines(Path.Combine(folderPath, @"DbForPrunedDbPruned.xml")).Length > 0);
        }
    }
}
