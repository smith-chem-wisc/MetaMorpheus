using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Nett;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Uses test cases found in EverythingRunnerEngineTestCase.cs
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class ReproducibilityTests
    {
        public static Array GetTestCases() => Enum.GetValues(typeof(EverythingRunnerEngineTestCases));

        //public static async Task DrainResources(CancellationToken cancellationToken)
        //{
        //await Task.Run(() =>
        //{
        //    while(!cancellationToken.IsCancellationRequested)
        //    {
        //        // Tie up the CPU
        //        int[] unsortedArray = new int[1000000];
        //        Random random = new Random();
        //        for (int i = 0; i < unsortedArray.Length; i++)
        //        {
        //            unsortedArray[i] = random.Next(0, 1000000);
        //        }
        //        Array.Sort(unsortedArray);
        //    }
        //});
        //}

        //[Test]
        //public static void ReproducibilityTest()
        //{
        //    EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpGPTMD, out var testCase);
        //    string outputFolder = testCase.OutputDirectory;
        //    var allResultsFile = Path.Combine(outputFolder, "allResults.txt");
        //}

        [Test]
        public static void PeptideIntersectTest()
        {
            string database1 = @"\\bison.chem.wisc.edu\share\Users\AlexanderS_Bison\Detritus\uniprot-mouse-reviewed-1-24-2018GPTMD_010525.xml";
            string database2 = @"\\bison.chem.wisc.edu\share\Users\AlexanderS_Bison\Detritus\uniprot-mouse-reviewed-1-24-2018GPTMD_010625.xml";
            string myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\Task1-SearchTaskconfig.toml");
            SearchTask searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip_2.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip.fasta");
            searchTaskLoaded.CommonParameters.QValueCutoffForPepCalculation = 0.01;

            searchTaskLoaded.LoadModifications("taskId", out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);

            var listOfProteins1 = searchTaskLoaded.LoadProteins("taskId", new List<DbForTask> { new DbForTask(database1, false) }, true, DecoyType.Reverse, localizeableModificationTypes,
                searchTaskLoaded.CommonParameters);

            var listOfProteins2 = searchTaskLoaded.LoadProteins("taskId", new List<DbForTask> { new DbForTask(database2, false) }, true, DecoyType.Reverse, localizeableModificationTypes,
                searchTaskLoaded.CommonParameters);

            var uniquePeptides1 = new HashSet<PeptideWithSetModifications>();
            var uniquePeptides2 = new HashSet<PeptideWithSetModifications>();
            for (int i = 0; i < listOfProteins1.Count; i++)
            {
                if (listOfProteins1[i].Accession != listOfProteins2[i].Accession)
                    throw new Exception("Proteins are not in the same order in the two databases.");
                var pep1 = listOfProteins1[i].Digest(searchTaskLoaded.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                var pep2 = listOfProteins2[i].Digest(searchTaskLoaded.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                uniquePeptides1.UnionWith(pep1.Except(pep2));
                uniquePeptides2.UnionWith(pep2.Except(pep1));
            }

            using StreamWriter file = new StreamWriter(Path.Combine(@"D:\MetaMorpheusVignette", "uniquePeptides1.txt"));
            foreach (var pep in uniquePeptides1)
            {
                file.WriteLine(pep.FullSequence + " - " + pep.Parent.Accession);
            }
            using StreamWriter file2 = new StreamWriter(Path.Combine(@"D:\MetaMorpheusVignette", @"uniquePeptides2.txt"));
            foreach (var pep in uniquePeptides2)
            {
                file2.WriteLine(pep.FullSequence + " - " + pep.Parent.Accession);
            }

            Console.WriteLine("Unique peptides in database 1: " + uniquePeptides1.Count);
            Console.WriteLine("Unique peptides in database 2: " + uniquePeptides2.Count);

        }

    }
}
