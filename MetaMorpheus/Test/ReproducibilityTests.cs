using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using MathNet.Numerics.Statistics;
using Nett;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using Readers;
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

        public static async Task DrainResources(CancellationToken cancellationToken)
        {
            Task.Run(() =>
            {
                while (!cancellationToken.IsCancellationRequested)
                {
                    // Tie up the CPU
                    double[] unsortedArray = new double[1000000];
                    Random random = new Random();
                    for (int i = 0; i < unsortedArray.Length; i++)
                    {
                        unsortedArray[i] = random.NextDouble();
                    }
                    Thread.Sleep((int)(random.NextDouble()*1000));
                    if (random.NextDouble() > 0.99)
                        Console.WriteLine(unsortedArray[0]);
                    Array.Sort(unsortedArray);
                }
            }, cancellationToken);
        }

        //[Test]
        //public static void ReproducibilityTest()
        //{
        //    EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpGPTMD, out var testCase);
        //    string outputFolder = testCase.OutputDirectory;
        //    var allResultsFile = Path.Combine(outputFolder, "allResults.txt");
        //}

        [Test]
        public static void TinyTest()
        {
            string vignetteFilePath = @"D:\MetaMorpheusVignette\04-30-13_CAST_Frac5_4uL.raw";
            string outputFolder = @"D:\MetaMorpheusVignette\CalibrationUnitTest";
            string dbPath = @"D:\MetaMorpheusVignette\uniprot-mouse-reviewed-1-24-2018.xml.gz";
            string contamPath = @"D:\MetaMorpheusVignette\uniprot-cRAP-1-24-2018.xml.gz";
            string myTomlPath = @"D:\MetaMorpheusVignette\Task2CalibrationTaskconfig.toml";
            CalibrationTask calibrationTask = Toml.ReadFile<CalibrationTask>(myTomlPath, MetaMorpheusTask.tomlConfig);

            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
            Directory.CreateDirectory(outputFolder);
            var everythingRunner = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("calibration", calibrationTask) },
                new List<string> { vignetteFilePath },
                new List<DbForTask> { new DbForTask(dbPath, false), new DbForTask(contamPath, true) },
                outputFolder
            );
            everythingRunner.Run();
        }

        [Test]
        public static void LocalCalTest()
        {

            string vingetteDataPath = @"D:\MetaMorpheusVignette\04-30-13_CAST_Frac5_4uL.raw";
            string vingetteDataPath2 = @"D:\MetaMorpheusVignette\04-30-13_CAST_Frac4_6uL.raw";
            string outputFolder = @"D:\MetaMorpheusVignette\CalibrationUnitTest";
            string dbPath = @"D:\MetaMorpheusVignette\uniprot-mouse-reviewed-1-24-2018.xml.gz";
            string contamPath = @"D:\MetaMorpheusVignette\uniprot-cRAP-1-24-2018.xml.gz";
            string myTomlPath = @"D:\MetaMorpheusVignette\Task2CalibrationTaskconfig.toml";
            CalibrationTask calibrationTask = Toml.ReadFile<CalibrationTask>(myTomlPath, MetaMorpheusTask.tomlConfig);

            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
            Directory.CreateDirectory(outputFolder);
            var everythingRunner = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("calibration", calibrationTask) },
                new List<string> { vingetteDataPath, vingetteDataPath2 },
                new List<DbForTask> { new DbForTask(dbPath, false), new DbForTask(contamPath, true) },
                outputFolder
            );
            everythingRunner.Run();


            List<double> mzDiscrepancies;
            do
            {
                // Start resource drainers
                CancellationTokenSource cts = new();
                var token = cts.Token;
                List<Task> resourceDrainers = new();
                Random rand = new Random();
                int numberOfDrainers = rand.Next(1, 22);
                for (int i = 0; i < numberOfDrainers; i++)
                {
                    Task drain = DrainResources(token);
                    resourceDrainers.Add(drain);
                }

                outputFolder = @"D:\MetaMorpheusVignette\CalibrationUnitTestAsync";
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                Directory.CreateDirectory(outputFolder);
                CalibrationTask calTask2 = Toml.ReadFile<CalibrationTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
                everythingRunner = new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)> { ("calibration", calTask2) },
                    new List<string> { vingetteDataPath, vingetteDataPath2 },
                    new List<DbForTask> { new DbForTask(dbPath, false), new DbForTask(contamPath, true) },
                    outputFolder
                );
                everythingRunner.Run();
            }
            while (TestIfCalibratedFilesAreEquivalent(out mzDiscrepancies));

           
            Console.WriteLine("No. unequal mz values: " + mzDiscrepancies.Count);

            Console.WriteLine("Average mz discrepancy: " + mzDiscrepancies.Mean());

            Console.WriteLine("Max mz Discrepancy: " + mzDiscrepancies.Max());
            // Check the calibrated outputs for equality

        }

        [Test]
        public static void ReadCalTest()
        {
            TestIfCalibratedFilesAreEquivalent(out var mzDiscrepancies);
            Console.WriteLine("No. unequal mz values: " + mzDiscrepancies.Count);

            Console.WriteLine("Average mz discrepancy: " + mzDiscrepancies.Mean());

            Console.WriteLine("Max mz Discrepancy: " + mzDiscrepancies.Max());
        }

        /// <summary>
        /// Returns true if the files are equivalent, false otherwise
        /// </summary>
        /// <param name="mzDiscrepancies"></param>
        /// <returns></returns>
        public static bool TestIfCalibratedFilesAreEquivalent(out List<double> mzDiscrepancies)
        {
            // Read in two mzmls, and go through every spectra point by point to make sure they are the same
            var reader = MsDataFileReader.GetDataFile(@"D:\MetaMorpheusVignette\CalibrationUnitTest\calibration\04-30-13_CAST_Frac5_4uL-calib.mzML");
            reader.LoadAllStaticData();

            var reader2 = MsDataFileReader.GetDataFile(@"D:\MetaMorpheusVignette\CalibrationUnitTestAsync\calibration\04-30-13_CAST_Frac5_4uL-calib.mzML");
            reader2.LoadAllStaticData();

            mzDiscrepancies = new List<double>();

            for (int i = 0; i < reader.Scans.Length; i++)
            {
                var spectrum1 = reader.Scans[i].MassSpectrum;
                var spectrum2 = reader2.Scans[i].MassSpectrum;

                if (spectrum1.Size != spectrum2.Size) Console.WriteLine("Spectra are unequal length. One based scan index: " + i);

                for (int j = 0; j < spectrum1.Size; j++)
                {
                    if (spectrum1.XArray[j] != spectrum2.XArray[j])
                    {
                        mzDiscrepancies.Add(Math.Abs(spectrum1.XArray[j] - spectrum2.XArray[j]));
                    }
                }
            }

            //Console.WriteLine("No. unequal mz values: " + mzDiscrepancies.Count);
            
            //Console.WriteLine("Average mz discrepancy" + mzDiscrepancies.Mean());

            //if(mzDiscrepancies.Any())
            //    Console.WriteLine("Maximum mz discrepancy" + mzDiscrepancies.Max());
            //if(intensityDiscrepancies.Any())
            //    Console.WriteLine("Maximum intensity discrepancy" + intensityDiscrepancies.Max());

            return !mzDiscrepancies.Any();
        }

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
