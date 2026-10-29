using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using GuiFunctions;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Omics.Fragmentation;
using Readers;
using TaskLayer;

namespace Test.MetaDraw
{
    [ExcludeFromCodeCoverage]
    [NonParallelizable]
    internal class FragmentReanalysisRaceConditionTest
    {
        [Test]
        [NonParallelizable]
        public static void MatchIonsWithNewTypes_ConcurrentCalls_DoNotThrow()
        {
            // Load test data similar to TestFragmentationReanalysisViewModel_RematchIons
            var viewModel = new FragmentationReanalysisViewModel();

            // run a quick search 
            var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml");
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TestRaceCondition");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta");
            var engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("SearchTOML", searchTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();
            string psmFile = Path.Combine(outputFolder, @"SearchTOML\AllPSMs.psmtsv");
            var dataFile = MsDataFileReader.GetDataFile(myFile);

            // parse out psm and its respective scan
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            var psmToResearch = parsedPsms.First();
            var scan = dataFile.GetOneBasedScan(psmToResearch.Ms2ScanNumber);

            var startSignal = new ManualResetEventSlim(false);
            var exceptions = new ConcurrentQueue<Exception>();
            var tasks = Enumerable.Range(0, Environment.ProcessorCount)
                .Select(_ => Task.Run(() =>
                {
                    startSignal.Wait();
                    for (int i = 0; i < 25; i++)
                    {
                        try
                        {
                            viewModel.MatchIonsWithNewTypes(scan, psmToResearch, true);
                        }
                        catch (Exception ex)
                        {
                            exceptions.Enqueue(ex);
                            break;
                        }
                    }
                }))
                .ToArray();

            try
            {
                startSignal.Set();
                var completed = Task.WaitAll(tasks, TimeSpan.FromSeconds(30));
                Assert.That(completed, Is.True, "Concurrent rematching did not complete within the timeout.");
                Assert.That(exceptions, Is.Empty,
                    $"Concurrent rematching threw: {string.Join(Environment.NewLine, exceptions.Select(p => p.ToString()))}");
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }
    }
}
