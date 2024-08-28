using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    public enum EverythingRunnerEngineTestCases
    {
        BottomUpQValue,
        BottomUpQValueNoIndividualFilesWriteMzId,
        BottomUpQValueNoIndividualFilesWritePepXml,
        BottomUpQValueSingle,
        BottomUpPepQValue,
        TopDownQValue,
        TopDownQValueSingle
    }

    /// <summary>
    /// Test cases for the post search analysis task. These test cases are used to verify that the post search analysis task is functioning correctly.
    /// This structure ensures that the database search is only ran once, and only ran once called.
    /// These directories are cleaned up in the Global Cleanup found in SetUpTests.GlobalTearDown
    /// </summary>
    [ExcludeFromCodeCoverage]
    internal class EverythingRunnerEngineTestCase : IDisposable
    {
        internal EverythingRunnerEngineTestCases TestCase { get; init; }
        internal List<(string, MetaMorpheusTask)> TaskList { get; init; }
        internal List<DbForTask> DatabaseList { get; init; }
        internal List<string> DataFileList { get; init; }
        internal string OutputDirectory => Path.Combine(ResultDirectory, TestCase.ToString());
        internal bool IsTopDown { get; init; }
        internal bool HasRun { get; private set; }
        internal bool WriteIndividualResults { get; init; }
        internal bool WritePepXml { get; init; }
        internal bool WriteMzId { get; init; }

        internal EverythingRunnerEngineTestCase(EverythingRunnerEngineTestCases testCase,
            List<(string, MetaMorpheusTask)> taskList, List<string> dataFileList,
            List<DbForTask> databaseList, bool isTopDown)
        {
            TestCase = testCase;
            TaskList = taskList;
            DatabaseList = databaseList;
            DataFileList = dataFileList;
            IsTopDown = isTopDown;
            HasRun = false;

            var firstSearchTask = taskList.Select(p => p.Item2).FirstOrDefault(p => p.TaskType == MyTask.Search);
            if (firstSearchTask is null) return;

            var searchTask = (SearchTask)firstSearchTask;
            WriteIndividualResults = searchTask.SearchParameters.WriteIndividualFiles;
            WritePepXml = searchTask.SearchParameters.WritePepXml;
            WriteMzId = searchTask.SearchParameters.WriteMzId;
        }

        internal void Run()
        {
            if (Directory.Exists(OutputDirectory))
                Directory.Delete(OutputDirectory, true);

            var runner = new EverythingRunnerEngine(TaskList, DataFileList, DatabaseList, OutputDirectory);
            runner.Run();
            HasRun = true;
        }

        public void Dispose()
        {
            if (Directory.Exists(OutputDirectory))
                Directory.Delete(OutputDirectory, true);
        }

        #region Case Setup

        internal static string ResultDirectory =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PostSearchAnalysisTaskTest");

        private static Dictionary<EverythingRunnerEngineTestCases, EverythingRunnerEngineTestCase> _cases;

        internal static bool TryGetTestCase(EverythingRunnerEngineTestCases testCase,
            out EverythingRunnerEngineTestCase outCase)
        {
            if (!_cases.TryGetValue(testCase, out outCase)) return false;

            if (!outCase.HasRun)
                outCase.Run();
            return true;
        }

        internal static EverythingRunnerEngineTestCase GetTestCase(EverythingRunnerEngineTestCases testCase)
        {
            if (!TryGetTestCase(testCase, out var outCase))
                throw new KeyNotFoundException($"Test case {testCase} not found");
            return outCase;
        }

        internal static void DisposeAll()
        {
            foreach (var testCase in _cases.Values)
                testCase.Dispose();
        }

        static EverythingRunnerEngineTestCase()
        {
            // Directory GlobalSetup
            if (Directory.Exists(ResultDirectory))
                Directory.Delete(ResultDirectory, true);

            if (!Directory.Exists(ResultDirectory))
                Directory.CreateDirectory(ResultDirectory);

            // Test Case Instantiation
            _cases = new();

            string myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\Task1-SearchTaskconfig.toml");
            SearchTask searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip_2.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\TaGe_SA_A549_3_snip.fasta");
            _cases.Add(EverythingRunnerEngineTestCases.BottomUpQValue,
                new EverythingRunnerEngineTestCase(EverythingRunnerEngineTestCases.BottomUpQValue,
                    new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                    false));
            _cases.Add(EverythingRunnerEngineTestCases.BottomUpQValueSingle,
                new EverythingRunnerEngineTestCase(EverythingRunnerEngineTestCases.BottomUpQValueSingle,
                    new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile2 },
                    new List<DbForTask> { new DbForTask(myDatabase, false) }, false));

            searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            searchTaskLoaded.SearchParameters.WriteIndividualFiles = false;
            searchTaskLoaded.SearchParameters.WriteMzId = true;
            searchTaskLoaded.SearchParameters.WritePepXml = false;
            _cases.Add(EverythingRunnerEngineTestCases.BottomUpQValueNoIndividualFilesWriteMzId,
                new EverythingRunnerEngineTestCase(EverythingRunnerEngineTestCases.BottomUpQValueNoIndividualFilesWriteMzId,
                    new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                    false));

            searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            searchTaskLoaded.SearchParameters.WriteIndividualFiles = false;
            searchTaskLoaded.SearchParameters.WriteMzId = false;
            searchTaskLoaded.SearchParameters.WritePepXml = true;
            _cases.Add(EverythingRunnerEngineTestCases.BottomUpQValueNoIndividualFilesWritePepXml,
                new EverythingRunnerEngineTestCase(EverythingRunnerEngineTestCases.BottomUpQValueNoIndividualFilesWritePepXml,
                    new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                    false));

            myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\Task2-SearchTaskconfig.toml");
            searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            searchTaskLoaded.CommonParameters.QValueCutoffForPepCalculation = 0.01;
            _cases.Add(EverythingRunnerEngineTestCases.BottomUpPepQValue,
                new EverythingRunnerEngineTestCase(EverythingRunnerEngineTestCases.BottomUpPepQValue,
                    new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                    false));

            myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TopDownTestData\TopDownSearchToml.toml");
            searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\slicedTDYeast.mzML");
            myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            _cases.Add(EverythingRunnerEngineTestCases.TopDownQValue,
                new EverythingRunnerEngineTestCase(EverythingRunnerEngineTestCases.TopDownQValue,
                    new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                    true));
            _cases.Add(EverythingRunnerEngineTestCases.TopDownQValueSingle,
                new EverythingRunnerEngineTestCase(EverythingRunnerEngineTestCases.TopDownQValueSingle,
                    new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) }, true));
        }

        #endregion
    }
}


