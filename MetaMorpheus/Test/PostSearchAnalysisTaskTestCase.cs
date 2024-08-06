using System;
using System.Collections.Generic;
using System.IO;
using Easy.Common.Extensions;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test;

/// <summary>
/// Test cases for the post search analysis task. These test cases are used to verify that the post search analysis task is functioning correctly.
/// <remarks>
/// This structure ensures that the database search is only ran once, and only ran once called.
/// 
/// Two ways to access the test cases.
///     1. To use a single specific result: Call teh Test case directly from the static implementations
///         e.x. PostSearchAnalysisTaskTests.AllResultsAndResultTxtContainsCorrectValues_QValue_BottomUp();
///     2. To use all test case results: Call AllTestCaseIdentifiers in TestCaseSource to get all test case identifiers
///         and then call TestCaseLocator to get the test case.
///         e.x. PostSearchAnalysisTaskTests.AllResultAndResultTxtContainsCorrectNumberOfLines();
///
/// To add a new test case
///     1. make a private and public static test case
///     2. Implement the get method on public test case to run the search once
///     3. Add identifier to AllTestCaseIdentifiers
///     4. Add locator to TestCaseLocator
///
/// TODO: Move the instantiation and deconstruction out of the set-up and tear-down in PostSearchAnalysisTaskTests to enable these to be used globally. Eliminating the need to perform redundant searches
/// </remarks>
/// </summary>
public class PostSearchAnalysisTaskTestCase
{
    internal string Identifier { get; set; }
    internal string OutputDirectory { get; set; }
    internal int SpectraFileCount { get; set; }
    internal bool IsTopDown { get; init; }
    internal bool WriteIndividualResults { get; init; }
    internal bool WritePepXml { get; init; }

    #region TestCases

    internal static string ResultDirectory =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PostSearchAnalysisTaskTest");

    internal static PostSearchAnalysisTaskTestCase TestCaseLocator(string identifier)
    {
        return identifier switch
        {
            "BottomUpQValue" => BottomUpQValue,
            "BottomUpQValueSingle" => BottomUpQValueSingle,
            "BottomUpPepQValue" => BottomUpPepQValue,
            "TopDownQValue" => TopDownQValue,
            "TopDownQValueSingle" => TopDownQValueSingle,
            _ => throw new ArgumentOutOfRangeException(nameof(identifier), identifier, null)
        };
    }

    internal static string[] AllTestCaseIdentifiers => new[]
    {
        "BottomUpQValue",
        "BottomUpQValueSingle",
        "BottomUpPepQValue",
        "TopDownQValue",
        "TopDownQValueSingle"
    };

    private static PostSearchAnalysisTaskTestCase _bottomUpQValue;
    internal static PostSearchAnalysisTaskTestCase BottomUpQValue
    {
        get
        {
            if (!_bottomUpQValue.IsDefault()) return _bottomUpQValue;

            string outputFolder = Path.Combine(ResultDirectory, @"BottomUpQValue");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            string myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml");
            SearchTask searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip_2.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta");

            EverythingRunnerEngine engineToml =
                new(new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                    new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                    outputFolder);
            engineToml.Run();

            _bottomUpQValue = new PostSearchAnalysisTaskTestCase()
            {
                Identifier = "BottomUpQValue",
                OutputDirectory = outputFolder,
                SpectraFileCount = 2,
                IsTopDown = false,
                WriteIndividualResults = searchTaskLoaded.SearchParameters.WriteIndividualFiles,
                WritePepXml = searchTaskLoaded.SearchParameters.WritePepXml,
            };

            return _bottomUpQValue;
        }
    }

    private static PostSearchAnalysisTaskTestCase _bottomUpQValueSingle;
    internal static PostSearchAnalysisTaskTestCase BottomUpQValueSingle
    {
        get
        {
            if (!_bottomUpQValueSingle.IsDefault()) return _bottomUpQValueSingle;

            var outputFolder = Path.Combine(ResultDirectory, @"BottomUpQValueSingle");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            string myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml");
            SearchTask searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip_2.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta");

            var engineToml = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                new List<string> { myFile2 },
                new List<DbForTask> { new DbForTask(myDatabase, false) },
                outputFolder);
            engineToml.Run();

            _bottomUpQValueSingle = new PostSearchAnalysisTaskTestCase()
            {
                Identifier = "BottomUpQValueSingle",
                OutputDirectory = outputFolder,
                SpectraFileCount = 1,
                IsTopDown = false,
                WriteIndividualResults = searchTaskLoaded.SearchParameters.WriteIndividualFiles,
                WritePepXml = searchTaskLoaded.SearchParameters.WritePepXml,
            };

            return _bottomUpQValueSingle;
        }
    }

    private static PostSearchAnalysisTaskTestCase _bottomUpPepQValue;
    internal static PostSearchAnalysisTaskTestCase BottomUpPepQValue
    {
        get
        {
            if (!_bottomUpPepQValue.IsDefault()) return _bottomUpPepQValue;

            var outputFolder = Path.Combine(ResultDirectory, @"BottomUpPepQValue");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task2-SearchTaskconfig.toml");
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip_2.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta");

            // TODO: Uncomment this line and change values for PR 2394
            //searchTaskLoaded.CommonParameters.QValueCutoffForPepCalculation = 0.01;
            var engineToml = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                outputFolder);
            engineToml.Run();

            _bottomUpPepQValue = new PostSearchAnalysisTaskTestCase()
            {
                Identifier = "BottomUpPepQValue",
                OutputDirectory = outputFolder,
                SpectraFileCount = 2,
                IsTopDown = false,
                WriteIndividualResults = searchTaskLoaded.SearchParameters.WriteIndividualFiles,
                WritePepXml = searchTaskLoaded.SearchParameters.WritePepXml,
            };

            return _bottomUpPepQValue;
        }
    }

    private static PostSearchAnalysisTaskTestCase _topDownQValue;
    internal static PostSearchAnalysisTaskTestCase TopDownQValue
    {
        get
        {
            if (!_topDownQValue.IsDefault()) return _topDownQValue;

            var outputFolder = Path.Combine(ResultDirectory, @"TopDownQValue");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\TopDownSearchToml.toml");
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            var myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\slicedTDYeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            var engineToml = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                outputFolder);
            engineToml.Run();

            _topDownQValue = new PostSearchAnalysisTaskTestCase()
            {
                Identifier = "TopDownQValue",
                OutputDirectory = outputFolder,
                SpectraFileCount = 2,
                IsTopDown = true,
                WriteIndividualResults = searchTaskLoaded.SearchParameters.WriteIndividualFiles,
                WritePepXml = searchTaskLoaded.SearchParameters.WritePepXml,
            };

            return _topDownQValue;
        }
    }

    private static PostSearchAnalysisTaskTestCase _topDownQValueSingle;
    internal static PostSearchAnalysisTaskTestCase TopDownQValueSingle
    {
        get
        {
            if (!_topDownQValueSingle.IsDefault()) return _topDownQValueSingle;

            var outputFolder = Path.Combine(ResultDirectory, "TopDownQValueSingle");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\TopDownSearchToml.toml");
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string myFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            var myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\slicedTDYeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            var engineToml = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("postSearchAnalysisTaskTestOutput", searchTaskLoaded) },
                new List<string> { myFile1, myFile2 }, new List<DbForTask> { new DbForTask(myDatabase, false) },
                outputFolder);
            engineToml.Run();

            _topDownQValueSingle = new PostSearchAnalysisTaskTestCase()
            {
                Identifier = "TopDownQValueSingle",
                OutputDirectory = outputFolder,
                SpectraFileCount = 2,
                IsTopDown = true,
                WriteIndividualResults = searchTaskLoaded.SearchParameters.WriteIndividualFiles,
                WritePepXml = searchTaskLoaded.SearchParameters.WritePepXml,
            };

            return _topDownQValueSingle;
        }
    }

    #endregion
}