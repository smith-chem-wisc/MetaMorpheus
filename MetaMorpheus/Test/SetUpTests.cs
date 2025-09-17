// Copyright 2016 Stefan Solntsev
using EngineLayer;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.IO;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [SetUpFixture]
    public class MySetUpClass
    {
        public static string outputFolder = "";

        private const string elementsLocation = @"elements.dat";

        [OneTimeSetUp]
        public static void GlobalSetup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
            Loaders.LoadElements();
            GlobalVariables.SetUpGlobalVariables();

            MetaMorpheusEngine.WarnHandler += WarnStatusHandler;
            MetaMorpheusTask.WarnHandler += WarnStatusHandler;

            EverythingRunnerEngine.FinishedAllTasksEngineHandler += SuccessfullyFinishedAllTasks;
        }

        [OneTimeTearDown]
        public static void GlobalTearDown()
        {
            EverythingRunnerEngineTestCase.DisposeAll();

            // Delete all "DatabaseIndex" folders in the test directory
            string testDir = TestContext.CurrentContext.TestDirectory;
            if (Directory.Exists(testDir))
            {
                foreach (var dir in Directory.GetDirectories(testDir, "DatabaseIndex", SearchOption.AllDirectories))
                {
                    try
                    {
                        Directory.Delete(dir, true);
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine($"Failed to delete directory '{dir}': {ex.Message}");
                    }
                }
            }
        }

        private static void SuccessfullyFinishedAllTasks(object sender, StringEventArgs rootOutputFolderPath)
        {
            outputFolder = rootOutputFolderPath.S;
        }

        private static void WarnStatusHandler(object sender, StringEventArgs e)
        {
            Console.WriteLine(e.S);
        }
    }
}