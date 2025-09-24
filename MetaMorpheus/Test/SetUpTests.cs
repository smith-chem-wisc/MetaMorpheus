// Copyright 2016 Stefan Solntsev
using EngineLayer;
using NUnit.Framework; 
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.IO;
using TaskLayer;
using UsefulProteomicsDatabases;
using GuiFunctions;

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
            MessageBoxHelper.SuppressMessageBoxes = true;
            EverythingRunnerEngine.FinishedAllTasksEngineHandler += SuccessfullyFinishedAllTasks;
        }

        [OneTimeTearDown]
        public static void GlobalTearDown()
        {
            EverythingRunnerEngineTestCase.DisposeAll();
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