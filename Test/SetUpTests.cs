// Copyright 2016 Stefan Solntsev
using EngineLayer;
using NUnit.Framework;
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
        public static void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
            Loaders.LoadElements();

            MetaMorpheusEngine.WarnHandler += WarnStatusHandler;
            MetaMorpheusTask.WarnHandler += WarnStatusHandler;

            EverythingRunnerEngine.FinishedAllTasksEngineHandler += SuccessfullyFinishedAllTasks;
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