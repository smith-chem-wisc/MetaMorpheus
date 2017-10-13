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
        #region Public Fields

        public static string outputFolder = "";

        #endregion Public Fields

        #region Private Fields

        private const string elementsLocation = @"elements.dat";
        private const string unimodLocation = @"unimod_tables.xml";
        private const string uniprotLocation = @"ptmlist.txt";

        #endregion Private Fields

        #region Public Methods

        [OneTimeSetUp]
        public static void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, elementsLocation));

            MetaMorpheusEngine.WarnHandler += WarnStatusHandler;
            MetaMorpheusTask.WarnHandler += WarnStatusHandler;

            EverythingRunnerEngine.FinishedAllTasksEngineHandler += SuccessfullyFinishedAllTasks;
        }

        #endregion Public Methods

        #region Private Methods

        private static void SuccessfullyFinishedAllTasks(object sender, StringEventArgs rootOutputFolderPath)
        {
            outputFolder = rootOutputFolderPath.s;
        }

        private static void WarnStatusHandler(object sender, StringEventArgs e)
        {
            Console.WriteLine(e.s);
        }

        #endregion Private Methods
    }
}