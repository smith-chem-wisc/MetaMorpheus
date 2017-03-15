// Copyright 2016 Stefan Solntsev
using EngineLayer;
using NUnit.Framework;
using System;
using System.IO;
using UsefulProteomicsDatabases;

namespace Test
{
    [SetUpFixture]
    public class MySetUpClass
    {

        #region Private Fields

        private const string elementsLocation = @"elements.dat";
        private const string unimodLocation = @"unimod_tables.xml";
        private const string uniprotLocation = @"ptmlist.txt";

        #endregion Private Fields

        #region Public Methods

        [OneTimeSetUp]
        public static void Setup()
        {
            //Console.WriteLine("Setting up tests...");
            //Console.WriteLine("Environment.CurrentDirectory is " + Environment.CurrentDirectory);
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
            //Console.WriteLine("Now Environment.CurrentDirectory is " + Environment.CurrentDirectory);
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, elementsLocation));

            //MyEngine.unimodDeserialized = Loaders.LoadUnimod(Path.Combine(TestContext.CurrentContext.TestDirectory, unimodLocation));
            //MyEngine.uniprotDeseralized = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, uniprotLocation));

            MyEngine.OutLabelStatusHandler += MyEngine_outLabelStatusHandler;
            MyEngine.FinishedSingleEngineHandler += MyEngine_FinishedSingleEngineHandler;
        }

        #endregion Public Methods

        #region Private Methods

        private static void MyEngine_FinishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
        {
            //Console.WriteLine(e.ToString());
        }

        private static void MyEngine_outLabelStatusHandler(object sender, StringEventArgs e)
        {
            // Console.WriteLine(e.s);
        }

        #endregion Private Methods

    }
}