// Copyright 2016 Stefan Solntsev
using InternalLogicEngineLayer;
using NUnit.Framework;
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
        public void setup()
        {
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, elementsLocation));
            //MyEngine.unimodDeserialized = Loaders.LoadUnimod(Path.Combine(TestContext.CurrentContext.TestDirectory, unimodLocation));
            //MyEngine.uniprotDeseralized = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, uniprotLocation));
        }

        #endregion Public Methods

    }
}