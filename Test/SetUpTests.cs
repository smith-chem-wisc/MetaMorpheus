// Copyright 2016 Stefan Solntsev
using NUnit.Framework;
using System.IO;
using UsefulProteomicsDatabases;

namespace Test
{
    [SetUpFixture]
    public class MySetUpClass
    {

        #region Public Methods

        [OneTimeSetUp]
        public void setup()
        {
            var elementLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "lal.dat");
            Loaders.LoadElements(elementLocation);
        }

        #endregion Public Methods

    }
}