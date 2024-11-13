using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using Nett;
using NUnit.Framework;
using TaskLayer;


namespace Test
{
    [TestFixture]
    public static class ProteaseTests
    {
        [Test]
        public static void TestChymotrypsinCleaveAfterLeucine()
        {
            var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\Task1-SearchTaskconfig.toml");
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\TestProtease");
            Directory.CreateDirectory(outputFolder);
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\Q9UHB6_Chym_snip.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\Q9UHB6.fasta");

            var engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("SearchTOML", searchTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();

            string psmFile = Path.Combine(outputFolder, @"SearchTOML\AllPSMs.psmtsv");

            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            PsmFromTsv psm = parsedPsms.First();
            Assert.That(psm.BaseSeq, Is.EqualTo("TTQNQKSQDVELWEGEVVKEL")); //base sequence ends in leucine as expected
            Directory.Delete(outputFolder,true);
        }
    }
}
