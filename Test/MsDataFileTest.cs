using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class MsDataFileTest
    {
        [OneTimeSetUp]
        public static void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        [Test]
        public static void TestLoadAndRunMgf()
        {
            //The purpose of this test is to ensure that mgfs can be run without crashing. 
            //Whenever a new feature is added that may require things an mgf does not have, 
            //there should be a check that prevents mgfs from using that feature.
            string mgfName = @"TestData\ok.mgf";
            string xmlName = @"TestData\okk.xml";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestLoadAndRunMgf");

            SearchTask task1 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DoQuantification = true
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>
            {
                ("task1", task1),
            };
            //run!

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mgfName }, new List<DbForTask> { new DbForTask(xmlName, false) }, outputFolder);
            engine.Run();
            //Just don't crash! There should also be at least one psm at 1% FDR, but can't check for that.
            Directory.Delete(outputFolder, true);
        }
    }
}
