using EngineLayer;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Threading;

namespace Test
{
    [TestFixture]
    public static class GlobalVariablesTest
    {
        [Test]
        public static void TestCustomDataDirectory()
        {
            string originalDataDir = GlobalVariables.DataDir.ToString();
            string customDataDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CustomDataDir");

            Assert.That(!Directory.Exists(customDataDir));

            // set the custom data directory + load it
            GlobalVariables.UserSpecifiedDataDir = customDataDir;
            GlobalVariables.SetUpGlobalVariables();

            // the DataDir should now be set to the custom data dir
            Assert.That(GlobalVariables.DataDir != originalDataDir);
            Assert.That(GlobalVariables.DataDir == customDataDir);

            // the custom data dir should have data in it
            Assert.That(!Directory.Exists(Path.Combine(customDataDir, @"Data\ptmlist.txt")));

            // reset the data dir to the default one so other unit tests will use that one 
            // and not the soon-to-be-deleted custom data dir
            GlobalVariables.UserSpecifiedDataDir = null;
            GlobalVariables.SetUpGlobalVariables();

            // clean up from unit test
            // these 2 lines need to be here to delete the CustomDataDir folder for some reason
            // https://stackoverflow.com/questions/3420737/file-delete-error-the-process-cannot-access-the-file-because-it-is-being-used-b/3421084
            GC.Collect();
            GC.WaitForPendingFinalizers();

            Directory.Delete(customDataDir, true);
        }
    }
}
