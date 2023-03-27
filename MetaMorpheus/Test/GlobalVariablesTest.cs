using EngineLayer;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Threading;

namespace Test
{
    [SetUpFixture]
    public static class SetUpGlobalVariables
    {
        [OneTimeSetUp]
        public static void LoadGlobalVariables()
        {
            // this loads the global variables once for all unit tests
            GlobalVariables.SetUpGlobalVariables();
        }
    }

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
            Assert.That(File.Exists(Path.Combine(customDataDir, @"Data\ptmlist.txt")));

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

        [Test]
        public static void TestCustomFileExtensionGetter()
        {
            string test1 = @"C:\myFile.fasta";
            Assert.That(GlobalVariables.GetFileExtension(test1, getUncompressedExtension: false) == ".fasta");
            Assert.That(GlobalVariables.GetFileExtension(test1, getUncompressedExtension: true) == ".fasta");
            Assert.That(GlobalVariables.GetFilenameWithoutExtension(test1) == "myFile");

            string test2 = @"C:\myFile.fasta.gz";
            Assert.That(GlobalVariables.GetFileExtension(test2, getUncompressedExtension: false) == ".fasta.gz");
            Assert.That(GlobalVariables.GetFileExtension(test2, getUncompressedExtension: true) == ".fasta");
            Assert.That(GlobalVariables.GetFilenameWithoutExtension(test2) == "myFile");

            string test3 = @"C:\myFile.11.1.mzML";
            Assert.That(GlobalVariables.GetFileExtension(test3, getUncompressedExtension: false) == ".mzML");
            Assert.That(GlobalVariables.GetFileExtension(test3, getUncompressedExtension: true) == ".mzML");
            Assert.That(GlobalVariables.GetFilenameWithoutExtension(test3) == "myFile.11.1");

            string test4 = @"C:\myFile.gz";
            Assert.That(GlobalVariables.GetFileExtension(test4, getUncompressedExtension: false) == ".gz");
            Assert.That(GlobalVariables.GetFileExtension(test4, getUncompressedExtension: true) == string.Empty);
            Assert.That(GlobalVariables.GetFilenameWithoutExtension(test4) == "myFile");

            string test5 = @"C:\myFile";
            Assert.That(GlobalVariables.GetFileExtension(test5, getUncompressedExtension: false) == string.Empty);
            Assert.That(GlobalVariables.GetFileExtension(test5, getUncompressedExtension: true) == string.Empty);
            Assert.That(GlobalVariables.GetFilenameWithoutExtension(test5) == "myFile");

            string test6 = @"C:\my.Fi.le.fasta.gz";
            Assert.That(GlobalVariables.GetFileExtension(test6, getUncompressedExtension: false) == ".fasta.gz");
            Assert.That(GlobalVariables.GetFileExtension(test6, getUncompressedExtension: true) == ".fasta");
            Assert.That(GlobalVariables.GetFilenameWithoutExtension(test6) == "my.Fi.le");

            string test7 = @"C:\my.Fi.le.fasta.bGZ";
            Assert.That(GlobalVariables.GetFileExtension(test7, getUncompressedExtension: false) == ".fasta.bGZ");
            Assert.That(GlobalVariables.GetFileExtension(test7, getUncompressedExtension: true) == ".fasta");
            Assert.That(GlobalVariables.GetFilenameWithoutExtension(test7) == "my.Fi.le");
        }
    }
}
