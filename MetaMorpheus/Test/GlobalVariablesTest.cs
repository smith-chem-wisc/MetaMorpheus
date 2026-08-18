using EngineLayer;
using NUnit.Framework;
using System;
using System.IO;

namespace Test
{
    [TestFixture]
    public static class GlobalVariablesTest
    {
        /// <summary>
        /// The UniProt proteome catalogue is a convenience for one GUI window, not a prerequisite for
        /// searching, so MetaMorpheus must still start without it. SetUpGlobalVariables runs in the
        /// MainWindow constructor and in Program.Main outside any try/catch, so an exception escaping
        /// LoadAvailableProteomes is an unhandled crash at launch — and mzLib now reports a missing
        /// catalogue by throwing rather than returning null (smith-chem-wisc/mzLib#1126).
        ///
        /// The file is shipped beside the executable, but DataDir can resolve to %LOCALAPPDATA%\MetaMorpheus
        /// for a Program Files install, or to a --customDataDir the user already created, and in those cases
        /// it may simply be absent. This hides it to reproduce that.
        /// </summary>
        [Test]
        [NonParallelizable]
        public static void TestStartUpSurvivesAMissingProteomeCatalogue()
        {
            string proteomesDirectory = Path.Combine(GlobalVariables.DataDir, "Proteomes");
            string hiddenDirectory = proteomesDirectory + ".hidden-for-test";

            // A run that died between the Move below and its finally block — a crash, a cancellation, a
            // killed test host — leaves the only copy of the real directory here. Put it back rather than
            // deleting it, or this run destroys the thing the test is meant to be protecting.
            if (Directory.Exists(hiddenDirectory))
            {
                if (Directory.Exists(proteomesDirectory))
                    Directory.Delete(hiddenDirectory, true);
                else
                    Directory.Move(hiddenDirectory, proteomesDirectory);
            }

            bool wasPresent = Directory.Exists(proteomesDirectory);
            if (wasPresent)
                Directory.Move(proteomesDirectory, hiddenDirectory);

            try
            {
                Assert.DoesNotThrow(() => GlobalVariables.SetUpGlobalVariables(),
                    "a missing proteome catalogue must not stop MetaMorpheus starting");

                // Empty, never null: DownloadUniProtDatabaseWindow enumerates this property and calls
                // FirstOrDefault on it with no null check, so null would be a latent NullReferenceException.
                Assert.That(GlobalVariables.AvailableUniProtProteomes, Is.Not.Null);
                Assert.That(GlobalVariables.AvailableUniProtProteomes, Is.Empty);
            }
            finally
            {
                if (wasPresent)
                {
                    if (Directory.Exists(proteomesDirectory))
                        Directory.Delete(proteomesDirectory, true);
                    Directory.Move(hiddenDirectory, proteomesDirectory);
                }

                // Leave the globals as the rest of the suite expects to find them.
                GlobalVariables.SetUpGlobalVariables();
            }
        }

        /// <summary>
        /// With the catalogue present it is read, and the human proteome is in it. The directory is
        /// asserted rather than tested for: the catalogue ships beside the executable, so if it is not
        /// there the fixture is broken (or the sibling test above left it hidden), and this test passing
        /// vacuously would hide that.
        /// </summary>
        [Test]
        [NonParallelizable]
        public static void TestProteomeCatalogueIsReadWhenPresent()
        {
            Assert.That(Directory.Exists(Path.Combine(GlobalVariables.DataDir, "Proteomes")), Is.True,
                "the shipped Proteomes directory is missing from the data directory");

            Assert.That(GlobalVariables.AvailableUniProtProteomes, Is.Not.Null);
            Assert.That(GlobalVariables.AvailableUniProtProteomes.ContainsKey("UP000005640"));
            Assert.That(GlobalVariables.AvailableUniProtProteomes["UP000005640"], Is.EqualTo("Homo sapiens (Human)"));
        }

        [Test]
        public static void TestCustomDataDirectory()
        {
            string originalDataDir = GlobalVariables.DataDir.ToString();
            string customDataDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CustomDataDir");

            if (Directory.Exists(customDataDir))
                Directory.Delete(customDataDir, true);

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

        [Test]
        public static void TestGetSpectralMatchLabel()
        {
            Assert.That(AnalyteType.Peptide.GetSpectralMatchLabel(), Is.EqualTo("PSM"));
            Assert.That(AnalyteType.Proteoform.GetSpectralMatchLabel(), Is.EqualTo("PSM"));
            Assert.That(AnalyteType.Oligo.GetSpectralMatchLabel(), Is.EqualTo("OSM"));
        }

        [Test]
        public static void TestGetSpectralMatchExtension()
        {
            Assert.That(AnalyteType.Peptide.GetSpectralMatchExtension(), Is.EqualTo("psmtsv"));
            Assert.That(AnalyteType.Proteoform.GetSpectralMatchExtension(), Is.EqualTo("psmtsv"));
            Assert.That(AnalyteType.Oligo.GetSpectralMatchExtension(), Is.EqualTo("osmtsv"));
        }

        [Test]
        public static void TestGetUniqueFormLabel()
        {
            Assert.That(AnalyteType.Peptide.GetUniqueFormLabel(), Is.EqualTo("Peptide"));
            Assert.That(AnalyteType.Peptide.ToString(), Is.EqualTo("Peptide"));
            Assert.That(AnalyteType.Proteoform.GetUniqueFormLabel(), Is.EqualTo("Proteoform"));
            Assert.That(AnalyteType.Proteoform.ToString(), Is.EqualTo("Proteoform"));
            Assert.That(AnalyteType.Oligo.GetUniqueFormLabel(), Is.EqualTo("Oligo"));
            Assert.That(AnalyteType.Oligo.ToString(), Is.EqualTo("Oligo"));
        }

        [Test]
        public static void TestGetBioPolymerLabel()
        {
            Assert.That(AnalyteType.Peptide.GetBioPolymerLabel(), Is.EqualTo("Protein"));
            Assert.That(AnalyteType.Proteoform.GetBioPolymerLabel(), Is.EqualTo("Protein"));
            Assert.That(AnalyteType.Oligo.GetBioPolymerLabel(), Is.EqualTo("Transcript"));
        }
        
        [Test]
        public static void TestDigestionAgentLabel()
        {
            Assert.That(AnalyteType.Peptide.GetDigestionAgentLabel(), Is.EqualTo("Protease"));
            Assert.That(AnalyteType.Proteoform.GetDigestionAgentLabel(), Is.EqualTo("Protease"));
            Assert.That(AnalyteType.Oligo.GetDigestionAgentLabel(), Is.EqualTo("RNase"));
        }
    }
}
