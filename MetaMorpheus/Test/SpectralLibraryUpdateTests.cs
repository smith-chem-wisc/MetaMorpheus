using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using TaskLayer;
using EngineLayer.DatabaseLoading;

namespace Test
{
    [TestFixture]
    public static class SpectralLibraryUpdateTests
    {
        /// <summary>
        /// Issue #2291. Asking to update a spectral library without giving one used to run the whole search
        /// and then throw NullReferenceException out of UpdateSpectralLibrary, which surfaced as the task
        /// hanging on "Writing PSM results" with the exception only visible in results.txt afterwards.
        ///
        /// The combination is refused before searching now, so the wasted run is what this pins: reaching
        /// the check costs a database load, not a search.
        /// </summary>
        [Test]
        public static void UpdatingASpectralLibraryWithoutOneIsRefusedBeforeSearching()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "SpectralLibraryUpdateTests", "NoLibrary");
            Directory.CreateDirectory(outputFolder);

            var task = new SearchTask();
            task.SearchParameters.UpdateSpectralLibrary = true;

            string database = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "TestData", "hela_snip_for_unitTest.fasta");
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "TestData", "TaGe_SA_A549_3_snip.mzML");

            var thrown = Assert.Throws<MetaMorpheusException>(() => task.RunTask(
                outputFolder,
                new List<DbForTask> { new DbForTask(database, false) },
                new List<string> { spectra },
                "TestUpdateWithoutLibrary"));

            Assert.That(thrown.Message, Does.Contain("spectral library"));
            Assert.That(thrown.Message, Does.Contain("no spectral library was given").IgnoreCase,
                "say what is missing, not just that something went wrong");

            Directory.Delete(outputFolder, recursive: true);
        }
    }
}
