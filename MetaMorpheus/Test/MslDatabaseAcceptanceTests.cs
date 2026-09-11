using EngineLayer;
using EngineLayer.DatabaseLoading;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Pins the contract added by PR #2674: a .msl file is admitted as an accepted database
    /// format and flagged as a spectral library, on par with .msp. These fast unit tests guard
    /// the two-line change against a future refactor silently dropping .msl recognition.
    /// </summary>
    [TestFixture]
    public static class MslDatabaseAcceptanceTests
    {
        [OneTimeSetUp]
        public static void SetUp() => GlobalVariables.SetUpGlobalVariables();

        [Test]
        public static void AcceptedDatabaseFormats_ContainsMsl()
        {
            Assert.That(GlobalVariables.AcceptedDatabaseFormats, Does.Contain(".msl"));
            // Regression guard: existing .msp acceptance is preserved.
            Assert.That(GlobalVariables.AcceptedDatabaseFormats, Does.Contain(".msp"));
        }

        [Test]
        [TestCase("library.msl", true)]   // the new behavior
        [TestCase("library.MSL", true)]   // case-insensitive (ctor lower-cases the extension)
        [TestCase("library.msp", true)]   // unchanged behavior
        [TestCase("proteins.fasta", false)]
        [TestCase("proteins.xml", false)]
        public static void DbForTask_IsSpectralLibrary_RecognizesMsl(string path, bool expected)
        {
            var db = new DbForTask(path, false);
            Assert.That(db.IsSpectralLibrary, Is.EqualTo(expected));
        }
    }
}
