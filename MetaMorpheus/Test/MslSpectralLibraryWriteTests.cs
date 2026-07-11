using EngineLayer.DatabaseLoading;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Covers the .msl co-write paths added to MetaMorpheusTask (WriteSpectrumLibrary /
    /// UpdateSpectralLibrary): both the legacy .msp and the binary .msl are emitted, and
    /// entries that fail conversion are skipped from the .msl (with a warning) rather than
    /// crashing the write, while the .msp still contains every spectrum.
    /// </summary>
    [TestFixture]
    public static class MslSpectralLibraryWriteTests
    {
        /// <summary>
        /// Minimal concrete task that exposes the protected spectral-library writers.
        /// </summary>
        private sealed class TestableMetaMorpheusTask : MetaMorpheusTask
        {
            public TestableMetaMorpheusTask() : base(MyTask.Search) { }

            public static void WriteForTest(List<LibrarySpectrum> library, string outputFolder)
                => WriteSpectrumLibrary(library, outputFolder);

            public string UpdateForTest(List<LibrarySpectrum> library, string outputFolder)
                => UpdateSpectralLibrary(library, outputFolder);

            protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList,
                List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
                => null;
        }

        private static List<LibrarySpectrum> BuildMixedLibrary()
        {
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 100, 1, 1, 0), 100, 1000, 1)
            };

            return new List<LibrarySpectrum>
            {
                new LibrarySpectrum("PEPTIDEK", 500.0, 2, peaks, 10.0),
                new LibrarySpectrum("PEPTIDER", 600.0, 3, peaks, 12.0),
                // Non-positive charge fails MslLibraryEntry.TryFromLibrarySpectrum -> exercises the
                // skip-and-warn branch; must NOT be present in the resulting .msl.
                new LibrarySpectrum("BADENTRY", 700.0, 0, peaks, 14.0),
            };
        }

        [Test]
        public static void WriteSpectrumLibrary_CoWritesMspAndMsl_SkippingInvalidEntries()
        {
            string outDir = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "MslWriteTest_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(outDir);
            try
            {
                TestableMetaMorpheusTask.WriteForTest(BuildMixedLibrary(), outDir);

                var files = Directory.GetFiles(outDir).Select(Path.GetFileName).ToList();
                Assert.That(files.Any(f => f.Contains("SpectralLibrary") && f.EndsWith(".msp")), "Expected a co-written .msp");
                string msl = files.Single(f => f.Contains("SpectralLibrary") && f.EndsWith(".msl"));

                var mslEntries = new SpectralLibrary(new List<string> { Path.Combine(outDir, msl) })
                    .GetAllLibrarySpectra().ToList();

                // Both valid entries survive; the invalid one is dropped from the binary library.
                Assert.That(mslEntries.Count, Is.EqualTo(2));
                Assert.That(mslEntries.All(e => e.Sequence != "BADENTRY"));
            }
            finally
            {
                Directory.Delete(outDir, true);
            }
        }

        [Test]
        public static void UpdateSpectralLibrary_ReturnsMslPath_SkippingInvalidEntries()
        {
            string outDir = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "MslUpdateTest_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(outDir);
            try
            {
                string mslPath = new TestableMetaMorpheusTask().UpdateForTest(BuildMixedLibrary(), outDir);

                Assert.That(mslPath, Does.EndWith(".msl"));
                Assert.That(File.Exists(mslPath));

                var mslEntries = new SpectralLibrary(new List<string> { mslPath })
                    .GetAllLibrarySpectra().ToList();
                Assert.That(mslEntries.Count, Is.EqualTo(2));

                // The legacy .msp is co-written alongside the returned .msl.
                Assert.That(Directory.GetFiles(outDir, "*.msp").Any());
            }
            finally
            {
                Directory.Delete(outDir, true);
            }
        }
    }
}
