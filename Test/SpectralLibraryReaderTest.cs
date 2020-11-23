using NUnit.Framework;
using System.IO;
using System;
using System.Linq;
using EngineLayer;
using MzLibUtil;
using TaskLayer;
using System.Collections.Generic;
using Proteomics.Fragmentation;

namespace Test
{
    [TestFixture]
    public static class SpectralLibraryTest
    {
        [Test]
        public static void SpectralLibraryReaderTest()
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\myPrositLib.msp");

            var testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { path });
            var librarySpectra = testLibraryWithoutDecoy.GetAllLibrarySpectra().ToList();

            Assert.That(librarySpectra.Count == 5);
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("QSQHM[Common Variable:Oxidation on M]TEVVR", 5, out var spectrum));
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("M[Common Variable:Oxidation on M]C[Common Fixed:Carbamidomethyl on C]SDSDGLAPPQHLIR", 2, out spectrum));

            testLibraryWithoutDecoy.TryGetSpectrum("ALAVDGAGKPGAEE", 2, out var test1);

            Assert.AreEqual(test1.ChargeState, 2);

            var frags = new List<(double mz, double intensity, ProductType ProductType, int fragmentNumber, int charge, double ppm)>
            {
                (148.06044, 0.03711248, ProductType.y, 1, 1, 0.0),
                (277.10303, 0.025135221, ProductType.y, 2, 1, 0.0),
                (185.12845, 0.8128169, ProductType.b, 2, 1, 0.0),
                (348.14014, 0.008474186, ProductType.y, 3, 1, 0.0),
                (256.16556, 1.0, ProductType.b, 3, 1, 0.0),
                (405.1616, 0.006187536, ProductType.y, 4, 1, 0.0),
                (203.08444, 0.00014141058, ProductType.y, 4, 2, 0.0),
                (355.23398, 0.1165214, ProductType.b, 4, 1, 0.0),
                (178.12064, 0.010349626, ProductType.b, 4, 2, 0.0),
                (502.21436, 0.7401104, ProductType.y, 5, 1, 0.0),
                (470.26093, 0.055574868, ProductType.b, 5, 1, 0.0),
                (235.6341, 0.005631463, ProductType.b, 5, 2, 0.0),
                (630.3093, 0.0679749, ProductType.y, 6, 1, 0.0),
                (527.2824, 0.02713329, ProductType.b, 6, 1, 0.0),
                (264.14484, 0.002669461, ProductType.b, 6, 2, 0.0),
                (687.3308, 0.38598263, ProductType.y, 7, 1, 0.0),
                (598.3195, 0.016116723, ProductType.b, 7, 1, 0.0),
                (758.3679, 0.1706151, ProductType.y, 8, 1, 0.0),
                (655.34094, 0.007904499, ProductType.b, 8, 1, 0.0),
                (815.3894, 0.5167622, ProductType.y, 9, 1, 0.0),
                (408.19833, 0.0012127026, ProductType.y, 9, 2, 0.0),
                (783.4359, 0.01972321, ProductType.b, 9, 1, 0.0),
                (930.4163, 0.8694488, ProductType.y, 10, 1, 0.0),
                (465.7118, 0.026939344, ProductType.y, 10, 2, 0.0),
                (1029.4847, 0.27091113, ProductType.y, 11, 1, 0.0),
                (515.24603, 0.020846518, ProductType.y, 11, 2, 0.0),
                (1100.5219, 0.22043262, ProductType.y, 12, 1, 0.0),
                (550.7646, 0.0036459658, ProductType.y, 12, 2, 0.0),
                (1008.5473, 0.0029647197, ProductType.b, 12, 1, 0.0),
                (1137.5898, 0.009047425, ProductType.b, 13, 1, 0.0),
                (569.2986, 7.061393e-05, ProductType.b, 13, 2, 0.0)
            };

            for (int i = 0; i < frags.Count; i++)
            {
                var frag = frags[i];
                var readFrag = test1.MatchedFragmentIons[i];

                Assert.That(frag.mz == readFrag.Mz);
                Assert.That(frag.intensity == readFrag.Intensity);
                Assert.That(frag.ProductType == readFrag.NeutralTheoreticalProduct.ProductType);
                Assert.That(frag.fragmentNumber == readFrag.NeutralTheoreticalProduct.FragmentNumber);
                Assert.That(frag.charge == readFrag.Charge);
                //Assert.That(frag.ppm == readFrag.MassErrorPpm);
            }

            // write the library w/ the ToString method
            var writtenPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\testLibraryToString.msp");
            var str = librarySpectra.SelectMany(p => p.ToString().Split(new char[] { '\n' }));
            File.WriteAllLines(writtenPath, str);

            testLibraryWithoutDecoy.CloseConnections();

            // read the written library and make sure the results are readable
            testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { writtenPath });
            librarySpectra = testLibraryWithoutDecoy.GetAllLibrarySpectra().ToList();

            Assert.That(librarySpectra.Count == 5);
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("QSQHM[Common Variable:Oxidation on M]TEVVR", 5, out spectrum));
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("M[Common Variable:Oxidation on M]C[Common Fixed:Carbamidomethyl on C]SDSDGLAPPQHLIR", 2, out spectrum));

            testLibraryWithoutDecoy.TryGetSpectrum("ALAVDGAGKPGAEE", 2, out test1);

            Assert.AreEqual(test1.ChargeState, 2);

            for (int i = 0; i < frags.Count; i++)
            {
                var frag = frags[i];
                var readFrag = test1.MatchedFragmentIons[i];

                Assert.That(frag.mz == readFrag.Mz);
                Assert.That(frag.intensity == readFrag.Intensity);
                Assert.That(frag.ProductType == readFrag.NeutralTheoreticalProduct.ProductType);
                Assert.That(frag.fragmentNumber == readFrag.NeutralTheoreticalProduct.FragmentNumber);
                Assert.That(frag.charge == readFrag.Charge);
                //Assert.That(frag.ppm == readFrag.MassErrorPpm);
            }

            testLibraryWithoutDecoy.CloseConnections();
            File.Delete(writtenPath);
        }

        [Test]
        public static void SpectralLibrarySearchTest()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            var outputDir = Path.Combine(testDir, @"SpectralLibrarySearchTest");

            string library1 = Path.Combine(testDir, @"P16858_target.msp");
            string library2 = Path.Combine(testDir, @"P16858_decoy.msp");
            string fastaDb = Path.Combine(testDir, @"P16858.fasta");
            string spectraFile = Path.Combine(testDir, @"slicedMouse.raw");

            Directory.CreateDirectory(outputDir);

            var searchTask = new SearchTask();

            searchTask.RunTask(outputDir,
                new List<DbForTask>
                {
                    new DbForTask(library1, false),
                    new DbForTask(library2, false),
                    new DbForTask(fastaDb, false)
                },
                new List<string> { spectraFile },
                "");

            var results = File.ReadAllLines(Path.Combine(outputDir, @"AllPSMs.psmtsv"));
            var split = results[0].Split('\t');
            int ind = Array.IndexOf(split, "Normalized Spectral Angle");
            Assert.That(ind >= 0);

            var spectralAngle = double.Parse(results[1].Split('\t')[ind]);
            Assert.That(Math.Round(spectralAngle, 2) == 0.68);

            Directory.Delete(outputDir, true);
        }
    }
}
