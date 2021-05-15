using NUnit.Framework;
using System.IO;
using System;
using System.Linq;
using EngineLayer;
using MzLibUtil;
using TaskLayer;
using System.Collections.Generic;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using MassSpectrometry;

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
        [Test]
        public static void TestDecoyOnTheFlyBottomUp()
        {
            //test decoy generation without modification
            List<Modification> mods1 = new List<Modification>();
            CommonParameters commonParameters1 = new CommonParameters();
            PeptideWithSetModifications target1 = new Protein("PEPTIDE", "TARGET").Digest(commonParameters1.DigestionParams, mods1, mods1).First();
            int[] newAAlocations1 = new int[target1.BaseSequence.Length];
            PeptideWithSetModifications decoy1 = DecoyOnTheFly.GetReverseDecoyFromTarget(target1, newAAlocations1);
            Assert.That(decoy1.FullSequence == "EDITPEP");
            Assert.That(decoy1.Protein.IsDecoy);

            //test decoy generation with modification 
            List<Modification> mods2 = new List<Modification>();
            ModificationMotif.TryGetMotif("D", out var motif);
            Modification m = new Modification("TEST", "", "OK", null, motif, "Anywhere.", null, 20);
            mods2.Add(m);
            PeptideWithSetModifications target2 = new Protein("PEPTIDE", "TARGET").Digest(commonParameters1.DigestionParams, mods2, mods1).First();
            int[] newAAlocations2 = new int[target2.BaseSequence.Length];
            PeptideWithSetModifications decoy2 = DecoyOnTheFly.GetReverseDecoyFromTarget(target2, newAAlocations2);
            Assert.That(decoy2.FullSequence == "ED[OK:TEST on D]ITPEP");
            Assert.That(decoy2.Protein.IsDecoy);

            //test decoy generation with modification in the digestion motif
            List<Modification> mods3 = new List<Modification>();
            ModificationMotif.TryGetMotif("K", out var motif3);
            Modification m3 = new Modification("TEST", "", "OK", null, motif3, "Anywhere.", null, 20);
            mods3.Add(m3);
            PeptideWithSetModifications target3 = new Protein("PEPTIDEK", "TARGET").Digest(commonParameters1.DigestionParams, mods3, mods1).First();
            int[] newAAlocations3 = new int[target3.BaseSequence.Length];
            PeptideWithSetModifications decoy3 = DecoyOnTheFly.GetReverseDecoyFromTarget(target3, newAAlocations3);
            Assert.That(decoy3.FullSequence == "EDITPEPK[OK:TEST on K]");
            Assert.That(decoy3.Protein.IsDecoy);

        }

        [Test]
        public static void TestDecoyOnTheFlyTopDown()
        {

            CommonParameters CommonParameters = new CommonParameters(
              digestionParams: new DigestionParams(protease: "top-down"),
              scoreCutoff: 1,
              assumeOrphanPeaksAreZ1Fragments: false);

            MetaMorpheusTask.DetermineAnalyteType(CommonParameters);

            // test output file name (should be proteoform and not peptide)
            Assert.That(GlobalVariables.AnalyteType == "Proteoform");

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein>
            {
                new Protein("MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE", "P40312")
            };

            //test decoy generation without modification
            List<Modification> mods1 = new List<Modification>();
            PeptideWithSetModifications target1 = proteinList.First().Digest(CommonParameters.DigestionParams, mods1, mods1).First();
            int[] newAAlocations1 = new int[target1.BaseSequence.Length];
            PeptideWithSetModifications decoy1 = DecoyOnTheFly.GetReverseDecoyFromTarget(target1, newAAlocations1);
            Assert.That(decoy1.FullSequence == "ENLLYYAVGLMLIALIVVLTGSGKSQNESTSVKEVSVRESTKDVDGIYLGKLLRLAEDSHGIDVFSETADQGGLDMIIEDGGPHEDKFQSVDYVKDDIIIWFNEPGNHEAVEQYSYVKPM");
            Assert.That(decoy1.Protein.IsDecoy);

            //test decoy generation with modification 
            List<Modification> mods2 = new List<Modification>();
            ModificationMotif.TryGetMotif("A", out var motif);
            Modification m = new Modification("TEST", "", "OK", null, motif, "Anywhere.", null, 20);
            mods2.Add(m);
            PeptideWithSetModifications target2 = proteinList.First().Digest(CommonParameters.DigestionParams, mods2, mods1).First();
            int[] newAAlocations2 = new int[target2.BaseSequence.Length];
            PeptideWithSetModifications decoy2 = DecoyOnTheFly.GetReverseDecoyFromTarget(target2, newAAlocations2);
            Assert.That(decoy2.FullSequence == "ENLLYYA[OK:TEST on A]VGLMLIA[OK:TEST on A]LIVVLTGSGKSQNESTSVKEVSVRESTKDVDGIYLGKLLRLA[OK:TEST on A]EDSHGIDVFSETA[OK:TEST on A]DQGGLDMIIEDGGPHEDKFQSVDYVKDDIIIWFNEPGNHEA[OK:TEST on A]VEQYSYVKPM");
            Assert.That(decoy2.Protein.IsDecoy);


        }

    }
}
