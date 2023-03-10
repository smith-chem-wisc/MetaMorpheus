using NUnit.Framework;
using System.IO;
using System;
using System.Linq;
using EngineLayer;
using TaskLayer;
using System.Collections.Generic;
using Proteomics.Fragmentation;
using System.Globalization;
using Proteomics;
using MassSpectrometry;
using EngineLayer.ClassicSearch;

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
            Assert.That(Math.Round(spectralAngle, 2) == 0.82);

            Directory.Delete(outputDir, true);
        }

        /// <summary>
        /// Test ensures peptide FDR is calculated and that it doesn't output PSM FDR results
        /// </summary>
        [Test]
        public static void AnotherSpectralLibrarySearchTest()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            string myFile = Path.Combine(testDir, @"slicedMouse.raw");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters(maxThreadsToUsePerFile: 1, scoreCutoff: 1);
            MsDataFile myMsDataFile = myFileManager.LoadFile(myFile, commonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("VIHDNFGIVEGLMTTVHAITATQK", "P16858") };

            string targetSpectralLibrary = Path.Combine(testDir, @"P16858_target.msp");
            string decoySpectralLibrary = Path.Combine(testDir, @"P16858_decoy.msp");

            List<string> specLibs = new List<string> { targetSpectralLibrary, decoySpectralLibrary };

            SpectralLibrary sl = new SpectralLibrary(specLibs);

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters,null, sl, new List<string>(), writeSpectralLibrary).Run();
             
            // Single search mode
            Assert.AreEqual(7, allPsmsArray.Length);
            Assert.IsTrue(allPsmsArray[5].Score > 38);
            Assert.AreEqual("VIHDNFGIVEGLMTTVHAITATQK", allPsmsArray[5].BaseSequence);

            SpectralLibrarySearchFunction.CalculateSpectralAngles(sl, allPsmsArray, listOfSortedms2Scans, commonParameters);
            Assert.That(allPsmsArray[5].SpectralAngle, Is.EqualTo(0.82).Within(0.01));
        }


        /// <summary>
        /// Test ensures peptide FDR is calculated and that it doesn't output PSM FDR results
        /// </summary>
        [Test]
        public static void AnotherSpectralLibrarySearchTestDecoy()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            string myFile = Path.Combine(testDir, @"slicedMouse.raw");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters(maxThreadsToUsePerFile: 1, scoreCutoff: 1);
            MsDataFile myMsDataFile = myFileManager.LoadFile(myFile, commonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("QTATIAHVTTMLGEVIGFNDHIVK", "P16858") };

            string targetSpectralLibrary = Path.Combine(testDir, @"P16858_target.msp");
            string decoySpectralLibrary = Path.Combine(testDir, @"P16858_decoy.msp");

            List<string> specLibs = new List<string> { targetSpectralLibrary, decoySpectralLibrary };

            SpectralLibrary sl = new SpectralLibrary(specLibs);

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters, null, sl, new List<string>(), writeSpectralLibrary).Run();

            // Single search mode
            Assert.AreEqual(7, allPsmsArray.Length);
            Assert.IsTrue(allPsmsArray[5].Score > 38);
            Assert.AreEqual("VIHDNFGIVEGLMTTVHAITATQK", allPsmsArray[5].BaseSequence);
            Assert.IsTrue(allPsmsArray[5].IsDecoy);

            SpectralLibrarySearchFunction.CalculateSpectralAngles(sl, allPsmsArray, listOfSortedms2Scans, commonParameters);
            Assert.That(allPsmsArray[5].SpectralAngle, Is.EqualTo(0.69).Within(0.01));

        }

        [Test]
        public static void SpectralLibraryReaderTestNeutralLoss()
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\spectralLibraryNeutralLossTest.msp");

            var testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { path });
            var librarySpectra = testLibraryWithoutDecoy.GetAllLibrarySpectra().ToList();
            Assert.That(librarySpectra[0].MatchedFragmentIons[9].NeutralTheoreticalProduct.NeutralLoss == 97.976895573);
            Assert.That(librarySpectra[0].MatchedFragmentIons[11].NeutralTheoreticalProduct.NeutralLoss == 97.976895573);
            Assert.That(librarySpectra[0].MatchedFragmentIons[13].NeutralTheoreticalProduct.NeutralLoss == 97.976895573);
            Assert.That(librarySpectra[0].MatchedFragmentIons[23].NeutralTheoreticalProduct.NeutralLoss == 97.976895573);

            testLibraryWithoutDecoy.TryGetSpectrum("ASVSELAC[Common Fixed:Carbamidomethyl on C]IYSALILHDDEVTVTEDKINALIKAAGVNVEPFWPGLFAKALANVNIGSLIC[Common Fixed:Carbamidomethyl on C]NVGAGGPAPAAGAAPAGGPAPSTAAAPAEEKKVEAKKEES[Common Biological:Phosphorylation on S]EES[Common Biological:Phosphorylation on S]DDDMGFGLFD", 11, out var test1);
            Assert.AreEqual(test1.ChargeState, 11);
            var frags = new List<(double mz, double intensity, ProductType ProductType, int fragmentNumber, int charge, double neutralLoss)>
            {
                 (474.22253552109225, 0.12327031966337068, ProductType.b, 5, 1, 0.0),
                (598.2911881692844, 0.17637685313108434, ProductType.y, 5, 1, 0.0),
                (587.3078660257297, 0.26448833020105406, ProductType.b, 6, 1, 0.0),
                (655.3122852061485, 0.2573965224084584, ProductType.y, 6, 1, 0.0),
                (786.3547054957799, 0.5973180470179946, ProductType.y, 7, 1, 0.0),
                (901.3828669653562, 0.1394747476004537, ProductType.y, 8, 1, 0.0),
                (931.4665832519531, 0.1434341615825174, ProductType.b, 9, 1, 0.0),
                (1298.4420166015627, 0.04799051344737827, ProductType.y, 11, 1, 0.0),
                (1194.1036493343909, 0.16497897552735782, ProductType.b, 22, 2, 0.0),
                (1412.6153999317123, 0.059408140384558945, ProductType.y, 24, 2, 97.976895573),

                (1610.16612195462, 0.06598126929193648, ProductType.y, 27, 2, 0.0),
                (1561.176395469285, 0.3655719849041352, ProductType.y, 27, 2, 97.976895573),
                (1596.6949731527957, 0.1196248007420883, ProductType.y, 28, 2,97.976895573),
                (1632.2140387967154, 0.1134099785167373, ProductType.y, 29, 2, 97.976895573),
                (1239.8634635821659, 0.19969145903784855, ProductType.y, 33, 3, 0.0),
                (1207.2033996967743, 0.4866798941848081, ProductType.y, 33, 3, 97.976895573),
                (1810.301461311722, 0.2822226465462767, ProductType.y, 33, 2, 97.976895573),
                (1282.2420159220214, 0.06365769264631967, ProductType.y, 36, 3, 97.976895573),
                (1333.9057047796402, 0.04957330666530949, ProductType.y, 37, 3, 0.0),
                (1042.7038942258248, 0.10714201465676373, ProductType.y, 39, 4, 0.0),
                (1389.9361001454733, 0.20027716863748268, ProductType.y, 39, 3, 0.0),
                (1018.2100179729548, 0.23609364592477405, ProductType.y, 39, 4, 97.976895573),
                (1357.2775984749799, 1, ProductType.y, 39, 3, 97.976895573),
                (1503.354187325224, 0.060227438209151066, ProductType.y, 45, 3, 97.976895573)
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
                Assert.That(frag.neutralLoss == readFrag.NeutralTheoreticalProduct.NeutralLoss);
            }

            // write the library w/ the ToString method
            var writtenPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\testLibraryToString.msp");
            var str = librarySpectra.SelectMany(p => p.ToString().Split(new char[] { '\n' }));
            File.WriteAllLines(writtenPath, str);

            testLibraryWithoutDecoy.CloseConnections();

            // read the written library and make sure the results are readable
            testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { writtenPath });
            librarySpectra = testLibraryWithoutDecoy.GetAllLibrarySpectra().ToList();

            Assert.That(librarySpectra.Count == 1);

            testLibraryWithoutDecoy.TryGetSpectrum("ASVSELAC[Common Fixed:Carbamidomethyl on C]IYSALILHDDEVTVTEDKINALIKAAGVNVEPFWPGLFAKALANVNIGSLIC[Common Fixed:Carbamidomethyl on C]NVGAGGPAPAAGAAPAGGPAPSTAAAPAEEKKVEAKKEES[Common Biological:Phosphorylation on S]EES[Common Biological:Phosphorylation on S]DDDMGFGLFD", 11, out test1);

            Assert.AreEqual(test1.ChargeState, 11);
            double maxOfIntensity = frags.Select(p => p.intensity).ToList().Max();
            for (int i = 0; i < frags.Count; i++)
            {
                var frag = frags[i];
                var readFrag = test1.MatchedFragmentIons[i];

                Assert.That(frag.mz == readFrag.Mz);
                Assert.That((frag.intensity / maxOfIntensity) == readFrag.Intensity);
                Assert.That(frag.ProductType == readFrag.NeutralTheoreticalProduct.ProductType);
                Assert.That(frag.fragmentNumber == readFrag.NeutralTheoreticalProduct.FragmentNumber);
                Assert.That(frag.charge == readFrag.Charge);
                Assert.That(frag.neutralLoss == readFrag.NeutralTheoreticalProduct.NeutralLoss);
            }

            testLibraryWithoutDecoy.CloseConnections();
            File.Delete(writtenPath);
        }

        [Test]
        public static void SpectralLibraryReaderTest_pDeep()
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\yeast2fake_pdeep_lib.msp");

            var testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { path });
            var librarySpectra = testLibraryWithoutDecoy.GetAllLibrarySpectra().ToList();

            Assert.That(librarySpectra.Count == 5);
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("VGIVPGEVIAPGM[Common Variable:Oxidation on M]R", 3, out var spectrum1));
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("C[Common Fixed:Carbamidomethyl on C]TSC[Common Fixed:Carbamidomethyl on C]NGQGIKFVTR", 3, out var spectrum2));
            Assert.AreEqual(spectrum2.PrecursorMz, 543.2608252287667);
            Assert.AreEqual(spectrum2.RetentionTime, 2789.812255859375);

            testLibraryWithoutDecoy.TryGetSpectrum("YHPDKNPSEEAAEK", 3, out var test1);

            Assert.AreEqual(test1.PrecursorMz, 538.9179945621);
            Assert.AreEqual(test1.RetentionTime, 1361.375244140625);
            Assert.AreEqual(test1.ChargeState, 3);

            var frags = new List<(double mz, double intensity, ProductType ProductType, int fragmentNumber, int charge, double ppm)>
            {
                (301.1295080000, 10000.0, ProductType.b, 2, 1, 0.0),
                (657.8122378432, 1102.3, ProductType.y, 12, 2, 0.0),
                (974.4425296863, 1476.0, ProductType.y, 9, 1, 0.0),
                (860.3996026863, 7228.2, ProductType.y, 8, 1, 0.0),
                (763.3468386863, 1201.9, ProductType.y, 7, 1, 0.0),
                (418.2296246863, 1178.9, ProductType.y, 4, 1, 0.0),
                (347.1925106863, 1042.8, ProductType.y, 3, 1, 0.0)
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
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("VGIVPGEVIAPGM[Common Variable:Oxidation on M]R", 3, out var spectrum3));
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("C[Common Fixed:Carbamidomethyl on C]TSC[Common Fixed:Carbamidomethyl on C]NGQGIKFVTR", 3, out var spectrum4));
            Assert.AreEqual(spectrum4.PrecursorMz, 543.2608252287667);
            Assert.AreEqual(spectrum4.RetentionTime, 2789.812255859375);


            testLibraryWithoutDecoy.TryGetSpectrum("YHPDKNPSEEAAEK", 3, out test1);

            Assert.AreEqual(test1.PrecursorMz, 538.9179945621);
            Assert.AreEqual(test1.RetentionTime, 1361.375244140625);
            Assert.AreEqual(test1.ChargeState, 3);
            double maxOfIntensity = frags.Select(p => p.intensity).ToList().Max();
            for (int i = 0; i < frags.Count; i++)
            {
                var frag = frags[i];
                var readFrag = test1.MatchedFragmentIons[i];

                Assert.That(frag.mz == readFrag.Mz);
                Assert.That((frag.intensity/ maxOfIntensity) == readFrag.Intensity);
                Assert.That(frag.ProductType == readFrag.NeutralTheoreticalProduct.ProductType);
                Assert.That(frag.fragmentNumber == readFrag.NeutralTheoreticalProduct.FragmentNumber);
                Assert.That(frag.charge == readFrag.Charge);
                //Assert.That(frag.ppm == readFrag.MassErrorPpm);
            }

            testLibraryWithoutDecoy.CloseConnections();
            File.Delete(writtenPath);
        }

    }
    
}
