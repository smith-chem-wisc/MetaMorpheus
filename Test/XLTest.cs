using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class XLTest
    {
        private static IndexingResults indexResults { get; set; }
        private static CommonParameters commonParameters { get; set; }
        private static XlSearchParameters xlSearchParameters { get; set; }
        private static List<Protein> proteinList { get; set; }
        private static List<Modification> variableModifications { get; set; }
        private static List<Modification> fixedModifications { get; set; }
        private static Crosslinker crosslinker { get; set; }
        private static List<PeptideWithSetModifications> digestedList { get; set; }

        [Test]
        public static void XlTestXlPosCal()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            List<DigestionMotif> motifs = new List<DigestionMotif> { new DigestionMotif("K", null, 1, null) };
            Protease protease = new Protease("New Custom Protease", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> variableModifications = new List<Modification>();

            var ye = prot.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();

            var pep = ye[0];
            Assert.AreEqual(pep.BaseSequence, "MNNNK");
            Crosslinker crosslinker = new Crosslinker();
            crosslinker.SelectCrosslinker(CrosslinkerType.DSS);
            Assert.AreEqual(crosslinker.CrosslinkerModSites, "K");
            Assert.AreEqual(Residue.GetResidue(crosslinker.CrosslinkerModSites).MonoisotopicMass, 128.09496301518999, 1e-9);
            var n = pep.Fragment(DissociationType.HCD, FragmentationTerminus.N);
            var c = pep.Fragment(DissociationType.HCD, FragmentationTerminus.C);
            Assert.AreEqual(n.Count(), 4);
            Assert.AreEqual(c.Count(), 4);
            Assert.AreEqual(c.First().NeutralMass, 146.10552769899999, 1e-6);
            var x = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep).ToArray();
            Assert.AreEqual(x[0], 5);

            var pep2 = ye[2];
            Assert.AreEqual("MNNNKQQQQ", pep2.BaseSequence);
            var n2 = pep2.Fragment(DissociationType.HCD, FragmentationTerminus.N);
            var c2 = pep2.Fragment(DissociationType.HCD, FragmentationTerminus.C);
            Assert.AreEqual(n2.Count(), 8);
            Assert.AreEqual(c2.Count(), 8);
            var x2 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep2).ToArray();
            Assert.AreEqual(x2[0], 5);

            //Test crosslinker with multiple types of mod
            var protSTC = new Protein("GASTACK", null);
            var peps = protSTC.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            var pepSTC = peps[0];
            Assert.AreEqual(pepSTC.BaseSequence, "GASTACK");
            Crosslinker crosslinker2 = new Crosslinker("ST", "C", "crosslinkerSTC", false, -18.01056, 0, 0, 0, 0, 0, 0);
            string crosslinkerModSitesAll = new string((crosslinker2.CrosslinkerModSites + crosslinker2.CrosslinkerModSites2).ToCharArray().Distinct().ToArray());
            Assert.AreEqual(crosslinkerModSitesAll, "STC");
        }

        [Test]
        public static void XlTestGenerateIntensityRanks()
        {
            double[] intensity = new double[] { 1.1, 1.1, 0.5, 3.2, 0.5, 6.0 };
            int[] rank = CrosslinkSpectralMatch.GenerateIntensityRanks(intensity);
            int[] Rank = new int[] { 4, 3, 6, 2, 5, 1 };
            Assert.AreEqual(rank, Rank);
        }

        [Test]
        public static void XlTest_BSA_DSSO()
        {
            //Generate parameters
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, dissociationType: DissociationType.EThcD,
                scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 5), precursorMassTolerance: new PpmTolerance(10));

            var xlSearchParameters = new XlSearchParameters();

            //Create databases contain two protein.
            var proteinList = new List<Protein> { new Protein("EKVLTSSAR", "Fake01"), new Protein("LSQKFPK", "Fake02") };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            var variableModifications = new List<Modification>() { mod1 };
            var fixedModifications = new List<Modification>() { mod2 };
            var localizeableModifications = new List<Modification>();

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, commonParameters, 30000, false, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            var indexedFragments = indexResults.FragmentIndex.Where(p => p != null).SelectMany(v => v).ToList();
            Assert.AreEqual(82, indexedFragments.Count);
            Assert.AreEqual(3, indexResults.PeptideIndex.Count);

            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFile();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //Generate crosslinker, which is DSSO here.
            Crosslinker crosslinker = new Crosslinker().SelectCrosslinker(CrosslinkerType.DSSO);

            CrosslinkSpectralMatch[] possiblePsms = new CrosslinkSpectralMatch[listOfSortedms2Scans.Length];
            new CrosslinkSearchEngine(possiblePsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, commonParameters, crosslinker, xlSearchParameters.RestrictToTopNHits, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, new List<string> { }).Run();

            var newPsms = possiblePsms.Where(p => p != null).ToList();
            foreach (var item in newPsms)
            {
                item.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            }

            //Test newPsms
            Assert.AreEqual(3, newPsms.Count);

            //Test Output
            var task = new XLSearchTask();
            task.WritePepXML_xl(newPsms, proteinList, null, variableModifications, fixedModifications, null, TestContext.CurrentContext.TestDirectory, "pep.XML", new List<string> { });

            //Test PsmCross.XlCalculateTotalProductMasses
            //var psmCrossAlpha = new CrosslinkSpectralMatch(digestedList[1], 0, 0, 0, listOfSortedms2Scans[0], commonParameters.DigestionParams, new List<MatchedFragmentIon>());
            //var psmCrossBeta = new CrosslinkSpectralMatch(digestedList[2], 0, 0, 0, listOfSortedms2Scans[0], commonParameters.DigestionParams, new List<MatchedFragmentIon>());
            //var linkPos = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), digestedList[1]);
            //var productMassesAlphaList = CrosslinkedPeptide.XlGetTheoreticalFragments(DissociationType.EThcD, false, crosslinker, linkPos, digestedList[2].MonoisotopicMass, digestedList[1]);
            //Assert.AreEqual(productMassesAlphaList.First().Value.Count, 50); //TO DO: The number here should be manually verified.
            File.Delete(@"singlePsms.tsv");
            File.Delete(@"pep.XML.pep.xml");
            File.Delete(@"allPsms.tsv");
        }

        [Test]
        public static void XlTest_BSA_DSS_file()
        {
            var task = Toml.ReadFile<XLSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/XLSearchTaskconfig_BSA_DSS_23747.toml"), MetaMorpheusTask.tomlConfig);
            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData"));
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/BSA.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/BSA_DSS_23747.mzML");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData")).Run();
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData"), true);
        }

        [Test]
        public static void XlTest_GenerateUserDefinedCrosslinker()
        {
            XlSearchParameters xlSearchParameters = new XlSearchParameters();
            xlSearchParameters.CrosslinkerType = CrosslinkerType.UserDefined;
            xlSearchParameters.CrosslinkerName = "CrossST-C";
            xlSearchParameters.CrosslinkerResidues = "ST";
            xlSearchParameters.CrosslinkerResidues2 = "C";
            xlSearchParameters.CrosslinkerTotalMass = -18.01056;
            var crosslinker = XLSearchTask.GenerateUserDefinedCrosslinker(xlSearchParameters);
        }

        [Test]
        public static void XlTest_DiffCrosslinkSites()
        {
            //Generate parameters
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 4));

            var xlSearchParameters = new XlSearchParameters
            {
                CrosslinkerType = CrosslinkerType.UserDefined,
                CrosslinkerName = "CrossST-C",
                CrosslinkerResidues = "ST",
                CrosslinkerResidues2 = "C",
                CrosslinkerTotalMass = -18.01056
            };

            //Create databases contain two protein.
            var proteinList = new List<Protein> { new Protein("VLTAR", "Fake01"), new Protein("LCQK", "Fake02") };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            var variableModifications = new List<Modification>() { mod1 };
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();

            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();

            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            //Generate digested peptide lists.
            List<PeptideWithSetModifications> digestedList = new List<PeptideWithSetModifications>();
            foreach (var item in proteinList)
            {
                var digested = item.Digest(commonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                digestedList.AddRange(digested);
            }

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, commonParameters, 30000, false, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFileDiffSite();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //Generate crosslinker, which is UserDefined here.
            var crosslinker = XLSearchTask.GenerateUserDefinedCrosslinker(xlSearchParameters);

            //TwoPassCrosslinkSearchEngine.Run().
            CrosslinkSpectralMatch[] possiblePsms = new CrosslinkSpectralMatch[listOfSortedms2Scans.Length];
            new CrosslinkSearchEngine(possiblePsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, commonParameters, crosslinker, xlSearchParameters.RestrictToTopNHits, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, new List<string> { }).Run();

            var newPsms = possiblePsms.Where(p => p != null).ToList();
            Assert.AreEqual(1, newPsms.Count);
        }

        /// <summary>
        /// Verifies that crosslinker is generated properly
        /// </summary>
        [Test]
        public static void CrosslinkCreateTest()
        {
            Assert.That((XLSearchTask.GenerateUserDefinedCrosslinker(new XlSearchParameters())).GetType().Equals(typeof(Crosslinker)));
        }

        [Test]
        public static void DeadendPeptideTest()
        {
            string myFileXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSSO_ETchD6010.mgf");
            string myDatabaseXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestDeadendPeptide");

            XLSearchTask xLSearchTask = new XLSearchTask()
            {
                CommonParameters = new CommonParameters(precursorMassTolerance: new PpmTolerance(51000))
            };

            XLSearchTask xLSearchTask2 = new XLSearchTask()
            {
                CommonParameters = new CommonParameters(precursorMassTolerance: new PpmTolerance(112000)),
                XlSearchParameters = new XlSearchParameters()
                {
                    XlQuench_Tris = false,
                    XlQuench_H2O = false,
                    XlQuench_NH2 = true
                }
            };
            Directory.CreateDirectory(outputFolder);

            xLSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabaseXl, false) }, new List<string> { myFileXl }, "test");
            xLSearchTask2.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabaseXl, false) }, new List<string> { myFileXl }, "test");
            Directory.Delete(outputFolder, true);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }

        /// <summary>
        /// Makes sure helper methods that generate indices function properly
        /// </summary>
        [Test]
        public static void XLSearchWithGeneratedIndices()
        {
            XLSearchTask xlSearchTask = new XLSearchTask();
            xlSearchTask.CommonParameters = new CommonParameters();
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSSO_ETchD6010.mgf");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestXLSearch");
            DbForTask db = new DbForTask(myDatabase, false);
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("TestXLSearch", xlSearchTask) };
            Directory.CreateDirectory(folderPath);

            //creates .params files if they do not exist
            xlSearchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");
            //tests .params files
            xlSearchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");

            var baseDir = Path.GetDirectoryName(db.FilePath);
            var directory = new DirectoryInfo(baseDir);
            DirectoryInfo[] directories = directory.GetDirectories();
            foreach (DirectoryInfo possibleFolder in directories)
            {
                if (File.Exists(Path.Combine(possibleFolder.FullName, "indexEngine.params")))
                {
                    File.Delete(possibleFolder.GetFiles().ElementAt(0).FullName);
                }
            }
            //tests without .params files
            xlSearchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");

            var lines = File.ReadAllLines(Path.Combine(folderPath, @"XL_Intralinks.tsv"));
            Assert.That(lines.Length == 2);
            Directory.Delete(folderPath, true);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }

        /// <summary>
        /// Tests that a crosslinker that links at proline ("P") generates the correct indices
        /// of potential crosslink sites in the sequence PEPTIDE. The indices should be positions 1 and 3
        /// </summary>
        [Test]
        public static void TestGetPossibleCrosslinkerSites()
        {
            PeptideWithSetModifications peptide = new PeptideWithSetModifications("PEPTIDE", null);
            List<int> sites = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(new char[] { 'P' }, peptide);
            Assert.That(sites.SequenceEqual(new int[] { 1, 3 }));
        }

        /// <summary>
        /// Generate and test fragments for this loop-linked peptide:
        ///    _
        ///   | |
        /// PEPTIDE
        /// </summary>
        [Test]
        public static void TestTheoreticalFragmentsLoop()
        {
            Protein p = new Protein("PEPTIDE", "");
            var peptide = p.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var loopMod = new Modification("Loop", _modificationType: "XLTest", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10000);

            var loopLocationsWithFragments = CrosslinkedPeptide.XlLoopGetTheoreticalFragments(DissociationType.HCD, loopMod, new List<int> { 3, 5 },
                peptide);

            Assert.That(loopLocationsWithFragments.Count == 1);
            var loopLocationWithFragments = loopLocationsWithFragments.First();

            Assert.That(loopLocationWithFragments.Key.Item1 == 3);
            Assert.That(loopLocationWithFragments.Key.Item2 == 5);

            var fragments = loopLocationWithFragments.Value;

            var bIons = fragments.Where(v => v.ProductType == ProductType.b).ToList();
            Assert.That(bIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 97, 226, 10537, 10652 }));

            var yIons = fragments.Where(v => v.ProductType == ProductType.y).ToList();
            Assert.That(yIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 147, 262, 10573, 10702 }));
        }

        /// <summary>
        /// Generate and test fragments for loop-linked peptides with a mod placed at several different locations
        /// </summary>
        [Test]
        public static void TestTheoreticalLoopFragmentsWithMod()
        {
            ModificationMotif.TryGetMotif("T", out var tMotif);
            Modification phospho = new Modification(_originalId: "Phospho", _modificationType: "Mod", _locationRestriction: "Anywhere.", _monoisotopicMass: 79.98, _target: tMotif);
            int[] modPositions = new int[] { 2, 3, 4, 5, 6 };

            foreach (var modPosition in modPositions)
            {
                Dictionary<int, List<Modification>> oneBasedMods = new Dictionary<int, List<Modification>>
                    { { modPosition, new List<Modification> { phospho } } };
                Protein p = new Protein("PTTTTTE", "", oneBasedModifications: oneBasedMods);
                var peptide = p.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>())
                    .Where(v => v.AllModsOneIsNterminus.Count == 1).First();

                ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
                var loopMod = new Modification("Loop", _modificationType: "XLTest", _target: motif,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: 10000);

                var loopLocationsWithFragments = CrosslinkedPeptide.XlLoopGetTheoreticalFragments(DissociationType.HCD,
                    loopMod, new List<int> { 3, 5 },
                    peptide);

                Assert.That(loopLocationsWithFragments.Count == 1);
                var loopLocationWithFragments = loopLocationsWithFragments.First();

                Assert.That(loopLocationWithFragments.Key.Item1 == 3);
                Assert.That(loopLocationWithFragments.Key.Item2 == 5);

                var fragments = loopLocationWithFragments.Value;

                var bIons = fragments.Where(v => v.ProductType == ProductType.b).ToList();
                var yIons = fragments.Where(v => v.ProductType == ProductType.y).ToList();

                if (modPosition == 2)
                {
                    //             _
                    //            | |
                    // PT[Phospho]TTTTE
                    Assert.That(bIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 97, 278, 10581, 10682 }));
                    Assert.That(yIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 147, 248, 10551, 10732 }));
                }
                else if (modPosition == 3)
                {
                    //    __________
                    //   |          |
                    // PTT[Phospho]TTTE
                    Assert.That(bIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 97, 198, 10581, 10682 }));
                    Assert.That(yIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 147, 248, 10631, 10732 }));
                }
                else if (modPosition == 4)
                {
                    //    __________
                    //   |          |
                    // PTTT[Phospho]TTE
                    Assert.That(bIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 97, 198, 10581, 10682 }));
                    Assert.That(yIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 147, 248, 10631, 10732 }));
                }
                else if (modPosition == 5)
                {
                    //    _
                    //   | |
                    // PTTTT[Phospho]TE
                    Assert.That(bIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 97, 198, 10581, 10682 }));
                    Assert.That(yIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 147, 248, 10631, 10732 }));
                }
                else if (modPosition == 6)
                {
                    //    _
                    //   | |
                    // PTTTTT[Phospho]E
                    Assert.That(bIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 97, 198, 10501, 10682 }));
                    Assert.That(yIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 147, 328, 10631, 10732 }));
                }
            }
        }

        /// <summary>
        /// Generate and test fragments for this dead-end peptide (quenched with tris):
        ///   
        ///   |
        /// PEPTIDE[Methyl]
        /// </summary>
        [Test]
        public static void TestDeadendTris()
        {
            Protein protein = new Protein("PEPTIDE", "");
            var csms = new CrosslinkSpectralMatch[1];

            // generate the scan with the deadend mod peptide's fragments
            var scans = new Ms2ScanWithSpecificMass[1];
            ModificationMotif.TryGetMotif("T", out var motif);
            var crosslinker = new Crosslinker("T", "T", "test", false, 100, 0, 0, 0, 0, 0, 50);
            Modification deadend = new Modification("TestId", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: crosslinker.DeadendMassTris, _modificationType: "Test");

            var deadendPeptide = protein.Digest(new DigestionParams(), new List<Modification> { deadend }, new List<Modification>()).First();

            double[] mz = deadendPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both).Select(p => p.NeutralMass.ToMz(1)).OrderBy(v => v).ToArray();
            double[] intensities = new[] { 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

            MzSpectrum spectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan sc = new MsDataScan(spectrum, 1, 2, true, Polarity.Positive, 1, spectrum.Range, "",
                MZAnalyzerType.Orbitrap, 12, 1.0, null, null);
            scans[0] = new Ms2ScanWithSpecificMass(sc, deadendPeptide.MonoisotopicMass.ToMz(2), 2, "", new CommonParameters());

            // search the data with the peptide WITHOUT the deadend mod annotated in the search database.
            // the search engine should be able to correctly identify the deadend mod on T
            var indexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { protein }, new List<Modification>(), new List<Modification>(), 0, DecoyType.None,
                new CommonParameters(), 1000, false, new List<FileInfo>(), new List<string>()).Run();

            new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, 0, new CommonParameters(), crosslinker,
                false, 0, false, false, true, new List<string>()).Run();

            CrosslinkSpectralMatch csm = csms.First();
            Assert.That(csm.CrossType == PsmCrossType.DeadEndTris);
            Assert.That(csm.MatchedFragmentIons.Count == 12);
        }

        /// <summary>
        /// Generate and test fragments for this crosslinked peptide pair that was crosslinked
        /// with a non-cleavable crosslinker:
        ///  PRLTEIN
        ///   |
        /// PEPTIDE
        /// </summary>
        [Test]
        public static void TestTheoreticalFragmentsNonCleavableCrosslink()
        {
            Crosslinker c = new Crosslinker("P", "R", "Test", false, 1000, 0, 0, 1000, 5, 5, 5);

            Protein p1 = new Protein("PEPTIDE", "");
            var alphaPeptide = p1.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            Protein p2 = new Protein("PRLTEIN", "");
            var betaPeptide = p2.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(v => v.MissedCleavages == 1).First();

            var theoreticalCrosslinkFragments = CrosslinkedPeptide.XlGetTheoreticalFragments(DissociationType.HCD,
                c, new List<int> { 3 }, betaPeptide.MonoisotopicMass, alphaPeptide).ToList();

            Assert.That(theoreticalCrosslinkFragments.Count == 1);
            var loopLocationWithFragments = theoreticalCrosslinkFragments.First();

            Assert.That(loopLocationWithFragments.Item1 == 3);

            var fragments = loopLocationWithFragments.Item2;

            var bIons = fragments.Where(v => v.ProductType == ProductType.b).ToList();
            Assert.That(bIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 97, 226, 2164, 2265, 2378, 2493 }));

            var yIons = fragments.Where(v => v.ProductType == ProductType.y).ToList();
            Assert.That(yIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 147, 262, 375, 476, 2414, 2543 }));
        }

        /// <summary>
        /// Generate and test fragments for this crosslinked peptide pair that was crosslinked
        /// with a cleavable crosslinker:
        ///   PROTEIN
        ///   |
        /// PEPTIDE
        /// </summary>
        [Test]
        public static void TestTheoreticalFragmentsCleavableCrosslink()
        {
            Crosslinker c = new Crosslinker("P", "R", "Test", true, 1000, 15, 25, 1000, 5, 5, 5);

            Protein p1 = new Protein("PEPTIDE", "");
            var alphaPeptide = p1.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();

            Protein p2 = new Protein("PRLTEIN", "");
            var betaPeptide = p2.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(v => v.MissedCleavages == 1).First();

            var theoreticalCrosslinkFragments = CrosslinkedPeptide.XlGetTheoreticalFragments(DissociationType.HCD,
                c, new List<int> { 3 }, 10000, alphaPeptide).ToList();

            Assert.That(theoreticalCrosslinkFragments.Count == 1);

            // cleaved fragments
            var linkLocationWithFragments = theoreticalCrosslinkFragments[0];

            Assert.That(linkLocationWithFragments.Item1 == 3);
            var fragmentsWithCleavedXlPieces = linkLocationWithFragments.Item2;

            var bIons = fragmentsWithCleavedXlPieces.Where(v => v.ProductType == ProductType.b).ToList();
            Assert.That(bIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 97, 226, 338, 439, 552, 667, 348, 449, 562, 677 }));

            var yIons = fragmentsWithCleavedXlPieces.Where(v => v.ProductType == ProductType.y).ToList();
            Assert.That(yIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 147, 262, 375, 476, 588, 717, 598, 727 }));

            var signatureIons = fragmentsWithCleavedXlPieces.Where(v => v.ProductType == ProductType.M).ToList();
            Assert.That(signatureIons.Count == 2);
            Assert.That(signatureIons.Select(v => (int)v.NeutralMass).SequenceEqual(new int[] { 814, 824 }));
        }

        [Test]
        public static void TestWriteToPercolator()
        {
            XLSearchTask xlst = new XLSearchTask()
            {
                XlSearchParameters = new XlSearchParameters
                {
                    WriteOutputForPercolator = true
                }
            };

            string myFileXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSSO_ETchD6010.mgf");
            string myDatabaseXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestXLSearch");
            DbForTask db = new DbForTask(myDatabaseXl, false);

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("TestPercolator", xlst) };

            var engine = new EverythingRunnerEngine(taskList, new List<string> { myFileXl }, new List<DbForTask> { db }, outputFolder);
            engine.Run();

            var results = Path.Combine(outputFolder, @"TestPercolator\XL_Intralinks_Percolator.txt");
            var lines = File.ReadAllLines(results);
            Assert.That(lines[0].Equals("SpecId\tLabel\tScannr\tScore\tdScore\tNormRank\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProtein"));
            Assert.That(lines[1].Equals("T-1-30.6190992666667\t1\t1\t20.0641008915522\t0\t7\t3\t1994.05202313843\t0.664979354397676\t7\t9\t16\t-.EKVLTSSAR2--LSQKFPK4.-\t3336842(211)\t3336842(245)"));
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestWriteNonSingleCross()
        {
            XLSearchTask xlst = new XLSearchTask();
            Protein protein = new Protein("PEPTIDE", "");
            var csms = new CrosslinkSpectralMatch[1];

            // generate the scan with the deadend mod peptide's fragments
            var scans = new Ms2ScanWithSpecificMass[1];
            ModificationMotif.TryGetMotif("T", out var motif);
            var crosslinker = new Crosslinker("T", "T", "test", false, 100, 0, 0, 0, 0, 0, 50);
            Modification deadend = new Modification("TestId", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: crosslinker.DeadendMassTris, _modificationType: "Test");

            var deadendPeptide = protein.Digest(new DigestionParams(), new List<Modification> { deadend }, new List<Modification>()).First();

            double[] mz = deadendPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both).Select(p => p.NeutralMass.ToMz(1)).OrderBy(v => v).ToArray();
            double[] intensities = new[] { 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

            MzSpectrum spectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan sc = new MsDataScan(spectrum, 1, 2, true, Polarity.Positive, 1, spectrum.Range, "",
                MZAnalyzerType.Orbitrap, 12, 1.0, null, null);
            scans[0] = new Ms2ScanWithSpecificMass(sc, deadendPeptide.MonoisotopicMass.ToMz(2), 2, "", new CommonParameters());

            var indexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { protein }, new List<Modification>(), new List<Modification>(), 0, DecoyType.None,
                new CommonParameters(), 1000, false, new List<FileInfo>(), new List<string>()).Run();

            new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, 0, new CommonParameters(), crosslinker, false, 0, false, false, true, new List<string>()).Run();

            csms[0].SetFdrValues(0, 0, 0.1, 0, 0, 0, 0, 0, 0, false);

            xlst.WritePepXML_xl(csms.ToList(), new List<Protein>(), "", new List<Modification> { deadend }, new List<Modification> { deadend }, new List<string>(), TestContext.CurrentContext.TestDirectory, "test", new List<string>());
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test.pep.XML"));
        }
    }

    internal class XLTestDataFile : MsDataFile
    {
        //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
        public XLTestDataFile() : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 1994.05.ToMz(3), 846.4963.ToMz(1), 1004.495.ToMz(1), 1093.544.ToMz(1), 1043.561.ToMz(1) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };

            var mz2 = new double[] { 100, 201.1234, 244.1656, 391.2340, 420.2201, 521.2678, 634.3519, 889.965, 1044.568, 1094.551, 1279.671, 1378.74, 1491.824 };
            var intensities2 = new double[] { 100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1.0, null, "scan=2", 1994.05.ToMz(3),
                3, 1, 1994.05.ToMz(3), 2, DissociationType.HCD, 1, 1994.05.ToMz(3)));

            var mz3 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
            var intensities3 = new double[] { 100, 1, 1, 1 };
            var MassSpectrum3 = new MzSpectrum(mz3, intensities3, false);
            ScansHere.Add(new MsDataScan(MassSpectrum3, 3, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=3", 846.4963.ToMz(1),
                1, 1, 846.4963.ToMz(1), 2, DissociationType.HCD, 1, 846.4963.ToMz(1)));

            var mz4 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
            var intensities4 = new double[] { 100, 1, 1, 1 };
            var MassSpectrum4 = new MzSpectrum(mz4, intensities4, false);
            ScansHere.Add(new MsDataScan(MassSpectrum4, 4, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=4", 1004.491.ToMz(1),
                1, 1, 1004.491.ToMz(1), 2, DissociationType.HCD, 1, 1004.491.ToMz(1)));

            Scans = ScansHere.ToArray();
        }

        public string FilePath
        {
            get
            {
                return "XLTestDataFile";
            }
        }

        public string Name
        {
            get
            {
                return "XLTestDataFile";
            }
        }

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            Scans[0] = new MsDataScan(massSpectrum, Scans[0].OneBasedScanNumber, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, null, Scans[0].NativeId);
        }
    }

    internal class XLTestDataFileDiffSite : MsDataFile
    {
        //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
        public XLTestDataFileDiffSite() : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 100, 1030.5956.ToMz(1) };
            var intensities1 = new double[] { 100, 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };

            var mz2 = new double[] { 100, 147.1128, 175.119, 213.1598, 246.1561, 275.1714, 757.4388, 786.4541, 819.4504, 857.4912, 885.4974, 918.5189, 932.5345 };
            var intensities2 = new double[] { 100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1.0, null, "scan=2", 1030.5956.ToMz(1),
                1, 1, 1030.5956.ToMz(1), 2, DissociationType.HCD, 1, 1030.5956.ToMz(1)));

            Scans = ScansHere.ToArray();
        }

        public string FilePath
        {
            get
            {
                return "XLTestDataFileDiffSite";
            }
        }

        public string Name
        {
            get
            {
                return "XLTestDataFileDiffSite";
            }
        }

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            Scans[0] = new MsDataScan(massSpectrum, Scans[0].OneBasedScanNumber, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, null, Scans[0].NativeId);
        }
    }
}
