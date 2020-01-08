using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.FdrAnalysis;
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
using System.IO.Compression;
using System.Linq;
using System.Text;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class XLTest
    {
        [Test]
        public static void XlTestCrosslinker()
        {
            var crosslinker = new Crosslinker("K", "K", "testCrosslinker", true, "CID|HCD", 100, 25, 60, 100, 0, 0, 50);
            var crosslinkerName = crosslinker.ToString();
            Assert.That(crosslinkerName == "testCrosslinker");
            var crosslinkerString = crosslinker.ToString(true);
            Assert.That(crosslinkerString == "testCrosslinker\tK\tK\tTrue\tCID|HCD|\t100\t25\t60\t0\t0\t50");
        }

        [Test]
        public static void XlTestXlPosCal()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            List<DigestionMotif> motifs = new List<DigestionMotif> { new DigestionMotif("K", null, 1, null), new DigestionMotif("R", null, 1, null) };
            Protease protease = new Protease("New Custom Protease", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, maxMissedCleavages: 4, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> variableModifications = new List<Modification>();

            var ye = prot.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();

            var pep = ye[0];
            Assert.AreEqual(pep.BaseSequence, "MNNNK");
            Crosslinker crosslinker = GlobalVariables.Crosslinkers.Where(p => p.CrosslinkerName == "DSS").First();
            Assert.AreEqual(crosslinker.CrosslinkerModSites, "K");
            Assert.AreEqual(Residue.GetResidue(crosslinker.CrosslinkerModSites).MonoisotopicMass, 128.09496301518999, 1e-9);
            var n = pep.Fragment(DissociationType.HCD, FragmentationTerminus.N);
            var c = pep.Fragment(DissociationType.HCD, FragmentationTerminus.C);
            Assert.AreEqual(n.Count(), 4);
            Assert.AreEqual(c.Count(), 4);
            Assert.AreEqual(c.First().NeutralMass, 146.10552769899999, 1e-6);
            var x = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep, digestionParams.InitiatorMethionineBehavior, false);
            Assert.That(x == null);

            var pep2 = ye[2];
            Assert.AreEqual("MNNNKQQQQ", pep2.BaseSequence);
            var n2 = pep2.Fragment(DissociationType.HCD, FragmentationTerminus.N);
            var c2 = pep2.Fragment(DissociationType.HCD, FragmentationTerminus.C);
            Assert.AreEqual(n2.Count(), 8);
            Assert.AreEqual(c2.Count(), 8);
            var x2 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep2, digestionParams.InitiatorMethionineBehavior, false);
            Assert.AreEqual(x2[0], 5);

            //TestXLPosCal on peptide with modification.
            var prot_mod = new Protein("KNNNKQRKQQK", null);
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Oxidation of K", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            var pep_mod = prot_mod.Digest(digestionParams, new List<Modification>(), new List<Modification> { mod1 }).ToList();

            var pep3 = pep_mod.Where(p => p.FullSequence == "KNNNK[Common Variable:Oxidation on K]").First();
            Assert.That(pep3.FullSequence == "KNNNK[Common Variable:Oxidation on K]");
            var x3_f = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep3, digestionParams.InitiatorMethionineBehavior, true);
            var x3_t = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep3, digestionParams.InitiatorMethionineBehavior, false);
            Assert.That(x3_f.Count() == 1 && x3_f[0] == 1);
            Assert.That(x3_t.Count() == 1);

            var pep4 = pep_mod.Where(p => p.FullSequence == "NNNK[Common Variable:Oxidation on K]QRKQQK").First();
            Assert.That(pep4.FullSequence == "NNNK[Common Variable:Oxidation on K]QRKQQK");
            var x4_f = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep4, digestionParams.InitiatorMethionineBehavior, true);
            var x4_t = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep4, digestionParams.InitiatorMethionineBehavior, false);
            //Both 'K' are crosslinked becuase the last 'K' is at protein C terminal
            Assert.That(x4_f[0] == 7 && x4_f[1] == 10 && x4_f.Count() == 2);
            Assert.That(x4_t[0] == 7 && x4_t[1] == 10 && x4_t.Count() == 2);

            var pep5 = pep_mod.Where(p => p.FullSequence == "KQQK").First();
            Assert.That(pep5.FullSequence == "KQQK");
            var x5_f = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep5, digestionParams.InitiatorMethionineBehavior, true);
            var x5_t = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep5, digestionParams.InitiatorMethionineBehavior, false);
            //Both 'K' are crosslinked becuase the last 'K' is at protein C terminal
            Assert.That(x5_f[0] == 1 && x5_f.Count() == 2);
            Assert.That(x5_t[0] == 1 && x5_t.Count() == 2);

            var pep6 = pep_mod.Where(p => p.FullSequence == "KNNNK").First();
            Assert.That(pep6.FullSequence == "KNNNK");
            var x6_f = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep6, digestionParams.InitiatorMethionineBehavior, true);
            var x6_t = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep6, digestionParams.InitiatorMethionineBehavior, false);
            //Both 'K' are crosslinked becuase the last 'K' is at protein C terminal
            Assert.That(x6_f[0] == 1 && x6_f.Count() == 2);
            Assert.That(x6_t[0] == 1 && x6_t.Count() == 1);

            //Test crosslinker with multiple types of mod
            var protSTC = new Protein("GASTACK", null);
            var peps = protSTC.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            var pepSTC = peps[0];
            Assert.AreEqual(pepSTC.BaseSequence, "GASTACK");
            Crosslinker crosslinker2 = new Crosslinker("ST", "C", "crosslinkerSTC", false, "", -18.01056, 0, 0, 0, 0, 0, 0);
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

            var xlSearchParameters = new XlSearchParameters
            {
                CrosslinkAtCleavageSite = true
            };

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
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse, commonParameters, null, 30000, false, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            var indexedFragments = indexResults.FragmentIndex.Where(p => p != null).SelectMany(v => v).ToList();
            Assert.AreEqual(82, indexedFragments.Count);
            Assert.AreEqual(3, indexResults.PeptideIndex.Count);

            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFile();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //Generate crosslinker, which is DSSO here.
            Crosslinker crosslinker = GlobalVariables.Crosslinkers.Where(p => p.CrosslinkerName == "DSSO").First();

            List<CrosslinkSpectralMatch>[] possiblePsms = new List<CrosslinkSpectralMatch>[listOfSortedms2Scans.Length];
            new CrosslinkSearchEngine(possiblePsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, null, 0, commonParameters, null, crosslinker, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.CrosslinkAtCleavageSite, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, new List<string> { }).Run();

            var newPsms = possiblePsms.Where(p => p != null).Select(p => p.First()).ToList();
            foreach (var item in newPsms)
            {
                item.ResolveAllAmbiguities();
                if (item.BetaPeptide != null)
                {
                    item.BetaPeptide.ResolveAllAmbiguities();
                }
                item.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                item.ResolveProteinPosAmbiguitiesForXl();
            }
            FdrAnalysisEngine fdrAnalysisEngine = new FdrAnalysisEngine(newPsms.ToList<PeptideSpectralMatch>(), 0, commonParameters, null, new List<string>(), "");

            Assert.AreEqual(4, newPsms.Count);
            Assert.That(newPsms[0].XlProteinPos == null); //single
            Assert.That(newPsms[1].XlProteinPos == 4 && newPsms[1].XlProteinPosLoop == 7); //loop
            Assert.That(newPsms[2].XlProteinPos == 4); //deadend
            Assert.That(newPsms[3].XlProteinPos == 2 && newPsms[3].BetaPeptide.XlProteinPos == 4); //cross

            //Test Output
            var task = new XLSearchTask();
            WriteFile.WritePepXML_xl(newPsms, proteinList, null, variableModifications, fixedModifications, null, TestContext.CurrentContext.TestDirectory, "pep.XML", commonParameters, xlSearchParameters);

            File.Delete(@"singlePsms.tsv");
            File.Delete(@"pep.XML.pep.xml");
            File.Delete(@"allPsms.tsv");

            // write percolator result
            WriteFile.WriteCrosslinkToTxtForPercolator(newPsms.Where(q => q.CrossType == PsmCrossType.Cross).ToList(), TestContext.CurrentContext.TestDirectory, "perc", new Crosslinker());
            var percOut = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, @"perc.txt"), Encoding.UTF8);
            string header = "SpecId\tLabel\tScannr\tScore\tdScore\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProtein";
            string dataRow = "T-2-1\t1\t2\t9.08035714285714\t9.08035714285714\t3\t1994.05\t79237.2823474838\t7\t9\t16\t-.EKVLTSSAR2--LSQKFPK4.-\tFake01(2)\tFake02(4)";
            Assert.AreEqual(header, percOut[0]);
            Assert.AreEqual(dataRow, percOut[1]);
        }

        [Test]
        public static void XlTest_BSA_DSS_file()
        {
            var task = Toml.ReadFile<XLSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/Rappsilber3-XLSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData"));
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/YeastPol2.fasta"), false);

            var origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/Rappsilber_CLMS_PolII_3-calib_slice.mzML.gz");
            string raw;

            FileInfo fileToDecompress = new FileInfo(origDataFile);

            using (FileStream originalFileStream = fileToDecompress.OpenRead())
            {
                string currentFileName = fileToDecompress.FullName;
                raw = currentFileName.Remove(currentFileName.Length - fileToDecompress.Extension.Length);

                using (FileStream decompressedFileStream = File.Create(raw))
                {
                    using (GZipStream decompressionStream = new GZipStream(originalFileStream, CompressionMode.Decompress))
                    {
                        decompressionStream.CopyTo(decompressedFileStream);
                        Console.WriteLine($"Decompressed: {fileToDecompress.Name}");
                    }
                }
            }

            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData")).Run();
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData"), true);
        }

        [Test]
        public static void XlTest_MoreComprehensive()
        {
            //Generate parameters
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, dissociationType: DissociationType.HCD,
                scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 5), precursorMassTolerance: new PpmTolerance(10));

            var xlSearchParameters = new XlSearchParameters
            {
                CrosslinkAtCleavageSite = true
            };

            //Create databases contain two protein.
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/YeastPol2.fasta"), true, DecoyType.Reverse, false, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, out var dbErrors, -1);

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            var variableModifications = new List<Modification>() { mod1 };
            var fixedModifications = new List<Modification>() { mod2 };
            var localizeableModifications = new List<Modification>();

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse, commonParameters, null, 30000, false, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            var origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/Rappsilber_CLMS_PolII_3-calib_slice.mzML.gz");



            string newFileName;

            FileInfo fileToDecompress = new FileInfo(origDataFile);

            using (FileStream originalFileStream = fileToDecompress.OpenRead())
            {
                string currentFileName = fileToDecompress.FullName;
                newFileName = currentFileName.Remove(currentFileName.Length - fileToDecompress.Extension.Length);

                using (FileStream decompressedFileStream = File.Create(newFileName))
                {
                    using (GZipStream decompressionStream = new GZipStream(originalFileStream, CompressionMode.Decompress))
                    {
                        decompressionStream.CopyTo(decompressedFileStream);
                    }
                }
            }

            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams());


            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add((Path.GetFileName(newFileName), CommonParameters));

            var myMsDataFile = myFileManager.LoadFile(newFileName, CommonParameters);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, newFileName, CommonParameters).OrderBy(b => b.PrecursorMass).ToArray();

            //Generate crosslinker, which is DSS here.
            Crosslinker crosslinker = GlobalVariables.Crosslinkers.Where(p => p.CrosslinkerName == "DSS").First();

            List<CrosslinkSpectralMatch>[] possiblePsms = new List<CrosslinkSpectralMatch>[listOfSortedms2Scans.Length];
            new CrosslinkSearchEngine(possiblePsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, null, 0, commonParameters, null,
                crosslinker, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.CrosslinkAtCleavageSite,
                xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, new List<string> { }).Run();

            var nonNullCsmsStillLists = possiblePsms.Where(p => p != null).ToList();

            #region Parsimony and assign crosslink

            XLSearchTask t = new XLSearchTask();

            List<List<CrosslinkSpectralMatch>> ListOfCsmsPerMS2ScanParsimony = new List<List<CrosslinkSpectralMatch>>();
            foreach (var csmsPerScan in nonNullCsmsStillLists)
            {
                foreach (var csm in csmsPerScan)
                {
                    csm.ResolveAllAmbiguities();
                    if (csm.BetaPeptide != null)
                    {
                        csm.BetaPeptide.ResolveAllAmbiguities();
                    }
                    csm.ResolveProteinPosAmbiguitiesForXl();
                }

                var orderedCsmsPerScan = XLSearchTask.RemoveDuplicateFromCsmsPerScan(csmsPerScan).OrderByDescending(p => p.XLTotalScore).ThenBy(p => p.FullSequence + ((p.BetaPeptide == null) ? "" : p.BetaPeptide.FullSequence)).ToList();

                ListOfCsmsPerMS2ScanParsimony.Add(orderedCsmsPerScan);
            }

            nonNullCsmsStillLists = ListOfCsmsPerMS2ScanParsimony;

            t.AssignCrossType(nonNullCsmsStillLists);

            #endregion Parsimony and assign crosslink

            t.SortCsmsAndCalculateDeltaScores(nonNullCsmsStillLists);
            foreach (List<CrosslinkSpectralMatch> xlinkCsmList in nonNullCsmsStillLists)
            {
                if (xlinkCsmList != null && xlinkCsmList.Any() && xlinkCsmList.Count > 1)
                {
                    double deltaScore = xlinkCsmList[0].XLTotalScore - xlinkCsmList[1].XLTotalScore;
                    Assert.That(deltaScore, Is.EqualTo(xlinkCsmList[0].DeltaScore).Within(0.05));
                }
            }

            int unnasignedCrossType = 0;
            int inter = 0;
            int intra = 0;
            int single = 0;
            int loop = 0;
            int deadend = 0;
            int deadendH2O = 0;
            int deadendNH2 = 0;
            int deadendTris = 0;

            List<CrosslinkSpectralMatch> firstCsmsFromListsOfCsms = nonNullCsmsStillLists.Select(c => c.First()).OrderBy(d => -d.XLTotalScore).ToList();

            foreach (CrosslinkSpectralMatch csm in firstCsmsFromListsOfCsms)
            {
                switch (csm.CrossType)
                {
                    case PsmCrossType.Inter:
                        inter++;
                        break;

                    case PsmCrossType.Intra:
                        intra++;
                        break;

                    case PsmCrossType.Single:
                        single++;
                        break;

                    case PsmCrossType.Loop:
                        loop++;
                        break;

                    case PsmCrossType.DeadEnd:
                        deadend++;
                        break;

                    case PsmCrossType.DeadEndH2O:
                        deadendH2O++;
                        break;

                    case PsmCrossType.DeadEndNH2:
                        deadendNH2++;
                        break;

                    case PsmCrossType.DeadEndTris:
                        deadendTris++;
                        break;

                    default:
                        unnasignedCrossType++;
                        break;
                }
            }

            Assert.AreEqual(521, inter);
            Assert.AreEqual(217, intra);
            Assert.AreEqual(287, single);
            Assert.AreEqual(13, loop);
            Assert.AreEqual(0, deadend);
            Assert.AreEqual(65, deadendH2O);
            Assert.AreEqual(0, deadendNH2);
            Assert.AreEqual(0, deadendTris);
            Assert.AreEqual(0, unnasignedCrossType);

            var fdrResultsXLink = new FdrAnalysisEngine(firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Inter || c.CrossType == PsmCrossType.Intra).ToList<PeptideSpectralMatch>(), 1, CommonParameters, fsp, new List<string>(), "crosslink").Run();


            fdrResultsXLink = new FdrAnalysisEngine(firstCsmsFromListsOfCsms.Where(c => c.CrossType != PsmCrossType.Inter && c.CrossType != PsmCrossType.Intra).ToList<PeptideSpectralMatch>(), 1, CommonParameters, fsp, new List<string>(), "standard").Run();

            unnasignedCrossType = 0;
            inter = 0;
            intra = 0;
            single = 0;
            loop = 0;
            deadend = 0;
            deadendH2O = 0;
            deadendNH2 = 0;
            deadendTris = 0;

            foreach (CrosslinkSpectralMatch csm in firstCsmsFromListsOfCsms.Where(c => c.FdrInfo.QValue <= 0.01).ToList())
            {
                switch (csm.CrossType)
                {
                    case PsmCrossType.Inter:
                        inter++;
                        break;

                    case PsmCrossType.Intra:
                        intra++;
                        break;

                    case PsmCrossType.Single:
                        single++;
                        break;

                    case PsmCrossType.Loop:
                        loop++;
                        break;

                    case PsmCrossType.DeadEnd:
                        deadend++;
                        break;

                    case PsmCrossType.DeadEndH2O:
                        deadendH2O++;
                        break;

                    case PsmCrossType.DeadEndNH2:
                        deadendNH2++;
                        break;

                    case PsmCrossType.DeadEndTris:
                        deadendTris++;
                        break;

                    default:
                        unnasignedCrossType++;
                        break;
                }
            }

            Assert.AreEqual(63, inter);
            Assert.AreEqual(83, intra);
            Assert.AreEqual(212, single);
            Assert.AreEqual(6, loop);
            Assert.AreEqual(0, deadend);
            Assert.AreEqual(42, deadendH2O);
            Assert.AreEqual(0, deadendNH2);
            Assert.AreEqual(0, deadendTris);
            Assert.AreEqual(0, unnasignedCrossType);

            var task = new PostXLSearchAnalysisTask();
            task.FileSpecificParameters = new List<(string FileName, CommonParameters Parameters)>();
            task.ComputeXlinkQandPValues(firstCsmsFromListsOfCsms, firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Intra).ToList(), firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Inter).ToList(), commonParameters, "");

            //check that alpha peptides have greater score than beta peptides
            foreach (CrosslinkSpectralMatch csm in firstCsmsFromListsOfCsms)
            {
                if (csm.CrossType == PsmCrossType.Intra || csm.CrossType == PsmCrossType.Inter)
                {
                    Assert.GreaterOrEqual(csm.Score, csm.BetaPeptide.Score);
                }
            }

            Dictionary<string, int> sequenceToPsmCount = new Dictionary<string, int>();
            List<string> sequences = new List<string>();
            foreach (CrosslinkSpectralMatch psm in firstCsmsFromListsOfCsms)
            {
                sequences.Add(psm.FullSequence);
            }
            var s = sequences.GroupBy(i => i);

            foreach (var grp in s)
            {
                sequenceToPsmCount.Add(grp.Key, grp.Count());
            }
            Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();
            Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();
            int chargeStateMode = 4;
            string[] trainingVariables = new[] { "TotalMatchingFragmentCount", "DeltaScore", "AlphaIntensity", "BetaIntensity", "LongestFragmentIonSeries_Alpha", "LongestFragmentIonSeries_Beta", "IsInter", "IsIntra" };
            //test PsmData for intra crosslink
            CrosslinkSpectralMatch intraCsm = firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Intra).First();

            Dictionary<string, float> medianFragmentMassError = new Dictionary<string, float>
            {
                { Path.GetFileName(intraCsm.FullFilePath), 0 }
            };

            var intraPsmData = PEP_Analysis.CreateOnePsmDataEntry("crosslink", fsp, intraCsm, sequenceToPsmCount, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, medianFragmentMassError, chargeStateMode, intraCsm.BestMatchingPeptides.First().Peptide, trainingVariables, intraCsm.BestMatchingPeptides.First().Notch, !intraCsm.BestMatchingPeptides.First().Peptide.Protein.IsDecoy);
            Assert.That(intraPsmData.AbsoluteAverageFragmentMassErrorFromMedian, Is.EqualTo(1.061).Within(0.001));
            Assert.That(intraPsmData.AlphaIntensity, Is.EqualTo(0.0098).Within(0.001));
            Assert.AreEqual(intraPsmData.Ambiguity, 0);
            Assert.That(intraPsmData.BetaIntensity, Is.EqualTo(0.01256).Within(0.0001));
            Assert.That(intraPsmData.DeltaScore, Is.EqualTo(0.467).Within(0.01));
            Assert.AreEqual(intraPsmData.HydrophobicityZScore, Double.NaN);
            Assert.AreEqual(intraPsmData.Intensity, 0);
            Assert.AreEqual(intraPsmData.IsDeadEnd, 0);
            Assert.AreEqual(intraPsmData.IsInter, 0);
            Assert.AreEqual(intraPsmData.IsIntra, 1);
            Assert.AreEqual(intraPsmData.IsLoop, 0);
            Assert.AreEqual(intraPsmData.IsVariantPeptide, 0);
            Assert.AreEqual(intraPsmData.Label, true);
            Assert.AreEqual(intraPsmData.LongestFragmentIonSeries, 0);
            Assert.That(intraPsmData.LongestFragmentIonSeries_Alpha, Is.EqualTo(0.857).Within(0.01));
            Assert.That(intraPsmData.LongestFragmentIonSeries_Beta, Is.EqualTo(0.75).Within(0.001));
            Assert.AreEqual(intraPsmData.MissedCleavagesCount, 0);
            Assert.AreEqual(intraPsmData.ModsCount, 0);
            Assert.AreEqual(intraPsmData.Notch, 0);
            Assert.AreEqual(intraPsmData.PrecursorChargeDiffToMode, -1);
            Assert.AreEqual(intraPsmData.PsmCount, 1);
            Assert.That(intraPsmData.TotalMatchingFragmentCount, Is.EqualTo(1.115).Within(0.001));

            CrosslinkSpectralMatch singleCsm = firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Single).OrderBy(c => -c.Score).First();

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();
            psms.AddRange(firstCsmsFromListsOfCsms);

            fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = PEP_Analysis.ComputeHydrophobicityValues(psms, fsp, false);
            fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = PEP_Analysis.ComputeHydrophobicityValues(psms, fsp, true);

            var singleCsmPsmData = PEP_Analysis.CreateOnePsmDataEntry("standard", fsp, singleCsm, sequenceToPsmCount, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, medianFragmentMassError, chargeStateMode, singleCsm.BestMatchingPeptides.First().Peptide, trainingVariables, singleCsm.BestMatchingPeptides.First().Notch, !singleCsm.BestMatchingPeptides.First().Peptide.Protein.IsDecoy);
            Assert.That(singleCsmPsmData.AbsoluteAverageFragmentMassErrorFromMedian, Is.EqualTo(0.4017).Within(0.001));
            Assert.AreEqual(singleCsmPsmData.AlphaIntensity, 0);
            Assert.AreEqual(singleCsmPsmData.Ambiguity, 0);
            Assert.AreEqual(singleCsmPsmData.BetaIntensity, 0);
            Assert.That(singleCsmPsmData.ComplementaryIonCount, Is.EqualTo(0.4).Within(0.01));
            Assert.That(singleCsmPsmData.DeltaScore, Is.EqualTo(1.05233).Within(0.001));
            Assert.That(singleCsmPsmData.HydrophobicityZScore, Is.EqualTo(0.8346).Within(0.001));
            Assert.That(singleCsmPsmData.Intensity, Is.EqualTo(0.01233).Within(0.0001));
            Assert.AreEqual(singleCsmPsmData.IsDeadEnd, 0);
            Assert.AreEqual(singleCsmPsmData.IsInter, 0);
            Assert.AreEqual(singleCsmPsmData.IsIntra, 0);
            Assert.AreEqual(singleCsmPsmData.IsLoop, 0);
            Assert.AreEqual(singleCsmPsmData.IsVariantPeptide, 0);
            Assert.AreEqual(singleCsmPsmData.Label, true);
            Assert.That(singleCsmPsmData.LongestFragmentIonSeries, Is.EqualTo(0.64).Within(0.01));
            Assert.AreEqual(singleCsmPsmData.LongestFragmentIonSeries_Alpha, 0);
            Assert.AreEqual(singleCsmPsmData.LongestFragmentIonSeries_Beta, 0);
            Assert.AreEqual(singleCsmPsmData.MissedCleavagesCount, 1);
            Assert.AreEqual(singleCsmPsmData.ModsCount, 1);
            Assert.AreEqual(singleCsmPsmData.Notch, 0);
            Assert.AreEqual(singleCsmPsmData.PrecursorChargeDiffToMode, -1);
            Assert.AreEqual(singleCsmPsmData.PsmCount, 4);
            Assert.That(singleCsmPsmData.TotalMatchingFragmentCount, Is.EqualTo(1.08).Within(0.01));

            CrosslinkSpectralMatch loopCsm = firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Loop).OrderBy(c => -c.Score).First();
            var loopCsmPsmData = PEP_Analysis.CreateOnePsmDataEntry("standard", fsp, loopCsm, sequenceToPsmCount, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, medianFragmentMassError, chargeStateMode, loopCsm.BestMatchingPeptides.First().Peptide, trainingVariables, loopCsm.BestMatchingPeptides.First().Notch, !loopCsm.BestMatchingPeptides.First().Peptide.Protein.IsDecoy);
            Assert.That(loopCsmPsmData.AbsoluteAverageFragmentMassErrorFromMedian, Is.EqualTo(0.511).Within(0.001));
            Assert.AreEqual(loopCsmPsmData.AlphaIntensity, 0);
            Assert.AreEqual(loopCsmPsmData.Ambiguity, 0);
            Assert.AreEqual(loopCsmPsmData.BetaIntensity, 0);
            Assert.That(loopCsmPsmData.ComplementaryIonCount, Is.EqualTo(0.3).Within(0.01));
            Assert.That(loopCsmPsmData.DeltaScore, Is.EqualTo(0.66497).Within(0.001));
            Assert.That(loopCsmPsmData.HydrophobicityZScore, Is.EqualTo(0.5264).Within(0.001));
            Assert.That(loopCsmPsmData.Intensity, Is.EqualTo(0.015).Within(0.0001));
            Assert.AreEqual(loopCsmPsmData.IsDeadEnd, 0);
            Assert.AreEqual(loopCsmPsmData.IsInter, 0);
            Assert.AreEqual(loopCsmPsmData.IsIntra, 0);
            Assert.AreEqual(loopCsmPsmData.IsLoop, 1);
            Assert.AreEqual(loopCsmPsmData.IsVariantPeptide, 0);
            Assert.AreEqual(loopCsmPsmData.Label, true);
            Assert.That(loopCsmPsmData.LongestFragmentIonSeries, Is.EqualTo(0.25).Within(0.01));
            Assert.AreEqual(loopCsmPsmData.LongestFragmentIonSeries_Alpha, 0);
            Assert.AreEqual(loopCsmPsmData.LongestFragmentIonSeries_Beta, 0);
            Assert.AreEqual(loopCsmPsmData.MissedCleavagesCount, 2);
            Assert.AreEqual(loopCsmPsmData.ModsCount, 0);
            Assert.AreEqual(loopCsmPsmData.Notch, 0);
            Assert.AreEqual(loopCsmPsmData.PrecursorChargeDiffToMode, -1);
            Assert.AreEqual(loopCsmPsmData.PsmCount, 8);
            Assert.That(loopCsmPsmData.TotalMatchingFragmentCount, Is.EqualTo(0.7).Within(0.01));

            unnasignedCrossType = 0;
            inter = 0;
            intra = 0;
            single = 0;
            loop = 0;
            deadend = 0;
            deadendH2O = 0;
            deadendNH2 = 0;
            deadendTris = 0;

            foreach (CrosslinkSpectralMatch csm in firstCsmsFromListsOfCsms.Where(c => c.FdrInfo.PEP_QValue <= 0.01).ToList())
            {
                switch (csm.CrossType)
                {
                    case PsmCrossType.Inter:
                        inter++;
                        break;

                    case PsmCrossType.Intra:
                        intra++;
                        break;

                    case PsmCrossType.Single:
                        single++;
                        break;

                    case PsmCrossType.Loop:
                        loop++;
                        break;

                    case PsmCrossType.DeadEnd:
                        deadend++;
                        break;

                    case PsmCrossType.DeadEndH2O:
                        deadendH2O++;
                        break;

                    case PsmCrossType.DeadEndNH2:
                        deadendNH2++;
                        break;

                    case PsmCrossType.DeadEndTris:
                        deadendTris++;
                        break;

                    default:
                        unnasignedCrossType++;
                        break;
                }
            }

            Assert.AreEqual(0, unnasignedCrossType);
            Assert.AreEqual(52, inter);
            Assert.AreEqual(77, intra);
            Assert.AreEqual(216, single);
            Assert.AreEqual(6, loop);
            Assert.AreEqual(0, deadend);
            Assert.AreEqual(44, deadendH2O);
            Assert.AreEqual(0, deadendNH2);
            Assert.AreEqual(0, deadendTris);
        }

        [Test]
        public static void XlTest_DiffCrosslinkSites()
        {
            //Generate parameters
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 4));

            var xlSearchParameters = new XlSearchParameters
            {
                Crosslinker = new Crosslinker(crosslinkerName: "CrossST-C", crosslinkerModSites: "ST", crosslinkerModSites2: "C", totalMass: -18.01056,
                cleavable: true, dissociationTypes: "HCD", deadendMassH2O: 0, deadendMassNH2: 0, deadendMassTris: 0, cleaveMassShort: 0, cleaveMassLong: 0, loopMass: 0)
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
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse, commonParameters, null, 30000, false, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFileDiffSite();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //Generate crosslinker, which is UserDefined here.
            var crosslinker = xlSearchParameters.Crosslinker;

            //TwoPassCrosslinkSearchEngine.Run().
            List<CrosslinkSpectralMatch>[] possiblePsms = new List<CrosslinkSpectralMatch>[listOfSortedms2Scans.Length];
            new CrosslinkSearchEngine(possiblePsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, null, 0, commonParameters, null,
                crosslinker, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.CrosslinkAtCleavageSite,
                xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, new List<string> { }).Run();

            var newPsms = possiblePsms.Where(p => p != null).ToList();
            Assert.AreEqual(1, newPsms.Count);
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

        [Test]
        public static void AmbiguiousDeadendPeptideTest() //The scan was previously identified as crosslink, but it is actually deadend
        {
            string myFileXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSSO_29061.mgf");
            string myDatabaseXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestDeadendPeptide");

            XLSearchTask xLSearchTask = new XLSearchTask() { };

            Directory.CreateDirectory(outputFolder);

            xLSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabaseXl, false) }, new List<string> { myFileXl }, "test");
            var lines = File.ReadAllLines(Path.Combine(outputFolder, @"Deadends.tsv"));
            Assert.That(lines.Length == 2);
            Assert.That(lines[1].Contains("Dead"));

            Directory.Delete(outputFolder, true);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }

        /// <summary>
        /// Makes sure helper methods that generate indices function properly
        /// </summary>
        [Test]
        public static void XLSearchWithGeneratedIndices()
        {
            XLSearchTask xlSearchTask = new XLSearchTask
            {
                CommonParameters = new CommonParameters()
            };
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
            var csms = new List<CrosslinkSpectralMatch>[1];

            // generate the scan with the deadend mod peptide's fragments
            var scans = new Ms2ScanWithSpecificMass[1];
            ModificationMotif.TryGetMotif("T", out var motif);
            var crosslinker = new Crosslinker("T", "T", "test", false, "", 100, 0, 0, 0, 0, 0, 50);
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
            var indexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { protein }, new List<Modification>(), new List<Modification>(), null, null, null,
                0, DecoyType.None, new CommonParameters(), null, 1000, false, new List<FileInfo>(), new List<string>()).Run();

            new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, null, 0, new CommonParameters(), null, crosslinker,
                50, true, false, false, true, new List<string>()).Run();

            CrosslinkSpectralMatch csm = csms.First().First();
            csm.ResolveAllAmbiguities();
            csm.ResolveProteinPosAmbiguitiesForXl();
            Assert.That(csm.XlProteinPos == 4);
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
            Crosslinker c = new Crosslinker("P", "R", "Test", false, "", 1000, 0, 0, 1000, 5, 5, 5);

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
            Crosslinker c = new Crosslinker("P", "R", "Test", true, "HCD", 1000, 15, 25, 1000, 5, 5, 5);

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
        public static void TestWriteNonSingleCross()
        {
            XLSearchTask xlst = new XLSearchTask();
            Protein protein = new Protein("PEPTIDE", "");
            var csms = new List<CrosslinkSpectralMatch>[1];

            // generate the scan with the deadend mod peptide's fragments
            var scans = new Ms2ScanWithSpecificMass[1];
            ModificationMotif.TryGetMotif("T", out var motif);
            var crosslinker = new Crosslinker("T", "T", "test", false, "", 100, 0, 0, 0, 0, 0, 50);
            Modification deadend = new Modification("TestId", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: crosslinker.DeadendMassTris, _modificationType: "Test");

            var deadendPeptide = protein.Digest(new DigestionParams(), new List<Modification> { deadend }, new List<Modification>()).First();

            double[] mz = deadendPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both).Select(p => p.NeutralMass.ToMz(1)).OrderBy(v => v).ToArray();
            double[] intensities = new[] { 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

            MzSpectrum spectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan sc = new MsDataScan(spectrum, 1, 2, true, Polarity.Positive, 1, spectrum.Range, "",
                MZAnalyzerType.Orbitrap, 12, 1.0, null, null);
            scans[0] = new Ms2ScanWithSpecificMass(sc, deadendPeptide.MonoisotopicMass.ToMz(2), 2, "", new CommonParameters());

            var indexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { protein }, new List<Modification>(), new List<Modification>(), null, null, null, 0, DecoyType.None,
                new CommonParameters(), null, 1000, false, new List<FileInfo>(), new List<string>()).Run();

            new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, null, 0, new CommonParameters(), null,
                crosslinker, 50, true, false, false, true, new List<string>()).Run();

            csms.First()[0].ResolveAllAmbiguities();
            csms.First()[0].SetFdrValues(0, 0, 0.1, 0, 0, 0, 0, 0);

            WriteFile.WritePepXML_xl(csms.SelectMany(p => p).ToList(), new List<Protein>(), "", new List<Modification> { deadend }, new List<Modification> { deadend }, new List<string>(), TestContext.CurrentContext.TestDirectory, "test", new CommonParameters(), new XlSearchParameters { Crosslinker = crosslinker });
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test.pep.XML"));
        }

        [Test]
        public static void TestMixedMs2Ms2()
        {
            string outputFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMixedMs2Ms2.tsv");
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.CID, childScanDissociationType: DissociationType.ETD,
                trimMsMsPeaks: false);

            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\ms2mixed_bsa_xlink.mzML");

            var fsp = new List<(string, CommonParameters)>();
            fsp.Add((spectraFile, commonParameters));
            
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);

            var scans = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).ToArray();

            Assert.That(scans.First().ChildScans.Count == 1);
            Assert.That(scans.Length == 1);

            Protein bsa = new Protein("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYL" +
                "QQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDS" +
                "PDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMRE" +
                "KVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDR" +
                "ADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFL" +
                "GSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEK" +
                "LGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKT" +
                "PVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKP" +
                "KATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA", "BSA");

            //Add the same seq to test protein ambiguious assignment.
            Protein bsa2 = new Protein("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYL" +
                "QQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDS" +
                "PDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMRE" +
                "KVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDR" +
                "ADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFL" +
                "GSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEK" +
                "LGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKT" +
                "PVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKP" +
                "KATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA", "BSA2");

            var indexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { bsa, bsa2 }, new List<Modification>(), new List<Modification>(),
                null, null, null, 0, DecoyType.None, commonParameters, fsp, 5000, false, new List<FileInfo>(), new List<string>()).Run();

            var secondCombinedParams = new CommonParameters(dissociationType: DissociationType.ETD, childScanDissociationType: DissociationType.ETD,
                trimMsMsPeaks: false);

            var fsp2 = new List<(string, CommonParameters)>();
            fsp2.Add((spectraFile, secondCombinedParams));


            var secondIndexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { bsa }, new List<Modification>(), new List<Modification>(),
                null, null, null, 0, DecoyType.None, secondCombinedParams, fsp2, 5000, false, new List<FileInfo>(), new List<string>()).Run();

            var csms = new List<CrosslinkSpectralMatch>[1];
            new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, secondIndexingResults.FragmentIndex, 0, commonParameters, fsp,
                GlobalVariables.Crosslinkers.First(p => p.CrosslinkerName == "DSSO"), 50, true, false, false, true, new List<string>()).Run();

            foreach (var c in csms.First())
            {
                c.ResolveAllAmbiguities();
                if (c.BetaPeptide != null)
                {
                    c.BetaPeptide.ResolveAllAmbiguities();
                }
            }

            //This function is important for crosslink protein ambiguious assignment.
            var csm = XLSearchTask.RemoveDuplicateFromCsmsPerScan(csms.First()).First();
            var isIntra = csm.IsIntraCsm();
            Assert.That(isIntra == true);
            csm.ResolveProteinPosAmbiguitiesForXl();
            Assert.That(csm.XlProteinPos == 455 && csm.BetaPeptide.XlProteinPos == 211);

            // test parent scan (CID)
            Assert.That(csm.MatchedFragmentIons.Count == 21);
            Assert.That(csm.BetaPeptide.MatchedFragmentIons.Count == 13);
            Assert.That(csm.ScanNumber == 2);

            // test child scan (ETD)
            Assert.That(csm.ChildMatchedFragmentIons.First().Key == 3);
            Assert.That(csm.ChildMatchedFragmentIons.First().Value.Count == 22);
            Assert.That(csm.BetaPeptide.ChildMatchedFragmentIons.First().Key == 3);
            Assert.That(csm.BetaPeptide.ChildMatchedFragmentIons.First().Value.Count == 25);

            // write results to TSV
            csm.SetFdrValues(1, 0, 0, 0, 0, 0, 0, 0);
            WriteFile.WritePsmCrossToTsv(new List<CrosslinkSpectralMatch> { csm }, outputFile, 2);

            // read results from TSV
            var psmFromTsv = PsmTsvReader.ReadTsv(outputFile, out var warnings).First();

            Assert.That(psmFromTsv.ChildScanMatchedIons.Count == 1
                && psmFromTsv.ChildScanMatchedIons.First().Key == 3
                && psmFromTsv.ChildScanMatchedIons.First().Value.Count == 22);

            Assert.That(psmFromTsv.BetaPeptideChildScanMatchedIons.Count == 1
                && psmFromTsv.BetaPeptideChildScanMatchedIons.First().Key == 3
                && psmFromTsv.BetaPeptideChildScanMatchedIons.First().Value.Count == 25);

            Assert.That(csm.ProteinAccession == null && csm.BetaPeptide.ProteinAccession == null);
            Assert.That(psmFromTsv.ProteinAccession == "BSA|BSA2");

            File.Delete(outputFile);
        }

        [Test]
        public static void TestMs2Ms3()
        {
            string outputFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMs2Ms3.tsv");
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.CID, childScanDissociationType: DissociationType.LowCID, precursorMassTolerance: new PpmTolerance(10));

            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\10226.mzML");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);

            var fsp = new List<(string, CommonParameters)>();
            fsp.Add((spectraFile, commonParameters));

            var scans = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).ToArray();

            Assert.That(scans.First().ChildScans.Count == 4);
            Assert.That(scans.Length == 2);

            Protein bsa = new Protein("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYL" +
                "QQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDS" +
                "PDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMRE" +
                "KVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDR" +
                "ADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFL" +
                "GSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEK" +
                "LGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKT" +
                "PVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKP" +
                "KATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA", "BSA");

            var fixedMods = GlobalVariables.AllModsKnown.Where(p => p.IdWithMotif == "Carbamidomethyl on C").ToList();

            var indexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { bsa }, new List<Modification>(), fixedMods,
               null, null, null, 0, DecoyType.None, commonParameters, fsp, 5000, false, new List<FileInfo>(), new List<string>()).Run();

            var csms = new List<CrosslinkSpectralMatch>[2];
            new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, null, 0, commonParameters, fsp,
                GlobalVariables.Crosslinkers.First(p => p.CrosslinkerName == "DSSO"), 50, true, false, false, true, new List<string>()).Run();

            var csm = csms.First()[0];
            csm.ResolveAllAmbiguities();
            if (csm.BetaPeptide != null)
            {
                csm.BetaPeptide.ResolveAllAmbiguities();
            }
            // test parent scan (CID)
            Assert.That(csm.MatchedFragmentIons.Count == 37);
            Assert.That(csm.ScanNumber == 2);

            // test child scan (low-resolution CID, alpha peptide signature ion)
            Assert.That(csm.ChildMatchedFragmentIons.First().Key == 6);
            Assert.That(csm.ChildMatchedFragmentIons.First().Value.Count == 43);

            // test child scan (low-resolution CID, beta peptide signature ion)
            Assert.That(csm.BetaPeptide.ChildMatchedFragmentIons.First().Key == 4);
            Assert.That(csm.BetaPeptide.ChildMatchedFragmentIons.First().Value.Count == 63);

            // write results to TSV
            csm.SetFdrValues(1, 0, 0, 0, 0, 0, 0, 0);
            WriteFile.WritePsmCrossToTsv(new List<CrosslinkSpectralMatch> { csm }, outputFile, 2);

            // read results from TSV
            var psmFromTsv = PsmTsvReader.ReadTsv(outputFile, out var warnings).First();

            Assert.That(psmFromTsv.ChildScanMatchedIons.Count == 2
                && psmFromTsv.ChildScanMatchedIons.First().Key == 6
                && psmFromTsv.ChildScanMatchedIons.First().Value.Count == 43);

            Assert.That(psmFromTsv.BetaPeptideChildScanMatchedIons.Count == 2
                && psmFromTsv.BetaPeptideChildScanMatchedIons.First().Key == 4
                && psmFromTsv.BetaPeptideChildScanMatchedIons.First().Value.Count == 63);

            File.Delete(outputFile);
        }

        internal class XLTestDataFile : MsDataFile
        {
            //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
            public XLTestDataFile() : base(2, new SourceFile(null, null, null, null, null))
            {
                var mz1 = new double[] { 1994.05.ToMz(3), 846.4963.ToMz(1), 1004.495.ToMz(1), 1022.511.ToMz(1), 1093.544.ToMz(1), 1500.00.ToMz(1) };
                var intensities1 = new double[] { 1, 1, 1, 1, 1, 1 };
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

                var mz5 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
                var intensities5 = new double[] { 100, 1, 1, 1 };
                var MassSpectrum5 = new MzSpectrum(mz5, intensities5, false);
                ScansHere.Add(new MsDataScan(MassSpectrum5, 5, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=5", 1022.511.ToMz(1),
                    1, 1, 1022.511.ToMz(1), 2, DissociationType.HCD, 1, 1022.511.ToMz(1)));

                Scans = ScansHere.ToArray();
            }

            public static string FilePath
            {
                get
                {
                    return "XLTestDataFile";
                }
            }

            public static string Name
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
}