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
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using Omics.Digestion;
using Omics.Modifications;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class XLTest
    {
        [Test]
        public static void TestDissociationTypeGenerateSameTypeOfIons()
        {
            Assert.That(CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(DissociationType.CID, DissociationType.CID));
            Assert.That(CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(DissociationType.CID, DissociationType.HCD));
            Assert.That(CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(DissociationType.HCD, DissociationType.CID));
            Assert.That(CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(DissociationType.ETD, DissociationType.ECD));
            Assert.That(CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(DissociationType.ECD, DissociationType.ETD));
            Assert.That(!CrosslinkSearchEngine.DissociationTypeGenerateSameTypeOfIons(DissociationType.CID, DissociationType.ETD));
        }

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
            Assert.That(pep.BaseSequence, Is.EqualTo("MNNNK"));
            Crosslinker crosslinker = GlobalVariables.Crosslinkers.Where(p => p.CrosslinkerName == "DSS").First();
            Assert.That(crosslinker.CrosslinkerModSites, Is.EqualTo("K"));
            Assert.That(Residue.GetResidue(crosslinker.CrosslinkerModSites).MonoisotopicMass, Is.EqualTo(128.09496301518999).Within(1e-9));
            List<Product> n = new List<Product>();
            pep.Fragment(DissociationType.HCD, FragmentationTerminus.N, n);
            List<Product> c = new List<Product>();
            pep.Fragment(DissociationType.HCD, FragmentationTerminus.C, c);
            Assert.That(n.Count(), Is.EqualTo(4));
            Assert.That(c.Count(), Is.EqualTo(4));
            Assert.That(c.First().NeutralMass, Is.EqualTo(146.10552769899999).Within(1e-6));
            var x = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep, digestionParams.InitiatorMethionineBehavior, false);
            Assert.That(x == null);

            var pep2 = ye[2];
            Assert.That(pep2.BaseSequence, Is.EqualTo("MNNNKQQQQ"));
            List<Product> n2 = new List<Product>();
            pep2.Fragment(DissociationType.HCD, FragmentationTerminus.N, n2);
            List<Product> c2 = new List<Product>();
            pep2.Fragment(DissociationType.HCD, FragmentationTerminus.C, c2);
            Assert.That(n2.Count(), Is.EqualTo(8));
            Assert.That(c2.Count(), Is.EqualTo(8));
            var x2 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites.ToCharArray(), pep2, digestionParams.InitiatorMethionineBehavior, false);
            Assert.That(x2[0], Is.EqualTo(5));

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
            Assert.That(pepSTC.BaseSequence, Is.EqualTo("GASTACK"));
            Crosslinker crosslinker2 = new Crosslinker("ST", "C", "crosslinkerSTC", false, "", -18.01056, 0, 0, 0, 0, 0, 0);
            string crosslinkerModSitesAll = new string((crosslinker2.CrosslinkerModSites + crosslinker2.CrosslinkerModSites2).ToCharArray().Distinct().ToArray());
            Assert.That(crosslinkerModSitesAll, Is.EqualTo("STC"));
        }

        [Test]
        public static void XlTestGenerateIntensityRanks()
        {
            double[] intensity = new double[] { 1.1, 1.1, 0.5, 3.2, 0.5, 6.0 };
            int[] rank = CrosslinkSpectralMatch.GenerateIntensityRanks(intensity);
            int[] Rank = new int[] { 4, 3, 6, 2, 5, 1 };
            Assert.That(rank, Is.EqualTo(Rank));
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
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            var variableModifications = new List<Modification>() { mod1 };
            var fixedModifications = new List<Modification>() { mod2 };
            var localizeableModifications = new List<Modification>();

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse,
                commonParameters, null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            var indexedFragments = indexResults.FragmentIndex.Where(p => p != null).SelectMany(v => v).ToList();
            Assert.That(indexedFragments.Count, Is.EqualTo(82));
            Assert.That(indexResults.PeptideIndex.Count, Is.EqualTo(3));

            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFile();
            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = MetaMorpheusTask.GetMs2ScansWrapByScanNum(myMsDataFile, null, new CommonParameters(), out List<List<(double, int, double)>> precursorss).ToArray();

            //Generate crosslinker, which is DSSO here.
            Crosslinker crosslinker = GlobalVariables.Crosslinkers.Where(p => p.CrosslinkerName == "DSSO").First();

            List<CrosslinkSpectralMatch>[] possiblePsms = new List<CrosslinkSpectralMatch>[listOfSortedms2Scans.Length];
            List<(int, int, int)>[] candidates = new List<(int, int, int)>[listOfSortedms2Scans.Length];

            var XLEngine = new CrosslinkSearchEngine(possiblePsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, null, 0, 
                commonParameters, null, crosslinker, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.CrosslinkAtCleavageSite, 
                xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, new List<string> { },
                candidates, 0, indexResults.PeptideIndex, precursorss);
            XLEngine.FirstRoundSearch();
            XLEngine.Run();

            var newPsms = possiblePsms.Where(p => p != null).Select(p => p.First()).ToList();
            foreach (var item in newPsms)
            {
                item.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                CrosslinkSpectralMatch.ResolveProteinPosAmbiguitiesForXl(item);
            }

            List<(string fileName, CommonParameters fileSpecificParameters)> fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("filename", commonParameters) };

            FdrAnalysisEngine fdrAnalysisEngine = new FdrAnalysisEngine(newPsms.ToList<SpectralMatch>(), 0, commonParameters, fsp, new List<string>(), "");

            Assert.That(newPsms.Count, Is.EqualTo(4));
            Assert.That(newPsms[0].XlProteinPos == 2 && newPsms[0].BetaPeptide.XlProteinPos == 4); //cross
            Assert.That(newPsms[1].XlProteinPos == null); //single
            Assert.That(newPsms[2].XlProteinPos == 4 && newPsms[2].XlProteinPosLoop == 7); //loop
            Assert.That(newPsms[3].XlProteinPos == 4); //deadend           

            //Test Output
            var task = new XLSearchTask();
            WriteXlFile.WritePepXML_xl(newPsms, proteinList, null, variableModifications, fixedModifications, null, TestContext.CurrentContext.TestDirectory, "pep.XML", commonParameters, xlSearchParameters);

            File.Delete(@"singlePsms.tsv");
            File.Delete(@"pep.XML.pep.xml");
            File.Delete(@"allPsms.tsv");

            // write percolator result
            WriteXlFile.WriteCrosslinkToTxtForPercolator(newPsms.Where(q => q.CrossType == PsmCrossType.Inter || q.CrossType == PsmCrossType.Intra).ToList(), TestContext.CurrentContext.TestDirectory, "perc", new Crosslinker());
            var percOut = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, @"perc.txt"), Encoding.UTF8);
            string header = "SpecId\tLabel\tScannr\tScore\tdScore\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProtein";
            string dataRow = "T-2-1\t1\t2\t9.080357142857142\t9.080357142857142\t3\t1994.05\t79237.2823474838\t7\t9\t16\t-.EKVLTSSAR2--LSQKFPK4.-\tFake01(2)\tFake02(4)";
            Assert.That(percOut[0], Is.EqualTo(header));
            Assert.That(percOut[1], Is.EqualTo(dataRow));
            File.Delete(@"perc.txt");
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
                    }
                }
            }

            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData")).Run();
            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTXlTestData"), true);
        }

        [Test]
        public static void TestCsmSort()
        {
            Protein protForward = new Protein(sequence: "VPEPTIDELPEPTIDEAPEPTIDE", accession: "", isDecoy: false);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(
                new MsDataScan(
                    new MzSpectrum(new double[] { }, new double[] { }, false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);
            Dictionary<int, Modification> mod = new Dictionary<int, Modification>();

            PeptideWithSetModifications pwsmV = new PeptideWithSetModifications(protForward, new DigestionParams(), 1, 8, CleavageSpecificity.Full, "VPEPTIDE", 0, mod, 0, null);
            //PeptideWithSetModifications pwsmL = new PeptideWithSetModifications(protForward, new DigestionParams(), 9, 16, CleavageSpecificity.Full, "LPEPTIDE", 0, mod, 0, null);
            PeptideWithSetModifications pwsmA = new PeptideWithSetModifications(protForward, new DigestionParams(), 17, 24, CleavageSpecificity.Full, "APEPTIDE", 0, mod, 0, null);

            PeptideWithSetModifications pwsmVbetaA = new PeptideWithSetModifications(protForward, new DigestionParams(), 17, 20, CleavageSpecificity.Full, "APEP", 0, mod, 0, null);
            PeptideWithSetModifications pwsmLbetaV = new PeptideWithSetModifications(protForward, new DigestionParams(), 1, 4, CleavageSpecificity.Full, "VPEP", 0, mod, 0, null);
            PeptideWithSetModifications pwsmAbetaL = new PeptideWithSetModifications(protForward, new DigestionParams(), 9, 12, CleavageSpecificity.Full, "LPEP", 0, mod, 0, null);

            CrosslinkSpectralMatch csmOne = new CrosslinkSpectralMatch(pwsmV, 0, 3, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            CrosslinkSpectralMatch csmThree = new CrosslinkSpectralMatch(pwsmV, 0, 3, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            CrosslinkSpectralMatch csmTwo = new CrosslinkSpectralMatch(pwsmA, 0, 2, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());

            CrosslinkSpectralMatch csmOneBetaA = new CrosslinkSpectralMatch(pwsmVbetaA, 0, 30, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            CrosslinkSpectralMatch csmThreeBetaV = new CrosslinkSpectralMatch(pwsmLbetaV, 0, 30, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            CrosslinkSpectralMatch csmTwoBetaL = new CrosslinkSpectralMatch(pwsmAbetaL, 0, 20, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());

            csmOne.BetaPeptide = csmOneBetaA;
            csmTwo.BetaPeptide = csmTwoBetaL;
            csmThree.BetaPeptide = csmThreeBetaV;

            csmOne.XLTotalScore = 33;
            csmTwo.XLTotalScore = 22;
            csmThree.XLTotalScore = 33;

            csmOne.ResolveAllAmbiguities();
            csmTwo.ResolveAllAmbiguities();
            csmThree.ResolveAllAmbiguities();

            csmOne.BetaPeptide.ResolveAllAmbiguities();
            csmTwo.BetaPeptide.ResolveAllAmbiguities();
            csmThree.BetaPeptide.ResolveAllAmbiguities();

            List<CrosslinkSpectralMatch> forward = new List<CrosslinkSpectralMatch> { csmOne, csmThree, csmTwo };

            Protein protReverse = new Protein(sequence: "AEDITPEPVEDITPEPLEDITPEP", accession: "", isDecoy: false);

            PeptideWithSetModifications mswpA = new PeptideWithSetModifications(protReverse, new DigestionParams(), 1, 8, CleavageSpecificity.Full, "AEDITPEP", 0, mod, 0, null);
            PeptideWithSetModifications mswpV = new PeptideWithSetModifications(protReverse, new DigestionParams(), 9, 16, CleavageSpecificity.Full, "VEDITPEP", 0, mod, 0, null);
            PeptideWithSetModifications mswpL = new PeptideWithSetModifications(protReverse, new DigestionParams(), 17, 24, CleavageSpecificity.Full, "LEDITPEP", 0, mod, 0, null);

            PeptideWithSetModifications mswpAbetaL = new PeptideWithSetModifications(protReverse, new DigestionParams(), 17, 20, CleavageSpecificity.Full, "LEDI", 0, mod, 0, null);
            PeptideWithSetModifications mswpVbetaA = new PeptideWithSetModifications(protReverse, new DigestionParams(), 1, 4, CleavageSpecificity.Full, "AEDI", 0, mod, 0, null);
            PeptideWithSetModifications mswpLbetaV = new PeptideWithSetModifications(protReverse, new DigestionParams(), 9, 12, CleavageSpecificity.Full, "VEDI", 0, mod, 0, null);

            CrosslinkSpectralMatch twoCsm = new CrosslinkSpectralMatch(mswpA, 0, 2, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            CrosslinkSpectralMatch oneCsm = new CrosslinkSpectralMatch(mswpV, 0, 1, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            CrosslinkSpectralMatch threeCsm = new CrosslinkSpectralMatch(mswpL, 0, 1, 3, scan, new CommonParameters(), new List<MatchedFragmentIon>());

            CrosslinkSpectralMatch twoCsmBetaL = new CrosslinkSpectralMatch(mswpAbetaL, 0, 20, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            CrosslinkSpectralMatch oneCsmBetaA = new CrosslinkSpectralMatch(mswpVbetaA, 0, 10, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            CrosslinkSpectralMatch threeCsmBetaV = new CrosslinkSpectralMatch(mswpLbetaV, 0, 40, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());

            oneCsm.BetaPeptide = oneCsmBetaA;
            twoCsm.BetaPeptide = twoCsmBetaL;
            threeCsm.BetaPeptide = threeCsmBetaV;

            twoCsm.XLTotalScore = 22;
            oneCsm.XLTotalScore = 11;
            threeCsm.XLTotalScore = 43;

            twoCsm.ResolveAllAmbiguities();
            oneCsm.ResolveAllAmbiguities();
            threeCsm.ResolveAllAmbiguities();

            twoCsm.BetaPeptide.ResolveAllAmbiguities();
            oneCsm.BetaPeptide.ResolveAllAmbiguities();
            threeCsm.BetaPeptide.ResolveAllAmbiguities();

            Assert.That(oneCsm.UniqueSequence, Is.EqualTo(oneCsm.FullSequence));
            twoCsm.LinkPositions = new List<int> { 1 };
            twoCsm.BetaPeptide.LinkPositions = new List<int> { 2 };
            twoCsm.CrossType = PsmCrossType.Cross;
            Assert.That(twoCsm.UniqueSequence, Is.EqualTo(twoCsm.FullSequence + "(1)" + twoCsm.BetaPeptide.FullSequence + "(2)"));
            threeCsm.CrossType = PsmCrossType.Loop;
            threeCsm.LinkPositions = new List<int> { 1, 3 };
            Assert.That(threeCsm.UniqueSequence, Is.EqualTo(threeCsm.FullSequence + "(1-3)")); //Because the beta peptide link positions wasn't set, the csm unique sequence is a loop link sequence, not a crosslink sequence

            List<CrosslinkSpectralMatch> reverse = new List<CrosslinkSpectralMatch> { twoCsm, oneCsm, threeCsm };

            CommonParameters commonParameters = new CommonParameters(scoreCutoff: 2);
            List<CrosslinkSpectralMatch> sortedForward = XLSearchTask.SortOneListCsmsSetSecondBestScore(forward, commonParameters);

            Assert.That(sortedForward.Select(s => s.XLTotalScore).ToList(), Is.EquivalentTo(new List<double> { 33d, 33d, 22d }));
            Assert.That(sortedForward.Select(s => s.DeltaScore).ToList(), Is.EquivalentTo(new List<double> { 0d, 0d, -11d }));
            Assert.That(sortedForward.Select(aSeq => aSeq.BaseSequence).ToList(), Is.EquivalentTo(new List<string> { "VPEPTIDE", "VPEPTIDE", "APEPTIDE" }));
            Assert.That(sortedForward.Select(bSeq => bSeq.BetaPeptide.BaseSequence).ToList(), Is.EquivalentTo(new List<string> { "APEP", "VPEP", "LPEP" }));

            List<List<CrosslinkSpectralMatch>> csmListList = new List<List<CrosslinkSpectralMatch>> { forward, reverse };
            csmListList = XLSearchTask.SortListsOfCsms(csmListList, commonParameters);
            Assert.That(csmListList.Select(c => c.First().XLTotalScore).ToList(), Is.EquivalentTo(new List<double> { 43d, 33d }));
            Assert.That(csmListList.Select(c => c.First().BaseSequence).ToList(), Is.EquivalentTo(new List<string> { "LEDITPEP", "VPEPTIDE" }));
            Assert.That(csmListList.Select(c => c.First().BetaPeptide.BaseSequence).ToList(), Is.EquivalentTo(new List<string> { "VEDI", "APEP" }));


            Protein protForwardAsDecoy = new Protein(sequence: "VPEPTIDELPEPTIDEAPEPTIDE", accession: "DECOY", isDecoy: true);
            PeptideWithSetModifications pwsmV_D = new PeptideWithSetModifications(protForwardAsDecoy, new DigestionParams(), 1, 8, CleavageSpecificity.Full, "VPEPTIDE", 0, mod, 0, null);

            csmOne.AddOrReplace(pwsmV_D, 3, 0, true, new List<MatchedFragmentIon>(), 0);
            csmOne.ResolveAllAmbiguities();
            Assert.That(csmOne.BestMatchingBioPolymersWithSetMods.Count() == 1);
            Assert.That(!csmOne.BestMatchingBioPolymersWithSetMods.First().WithSetMods.Parent.IsDecoy);
        }

        [Test]
        [NonParallelizable]
        public static void XlTest_MoreComprehensive()
        {
            //override to be only used for unit tests in non-parallelizable format
            //must set to false at the end of this method
            var type = typeof(FdrAnalysisEngine);
            var property = type.GetProperty("QvalueThresholdOverride");
            property.SetValue(null, true);

            //Generate parameters
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, dissociationType: DissociationType.HCD,
                scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 5), precursorMassTolerance: new PpmTolerance(10), maxThreadsToUsePerFile: 1);

            var xlSearchParameters = new XlSearchParameters
            {
                CrosslinkAtCleavageSite = true
            };

            //Create databases contain two protein.
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/YeastPol2.fasta"), true, DecoyType.Reverse, false, out var dbErrors,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, -1);

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            var variableModifications = new List<Modification>() { mod1 };
            var fixedModifications = new List<Modification>() { mod2 };
            var localizeableModifications = new List<Modification>();

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse,
                commonParameters, null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

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
            CommonParameters commonParameters2 = new CommonParameters(digestionParams: new DigestionParams(), maxThreadsToUsePerFile: 1);

            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add((Path.GetFileName(newFileName), commonParameters2));

            var myMsDataFile = myFileManager.LoadFile(newFileName, commonParameters2);

            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = MetaMorpheusTask.GetMs2ScansWrapByScanNum(myMsDataFile, newFileName, commonParameters2, out List<List<(double, int, double)>> precursorss).ToArray();

            //Generate crosslinker, which is DSS here.
            Crosslinker crosslinker = GlobalVariables.Crosslinkers.Where(p => p.CrosslinkerName == "DSS").First();

            List<CrosslinkSpectralMatch>[] possiblePsms = new List<CrosslinkSpectralMatch>[listOfSortedms2Scans.Length];

            List<(int, int, int)>[] candidates = new List<(int, int, int)>[listOfSortedms2Scans.Length];

            var XLEngine = new CrosslinkSearchEngine(possiblePsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, null, 0,
                commonParameters, null, crosslinker, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.CrosslinkAtCleavageSite,
                xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, new List<string> { },
                candidates, 0, indexResults.PeptideIndex, precursorss);
            XLEngine.FirstRoundSearch();
            XLEngine.Run();

            #region Parsimony and assign crosslink

            List<List<CrosslinkSpectralMatch>> ListOfCsmsPerMS2ScanParsimony = new List<List<CrosslinkSpectralMatch>>();
            foreach (var csmsPerScan in possiblePsms.Where(p => p != null && p.Count > 0))
            {
                foreach (var csm in csmsPerScan)
                {
                    CrosslinkSpectralMatch.ResolveProteinPosAmbiguitiesForXl(csm);
                }

                ListOfCsmsPerMS2ScanParsimony.Add(csmsPerScan.ToList());
            }

            #endregion Parsimony and assign crosslink

            var nonNullCsmsStillLists = XLSearchTask.SortListsOfCsms(ListOfCsmsPerMS2ScanParsimony, commonParameters);
            foreach (List<CrosslinkSpectralMatch> xlinkCsmList in nonNullCsmsStillLists.Where(x => x.First().XLTotalScore > 0))
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

            Assert.That(inter, Is.EqualTo(435));
            Assert.That(intra, Is.EqualTo(215));
            Assert.That(single, Is.EqualTo(318));
            Assert.That(loop, Is.EqualTo(18));
            Assert.That(deadend, Is.EqualTo(0));
            Assert.That(deadendH2O, Is.EqualTo(82));
            Assert.That(deadendNH2, Is.EqualTo(0));
            Assert.That(deadendTris, Is.EqualTo(0));
            Assert.That(unnasignedCrossType, Is.EqualTo(0));

            // We have pretty high peptide-level q values for crosslinks, so we need to up the cut-off is we want PEP to run
            commonParameters2.QValueCutoffForPepCalculation = 0.05;
            var fdrResultsXLink = new FdrAnalysisEngine(firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Inter || c.CrossType == PsmCrossType.Intra).ToList<SpectralMatch>(), 1, commonParameters2, fsp, new List<string>(), "crosslink").Run();

            unnasignedCrossType = 0;
            inter = 0;
            intra = 0;
            single = 0;
            loop = 0;
            deadend = 0;
            deadendH2O = 0;
            deadendNH2 = 0;
            deadendTris = 0;

            foreach (CrosslinkSpectralMatch csm in firstCsmsFromListsOfCsms.Where(c => (c.CrossType == PsmCrossType.Inter || c.CrossType == PsmCrossType.Intra) && c.FdrInfo.PEP_QValue <= 0.02).ToList())
            {
                switch (csm.CrossType)
                {
                    case PsmCrossType.Inter:
                        inter++;
                        break;

                    case PsmCrossType.Intra:
                        intra++;
                        break;

                    default:
                        unnasignedCrossType++;
                        break;
                }
            }

            Assert.That(inter, Is.EqualTo(53));
            Assert.That(intra, Is.EqualTo(81));
            Assert.That(unnasignedCrossType, Is.EqualTo(0));

            // We have pretty high peptide-level q values for crosslinks, so we need to up the cut-off is we want PEP to run
            fdrResultsXLink = new FdrAnalysisEngine(firstCsmsFromListsOfCsms.Where(c => c.CrossType != PsmCrossType.Inter && c.CrossType != PsmCrossType.Intra).ToList<SpectralMatch>(), 1, commonParameters2, fsp, new List<string>(), "standard").Run();

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

            Assert.That(inter, Is.EqualTo(55));
            Assert.That(intra, Is.EqualTo(83));
            Assert.That(single, Is.EqualTo(229));
            Assert.That(loop, Is.EqualTo(8));
            Assert.That(deadend, Is.EqualTo(0));
            Assert.That(deadendH2O, Is.EqualTo(62));
            Assert.That(deadendNH2, Is.EqualTo(0));
            Assert.That(deadendTris, Is.EqualTo(0));
            Assert.That(unnasignedCrossType, Is.EqualTo(0));

            var task = new PostXLSearchAnalysisTask();
            task.FileSpecificParameters = new List<(string FileName, CommonParameters commonParameters)> { ("filename", new CommonParameters(maxThreadsToUsePerFile: 1)) };
            task.ComputeXlinkQandPValues(firstCsmsFromListsOfCsms, firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Intra).ToList(), firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Inter).ToList(), commonParameters, "");

            //check that alpha peptides have greater score than beta peptides
            foreach (CrosslinkSpectralMatch csm in firstCsmsFromListsOfCsms)
            {
                if (csm.CrossType == PsmCrossType.Intra || csm.CrossType == PsmCrossType.Inter)
                {
                    Assert.That(csm.Score >= csm.BetaPeptide.Score);
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

            //test PsmData for intra crosslink
            CrosslinkSpectralMatch intraCsm = firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Intra).First();

            Dictionary<string, float> medianFragmentMassError = new Dictionary<string, float>
            {
                { Path.GetFileName(intraCsm.FullFilePath), 0 }
            };

            // Set values within PEP_Analysis through reflection
            
            PepAnalysisEngine pepEngine = new PepAnalysisEngine(new List<SpectralMatch>(firstCsmsFromListsOfCsms), "standard", fsp, Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\"));
            var pepEngineProperties = pepEngine.GetType().GetProperties();
            foreach (var p in pepEngineProperties)
            {
                switch (p.Name)
                {
                    case "ChargeStateMode":
                        p.SetValue(pepEngine, chargeStateMode);
                        break;
                    case "FileSpecificMedianFragmentMassErrors":
                        p.SetValue(pepEngine, medianFragmentMassError);
                        break;
                    default:
                        break;
                }
            }

            var intraPsmData = pepEngine.CreateOnePsmDataEntry("crosslink", intraCsm, intraCsm.BestMatchingBioPolymersWithSetMods.First(), !intraCsm.BestMatchingBioPolymersWithSetMods.First().WithSetMods.Parent.IsDecoy);
            Assert.That(intraPsmData.AbsoluteAverageFragmentMassErrorFromMedian, Is.EqualTo(1.0).Within(0.1));
            Assert.That(intraPsmData.AlphaIntensity, Is.EqualTo(1).Within(0.1));
            Assert.That(intraPsmData.Ambiguity, Is.EqualTo(0));
            Assert.That(intraPsmData.BetaIntensity, Is.EqualTo(1).Within(0.1));
            Assert.That(intraPsmData.DeltaScore, Is.EqualTo(0).Within(0.1));
            Assert.That(intraPsmData.HydrophobicityZScore, Is.EqualTo(Double.NaN));
            Assert.That(intraPsmData.Intensity, Is.EqualTo(0));
            Assert.That(intraPsmData.IsDeadEnd, Is.EqualTo(0));
            Assert.That(intraPsmData.IsInter, Is.EqualTo(0));
            Assert.That(intraPsmData.IsIntra, Is.EqualTo(1));
            Assert.That(intraPsmData.IsLoop, Is.EqualTo(0));
            Assert.That(intraPsmData.IsVariantPeptide, Is.EqualTo(0));
            Assert.That(intraPsmData.Label, Is.EqualTo(true));
            Assert.That(intraPsmData.LongestFragmentIonSeries, Is.EqualTo(0));
            Assert.That(intraPsmData.LongestFragmentIonSeries_Alpha, Is.EqualTo(9).Within(0.1));
            Assert.That(intraPsmData.LongestFragmentIonSeries_Beta, Is.EqualTo(8).Within(0.1));
            Assert.That(intraPsmData.MissedCleavagesCount, Is.EqualTo(0));
            Assert.That(intraPsmData.ModsCount, Is.EqualTo(0));
            Assert.That(intraPsmData.Notch, Is.EqualTo(0));
            Assert.That(intraPsmData.PrecursorChargeDiffToMode, Is.EqualTo(-1));
            Assert.That(intraPsmData.TotalMatchingFragmentCount, Is.EqualTo(11).Within(0.1));

            CrosslinkSpectralMatch singleCsm = firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Single).OrderBy(c => -c.Score).First();

            var singleCsmPsmData = pepEngine.CreateOnePsmDataEntry("standard",
                singleCsm,
                singleCsm.BestMatchingBioPolymersWithSetMods.FirstOrDefault(),
                !singleCsm.BestMatchingBioPolymersWithSetMods.FirstOrDefault().IsDecoy);
            Assert.That(singleCsmPsmData.AbsoluteAverageFragmentMassErrorFromMedian, Is.EqualTo(8).Within(0.1));
            Assert.That(singleCsmPsmData.AlphaIntensity, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.Ambiguity, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.BetaIntensity, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.ComplementaryIonCount, Is.EqualTo(2).Within(0.1));
            Assert.That(singleCsmPsmData.DeltaScore, Is.EqualTo(8).Within(0.1));
            Assert.That(singleCsmPsmData.HydrophobicityZScore, Is.EqualTo(5).Within(0.1));
            Assert.That(singleCsmPsmData.Intensity, Is.EqualTo(0).Within(0.1));
            Assert.That(singleCsmPsmData.IsDeadEnd, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.IsInter, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.IsIntra, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.IsLoop, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.IsVariantPeptide, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.Label, Is.EqualTo(true));
            Assert.That(singleCsmPsmData.LongestFragmentIonSeries, Is.EqualTo(4).Within(0.1));
            Assert.That(singleCsmPsmData.LongestFragmentIonSeries_Alpha, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.LongestFragmentIonSeries_Beta, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.MissedCleavagesCount, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.ModsCount, Is.EqualTo(1));
            Assert.That(singleCsmPsmData.Notch, Is.EqualTo(0));
            Assert.That(singleCsmPsmData.PrecursorChargeDiffToMode, Is.EqualTo(-1));
            Assert.That(singleCsmPsmData.TotalMatchingFragmentCount, Is.EqualTo(8).Within(0.1));


            CrosslinkSpectralMatch loopCsm = firstCsmsFromListsOfCsms.Where(c => c.CrossType == PsmCrossType.Loop).OrderBy(c => -c.Score).First();
            var loopCsmPsmData = pepEngine.CreateOnePsmDataEntry("standard", loopCsm, loopCsm.BestMatchingBioPolymersWithSetMods.First(), !loopCsm.BestMatchingBioPolymersWithSetMods.First().WithSetMods.Parent.IsDecoy); Assert.That(loopCsmPsmData.AbsoluteAverageFragmentMassErrorFromMedian, Is.EqualTo(6).Within(0.1));
            Assert.That(loopCsmPsmData.AlphaIntensity, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.Ambiguity, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.BetaIntensity, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.ComplementaryIonCount, Is.EqualTo(3).Within(0.1));
            Assert.That(loopCsmPsmData.DeltaScore, Is.EqualTo(8).Within(0.1));
            Assert.That(loopCsmPsmData.HydrophobicityZScore, Is.EqualTo(9).Within(0.1));
            Assert.That(loopCsmPsmData.Intensity, Is.EqualTo(1).Within(0.1));
            Assert.That(loopCsmPsmData.IsDeadEnd, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.IsInter, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.IsIntra, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.IsLoop, Is.EqualTo(1));
            Assert.That(loopCsmPsmData.IsVariantPeptide, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.Label, Is.EqualTo(true));
            Assert.That(loopCsmPsmData.LongestFragmentIonSeries, Is.EqualTo(3).Within(0.1));
            Assert.That(loopCsmPsmData.LongestFragmentIonSeries_Alpha, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.LongestFragmentIonSeries_Beta, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.MissedCleavagesCount, Is.EqualTo(2));
            Assert.That(loopCsmPsmData.ModsCount, Is.EqualTo(2));
            Assert.That(loopCsmPsmData.Notch, Is.EqualTo(0));
            Assert.That(loopCsmPsmData.PrecursorChargeDiffToMode, Is.EqualTo(-1));
            Assert.That(loopCsmPsmData.TotalMatchingFragmentCount, Is.EqualTo(8).Within(0.1));

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

            Assert.That(unnasignedCrossType, Is.EqualTo(0));
            Assert.That(inter, Is.EqualTo(40));
            Assert.That(intra, Is.EqualTo(49));
            Assert.That(single, Is.EqualTo(231));
            Assert.That(loop, Is.EqualTo(0));
            Assert.That(deadend, Is.EqualTo(0));
            Assert.That(deadendH2O, Is.EqualTo(0));
            Assert.That(deadendNH2, Is.EqualTo(0));
            Assert.That(deadendTris, Is.EqualTo(0));

            property.SetValue(null, false);
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
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse,
                commonParameters, null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFileDiffSite();
            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = MetaMorpheusTask.GetMs2ScansWrapByScanNum(myMsDataFile, null, new CommonParameters(), out List<List<(double, int, double)>> precursorss).ToArray();

            //Generate crosslinker, which is UserDefined here.
            var crosslinker = xlSearchParameters.Crosslinker;

            //TwoPassCrosslinkSearchEngine.Run().
            List<CrosslinkSpectralMatch>[] possiblePsms = new List<CrosslinkSpectralMatch>[listOfSortedms2Scans.Length];
            List<(int, int, int)>[] candidates = new List<(int, int, int)>[listOfSortedms2Scans.Length];

            var XLEngine = new CrosslinkSearchEngine(possiblePsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, null, 0,
                commonParameters, null, crosslinker, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.CrosslinkAtCleavageSite,
                xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, new List<string> { },
                candidates, 0, indexResults.PeptideIndex, precursorss);
            XLEngine.FirstRoundSearch();
            XLEngine.Run();

            var newPsms = possiblePsms.Where(p => p != null).ToList();
            Assert.That(newPsms.Count, Is.EqualTo(1));
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
                CommonParameters = new CommonParameters(trimMsMsPeaks: false)
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

            List<Product> products = new List<Product>();
            deadendPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            double[] mz = products.Select(p => p.NeutralMass.ToMz(1)).OrderBy(v => v).ToArray();
            double[] intensities = new[] { 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

            MzSpectrum spectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan sc = new MsDataScan(spectrum, 1, 2, true, Polarity.Positive, 1, spectrum.Range, "",
                MZAnalyzerType.Orbitrap, 12, 1.0, null, null);
            scans[0] = new Ms2ScanWithSpecificMass(sc, deadendPeptide.MonoisotopicMass.ToMz(2), 2, "", new CommonParameters());
            Assert.That(scans[0].PrecursorIntensity, Is.EqualTo(1));
            Assert.That(scans[0].PrecursorEnvelopePeakCount, Is.EqualTo(1));

            List<List<(double, int, double)>> precursorss = new List<List<(double, int, double)>> 
                { new List<(double, int, double)> { (scans[0].PrecursorMass, scans[0].PrecursorCharge, scans[0].PrecursorMonoisotopicPeakMz)} };


            // search the data with the peptide WITHOUT the deadend mod annotated in the search database.
            // the search engine should be able to correctly identify the deadend mod on T
            var indexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { protein }, new List<Modification>(), new List<Modification>(), null, null, null,
                0, DecoyType.None, new CommonParameters(), null, 1000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();

            List<(int, int, int)>[] candidates = new List<(int, int, int)>[scans.Length];
            var XLEngine = new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, null, 0,
                new CommonParameters(), null, crosslinker, 50, true, false, false, true, new List<string>(),
                candidates, 0, indexingResults.PeptideIndex, precursorss);
            XLEngine.FirstRoundSearch();
            XLEngine.Run();

            CrosslinkSpectralMatch csm = csms.First().First();
            csm.ResolveAllAmbiguities();
            CrosslinkSpectralMatch.ResolveProteinPosAmbiguitiesForXl(csm);
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
        public static void XLSearchTastWriteFileTest()
        {
            //sending zero CSMs to write doesn't error. Method Simply Returns.
            WriteXlFile.WritePsmCrossToTsv(new List<CrosslinkSpectralMatch>(), "", 0);

            //sending the wrong writeType doesn't error. Method Simply breaks.
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\TestXLWrite");
            Directory.CreateDirectory(outputFolder);
            
            Ms2ScanWithSpecificMass scan = new(new MsDataScan(
                    new MzSpectrum(new double[] { }, new double[] { }, false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);
            Dictionary<int, Modification> mod = new();

            string sequence = "PEPTIDE";
            Dictionary<string, Modification> allKnownMods = new();
            int numFixedMods = 0;
            DigestionParams digestionParams = null;
            Protein protForward = new Protein("PEPTIDE", "ACCESSION", isDecoy: true);
            int oneBasedStartResidueInProtein = int.MinValue;
            int oneBasedEndResidueInProtein = int.MinValue;
            int missedCleavages = int.MinValue;
            CleavageSpecificity cleavageSpecificity = CleavageSpecificity.Full;
            string peptideDescription = null;

            PeptideWithSetModifications alpha = new(sequence: sequence, allKnownMods, numFixedMods, digestionParams, protForward, oneBasedStartResidueInProtein,
                oneBasedEndResidueInProtein, missedCleavages, cleavageSpecificity, peptideDescription);
            PeptideWithSetModifications beta = new(sequence: sequence, allKnownMods, numFixedMods, digestionParams, protForward, oneBasedStartResidueInProtein,
                oneBasedEndResidueInProtein, missedCleavages, cleavageSpecificity, peptideDescription);

            CrosslinkSpectralMatch csmAlpha = new(alpha, 0, 3, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            CrosslinkSpectralMatch csmBeta = new(beta, 0, 3, 1, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            csmAlpha.ResolveAllAmbiguities();
            csmBeta.ResolveAllAmbiguities();

            csmAlpha.LinkPositions = new() { 1 };
            csmBeta.LinkPositions = new() { 1 };

            WriteXlFile.WritePsmCrossToTsv(new List<CrosslinkSpectralMatch>() { csmAlpha}, outputFolder + "csm.psmtsv", 0);

            //check decoy label
            csmAlpha.BetaPeptide = csmBeta;
            Crosslinker xlinker = Crosslinker.ParseCrosslinkerFromString("DSSO\tK\tK\tT\tCID|HCD\t158.0038\t54.01056\t85.982635\t176.0143\t175.0303\t279.0777");
            double someSortOfMass = csmAlpha.ScanPrecursorMass - csmAlpha.BetaPeptide.BioPolymerWithSetModsMonoisotopicMass.Value - csmAlpha.BioPolymerWithSetModsMonoisotopicMass.Value - xlinker.TotalMass;
            string someLongString = "-." + csmAlpha.FullSequence + csmAlpha.LinkPositions.First().ToString(CultureInfo.InvariantCulture) + "--" + csmAlpha.BetaPeptide.FullSequence + csmAlpha.BetaPeptide.LinkPositions.First().ToString(CultureInfo.InvariantCulture) + ".-";
            string someOtherString = csmAlpha.BestMatchingBioPolymersWithSetMods.First().WithSetMods.Parent.Accession.ToString(CultureInfo.InvariantCulture)
                                   + "(" + (csmAlpha.XlProteinPos.HasValue ? csmAlpha.XlProteinPos.Value.ToString(CultureInfo.InvariantCulture) : string.Empty) + ")";
            string lastRandomString = csmAlpha.BetaPeptide.BestMatchingBioPolymersWithSetMods.First().WithSetMods.Parent.Accession.ToString(CultureInfo.InvariantCulture)
                                   + "(" + (csmAlpha.BetaPeptide.XlProteinPos.HasValue ? csmAlpha.BetaPeptide.XlProteinPos.Value.ToString(CultureInfo.InvariantCulture) : string.Empty) + ")";

            Assert.That(csmAlpha.ScanRetentionTime.ToString(CultureInfo.InvariantCulture), Is.EqualTo("NaN"));
            Assert.That(csmAlpha.ScanNumber.ToString(CultureInfo.InvariantCulture), Is.EqualTo("2"));
            Assert.That(csmAlpha.XLTotalScore.ToString(CultureInfo.InvariantCulture), Is.EqualTo("3"));
            Assert.That(csmAlpha.DeltaScore.ToString(CultureInfo.InvariantCulture), Is.EqualTo("3"));
            Assert.That(csmAlpha.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture), Is.EqualTo("1"));
            Assert.That(csmAlpha.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture), Is.EqualTo("98.992723533121"));
            Assert.That(someSortOfMass, Is.EqualTo(-1657.7310045328791));
            Assert.That(csmAlpha.BetaPeptide.BaseSequence.Length.ToString(CultureInfo.InvariantCulture), Is.EqualTo("7"));
            Assert.That(csmAlpha.BaseSequence.Length.ToString(CultureInfo.InvariantCulture), Is.EqualTo("7"));
            Assert.That((csmAlpha.BetaPeptide.BaseSequence.Length + csmAlpha.BaseSequence.Length).ToString(CultureInfo.InvariantCulture), Is.EqualTo("14"));
            Assert.That(someLongString, Is.EqualTo("-.PEPTIDE1--PEPTIDE1.-"));
            Assert.That(someOtherString, Is.EqualTo("ACCESSION()"));
            Assert.That(lastRandomString, Is.EqualTo("ACCESSION()"));

            WriteXlFile.WriteCrosslinkToTxtForPercolator(new List<CrosslinkSpectralMatch>() { csmAlpha }, outputFolder, "percolator.tsv", xlinker);
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestWriteToPercolator()
        {
            XLSearchTask xlst = new XLSearchTask()
            {
                XlSearchParameters = new XlSearchParameters
                {
                    WriteOutputForPercolator = true
                },
                CommonParameters = new CommonParameters(trimMsMsPeaks: false, addCompIons:false)
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
            Assert.That(lines[0].Equals("SpecId\tLabel\tScannr\tScore\tdScore\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProtein"));

            Assert.That(lines[1].Equals("T-1-30.61909926666667\t1\t1\t26.06004534434461\t11.026813997502483\t3\t1994.0520231384269\t0.6649793543976755\t7\t9\t16\t-.EKVLTSSAR2--LSQKFPK4.-\tP02769(211)\tP02769(245)"));
            Directory.Delete(outputFolder, true);
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

            List<Product> products = new List<Product>();
            deadendPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            double[] mz = products.Select(p => p.NeutralMass.ToMz(1)).OrderBy(v => v).ToArray();
            double[] intensities = new[] { 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

            MzSpectrum spectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan sc = new MsDataScan(spectrum, 1, 2, true, Polarity.Positive, 1, spectrum.Range, "",
                MZAnalyzerType.Orbitrap, 12, 1.0, null, null);
            scans[0] = new Ms2ScanWithSpecificMass(sc, deadendPeptide.MonoisotopicMass.ToMz(2), 2, "", new CommonParameters());
            List<List<(double, int, double)>> precursorss = new List<List<(double, int, double)>>
                { new List<(double, int, double)> { (scans[0].PrecursorMass, scans[0].PrecursorCharge, scans[0].PrecursorMonoisotopicPeakMz)} };

            var indexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { protein }, new List<Modification>(), new List<Modification>(), null, null, null, 0, DecoyType.None,
                new CommonParameters(), null, 1000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();

            List<(int, int, int)>[] candidates = new List<(int, int, int)>[scans.Length];
            var XLEngine = new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, null, 0,
                new CommonParameters(), null, crosslinker, 50, true, false, false, true, new List<string>(),
                candidates, 0, indexingResults.PeptideIndex, precursorss);
            XLEngine.FirstRoundSearch();
            XLEngine.Run();

            csms[0].First().ResolveAllAmbiguities();
            csms[0].First().SetFdrValues(0, 0, 0.1, 0, 0, 0, 0, 0);

            WriteXlFile.WritePepXML_xl(csms.SelectMany(p => p).ToList(), new List<Protein>(), "", new List<Modification> { deadend }, new List<Modification> { deadend }, new List<string>(), TestContext.CurrentContext.TestDirectory, "test", new CommonParameters(), new XlSearchParameters { Crosslinker = crosslinker });
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test.pep.XML"));
        }

        [Test]
        public static void TestMixedMs2Ms2()
        {
            string outputFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMixedMs2Ms2.tsv");
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.CID, ms2childScanDissociationType: DissociationType.ETD,
                trimMsMsPeaks: false);

            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\ms2mixed_bsa_xlink.mzML");

            var fsp = new List<(string, CommonParameters)>();
            fsp.Add((spectraFile, commonParameters));

            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);

            Ms2ScanWithSpecificMass[] scans = MetaMorpheusTask.GetMs2ScansWrapByScanNum(file, spectraFile, commonParameters, out List<List<(double, int, double)>> precursorss).ToArray();

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
                null, null, null, 0, DecoyType.None, commonParameters, fsp, 5000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();

            var secondCombinedParams = new CommonParameters(dissociationType: DissociationType.ETD, ms2childScanDissociationType: DissociationType.ETD,
                trimMsMsPeaks: false);

            var fsp2 = new List<(string, CommonParameters)>();
            fsp2.Add((spectraFile, secondCombinedParams));

            var secondIndexingResults = (IndexingResults)new IndexingEngine(new List<Protein> { bsa }, new List<Modification>(), new List<Modification>(),
                null, null, null, 0, DecoyType.None, secondCombinedParams, fsp2, 5000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();

            var csms = new List<CrosslinkSpectralMatch>[1];

            List<(int, int, int)>[] candidates = new List<(int, int, int)>[scans.Length];
            var XLEngine = new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, secondIndexingResults.FragmentIndex, 0,
                commonParameters, fsp, GlobalVariables.Crosslinkers.First(p => p.CrosslinkerName == "DSSO"), 50, true, false, false, true, new List<string>(),
                candidates, 0, indexingResults.PeptideIndex, precursorss);
            XLEngine.FirstRoundSearch();
            XLEngine.Run();

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
            var isIntra = CrosslinkSpectralMatch.IsIntraCsm(csm);
            Assert.That(isIntra == true);
            csm.ResolveProteinPosAmbiguitiesForXl();
            Assert.That(csm.XlProteinPos == 455 && csm.BetaPeptide.XlProteinPos == 211);

            // test parent scan (CID)
            Assert.That(csm.MatchedFragmentIons.Count == 21);
            Assert.That(csm.BetaPeptide.MatchedFragmentIons.Count == 13);
            Assert.That(csm.ScanNumber == 2);

            // test child scan (ETD)
            Assert.That(csm.ChildMatchedFragmentIons.First().Key == 3);
            Assert.That(csm.ChildMatchedFragmentIons.First().Value.Count == 21);
            Assert.That(csm.BetaPeptide.ChildMatchedFragmentIons.First().Key == 3);
            Assert.That(csm.BetaPeptide.ChildMatchedFragmentIons.First().Value.Count == 25);

            // write results to TSV
            csm.SetFdrValues(1, 0, 0, 0, 0, 0, 0, 0);
            WriteXlFile.WritePsmCrossToTsv(new List<CrosslinkSpectralMatch> { csm }, outputFile, 2);

            // read results from TSV
            var psmFromTsv = PsmTsvReader.ReadTsv(outputFile, out var warnings).First();

            Assert.That(psmFromTsv.ChildScanMatchedIons.Count, Is.EqualTo(1));
            Assert.That(psmFromTsv.ChildScanMatchedIons.First().Key, Is.EqualTo(3));
            Assert.That(psmFromTsv.ChildScanMatchedIons.First().Value.Count, Is.EqualTo(21));

            Assert.That(psmFromTsv.BetaPeptideChildScanMatchedIons.Count, Is.EqualTo(1));
            Assert.That(psmFromTsv.BetaPeptideChildScanMatchedIons.First().Key, Is.EqualTo(3));
            Assert.That(psmFromTsv.BetaPeptideChildScanMatchedIons.First().Value.Count, Is.EqualTo(25));
            // Race condition causes the order of accessions to differ 
            //Assert.That(psmFromTsv.BetaPeptideProteinAccession, Is.EqualTo("BSA|BSA2"));
            Assert.That(psmFromTsv.BetaPeptideProteinAccession, Does.Contain("BSA"));
            Assert.That(psmFromTsv.BetaPeptideProteinAccession, Does.Contain("BSA2"));
            Assert.That(psmFromTsv.BetaPeptideProteinLinkSite, Is.EqualTo(211));
            Assert.That(psmFromTsv.BetaPeptideTheoreticalMass, Is.EqualTo("989.550558768"));

            Assert.That(psmFromTsv.CrossType, Is.EqualTo("Cross"));
            Assert.That(psmFromTsv.LinkResidues, Is.EqualTo("K"));
            Assert.That(psmFromTsv.ProteinLinkSite, Is.EqualTo(455));
            Assert.That(psmFromTsv.ParentIons, Is.EqualTo("2;2"));

            Assert.That(csm.Accession == null && csm.BetaPeptide.Accession == null);

            // Race condition causes the order of accessions to differ 
            //Assert.That(psmFromTsv.ProteinAccession == "BSA|BSA2");
            Assert.That(psmFromTsv.ProteinAccession, Does.Contain("BSA"));
            Assert.That(psmFromTsv.ProteinAccession, Does.Contain("BSA2"));
            Assert.That(psmFromTsv.UniqueSequence, Is.EqualTo(psmFromTsv.FullSequence + psmFromTsv.BetaPeptideFullSequence));

            File.Delete(outputFile);
        }

        [Test]
        public static void TestMs2Ms3()
        {
            string outputFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMs2Ms3.tsv");
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.CID, ms3childScanDissociationType: DissociationType.LowCID, precursorMassTolerance: new PpmTolerance(10));

            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\10226.mzML");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);

            var fsp = new List<(string, CommonParameters)>();
            fsp.Add((spectraFile, commonParameters));

            Ms2ScanWithSpecificMass[] scans = MetaMorpheusTask.GetMs2ScansWrapByScanNum(file, spectraFile, commonParameters, out List<List<(double, int, double)>> precursorss).ToArray();

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
               null, null, null, 0, DecoyType.None, commonParameters, fsp, 5000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();

            var csms = new List<CrosslinkSpectralMatch>[2];

            List<(int, int, int)>[] candidates = new List<(int, int, int)>[scans.Length];
            var XLEngine = new CrosslinkSearchEngine(csms, scans, indexingResults.PeptideIndex, indexingResults.FragmentIndex, null, 0,
                commonParameters, fsp, GlobalVariables.Crosslinkers.First(p => p.CrosslinkerName == "DSSO"), 50, true, false, false, true, new List<string>(),
                candidates, 0, indexingResults.PeptideIndex, precursorss);
            XLEngine.FirstRoundSearch();
            XLEngine.Run();

            var csm = csms[0].First();
            csm.ResolveAllAmbiguities();
            if (csm.BetaPeptide != null)
            {
                csm.BetaPeptide.ResolveAllAmbiguities();
            }
            // test parent scan (CID)
            Assert.That(csm.MatchedFragmentIons.Count, Is.EqualTo(36));
            Assert.That(csm.ScanNumber == 2);

            // test child scan (low-resolution CID, alpha peptide signature ion)
            Assert.That(csm.ChildMatchedFragmentIons.First().Key == 4);
            Assert.That(csm.ChildMatchedFragmentIons.First().Value.Count == 1);

            // test child scan (low-resolution CID, beta peptide signature ion)
            Assert.That(csm.BetaPeptide.ChildMatchedFragmentIons.First().Key == 4);
            Assert.That(csm.BetaPeptide.ChildMatchedFragmentIons.First().Value.Count == 5);

            // write results to TSV
            csm.SetFdrValues(1, 0, 0, 0, 0, 0, 0, 0);
            WriteXlFile.WritePsmCrossToTsv(new List<CrosslinkSpectralMatch> { csm }, outputFile, 2);

            // read results from TSV
            var psmFromTsv = PsmTsvReader.ReadTsv(outputFile, out var warnings).First();

            Assert.That(psmFromTsv.ChildScanMatchedIons.Count == 4
                && psmFromTsv.ChildScanMatchedIons.First().Key == 4
                && psmFromTsv.ChildScanMatchedIons.First().Value.Count == 1);

            Assert.That(psmFromTsv.BetaPeptideChildScanMatchedIons.Count == 4
                && psmFromTsv.BetaPeptideChildScanMatchedIons.First().Key == 4
                && psmFromTsv.BetaPeptideChildScanMatchedIons.First().Value.Count == 5);

            File.Delete(outputFile);
        }

        [Test]
        [TestCase(1.0d, 1)]//below smallest number and outside tolerance returns index 1 because value 5 is closest
        [TestCase(30.0d, 10)]//above highest number and outside tolerance returns index 10, which is not an acceptable index.
        [TestCase(4.9999999995d, 1)]//below smallest number and within tolerance returns index 1 because value 5 is closest. the fact that this is within tolerance is meaningless
        [TestCase(5.0000000005, 3)]//above smallest number and within tolerance returns 3 because its the closest number greater.
        [TestCase(7, 3)]//between numbers outside tolerance returns 3 because its the closest number greater.
        [TestCase(10.000000000, 4)]//between numbers exactly at tolerance but matches a number returns 4
        [TestCase(9.9999999999, 4)]//between numbers exactly at tolerance but matches a number returns 4 because 10 is the closest above.
        [TestCase(10.00000000005, 5)]//exactly between numbers and within tolerance of both returns 5 because its the larger.
        [TestCase(11, 6)]//when two numbers in the array are the same and the lookup matches them both you should get the index of the lowest
        public static void Test_BinarySearchGetIndex(double targetMass, int arrayIndex)
        {
            double[] massArray = new double[] { Double.NaN, 5, Double.NaN, 9.999999999, 10.000000000, 10.000000001, 11, 11, Double.NaN, Double.NaN };
            Assert.That(CrosslinkSearchEngine.BinarySearchGetIndex(massArray, targetMass), Is.EqualTo(arrayIndex));
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

                //Cross
                var mz2 = new double[] { 100, 201.1234, 244.1656, 391.2340, 420.2201, 521.2678, 634.3519, 889.965, 1044.568, 1094.551, 1279.671, 1378.74, 1491.824 };
                var intensities2 = new double[] { 100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
                var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
                ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1.0, null, "scan=2", 1994.05.ToMz(3),
                    3, 1, 1994.05.ToMz(3), 2, DissociationType.HCD, 1, 1994.05.ToMz(3)));

                //single
                var mz3 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
                var intensities3 = new double[] { 100, 1, 1, 1 };
                var MassSpectrum3 = new MzSpectrum(mz3, intensities3, false);
                ScansHere.Add(new MsDataScan(MassSpectrum3, 3, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=3", 846.4963.ToMz(1),
                    1, 1, 846.4963.ToMz(1), 2, DissociationType.HCD, 1, 846.4963.ToMz(1)));

                //loop
                var mz4 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
                var intensities4 = new double[] { 100, 1, 1, 1 };
                var MassSpectrum4 = new MzSpectrum(mz4, intensities4, false);
                ScansHere.Add(new MsDataScan(MassSpectrum4, 4, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=4", 1004.491.ToMz(1),
                    1, 1, 1004.491.ToMz(1), 2, DissociationType.HCD, 1, 1004.491.ToMz(1)));

                //Deadend
                var mz5 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
                var intensities5 = new double[] { 100, 1, 1, 1 };
                var MassSpectrum5 = new MzSpectrum(mz5, intensities5, false);
                ScansHere.Add(new MsDataScan(MassSpectrum5, 5, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=5", 1022.511.ToMz(1),
                    1, 1, 1022.511.ToMz(1), 2, DissociationType.HCD, 1, 1022.511.ToMz(1)));

                Scans = ScansHere.ToArray();
            }

            public new static string FilePath
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

            #region MsDataFile Abstract Methods

            public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
            {
                throw new NotImplementedException();
            }

            public override SourceFile GetSourceFile()
            {
                throw new NotImplementedException();
            }

            public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
            {
                throw new NotImplementedException();
            }

            public override void CloseDynamicConnection()
            {
                throw new NotImplementedException();
            }

            public override void InitiateDynamicConnection()
            {
                throw new NotImplementedException();
            }

            public override MsDataScan GetOneBasedScan(int scanNumber)
            {
                return Scans[scanNumber - 1];
            }

            public override IEnumerable<MsDataScan> GetMS1Scans()
            {
                throw new NotImplementedException();
            }

            #endregion
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

            public new string FilePath
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

            #region MsDataFile Abstract Methods

            public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
            {
                throw new NotImplementedException();
            }

            public override SourceFile GetSourceFile()
            {
                throw new NotImplementedException();
            }

            public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
            {
                throw new NotImplementedException();
            }

            public override void CloseDynamicConnection()
            {
                throw new NotImplementedException();
            }

            public override void InitiateDynamicConnection()
            {
                throw new NotImplementedException();
            }

            public override MsDataScan GetOneBasedScan(int scanNumber)
            {
                return Scans[scanNumber - 1];
            }

            public override IEnumerable<MsDataScan> GetMS1Scans()
            {
                throw new NotImplementedException();
            }

            #endregion
        }
    }
}