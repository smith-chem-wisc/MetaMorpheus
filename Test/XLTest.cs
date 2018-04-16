using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.CrosslinkAnalysis;
using EngineLayer.Indexing;
using IO.MzML;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using Nett;
using System;
using System.IO;


namespace Test
{
    [TestFixture]
    public static class XLTest
    {
        #region Public Methods

        [Test]
        public static void XlTestXlPosCal()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            DigestionParams digestionParams = new DigestionParams
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                MaxMissedCleavages = 2,
                Protease = protease,
                MinPeptideLength = 1
            };
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();

            var ye = prot.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();

            var pep = ye[0];
            Assert.AreEqual(pep.BaseSequence, "MNNNK");
            CrosslinkerTypeClass crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(CrosslinkerType.DSS);
            Assert.AreEqual(crosslinker.CrosslinkerModSites, "K");
            Assert.AreEqual(Residue.GetResidue(crosslinker.CrosslinkerModSites).MonoisotopicMass, 128.09496301518999, 1e-9);
            var n = pep.CompactPeptide(TerminusType.None).NTerminalMasses;
            var c = pep.CompactPeptide(TerminusType.None).CTerminalMasses;
            Assert.AreEqual(n.Count(), 4);
            Assert.AreEqual(c.Count(), 4);
            Assert.AreEqual(c[0], 128.09496301518999, 1e-6);
            var x = PsmCross.XlPosCal(pep.CompactPeptide(TerminusType.None), crosslinker.CrosslinkerModSites).ToArray();
            Assert.AreEqual(x[0], 4);

            var pep2 = ye[2];
            Assert.AreEqual("MNNNKQQQQ", pep2.BaseSequence);
            var n2 = pep2.CompactPeptide(TerminusType.None).NTerminalMasses;
            var c2 = pep2.CompactPeptide(TerminusType.None).CTerminalMasses;
            Assert.AreEqual(n2.Count(), 8);
            Assert.AreEqual(c2.Count(), 8);
            Assert.AreEqual(n2[4] - n2[3], 128.09496301518999, 1e-6);
            var x2 = PsmCross.XlPosCal(pep2.CompactPeptide(TerminusType.None), crosslinker.CrosslinkerModSites).ToArray();
            Assert.AreEqual(x2[0], 4);

            //Test crosslinker with multiple types of mod
            var protSTC = new Protein("GASTACK", null);
            var peps = protSTC.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();
            var pepSTC = peps[0];
            Assert.AreEqual(pepSTC.BaseSequence, "GASTACK");
            CrosslinkerTypeClass crosslinker2 = new CrosslinkerTypeClass("ST","C", "crosslinkerSTC", false, -18.01056, 0, 0, 0, 0, 0, 0);
            string crosslinkerModSitesAll = new string((crosslinker2.CrosslinkerModSites + crosslinker2.CrosslinkerModSites2).ToCharArray().Distinct().ToArray());
            Assert.AreEqual(crosslinkerModSitesAll, "STC");
        }

        [Test]
        public static void XlTestGenerateIntensityRanks()
        {
            double[] mz = new double[] { 1.0, 1.3, 1.5, 1.7, 1.9, 2.1 };
            double[] intensity = new double[] { 1.1, 1.1, 0.5, 3.2, 0.5, 6.0 };
            int[] rank = PsmCross.GenerateIntensityRanks(mz, intensity);
            int[] Rank = new int[] { 4, 3, 6, 2, 5, 1 };
            Assert.AreEqual(rank, Rank);
        }

        [Test]
        public static void XlTest_BSA_DSSO()
        {
            //Generate parameters
            var commonParameters = new CommonParameters
            {
                DoPrecursorDeconvolution = false,
                CIons = true,
                ZdotIons = true,
                ScoreCutoff = 2,
                DigestionParams = new DigestionParams
                {
                    MinPeptideLength = 5
                }
            };

            var xlSearchParameters = new XlSearchParameters { XlCharge_2_3_PrimeFragment = true };

            //Create databases contain two protein. 
            var proteinList = new List<Protein> { new Protein("EKVLTSSAR", "Fake01"), new Protein("LSQKFPK", "Fake02") };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Carbamidomethyl of C", "Common Fixed", motif2, TerminusLocalization.Any, 57.02146372068994);
            var variableModifications = new List<ModificationWithMass>() { mod1 };
            var fixedModifications = new List<ModificationWithMass>() { mod2 };
            var localizeableModifications = new List<ModificationWithMass>();

            var lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y, ProductType.C, ProductType.Zdot };
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
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
            
            foreach (var fdfd in digestedList)
            {
                fdfd.CompactPeptide(TerminusType.None);
            }

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, lp, 1, DecoyType.Reverse, new List<IDigestionParams> { commonParameters.DigestionParams }, commonParameters, 30000, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            var fragmentIndexCount = indexResults.FragmentIndex.Count(p => p != null);
            var fragmentIndexAll = indexResults.FragmentIndex.Select((s,j)=> new {j, s }).Where(p => p.s != null).Select(t=> t.j).ToList();
            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFile();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, commonParameters.DoPrecursorDeconvolution, commonParameters.UseProvidedPrecursorInfo, commonParameters.DeconvolutionIntensityRatio, commonParameters.DeconvolutionMaxAssumedChargeState, commonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            //Generate crosslinker, which is DSSO here.
            CrosslinkerTypeClass crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(xlSearchParameters.CrosslinkerType);

            //TwoPassCrosslinkSearchEngine.Run().
            List<PsmCross> newPsms = new List<PsmCross>();
            new TwoPassCrosslinkSearchEngine(newPsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, lp, 0, commonParameters, false, xlSearchParameters.XlPrecusorMsTl, crosslinker, xlSearchParameters.CrosslinkSearchTop, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, xlSearchParameters.XlCharge_2_3, xlSearchParameters.XlCharge_2_3_PrimeFragment, new List<string> { }).Run();

            var compactPeptideToProteinPeptideMatch = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            MetaMorpheusEngineResults allcrosslinkanalysisResults;
            allcrosslinkanalysisResults = new CrosslinkAnalysisEngine(newPsms, compactPeptideToProteinPeptideMatch, proteinList, variableModifications, fixedModifications, lp, null, crosslinker, TerminusType.None, commonParameters, new List<string> { }).Run();
            foreach (var item in newPsms)
            {
                item.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            }
            //Test newPsms
            Assert.AreEqual(newPsms.Count(), 3);

            //Test Output

            var task = new XLSearchTask();
            task.WriteAllToTsv(newPsms, TestContext.CurrentContext.TestDirectory, "allPsms", new List<string> { });
            task.WritePepXML_xl(newPsms, proteinList, null, variableModifications, fixedModifications, null, TestContext.CurrentContext.TestDirectory, "pep.XML", new List<string> { });
            task.WriteSingleToTsv(newPsms.Where(p => p.CrossType == PsmCrossType.Singe).ToList(), TestContext.CurrentContext.TestDirectory, "singlePsms", new List<string> { });


            //Test PsmCross.XlCalculateTotalProductMasses.
            var psmCrossAlpha = new PsmCross(digestedList[1].CompactPeptide(TerminusType.None), 0, 0, i, listOfSortedms2Scans[0]);
            var psmCrossBeta = new PsmCross(digestedList[2].CompactPeptide(TerminusType.None), 0, 0, i, listOfSortedms2Scans[0]);
            var linkPos = PsmCross.XlPosCal(psmCrossAlpha.compactPeptide, crosslinker.CrosslinkerModSites);
            var productMassesAlphaList = PsmCross.XlCalculateTotalProductMasses(psmCrossAlpha, psmCrossBeta.compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.TotalMass, crosslinker, lp, true, false, linkPos);
            Assert.AreEqual(productMassesAlphaList[0].ProductMz.Length, 99);

        }

        [Test]
        public static void XlTest_BSA_DSS_file()
        {
            var task = Toml.ReadFile<XLSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/XLSearchTaskconfig_BSA_DSS_23747.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/BSA.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/BSA_DSS_23747.mzML");
            EverythingRunnerEngine a = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"XlTestData"));

            a.Run();
        }

        [Test]
        public static void XlTest_GenerateUserDefinedCrosslinker()
        {
            XlSearchParameters xlSearchParameters = new XlSearchParameters();
            xlSearchParameters.CrosslinkerType = CrosslinkerType.UserDefined;
            xlSearchParameters.UdXLkerName = "CrossST-C";
            xlSearchParameters.UdXLkerResidues = "ST";
            xlSearchParameters.UdXLkerResidues2 = "C";
            xlSearchParameters.UdXLkerTotalMass = -18.01056; 
            var crosslinker = XLSearchTask.GenerateUserDefinedCrosslinker(xlSearchParameters);
            Assert.AreEqual(crosslinker.DeadendMassH2O, (double)9999);

        }

        [Test]
        public static void XlTest_DiffCrosslinkSites()
        {
            //Generate parameters
            var commonParameters = new CommonParameters
            {
                DoPrecursorDeconvolution = false,
                ScoreCutoff = 1,
                DigestionParams = new DigestionParams
                {
                    MinPeptideLength = 4
                }
            };

            var xlSearchParameters = new XlSearchParameters
            {
                CrosslinkerType = CrosslinkerType.UserDefined,
                UdXLkerName = "CrossST-C",
                UdXLkerResidues = "ST",
                UdXLkerResidues2 = "C",
                UdXLkerTotalMass = -18.01056
            };

            //Create databases contain two protein. 
            var proteinList = new List<Protein> { new Protein("VLTAR", "Fake01"), new Protein("LCQK", "Fake02") };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            var variableModifications = new List<ModificationWithMass>() { mod1 };
            var fixedModifications = new List<ModificationWithMass>() {  };
            var localizeableModifications = new List<ModificationWithMass>();

            var lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y};
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

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

            foreach (var fdfd in digestedList)
            {
                fdfd.CompactPeptide(TerminusType.None);
            }

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, lp, 1, DecoyType.Reverse, new List<IDigestionParams> { commonParameters.DigestionParams }, commonParameters, 30000, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            var fragmentIndexCount = indexResults.FragmentIndex.Count(p => p != null);
            var fragmentIndexAll = indexResults.FragmentIndex.Select((s, j) => new { j, s }).Where(p => p.s != null).Select(t => t.j).ToList();
            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFileDiffSite();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, commonParameters.DoPrecursorDeconvolution, commonParameters.UseProvidedPrecursorInfo, commonParameters.DeconvolutionIntensityRatio, commonParameters.DeconvolutionMaxAssumedChargeState, commonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            //Generate crosslinker, which is UserDefined here.
            var crosslinker = XLSearchTask.GenerateUserDefinedCrosslinker(xlSearchParameters);

            //TwoPassCrosslinkSearchEngine.Run().
            List<PsmCross> newPsms = new List<PsmCross>();
            new TwoPassCrosslinkSearchEngine(newPsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, lp, 0, commonParameters, false, xlSearchParameters.XlPrecusorMsTl, crosslinker, xlSearchParameters.CrosslinkSearchTop, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, xlSearchParameters.XlCharge_2_3, xlSearchParameters.XlCharge_2_3_PrimeFragment, new List<string> { }).Run();
            Assert.AreEqual(newPsms.Count(), 1);
        }

        #endregion Public Methods
    }

    internal class XLTestDataFile : MsDataFile<IMzmlScan>
    {
        #region Public Constructors
        //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match. 
        public XLTestDataFile() : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 1994.05.ToMz(3), 846.4963.ToMz(1), 1005.498.ToMz(1), 1093.544.ToMz(1), 1043.561.ToMz(1) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1};
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, "scan=1") };

            var mz2 = new double[] { 100, 201.1234, 244.1656, 391.2340, 420.2201, 521.2678, 634.3519, 889.965, 1044.568, 1094.551, 1279.671, 1378.74, 1491.824 };
            var intensities2 = new double[] {100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzmlMzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1994.05.ToMz(3), 
                3, 1, 1994.05.ToMz(3), 2, DissociationType.HCD, 1, 1994.05.ToMz(3), 1.0, "scan=2"));

            var mz3 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
            var intensities3 = new double[] { 100, 1, 1, 1 };
            var MassSpectrum3 = new MzmlMzSpectrum(mz3, intensities3, false);
            ScansHere.Add(new MzmlScanWithPrecursor(3, MassSpectrum3, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 846.4963.ToMz(1),
                1, 1, 846.4963.ToMz(1), 2, DissociationType.HCD, 1, 846.4963.ToMz(1), 1.0, "scan=3"));

            var mz4 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
            var intensities4 = new double[] { 100, 1, 1, 1 };
            var MassSpectrum4 = new MzmlMzSpectrum(mz3, intensities4, false);
            ScansHere.Add(new MzmlScanWithPrecursor(4, MassSpectrum3, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1005.498.ToMz(1),
                1, 1, 1005.498.ToMz(1), 2, DissociationType.HCD, 1, 1005.498.ToMz(1), 1.0, "scan=4"));

            Scans = ScansHere.ToArray();
        }

        #endregion Public Constructors

        #region Public Properties

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

        #endregion Public Properties

        #region Public Methods

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            Scans[0] = new MzmlScan(Scans[0].OneBasedScanNumber, massSpectrum, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, Scans[0].NativeId);
        }

        public override IMzmlScan GetOneBasedScan(int scanNumber)
        {
            return Scans[scanNumber - 1];
        }

        public override IEnumerable<IMzmlScan> GetMS1Scans()
        {
            throw new System.NotImplementedException();
        }

        #endregion Public Methods
    }

    internal class XLTestDataFileDiffSite : MsDataFile<IMzmlScan>
    {
        #region Public Constructors
        //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match. 
        public XLTestDataFileDiffSite() : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 100, 1030.5956.ToMz(1) };
            var intensities1 = new double[] { 100, 1 };
            var MassSpectrum1 = new MzmlMzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<IMzmlScan> { new MzmlScan(1, MassSpectrum1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, "scan=1") };

            var mz2 = new double[] { 100, 147.1128, 175.119, 213.1598, 246.1561, 275.1714, 757.4388, 786.4541, 819.4504, 857.4912, 885.4974, 918.5189, 932.5345 };
            var intensities2 = new double[] { 100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzmlMzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MzmlScanWithPrecursor(2, MassSpectrum2, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1030.5956.ToMz(1),
                1, 1, 1030.5956.ToMz(1), 2, DissociationType.HCD, 1, 1030.5956.ToMz(1), 1.0, "scan=2"));

            Scans = ScansHere.ToArray();
        }

        #endregion Public Constructors

        #region Public Properties

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

        #endregion Public Properties

        #region Public Methods

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            Scans[0] = new MzmlScan(Scans[0].OneBasedScanNumber, massSpectrum, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, Scans[0].NativeId);
        }

        public override IMzmlScan GetOneBasedScan(int scanNumber)
        {
            return Scans[scanNumber - 1];
        }

        public override IEnumerable<IMzmlScan> GetMS1Scans()
        {
            throw new System.NotImplementedException();
        }

        #endregion Public Methods
    }
}