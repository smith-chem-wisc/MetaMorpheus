using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using MzLibUtil;
using Nett;
using EngineLayer.GlycoSearch;
using MathNet.Numerics.LinearRegression;

namespace Test
{
    [TestFixture]
    public class XLTestNOGlyco
    {
        private static TestDataFile myMsDataFile { get; set; }

        private static IndexingResults indexResults { get; set; }

        [OneTimeSetUp]
        public static void Setup()
        {
            //Generate parameters
            var digestionParameters = new DigestionParams(protease: "semi-trypsin");

            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, dissociationType: DissociationType.HCD,
                ms2childScanDissociationType: DissociationType.EThcD, digestionParams: digestionParameters);

            //Create databases.
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/P16150.fasta"), true, DecoyType.Reverse, false, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
        ProteinDbLoader.UniprotOrganismRegex, out var dbErrors, -1);

            //Create modification.
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            var variableModifications = new List<Modification>() { mod1 };
            var fixedModifications = new List<Modification>() { mod2 };

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse,
                commonParameters, null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            indexResults = (IndexingResults)indexEngine.Run();

            //Get MS2 scans.
            myMsDataFile = new TestDataFile();

            {
                string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\11898_AIETD.mgf");
                var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
                var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();
                myMsDataFile.ReplaceScanArrays(scan.TheScan.MassSpectrum.XArray, scan.TheScan.MassSpectrum.YArray, 1);
            }
            {
                string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\11901_AIETD.mgf");
                var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
                var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();
                myMsDataFile.ReplaceScanArrays(scan.TheScan.MassSpectrum.XArray, scan.TheScan.MassSpectrum.YArray, 2);
            }
            {
                string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4562.mgf");
                var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
                var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();
                myMsDataFile.ReplaceScanArrays(scan.TheScan.MassSpectrum.XArray, scan.TheScan.MassSpectrum.YArray, 3);
            }
            {
                string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
                var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
                var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();
                myMsDataFile.ReplaceScanArrays(scan.TheScan.MassSpectrum.XArray, scan.TheScan.MassSpectrum.YArray, 4);
            }

        }

        [Test]
        public static void NOGlycoTest_Complicated()
        {
            //Generate parameters
            var digestionParameters = new DigestionParams(protease: "semi-trypsin");
            
            var commonParameters = new CommonParameters(doPrecursorDeconvolution:false, dissociationType: DissociationType.HCD, 
                ms2childScanDissociationType: DissociationType.EThcD, digestionParams: digestionParameters);

            var _glycoSearchParameters = new GlycoSearchParameters
            {
                GlycoSearchType = GlycoSearchType.N_O_GlycanSearch
            };

            //Get MS2 scans
            var listOfMs2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, commonParameters).ToArray();

            //Run GlycoSearchEngine
            List<GlycoSpectralMatch>[] possiblePsms = new List<GlycoSpectralMatch>[listOfMs2Scans.Length];
            new GlycoSearchEngine(possiblePsms, listOfMs2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, null, 0, 
                commonParameters, null, _glycoSearchParameters.OGlycanDatabasefile, _glycoSearchParameters.NGlycanDatabasefile, _glycoSearchParameters.GlycoSearchType, _glycoSearchParameters.GlycoSearchTopNum,
                        _glycoSearchParameters.MaximumOGlycanAllowed, _glycoSearchParameters.OxoniumIonFilt, _glycoSearchParameters.IndexingChildScan, new List<string> { }).Run();

            var newPsms = possiblePsms.Where(p => p != null).Select(p => p.First()).ToList();

            Assert.That(newPsms[0].FullSequence == "TN[N-Glycosylation:H9N2 on N]SSFIQGFVDHVKEDC[Common Fixed:Carbamidomethyl on C]DR");
            Assert.That(newPsms[1].FullSequence == "T[O-Glycosylation:H1N1 on X]T[O-Glycosylation:H1N1 on X]GSLEPSS[O-Glycosylation:N1 on X]GASGPQVSSVK");


            Assert.That(newPsms[0].R138to144 < 1.0);
            var kind_ngly = GlycoSpectralMatch.GetKind(newPsms[0]);
            Assert.That(Glycan.GetKindString(kind_ngly) == "H9N2");


            Assert.That(newPsms[1].NGlycanMotifExist == false);
            var kind_ogly = GlycoSpectralMatch.GetKind(newPsms[1]);
            Assert.That(Glycan.GetKindString(kind_ogly) == "H2N3");
            var output = newPsms[1].ToString();
            Assert.That(output.Contains("Leukosialin"));

        }

        [Test]
        public static void NOGlycoTest_Task()
        {
            Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"));
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData/mzMLfile.mzML");
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, spectraFile, false);

            var task = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/GlycoSearchTaskconfig_HCD-pd-EThcD.toml"), MetaMorpheusTask.tomlConfig);
         
            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData/P16150.fasta"), false);

            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { spectraFile }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData")).Run();

            var resultsPath_nglyco = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData/Task/nglyco.psmtsv"));
            Assert.That(resultsPath_nglyco.Length == 2);

            var resultsPath_oglyco = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TESTGlycoData/Task/oglyco.psmtsv"));
            Assert.That(resultsPath_oglyco.Length == 2);

            Directory.Delete(Path.Combine(Environment.CurrentDirectory, @"TESTGlycoData"), true);
        }

        internal class TestDataFile : MsDataFile
        {
            //Create fake N-OGlyco MS data. Including: Full Scan, HCD-pd-EThcd (N-glycopeptide), HCD-pd-EThcD (O-glycopeptide), HCD (single peptide). 
            public TestDataFile() : base(2, new SourceFile(null, null, null, null, null))
            {
                var mz = new double[] {  };
                var intensities = new double[] { };         
                var ScansHere = new List<MsDataScan>();

                var mz1 = new double[] { 937.0963, 937.4288, 937.7654, 938.0994, 1030.4214, 1030.6719, 1030.9227, 1031.1735, 1031.4242, 1031.6741 };
                var intensities1 = new double[] { 2191194.3, 1620103.4, 3470064.5, 807049.2, 828013.6, 1617673.9, 2045691.8, 1464751.9, 825144.2, 599464.8 };
                var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
                ScansHere.Add(new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1,
                    new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1"));

                var MassSpectrum2 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1.0, null, "scan=2", 1030.421508789063,
                    4, 2045691.75, 1030.421508789063, 2, DissociationType.HCD, null, 1030.421508789063));

                var MassSpectrum3 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum3, 3, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=3", 1030.421508789063,
                    4, 2045691.75, 1030.421508789063, 2, DissociationType.ETD, 2, 1030.421508789063));


                var MassSpectrum4 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum4, 4, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=4", 937.096923828125,
                    3, 3470064.5, 937.096923828125, 2, DissociationType.HCD, null, 937.096923828125));

                var MassSpectrum5 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum5, 5, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=5", 937.096923828125,
                    3, 3470064.5, 937.096923828125, 2, DissociationType.ETD, 4, 937.096923828125));

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

            public void ReplaceScanArrays(double[] mz, double[] intensities, int i)
            {
                MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
                Scans[i] = new MsDataScan(massSpectrum, Scans[i].OneBasedScanNumber, Scans[i].MsnOrder, Scans[i].IsCentroid, 
                    Scans[i].Polarity, Scans[i].RetentionTime, Scans[i].ScanWindowRange, Scans[i].ScanFilter, Scans[i].MzAnalyzer, 
                    massSpectrum.SumOfAllY, Scans[i].InjectionTime, null, Scans[i].NativeId, Scans[i].SelectedIonMZ, Scans[i].SelectedIonChargeStateGuess,
                    Scans[i].SelectedIonIntensity, Scans[i].IsolationMz, Scans[i].IsolationWidth, Scans[i].DissociationType, Scans[i].OneBasedPrecursorScanNumber,
                    Scans[i].SelectedIonMonoisotopicGuessMz);
            }
        }

    }
}
