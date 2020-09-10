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

        ////This is not a real test. It is for the local development. Please comment it before push.
        //    [Test]
        //    public static void NOGlycoTest_SVM()
        //    {
        //        string fileIn = @"E:\MassData\Glycan\Nick_2019_StcE\Rep1\_temp2\2020-09-04-23-06-46\Task1-GlycoSearchTask\all.psmtsv";
        //        string fileOut = @"E:\MassData\Glycan\Nick_2019_StcE\Rep1\_temp2\2020-09-04-23-06-46\Task1-GlycoSearchTask\svm_report.txt";

        //        SVM.RunSVM(fileIn, fileOut);     
        //    }


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
            var localizeableModifications = new List<Modification>();

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse,
                commonParameters, null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            var indexedFragments = indexResults.FragmentIndex.Where(p => p != null).SelectMany(v => v).ToList();
   
            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFile();

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



            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, commonParameters).ToArray();

            //Run GlycoSearchEngine
            List<GlycoSpectralMatch>[] possiblePsms = new List<GlycoSpectralMatch>[listOfSortedms2Scans.Length];
            new GlycoSearchEngine(possiblePsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, null, 0, 
                commonParameters, null, _glycoSearchParameters.OGlycanDatabasefile, _glycoSearchParameters.NGlycanDatabasefile, _glycoSearchParameters.GlycoSearchType, _glycoSearchParameters.GlycoSearchTopNum,
                        _glycoSearchParameters.MaximumOGlycanAllowed, _glycoSearchParameters.OxoniumIonFilt, _glycoSearchParameters.IndexingChildScan, new List<string> { }).Run();

            var newPsms = possiblePsms.Where(p => p != null).Select(p => p.First()).ToList();

            Assert.That(newPsms[0].FullSequence == "TN[N-Glycosylation:H9N2 on N]SSFIQGFVDHVKEDC[Common Fixed:Carbamidomethyl on C]DR");
            Assert.That(newPsms[1].FullSequence == "T[O-Glycosylation:H1N1 on X]T[O-Glycosylation:H1N1 on X]GSLEPSS[O-Glycosylation:N1 on X]GASGPQVSSVK");
        }


        internal class XLTestDataFile : MsDataFile
        {
            //Create fake N-OGlyco MS data. Including: Full Scan, HCD-pd-EThcd (N-glycopeptide), HCD-pd-EThcD (O-glycopeptide), HCD (single peptide). 
            public XLTestDataFile() : base(2, new SourceFile(null, null, null, null, null))
            {
                var mz = new double[] {  };
                var intensities = new double[] { };         
                var ScansHere = new List<MsDataScan>();

                var MassSpectrum1 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add( new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, 
                    new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1"));

                var MassSpectrum2 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1.0, null, "scan=2", 1030.421508789063,
                    4, 2045691.75, 1030.421508789063, 2, DissociationType.HCD, null, 1030.421508789063));

                var MassSpectrum3 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum3, 3, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=3", 1030.421508789063,
                    4, 2045691.75, 1030.421508789063, 2, DissociationType.EThcD, 2, 1030.421508789063));


                var MassSpectrum4 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum4, 4, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=4", 937.096923828125,
                    3, 3470064.5, 937.096923828125, 2, DissociationType.HCD, null, 937.096923828125));

                var MassSpectrum5 = new MzSpectrum(mz, intensities, false);
                ScansHere.Add(new MsDataScan(MassSpectrum5, 5, 2, true, Polarity.Positive, 1.0,
                    new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=5", 937.096923828125,
                    3, 3470064.5, 937.096923828125, 2, DissociationType.EThcD, 4, 937.096923828125));

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
