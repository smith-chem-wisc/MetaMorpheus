using Chemistry;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using MzLibUtil;
using ClassExtensions = Chemistry.ClassExtensions;
using Nett;
using TaskLayer;
using Omics.Modifications;
using Readers;
using Mzml = IO.MzML.Mzml;

namespace Test
{
    [TestFixture]
    internal static class Multiplex_Labeling_TMT_iTRAQ
    {
        [Test]
        [TestCase("C8 N1 H16", 126.128274520)]
        [TestCase("C8 H16 N{15}1", 127.125309415)]
        [TestCase("C7 H16 C{13}1 N{15}1", 128.128664250)]
        [TestCase("C7 N1 H16 C{13}1", 127.131629355)]
        [TestCase("C6 H16 C{13}2 N{15}1", 129.132019085)]
        [TestCase("C6 N1 H16 C{13}2", 128.134984190)]
        [TestCase("C5 N1 H16 C{13}3", 129.138339025)]
        [TestCase("C5 H16 C{13}3 N{15}1", 130.135373920)]
        [TestCase("C4 N1 H16 C{13}4", 130.141693860)]
        [TestCase("C4 H16 C{13}4 N{15}1", 131.138728755)]
        [TestCase("C3 N1 H16 C{13}5", 131.145048695)]
        [TestCase("C3 H16 C{13}5 N{15}1", 132.142083590)]
        [TestCase("C2 N1 H16 C{13}6", 132.148403531)]
        [TestCase("C2 H16 C{13}6 N{15}1", 133.145438425)]
        [TestCase("C1 N1 H16 C{13}7", 133.151758366)]
        [TestCase("C1 H16 C{13}7 N{15}1", 134.148793260)]
        [TestCase("N1 H16 C{13}8", 134.155113201)]
        [TestCase("H16 C{13}8 N{15}1", 135.152148095)]
        public static void TestChemicalFormulaWithIsotopesTMT(string formula, double mass)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            Assert.That(ClassExtensions.RoundedDouble(cf.MonoisotopicMass), Is.EqualTo(mass));
        }

        [Test]
        [TestCase("PEPTIDE", 1104.5743)]
        [TestCase("PEPTIDEK", 1536.8764)]
        public static void TestPeptideLabelledWithTMT18(string peptide, double totalMass)
        {
            List<Modification> fixedModifications = new List<Modification>();
            fixedModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> tmt18Mods = fixedModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("TMT18")).ToList();

            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, tmt18Mods, new List<Modification>()).First();
            var f = new List<Product>();
            p.Fragment(DissociationType.HCD, FragmentationTerminus.Both, f);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.That(ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4), Is.EqualTo(totalMass));
        }

        [Test]
        [TestCase("PEPTIDE", 1029.5302)]
        [TestCase("PEPTIDEK", 1386.7881)]
        public static void TestPeptideLabelledWithTMT(string peptide, double totalMass)
        {
            List<Modification> gptmdModifications = new List<Modification>();
            gptmdModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> tmt10Mods = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("TMT10")).ToList();
            List<Modification> tmtMods = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label").ToList();


            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, tmt10Mods, new List<Modification>()).First();
            var f = new List<Product>();
            p.Fragment(DissociationType.HCD, FragmentationTerminus.Both, f);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.That(ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4), Is.EqualTo(totalMass));
        }

        [Test]
        public static void TestAbilityToIDTMTDiagnosticIons()
        {
            var origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\tmt18test2.mzML");
            var scans = Mzml.LoadAllStaticData(origDataFile).GetAllScansList();
            var tolerance = new PpmTolerance(20);
            int roundTo = 5;
            var diagnosticIons = GlobalVariables.AllModsKnown
                .First(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("TMT18")).DiagnosticIons
                .First().Value.Select(p => Math.Round(p, roundTo)).ToList();

            // iterate through each scan
            foreach (var scan in scans)
            {
                if (!diagnosticIons.Any())
                    break;

                // find all possible matches to diagnostic ions
                var possibleMatches = scan.MassSpectrum.XArray
                    .Where(p => p <= diagnosticIons.Max() + tolerance.GetMaximumValue(diagnosticIons.Max()))
                    .Select(m => Math.Round(m, roundTo)).ToList();

                // check to see if the spectrum has each diagnostic ions
                for (int j = 0; j < diagnosticIons.Count; j++)
                {
                    // remove from diagnostic ion list if found within scan
                    if (possibleMatches.Any(p => tolerance.Within(p, diagnosticIons[j])))
                    {
                        diagnosticIons.Remove(diagnosticIons[j]);
                    }
                }
            }
            // will pass if all diagnostic ions are found in scans
            Assert.That(diagnosticIons.Count(), Is.EqualTo(0));
        }


        //[TestCase("PEPTIDE", 1104.5694)]
        //[TestCase("PEPTIDEK", 1536.8666)]
        public static void TestPeptideLabelledWiddth_iTRAQ_8plex(string peptide, double totalMass)
        {
            List<Modification> gptmdModifications = new List<Modification>();
            gptmdModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> itraq8plex = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("iTRAQ-8plex")).ToList();

            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, itraq8plex, new List<Modification>()).First();
            var f = new List<Product>();
            p.Fragment(DissociationType.HCD, FragmentationTerminus.Both, f);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.That(ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4), Is.EqualTo(totalMass));
        }

        [Test]
        public static void TestTmtIonsArentTreatedLikePeptideIsotopicEnvelopes()
        {
            // Previously, TMT reporter ions could sometime be mis-identified as isotopic envelopes
            // When MS2 Spectra are deconvoluted, if adjacent reporter ions have ratios that mimic an isotopic distribution,
            // Ex: 127C:128C 92:8, the 128C ion could be mis-identified as the M+1 peak of the 127C ion
            // and not reported as a distinct ion, resulting in loss of quantificative information.
            // 
            // In the new implementation, we ensure that TMT reporter ions are never treated as isotopic envelopes,
            // but instead are recovered directly from the raw MS2 spectra. 


            // Get the TMT 11 plex modifications (lysine and N-terminal)
            List<Modification> tmt11Mods = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("TMT11")).ToList();

            Protein protein = new Protein("PEPTIDEK", "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var peptide = protein.Digest(CommonParameters.DigestionParams, tmt11Mods, new List<Modification>()).First();
            var fragments = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);

            // We're going to create a fake MsDataScan that contains all fragment ions for PEPTIDEK and all TMT reporter ions
            // The TMT reporter ions will have intensities that mimic the isotopic distribution predicted by the Averagine Model
           
            List<double> allIonMzs = fragments.Select(m => m.NeutralMass.ToMz(1)).OrderBy(m => m).ToList();
            
            double[] mzArray = allIonMzs.ToArray();
            double[] intensityArray = new double[mzArray.Length];

            // Corresponding to                 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C
            var tmtIntensities = new double[] {  92,   92,    8,    12,  92,   92,    16,    20,    2,   3,   5 };

            for (int i = 0; i < intensityArray.Length; i++)
            {
                if(i < tmtIntensities.Length)
                {
                    intensityArray[i] = tmtIntensities[i];
                }
                else
                {
                    intensityArray[i] = 100; // arbitrary intensity for fragment ions
                }
            }

            // Create a MS1 scan
            var ms1Spectrum = new MzSpectrum(new double[] { peptide.MonoisotopicMass.ToMz(2), (peptide.MonoisotopicMass + Constants.C13MinusC12).ToMz(2) }, new double[] { 10000, 5000 }, false);
            var ms1Scan = new MsDataScan(ms1Spectrum, 1, 1, true, Polarity.Positive, 1.0, new MzRange(100, 2000), 
                scanFilter: "",
                MZAnalyzerType.Orbitrap, 
                totalIonCurrent: 1000,
                null,
                noiseData: null,
                nativeId: "scan=1",
                selectedIonMz: double.NaN,
                selectedIonChargeStateGuess: null,
                selectedIonIntensity: null,
                isolationMZ: double.NaN, null, null, null,
                selectedIonMonoisotopicGuessMz: null);

            // Create a MS2 scan
            var ms2Spectrum = new MzSpectrum(mzArray, intensityArray, false);
            var ms2Scan = new MsDataScan(ms2Spectrum, 2, 2, true, Polarity.Positive, 1.1, new MzRange(100, 2000),
                scanFilter: "", MZAnalyzerType.Orbitrap,
                totalIonCurrent: 1000, null, null, "scan=2", 
                selectedIonMz: peptide.MonoisotopicMass.ToMz(2), 
                selectedIonChargeStateGuess: 2,
                selectedIonIntensity: 10000, 
                isolationMZ: peptide.MonoisotopicMass.ToMz(2), 
                isolationWidth: 2.2,
                DissociationType.HCD, 1, 
                selectedIonMonoisotopicGuessMz: peptide.MonoisotopicMass.ToMz(2));

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestTmtOutput");
            var mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TMT_test", "SingleScan_TMT10.mzML");
            try
            {
                GenericMsDataFile singleScanDataFile = new GenericMsDataFile(new MsDataScan[] { ms1Scan, ms2Scan }, new SourceFile("no nativeID format", "mzML format", null, null, null));

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(singleScanDataFile, mzmlPath, true);

                // Now, run a task with this mzml and a simple fasta (contains only PEPTIDEK)
                var searchTask = Toml.ReadFile<SearchTask>(
                    Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\TMT-Task1-SearchTaskconfig.toml"),
                    MetaMorpheusTask.tomlConfig);

                List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) };
                string mzmlName = @"TMT_test\SingleScan_TMT10.mzML";
                string fastaName = @"TMT_test\PEPTIDEK.fasta";
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
                engine.Run();

                string[] peaksResults = File.ReadAllLines(Path.Combine(outputFolder, "search", "AllPeptides.psmtsv")).ToArray();


                string[] header = peaksResults[0].Trim().Split('\t');
                string[] ionLabelsInHeader = header[^11..]; // Last 11 columns should be the TMT labels
                Assert.That(ionLabelsInHeader, Is.EquivalentTo(new string[]
                    { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" }));

                var reportedIntensities = peaksResults[1].Split('\t')[^11..];
                for (int i = 0; i < reportedIntensities.Length; i++)
                {
                    Assert.That(reportedIntensities[i], Is.EqualTo(tmtIntensities[i].ToString("F1", CultureInfo.InvariantCulture)), $"The intensity for TMT channel {ionLabelsInHeader[i]} was not reported correctly. Expected {tmtIntensities[i]}, but was reported as {reportedIntensities[i]}.");
                }
                Assert.That(reportedIntensities.All(i => Double.Parse(i) > 0), Is.True, "All TMT channels should have intensities reported.");
            }
            finally
            {
                Directory.Delete(outputFolder, true);
                File.Delete(mzmlPath);
            }
        }

        [Test]
        public static void TestTmtQuantificationOutput()
        {
            var searchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\TMT-Task1-SearchTaskconfig.toml"),
                MetaMorpheusTask.tomlConfig);

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) };
            string mzmlName = @"TMT_test\VA084TQ_6.mzML";
            string fastaName =  @"TMT_test\mouseTMT.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestTmtOutput");
            if(Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            string[] peaksResults = File.ReadAllLines(Path.Combine(outputFolder, "search", "AllPeptides.psmtsv")).ToArray();
            Assert.That(peaksResults.Length == 5);

            string[] header = peaksResults[0].Trim().Split('\t');
            string[] ionLabelsInHeader = header[^11..]; // Last 11 columns should be the TMT labels
            Assert.That(ionLabelsInHeader, Is.EquivalentTo(new string[]
                { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" }));
            
            double channelSum127N = 0;
            for (int i = 1; i < peaksResults.Length; i++)
            {
                channelSum127N += Double.Parse(peaksResults[i].Trim().Split('\t')[^10]);
            }
            Assert.That(channelSum127N, Is.EqualTo(577226.34).Within(0.1));

            Directory.Delete(outputFolder, true);
        }


        [Test]
        public static void TestMs3TmtQuantificationWith()
        {
            var searchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\TMT11_LowCID_SearchTask.toml"),
                MetaMorpheusTask.tomlConfig);

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) };
            string mzmlName = @"TMT_test\MS3_TMT11_Mouse_snip.mzML";
            string fastaName = @"TMT_test\MUS_snip.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestTmtOutput");

            try
            {
                if(Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);

                var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
                engine.Run();

                var peptidePath = Path.Combine(outputFolder, "search", "AllPeptides.psmtsv");

                string[] peaksResults = File.ReadAllLines(peptidePath).ToArray();
                Assert.That(peaksResults.Length >= 5);

                string[] header = peaksResults[0].Trim().Split('\t');
                string[] ionLabelsInHeader = header[^11..]; // Last 11 columns should be the TMT labels
                Assert.That(ionLabelsInHeader, Is.EquivalentTo(new string[]
                    { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" }));

                var baseSeqIndex = Array.IndexOf(header, "Base Sequence");
                var scanNoIndex = Array.IndexOf(header, "Scan Number");

                var baseSeqsToCheck = new List<string>
                {
                    "TDTLCLNNTEISENGSDLSQK",
                    "LVQDVANNTNEEAGDGTTTATVLAR",
                    "TERPVNSAALSPNYDHVVLGGGQEAMDVTTTSTR",
                    "TDTLCLNNTEISENGSDLSQK",
                    "PMQFLGDEETVR",
                    "SLPGCQEIAEEFR"
                };

                bool checkedScan20 = false;
                foreach (var line in peaksResults)
                {
                    var columns = line.Trim().Split('\t');
                    if (baseSeqsToCheck.Contains(columns[baseSeqIndex]))
                    {
                        Assert.That(columns[^11..].All(i => Double.Parse(i) > 0), Is.True, $"Not all reporter ions have intensities > 0 for peptide {columns[baseSeqIndex]}");
                        if (columns[scanNoIndex] == "20") // Note, the MS3 scan is number 22
                        {
                            Assert.That(columns[^11..], Is.EquivalentTo(new string[]
                            {
                                "27274.3", "25929.9", "2745.8", "38435.7","4142.3",
                                "3555.4", "38433.3","31899.3", "34451.5", "27501.3", "25194.2"
                            }), $"The reporter ion intensities for scan 20 do not match the expected values.");
                        }
                    }
                }
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        [Test]
        public static void TestTMT10vs11()
        {
            var tmt10 = GlobalVariables.AllModsKnown
                .First(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("TMT10"));
            var tmt11 = GlobalVariables.AllModsKnown
                .First(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("TMT11"));

            Assert.That(tmt10.DiagnosticIons[DissociationType.HCD].Count, Is.EqualTo(10));
            Assert.That(tmt10.DiagnosticIons[DissociationType.HCD].Max(), Is.LessThan(tmt11.DiagnosticIons[DissociationType.HCD].Max()));
        }

        [Test]
        public static void TestDiLeuQuantificationOutput()
        {
            var searchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\DiLeu-12plex-Task1-SearchTaskconfig.toml"),
                MetaMorpheusTask.tomlConfig);

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) };
            string mzmlName = @"TMT_test\DiLeu_Slice_PXD029269.mzML";
            string fastaName = @"TMT_test\DiLeuSlice.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestDiLeuOutput");
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            string[] peaksResults = File.ReadAllLines(Path.Combine(outputFolder, "search", "AllPeptides.psmtsv")).ToArray();
            Assert.That(peaksResults.Length == 2);

            string[] header = peaksResults[0].Trim().Split('\t');
            string[] ionLabelsInHeader = header[^12..]; // Last 11 columns should be the TMT labels
            Assert.That(ionLabelsInHeader, Is.EquivalentTo(new string[]
                { "115a", "115b", "116a", "116b", "116c", "117a", "117b", "117c", "118a", "118b", "118c", "118d" }));

            double ionSum = peaksResults[1].Trim().Split('\t')[^12..].Select(s => double.Parse(s)).Sum();
            Assert.That(ionSum, Is.EqualTo(173357).Within(1));

            Directory.Delete(outputFolder, true);

            // Same data as above, but here we're pretending that DiLeu 4-plex was used
            searchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\DiLeu-4plex-Task1-SearchTaskconfig.toml"),
                MetaMorpheusTask.tomlConfig);
            taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) };
            engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            peaksResults = File.ReadAllLines(Path.Combine(outputFolder, "search", "AllPeptides.psmtsv")).ToArray();
            Assert.That(peaksResults.Length == 2);

            ionLabelsInHeader = peaksResults[0].Trim().Split('\t')[^4..];
            Assert.That(ionLabelsInHeader, Is.EquivalentTo(new string[] {"115", "116", "117", "118"}));

            ionSum = peaksResults[1].Trim().Split('\t')[^4..].Select(s => double.Parse(s)).Sum();
            Assert.That(ionSum, Is.EqualTo(115537).Within(1));

            Directory.Delete(outputFolder, true);
        }

        [Test]
        [TestCase("PEPTIDE", 944.4712)]
        [TestCase("PEPTIDEK", 1216.6702)]
        public static void TestPeptideLabelledWith_iTRAQ_4plex(string peptide, double totalMass)
        {
            List<Modification> gptmdModifications = new List<Modification>();
            gptmdModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> itraq4plex = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("iTRAQ-4plex")).ToList();

            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, itraq4plex, new List<Modification>()).First();
            var f = new List<Product>();
            p.Fragment(DissociationType.HCD, FragmentationTerminus.Both, f);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.That(ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4), Is.EqualTo(totalMass));
        }

        [Test]
        [TestCase("PEPTIDE", 945.4975)]
        [TestCase("PEPTIDEK", 1218.7227)]
        public static void TestPeptideLabelledWith_DiLeu_4plex(string peptide, double totalMass)
        {
            List<Modification> gptmdModifications = new List<Modification>();
            gptmdModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> itraq4plex = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("DiLeu-4plex")).ToList();

            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, itraq4plex, new List<Modification>()).First();
            var f = new List<Product>();
            p.Fragment(DissociationType.HCD, FragmentationTerminus.Both, f);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.That(ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4), Is.EqualTo(totalMass));
        }

        [Test]
        [TestCase("PEPTIDE", 1104.5694)]
        [TestCase("PEPTIDEK", 1536.8666)]
        public static void TestPeptideLabelledWith_iTRAQ_8plex(string peptide, double totalMass)
        {
            List<Modification> gptmdModifications = new List<Modification>();
            gptmdModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> itraq8plex = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("iTRAQ-8plex")).ToList();

            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, itraq8plex, new List<Modification>()).First();
            var f = new List<Product>();
            p.Fragment(DissociationType.HCD, FragmentationTerminus.Both, f);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.That(ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4), Is.EqualTo(totalMass));
        }

        [Test]
        [TestCase("C5 N2 H12 C{13}1", 114.110679698, true)]
        [TestCase("C5 C{13}1 N1 N{15}1 H12", 115.107714592, true)]
        [TestCase("C4 N1 H12 C{13}2 N{15}1", 116.111069427, true)]
        [TestCase("C3 N1 H12 C{13}3 N{15}1", 117.114424262, true)]
        [TestCase("C5 N2 H12 C{13}2 O{18}1", 144.105917679, false)]
        [TestCase("C4 N1O1 H12 C{13}3 N{15}1", 144.102062415, false)]
        [TestCase("C7 N3O3 H24 C{13}7 N{15}1", 304.205359390, false)]//for tmt 8-plex 113, 114, 116, 117
        [TestCase("C8 N2O3 H24 C{13}6 N{15}2", 304.199039449, false)]//for tmt 8-plex 115, 118, 119, 121
        //[TestCase("C4 N2 H12 C{13}2", 115.114034533, true)] this is old style (ABI or Thermo, don't know). no longer used
        public static void TestChemicalFormulaWithIsotopes_iTRAQ(string formula, double mass, bool mz)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            if (mz)
            {
                Assert.That(ClassExtensions.RoundedDouble(cf.MonoisotopicMass.ToMz(1)), Is.EqualTo(mass));
            }
            else
            {
                Assert.That(ClassExtensions.RoundedDouble(cf.MonoisotopicMass), Is.EqualTo(mass));
            }
        }

        [Test]
        [TestCase("C7 N{15}1 H15", 115.124760849, true)]
        [TestCase("C7 N1 H13 H{2}2", 116.140279447, true)]
        [TestCase("C7 N{15}1 H13 H{2}2", 117.137314341, true)]
        [TestCase("C7 N1 H11 H{2}4", 118.152832938, true)]
        [TestCase("C7 N{15}1 H15 C{13}1 O{18}1", 145.119998830, false)]
        [TestCase("C8 N1 H13 H{2}2 O{18}1", 145.132162593, false)]
        [TestCase("C7 N{15}1 H13 H{2}2 C{13}1 O1", 145.128307329, false)]
        [TestCase("C8 N1 H11 H{2}4 O1", 145.140471091, false)]
        public static void TestChemicalFormulaWithIsotopes_DiLeu4plex(string formula, double mass, bool mz)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            if (mz)
            {
                Assert.That(ClassExtensions.RoundedDouble(cf.MonoisotopicMass.ToMz(1)), Is.EqualTo(mass));
            }
            else
            {
                Assert.That(ClassExtensions.RoundedDouble(cf.MonoisotopicMass), Is.EqualTo(mass));
            }
        }
        [Test]
        [TestCase("C7 H15 N{15}1 ", 115.12476, true)]
        [TestCase("C6 C{13}1 H15 N{15}1 ", 116.12812, true)]
        [TestCase("C5 C{13}2 H15 N1 ", 116.13444, true)]
        [TestCase("C7 H13 H{2}2 N1 ", 116.14028, true)]
        [TestCase("C6 C{13}1 H15 N1 ", 115.13108, true)]
        [TestCase("C5 C{13}2 H15 N{15}1 ", 117.13147, true)]
        [TestCase("C7 H13 H{2}2 N{15}1 ", 117.13731, true)]
        [TestCase("C6 C{13}1 H13 H{2}2 N1 ", 117.14363, true)]
        [TestCase("C4 C{13}3 H15 N{15}1 ", 118.13483, true)]
        [TestCase("C6 C{13}1 H13 H{2}2 N{15}1 ", 118.14067, true)]
        [TestCase("C5 C{13}2 H13 H{2}2 N1 ", 118.14699, true)]
        [TestCase("C7 H11 H{2}4 N1 ", 118.15283, true)]

        [TestCase("C7 C{13}1 H15 N{15}1 O{18}1 ", 145.119998830, false)]
        [TestCase("C6 C{13}2 H15 N1 O{18}1 ", 145.126318771, false)]
        [TestCase("C7 C{13}1 H15 N{15}1 O{18}1 ", 145.119998830, false)]
        [TestCase("C6 C{13}2 H15 N1 O{18}1 ", 145.126318771, false)]
        [TestCase("C8 H13 H{2}2 N1 O{18}1 ", 145.132162593, false)]
        [TestCase("C5 C{13}3 H15 N{15}1 O1 ", 145.122463507, false)]
        [TestCase("C7 C{13}1 H13 H{2}2 N{15}1 O1 ", 145.128307329, false)]
        [TestCase("C6 C{13}2 H13 H{2}2 N1 O1 ", 145.134627269, false)]
        [TestCase("C5 C{13}3 H15 N{15}1 O1 ", 145.122463507, false)]
        [TestCase("C7 C{13}1 H13 H{2}2 N{15}1 O1 ", 145.128307329, false)]
        [TestCase("C6 C{13}2 H13 H{2}2 N1 O1 ", 145.134627269, false)]
        [TestCase("C8 H11 H{2}4 N1 O1 ", 145.140471091, false)]

        public static void TestChemicalFormulaWithIsotopes_DiLeu12plex(string formula, double mass, bool mz)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            if (mz)
            {
                Assert.That(ClassExtensions.RoundedDouble(cf.MonoisotopicMass.ToMz(1), 5), Is.EqualTo(mass));
            }
            else
            {
                Assert.That(ClassExtensions.RoundedDouble(cf.MonoisotopicMass), Is.EqualTo(mass));
            }
        }
        /// <summary>
        /// We choose not to count diagnostic ions in the score directly. But we still want to detect and annotate them. This test shows that we identify
        /// diagnostic ions but do not include them in the score.
        /// 
        /// NOTE: You can detect TMT diagnostic ions in LowCID because it lacks the resolution to isolate them.
        /// 
        /// </summary>
        [Test]
        public static void TestDoNotCountDiagnosticIonsInScore_HCD()
        {
            var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_Test\TMT-Task1-SearchTaskconfig.toml");
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_Test\TestOutput");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_Test\VA084TQ_6.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_Test\mouseTmt.fasta");

            var engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("SearchTOML", searchTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();

            string psmFile = Path.Combine(outputFolder, @"SearchTOML\AllPSMs.psmtsv");

            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            PsmFromTsv psm = parsedPsms.First();

            Assert.That(psm.MatchedIons.Count, Is.EqualTo(38)); //matched ions include b, y and D (diagnostic ions in TMT search)
            Assert.That(psm.MatchedIons.Where(i => i.NeutralTheoreticalProduct.ProductType == ProductType.D).Count(), Is.EqualTo(8)); //There are 8 discovered diagnostic ions
            Assert.That(psm.MatchedIons.Where(i => i.NeutralTheoreticalProduct.ProductType != ProductType.D).Count(), Is.EqualTo(30)); //There are 30 b and y ions (excluding diagnostic ions)
            Assert.That(psm.Score, Is.EqualTo(30.273)); //score should only use non-diagnostic ions (start with 30 in this case).
            Assert.That((int)psm.Score, Is.EqualTo(psm.MatchedIons.Where(i => i.NeutralTheoreticalProduct.ProductType != ProductType.D).Count())); //integer part of the MM score should match the count of non-diagnostic ions

            Assert.That(psm.BaseSeq, Is.EqualTo("VFNTTPDDLDLHVIYDVSHNIAK"));
            Assert.That(psm.DecoyContamTarget, Is.EqualTo("T"));
            Assert.That(psm.FullSequence, Is.EqualTo("[Multiplex Label:TMT11 on X]VFNTTPDDLDLHVIYDVSHNIAK[Multiplex Label:TMT11 on K]"));
            Assert.That(psm.OrganismName, Is.EqualTo("Mus musculus"));
            Assert.That(psm.PeptideDescription, Is.EqualTo("full"));
            Assert.That(psm.ProteinAccession, Is.EqualTo("Q99LF4"));

            Directory.Delete(outputFolder,true);
        }
        [Test]
        public static void TestDoNotCountDiagnosticIonsInScore_LowCID()
        {
            TestDataFile testDataFile = new();
            //Tolerance productMassTolerance = new AbsoluteTolerance(0.01);
            double precursorMass = 300;
            //The below theoretical does not accurately represent B-Y ions
            double[] sorted_theoretical_product_masses_for_this_peptide = new double[] { precursorMass + (2 * Constants.ProtonMass) - 275.1350, precursorMass + (2 * Constants.ProtonMass) - 258.127, precursorMass + (2 * Constants.ProtonMass) - 257.1244, 50, 60, 70, 147.0764, precursorMass + (2 * Constants.ProtonMass) - 147.0764, precursorMass + (2 * Constants.ProtonMass) - 70, precursorMass + (2 * Constants.ProtonMass) - 60, precursorMass + (2 * Constants.ProtonMass) - 50, 257.1244, 258.127, 275.1350 }; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
            List<Product> productsWithLocalizedMassDiff = new();
            
            //add one diagnostic ion
            productsWithLocalizedMassDiff.Add(new Product(ProductType.D, FragmentationTerminus.Both, sorted_theoretical_product_masses_for_this_peptide[11], 1, 1, 0));

            for (int i = 0; i < sorted_theoretical_product_masses_for_this_peptide.Length; i++)
            {
                if(i != 11)
                {
                    productsWithLocalizedMassDiff.Add(new Product(ProductType.b, FragmentationTerminus.Both, sorted_theoretical_product_masses_for_this_peptide[i], 1, 1, 0));
                }
            }

            //ensure there is only one diagnostic ion
            Assert.That(productsWithLocalizedMassDiff.Where(p => p.ProductType == ProductType.D).Count(), Is.EqualTo(1));

            //Check total ion count
            Assert.That(productsWithLocalizedMassDiff.Count, Is.EqualTo(14));

            MsDataScan scan = testDataFile.GetOneBasedScan(2);
            scan.MassSpectrum.XCorrPrePreprocessing(1.0, 500.0, 300.0);

            //check that the scan is noted as xcorr processed
            Assert.That(scan.MassSpectrum.XcorrProcessed);

            Tolerance tolerance = new AbsoluteTolerance(1.0);
            CommonParameters commonParams = new(productMassTolerance: tolerance);
            var scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(1), 1, "", commonParams);
            List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParams);

            // score when the mass-diff is on this residue
            double score = MetaMorpheusEngine.CalculatePeptideScore(scan, matchedIons);

            Assert.That((int)score, Is.EqualTo(0));
        }

        [Test]
        public static void TestGetIsobaricMassTagWithModificationId()
        {
            // Test valid modification IDs
            var tmt10Tag = IsobaricMassTag.GetIsobaricMassTag("TMT10 on K");
            Assert.That(tmt10Tag, Is.Not.Null);
            Assert.That(tmt10Tag.TagType, Is.EqualTo(IsobaricMassTagType.TMT10));
            Assert.That(tmt10Tag.ReporterIonMzs.Length, Is.EqualTo(10));
            Assert.That(tmt10Tag.ReporterIonMzRanges, Is.Not.Null);
            Assert.That(tmt10Tag.ReporterIonMzRanges, Is.TypeOf<DoubleRange[]>());

            var tmt11Tag = IsobaricMassTag.GetIsobaricMassTag("TMT11 on X");
            Assert.That(tmt11Tag, Is.Not.Null);
            Assert.That(tmt11Tag.TagType, Is.EqualTo(IsobaricMassTagType.TMT11));
            Assert.That(tmt11Tag.ReporterIonMzs.Length, Is.EqualTo(11));

            // Test null/invalid input
            var nullTag = IsobaricMassTag.GetIsobaricMassTag((string)null);
            Assert.That(nullTag, Is.Null);

            var invalidTag = IsobaricMassTag.GetIsobaricMassTag("InvalidModification");
            Assert.That(invalidTag, Is.Null);
        }

        [Test]
        public static void TestGetIsobaricMassTagWithTagType()
        {
            // Test all valid tag types
            var tmt6 = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT6);
            Assert.That(tmt6, Is.Not.Null);
            Assert.That(tmt6.ReporterIonMzs.Length, Is.EqualTo(6));
            Assert.That(tmt6.TagType, Is.EqualTo(IsobaricMassTagType.TMT6));

            var tmt10 = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);
            Assert.That(tmt10, Is.Not.Null);
            Assert.That(tmt10.ReporterIonMzs.Length, Is.EqualTo(10));

            var tmt11 = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT11);
            Assert.That(tmt11, Is.Not.Null);
            Assert.That(tmt11.ReporterIonMzs.Length, Is.EqualTo(11));

            var tmt18 = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT18);
            Assert.That(tmt18, Is.Not.Null);
            Assert.That(tmt18.ReporterIonMzs.Length, Is.EqualTo(18));

            var itraq4 = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.iTRAQ4);
            Assert.That(itraq4, Is.Not.Null);
            Assert.That(itraq4.ReporterIonMzs.Length, Is.EqualTo(4));

            var itraq8 = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.iTRAQ8);
            Assert.That(itraq8, Is.Not.Null);
            Assert.That(itraq8.ReporterIonMzs.Length, Is.EqualTo(8));

            var dileu4 = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.diLeu4);
            Assert.That(dileu4, Is.Not.Null);
            Assert.That(dileu4.ReporterIonMzs.Length, Is.EqualTo(4));

            var dileu12 = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.diLeu12);
            Assert.That(dileu12, Is.Not.Null);
            Assert.That(dileu12.ReporterIonMzs.Length, Is.EqualTo(12));

            // Test null input
            var nullTag = IsobaricMassTag.GetIsobaricMassTag((IsobaricMassTagType?)null);
            Assert.That(nullTag, Is.Null);
        }

        [Test]
        public static void TestGetReporterIonIntensitiesWithMatchingPeaks()
        {
            var tmt10Tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);

            // Create a mock spectrum with peaks at reporter ion m/z values
            double[] mzArray = tmt10Tag.ReporterIonMzs.ToArray();
            double[] intensityArray = new double[] { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };

            var spectrum = new MzSpectrum(mzArray, intensityArray, false);

            var intensities = tmt10Tag.GetReporterIonIntensities(spectrum);

            Assert.That(intensities, Is.Not.Null);
            Assert.That(intensities.Length, Is.EqualTo(10));
            Assert.That(intensities, Is.EqualTo(intensityArray));
        }

        [Test]
        public static void TestGetReporterIonIntensitiesWithMissingPeaks()
        {
            var tmt10Tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);

            // Create a spectrum with only some reporter ions
            double[] mzArray = new double[] { tmt10Tag.ReporterIonMzs[0], tmt10Tag.ReporterIonMzs[2], tmt10Tag.ReporterIonMzs[5] };
            double[] intensityArray = new double[] { 100, 300, 600 };

            var spectrum = new MzSpectrum(mzArray, intensityArray, false);

            var intensities = tmt10Tag.GetReporterIonIntensities(spectrum);

            Assert.That(intensities, Is.Not.Null);
            Assert.That(intensities.Length, Is.EqualTo(10));
            Assert.That(intensities[0], Is.EqualTo(100));
            Assert.That(intensities[1], Is.EqualTo(0)); // Missing
            Assert.That(intensities[2], Is.EqualTo(300));
            Assert.That(intensities[3], Is.EqualTo(0)); // Missing
            Assert.That(intensities[4], Is.EqualTo(0)); // Missing
            Assert.That(intensities[5], Is.EqualTo(600));
        }

        [Test]
        public static void TestGetReporterIonIntensitiesWithEmptySpectrum()
        {
            var tmt10Tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);

            var emptySpectrum = new MzSpectrum(new double[0], new double[0], false);

            var intensities = tmt10Tag.GetReporterIonIntensities(emptySpectrum);

            Assert.That(intensities, Is.Null);
        }

        [Test]
        public static void TestGetReporterIonIntensitiesWithSpectrumStartingAfterReporterIons()
        {
            var tmt10Tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);

            // Spectrum starts at m/z 200, well after reporter ions
            double[] mzArray = new double[] { 200, 300, 400, 500 };
            double[] intensityArray = new double[] { 100, 200, 300, 400 };

            var spectrum = new MzSpectrum(mzArray, intensityArray, false);

            var intensities = tmt10Tag.GetReporterIonIntensities(spectrum);

            Assert.That(intensities, Is.Not.Null);
            Assert.That(intensities.Length, Is.EqualTo(10));
            Assert.That(intensities.All(i => i == 0), Is.True);
        }

        [Test]
        public static void TestGetReporterIonIntensitiesWithTolerance()
        {
            var tmt10Tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);

            // Create peaks slightly offset from theoretical values (within tolerance)
            double[] mzArray = tmt10Tag.ReporterIonMzs.Select(mz => mz + 0.002).ToArray();
            double[] intensityArray = new double[] { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };

            var spectrum = new MzSpectrum(mzArray, intensityArray, false);

            var intensities = tmt10Tag.GetReporterIonIntensities(spectrum);

            Assert.That(intensities, Is.Not.Null);
            Assert.That(intensities, Is.EqualTo(intensityArray));
        }

        [Test]
        public static void TestGetReporterIonIntensitiesOutsideTolerance()
        {
            var tmt10Tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);

            // Create peaks outside tolerance (>0.003 Da)
            double[] mzArray = tmt10Tag.ReporterIonMzs.Select(mz => mz + 0.00301).ToArray();
            double[] intensityArray = new double[] { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };

            var spectrum = new MzSpectrum(mzArray, intensityArray, false);

            var intensities = tmt10Tag.GetReporterIonIntensities(spectrum);

            Assert.That(intensities, Is.Not.Null);
            Assert.That(intensities.All(i => i == 0), Is.True);
        }

        [Test]
        public static void TestGetReporterIonLabelsWithModificationId()
        {
            var tmt10Labels = IsobaricMassTag.GetReporterIonLabels("TMT10 on K");
            Assert.That(tmt10Labels, Is.Not.Null);
            Assert.That(tmt10Labels.Count, Is.EqualTo(10));
            Assert.That(tmt10Labels, Is.EqualTo(new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N" }));

            var tmt11Labels = IsobaricMassTag.GetReporterIonLabels("TMT11 on X");
            Assert.That(tmt11Labels.Count, Is.EqualTo(11));
            Assert.That(tmt11Labels[10], Is.EqualTo("131C"));

            var nullLabels = IsobaricMassTag.GetReporterIonLabels((string)null);
            Assert.That(nullLabels, Is.Null);

            var invalidLabels = IsobaricMassTag.GetReporterIonLabels("InvalidMod");
            Assert.That(invalidLabels, Is.Null);
        }

        [Test]
        public static void TestGetReporterIonLabelsWithTagType()
        {
            var tmt6Labels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.TMT6);
            Assert.That(tmt6Labels, Is.EqualTo(new List<string> { "126", "127", "128", "129", "130", "131" }));

            var tmt10Labels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.TMT10);
            Assert.That(tmt10Labels, Is.EqualTo(new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N" }));

            var tmt11Labels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.TMT11);
            Assert.That(tmt11Labels, Is.EqualTo(new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" }));

            var tmt18Labels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.TMT18);
            Assert.That(tmt18Labels.Count, Is.EqualTo(18));

            var itraq4Labels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.iTRAQ4);
            Assert.That(itraq4Labels, Is.EqualTo(new List<string> { "114", "115", "116", "117" }));

            var itraq8Labels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.iTRAQ8);
            Assert.That(itraq8Labels, Is.EqualTo(new List<string> { "113", "114", "115", "116", "117", "118", "119", "120" }));

            var dileu4Labels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.diLeu4);
            Assert.That(dileu4Labels, Is.EqualTo(new List<string> { "115", "116", "117", "118" }));

            var dileu12Labels = IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType.diLeu12);
            Assert.That(dileu12Labels, Is.EqualTo(new List<string> { "115a", "115b", "116a", "116b", "116c", "117a", "117b", "117c", "118a", "118b", "118c", "118d" }));
        }

        [Test]
        public static void TestGetTagTypeFromModificationId()
        {
            // Test TMT variants
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("TMT6 on K"), Is.EqualTo(IsobaricMassTagType.TMT6));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("TMT10 on K"), Is.EqualTo(IsobaricMassTagType.TMT10));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("TMT11 on X"), Is.EqualTo(IsobaricMassTagType.TMT11));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("TMT18 on X"), Is.EqualTo(IsobaricMassTagType.TMT18));

            // Test iTRAQ variants
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("iTRAQ-4plex on K"), Is.EqualTo(IsobaricMassTagType.iTRAQ4));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("iTRAQ-8plex on K"), Is.EqualTo(IsobaricMassTagType.iTRAQ8));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("iTRAQ4 on K"), Is.EqualTo(IsobaricMassTagType.iTRAQ4));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("iTRAQ8 on K"), Is.EqualTo(IsobaricMassTagType.iTRAQ8));

            // Test DiLeu variants
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("DiLeu-4plex on K"), Is.EqualTo(IsobaricMassTagType.diLeu4));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("DiLeu-12plex on X"), Is.EqualTo(IsobaricMassTagType.diLeu12));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("DiLeu4 on K"), Is.EqualTo(IsobaricMassTagType.diLeu4));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("DiLeu12 on X"), Is.EqualTo(IsobaricMassTagType.diLeu12));

            // Test case insensitivity
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("tmt10 on k"), Is.EqualTo(IsobaricMassTagType.TMT10));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("ITRAQ-4PLEX ON K"), Is.EqualTo(IsobaricMassTagType.iTRAQ4));
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("dileu-12plex on x"), Is.EqualTo(IsobaricMassTagType.diLeu12));

            // Test null and invalid inputs
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId(null), Is.Null);
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId(""), Is.Null);
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("   "), Is.Null);
            Assert.That(IsobaricMassTag.GetTagTypeFromModificationId("InvalidModification"), Is.Null);
        }

        [Test]
        public static void TestAbsoluteToleranceValue()
        {
            Assert.That(IsobaricMassTag.AbsoluteToleranceValue, Is.EqualTo(0.003));
        }

        [Test]
        public static void TestReporterIonMzsAreOrdered()
        {
            var tmt10Tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);

            for (int i = 0; i < tmt10Tag.ReporterIonMzs.Length - 1; i++)
            {
                Assert.That(tmt10Tag.ReporterIonMzs[i], Is.LessThan(tmt10Tag.ReporterIonMzs[i + 1]));
            }
        }

        [Test]
        public static void TestGetReporterIonIntensitiesWithComplexSpectrum()
        {
            var tmt10Tag = IsobaricMassTag.GetIsobaricMassTag(IsobaricMassTagType.TMT10);

            // Create a complex spectrum with reporter ions mixed with other peaks
            var mzList = new List<double>();
            var intensityList = new List<double>();

            // Add some noise peaks before reporter ions
            mzList.Add(tmt10Tag.ReporterIonMzs[0] - 0.002);
            intensityList.Add(50);

            // Add reporter ions with noise peaks between them
            for (int i = 0; i < tmt10Tag.ReporterIonMzs.Length; i++)
            {
                mzList.Add(tmt10Tag.ReporterIonMzs[i]);
                intensityList.Add((i + 1) * 100);

                // Add noise peak after each reporter ion
                mzList.Add(tmt10Tag.ReporterIonMzs[i] + 0.001);
                intensityList.Add(25);
            }

            // Add some noise peaks after reporter ions
            mzList.Add(tmt10Tag.ReporterIonMzs[^1] + 0.002);
            intensityList.Add(75);

            var spectrum = new MzSpectrum(mzList.ToArray(), intensityList.ToArray(), false);

            var intensities = tmt10Tag.GetReporterIonIntensities(spectrum);

            Assert.That(intensities, Is.Not.Null);
            Assert.That(intensities.Length, Is.EqualTo(10));
            for (int i = 0; i < 10; i++)
            {
                Assert.That(intensities[i], Is.EqualTo((i + 1) * 100));
            }
        }

        [Test]
        public static void TestSearchTaskExceptionOnNullMassTag()
        {
            // This test simulates what happens in SearchTask when IsobaricMassTag.GetIsobaricMassTag returns null
            // The actual SearchTask code throws MetaMorpheusException in this case
            
            string invalidModId = "InvalidModification";
            var massTag = IsobaricMassTag.GetIsobaricMassTag(invalidModId);

            // When massTag is null, SearchTask should throw MetaMorpheusException
            if (massTag == null)
            {
                Assert.Throws<MetaMorpheusException>(() => throw new MetaMorpheusException("Could not find isobaric mass tag with the name " + invalidModId));
            }
        }
    }
}