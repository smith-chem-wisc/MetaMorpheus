using Chemistry;
using EngineLayer;
using IO.MzML;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using MzLibUtil;
using ClassExtensions = Chemistry.ClassExtensions;
using Nett;
using TaskLayer;
using UsefulProteomicsDatabases;
using System.Threading.Tasks;

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
            Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass));
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
           
            Assert.AreEqual(totalMass, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4));
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

            Assert.AreEqual(totalMass, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4));
        }

        [Test]
        public static void TestMultiplexIonIntensityDetection()
        {
            double[] ionMzs = { 1, 2, 3, 4, 5 };
            Tolerance tol = new PpmTolerance(10);

            // 1-to-1 concordance between theoretical and observed diagnostic ions
            MzSpectrum fakeSpectrum = new MzSpectrum(
                mz: new double[] { 1, 2, 3, 4, 5 },
                intensities: new double[] { 2, 4, 6, 8, 10 },
                shouldCopy: false);
            Assert.AreEqual(
                PostSearchAnalysisTask.GetMultiplexIonIntensities(fakeSpectrum, ionMzs, tol),
                fakeSpectrum.YArray);

            // Every other diagnostic ion is present in spectrum
            fakeSpectrum = new MzSpectrum(
                mz: new double[] { 1, 2.5, 3, 4.5, 5 },
                intensities: new double[] { 2, 4, 6, 8, 10 },
                shouldCopy: false);
            Assert.AreEqual(
                PostSearchAnalysisTask.GetMultiplexIonIntensities(fakeSpectrum, ionMzs, tol),
                new double[] { 2, 0, 6, 0, 10 });

            // Last two diagnostic ions (highest m/z) are not observed
            fakeSpectrum = new MzSpectrum(
                mz: new double[] { 1, 2, 3 },
                intensities: new double[] { 2, 4, 6 },
                shouldCopy: false);
            Assert.AreEqual(
                PostSearchAnalysisTask.GetMultiplexIonIntensities(fakeSpectrum, ionMzs, tol),
                new double[] { 2, 4, 6, 0, 0 });
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
            Assert.AreEqual(0, diagnosticIons.Count());
        }

        [Test]
        public static void TestTmtQuantificationOutput()
        {
            var searchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\TMT-Task1-SearchTaskconfig.toml"),
                MetaMorpheusTask.tomlConfig);
            searchTask.SearchParameters.DoLabelFreeQuantification = false;

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) };
            string mzmlName = @"TMT_test\VA084TQ_6.mzML";
            string fastaName =  @"TMT_test\mouseTMT.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestTmtOutput");
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
            Assert.That(channelSum127N, Is.EqualTo(577226.336).Within(0.001));
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

            Assert.AreEqual(totalMass, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4));
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

            Assert.AreEqual(totalMass, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4));
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

            Assert.AreEqual(totalMass, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4));
        }

        [Test]
        [TestCase("C5 N2 H12 C{13}1", 114.110679698, true)]
        //[TestCase("C4 N2 H12 C{13}2", 115.114034533, true)] this is old style (ABI or Thermo, don't know). no longer used
        [TestCase("C5 C{13}1 N1 N{15}1 H12", 115.107714592, true)]
        [TestCase("C4 N1 H12 C{13}2 N{15}1", 116.111069427, true)]
        [TestCase("C3 N1 H12 C{13}3 N{15}1", 117.114424262, true)]
        [TestCase("C5 N2 H12 C{13}2 O{18}1", 144.105917679, false)]
        [TestCase("C4 N1O1 H12 C{13}3 N{15}1", 144.102062415, false)]
        [TestCase("C7 N3O3 H24 C{13}7 N{15}1", 304.205359390, false)]//for tmt 8-plex 113, 114, 116, 117
        [TestCase("C8 N2O3 H24 C{13}6 N{15}2", 304.199039449, false)]//for tmt 8-plex 115, 118, 119, 121
        public static void TestChemicalFormulaWithIsotopes_iTRAQ(string formula, double mass, bool mz)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            if (mz)
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass.ToMz(1)));
            }
            else
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass));
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
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass.ToMz(1)));
            }
            else
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass));
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
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass.ToMz(1), 5));
            }
            else
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass));
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

            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            PsmFromTsv psm = parsedPsms.First();
            

            Assert.AreEqual(38, psm.MatchedIons.Count); //matched ions include b, y and D (diagnostic ions in TMT search)
            Assert.AreEqual(8,psm.MatchedIons.Where(i=>i.NeutralTheoreticalProduct.ProductType == ProductType.D).Count()); //There are 8 discovered diagnostic ions
            Assert.AreEqual(30, psm.MatchedIons.Where(i => i.NeutralTheoreticalProduct.ProductType != ProductType.D).Count()); //There are 30 b and y ions (excluding diagnostic ions)
            Assert.AreEqual(30.273, psm.Score); //score should only use non-diagnostic ions (start with 30 in this case).
            Assert.AreEqual((int)psm.Score, psm.MatchedIons.Where(i => i.NeutralTheoreticalProduct.ProductType != ProductType.D).Count()); //integer part of the MM score should match the count of non-diagnostic ions

            Assert.AreEqual("VFNTTPDDLDLHVIYDVSHNIAK", psm.BaseSeq);
            Assert.AreEqual("T", psm.DecoyContamTarget);
            Assert.AreEqual("[Multiplex Label:TMT11 on X]VFNTTPDDLDLHVIYDVSHNIAK[Multiplex Label:TMT11 on K]", psm.FullSequence);
            Assert.AreEqual("Mus musculus", psm.OrganismName);
            Assert.AreEqual("full", psm.PeptideDescription);
            Assert.AreEqual("Q99LF4", psm.ProteinAccession);

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
            Assert.AreEqual(1, productsWithLocalizedMassDiff.Where(p => p.ProductType == ProductType.D).Count());

            //Check total ion count
            Assert.AreEqual(14, productsWithLocalizedMassDiff.Count);

            MsDataScan scan = testDataFile.GetOneBasedScan(2);
            scan.MassSpectrum.XCorrPrePreprocessing(1.0,500.0,300.0);

            //check that the scan is noted as xcorr processed
            Assert.IsTrue(scan.MassSpectrum.XcorrProcessed);

            Tolerance tolerance = new AbsoluteTolerance(1.0);
            CommonParameters commonParams = new(productMassTolerance: tolerance);
            var scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(1), 1, "", commonParams);
            List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParams);

            // score when the mass-diff is on this residue
            double score = MetaMorpheusEngine.CalculatePeptideScore(scan, matchedIons);

            Assert.AreEqual(0, (int)score);
        }
    }
}