using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using MzLibUtil;

namespace Test
{
    [TestFixture]
    public class BinGenerationTest
    {
        #region Public Methods

        [Test]
        public static void TestBinGeneration()
        {
            List<MassDiffAcceptor> massDiffAcceptors = new List<MassDiffAcceptor> { new OpenSearchMode() };
            SearchTask st = new SearchTask()
            {
                CommonParameters = new CommonParameters
                {
                    ScoreCutoff = 1,
                    Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"],
                    InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,                            
                    MinPeptideLength = 5,
                    ConserveMemory = false,
                    MaxMissedCleavages = 2,
                    MaxPeptideLength = null,
                    MaxModificationIsoforms = 4096,
                    BIons = true,
                    YIons = true, 
                    TotalPartitions = 1,
                    LocalizeAll = true,
                    MaxDegreeOfParallelism = 1,
                    ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") },
                    ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") },
                    ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList(),
                    ProductMassTolerance = new AbsoluteTolerance(0.01),
                    ZdotIons = false,
                    CIons = false,

                    Max_mods_for_peptide = 3,

                    // Deconvolution stuff
                    DoPrecursorDeconvolution = true,
                    UseProvidedPrecursorInfo = true,
                    DeconvolutionIntensityRatio = 4,
                    DeconvolutionMaxAssumedChargeState = 10,
                    DeconvolutionMassTolerance = new PpmTolerance(5),
                },
                SearchParameters = new SearchParameters
                {
                    DoHistogramAnalysis = true,
                    MassDiffAcceptors = massDiffAcceptors,                   
                    SearchTarget = true,
                    SearchDecoy = false,
                    DoParsimony = true,
                    DoQuantification = true,
                    DisposeOfFileWhenDone = true,
                    AddCompIons = false,
                    NoOneHitWonders = false,
                    ModPeptidesAreUnique = true,
                    QuantifyPpmTol = 5,             
                    DoLocalizationAnalysis = true,
                    WritePrunedDatabase = false,
                    KeepAllUniprotMods = true,

                },
        };

            string proteinDbFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.xml");
            string mzmlFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.mzml");

            Protein prot1 = new Protein("MEDEEK", "prot1");
            Protein prot2 = new Protein("MENEEK", "prot2");

            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            ModificationWithMass mod = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, 10);

            var possMod1 = prot1.Digest(st.CommonParameters.Protease, st.CommonParameters.MaxMissedCleavages, st.CommonParameters.MinPeptideLength, st.CommonParameters.MaxPeptideLength, st.CommonParameters.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep1_0 = possMod1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            var pep1_10 = possMod1.GetPeptidesWithSetModifications(new List<ModificationWithMass> { mod }, 4096, 3).Last();

            Protein prot3 = new Protein("MAAADAAAAAAAAAAAAAAA", "prot3");

            var possMod3 = prot3.Digest(st.CommonParameters.Protease, st.CommonParameters.MaxMissedCleavages, st.CommonParameters.MinPeptideLength, st.CommonParameters.MaxPeptideLength, st.CommonParameters.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep2_0 = possMod3.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            var pep2_10 = possMod3.GetPeptidesWithSetModifications(new List<ModificationWithMass>() { mod }, 4096, 3).Last();

            Protein prot4 = new Protein("MNNDNNNN", "prot4");
            var possMod4 = prot4.Digest(st.CommonParameters.Protease, st.CommonParameters.MaxMissedCleavages, st.CommonParameters.MinPeptideLength, st.CommonParameters.MaxPeptideLength, st.CommonParameters.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep3_10 = possMod4.GetPeptidesWithSetModifications(new List<ModificationWithMass>() { mod }, 4096, 3).Last();

            List<PeptideWithSetModifications> pepsWithSetMods = new List<PeptideWithSetModifications> { pep1_0, pep1_10, pep2_0, pep2_10, pep3_10 };
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepsWithSetMods);

            List<Protein> proteinList = new List<Protein> { prot1, prot2, prot3, prot4 };

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlFilePath, false);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList, proteinDbFilePath);

            string output_folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestBinGeneration");
            Directory.CreateDirectory(output_folder);
            st.RunTask(
                output_folder,
                new List<DbForTask> { new DbForTask(proteinDbFilePath, false) },
                new List<string> { mzmlFilePath },
                null);

            Assert.AreEqual(3, File.ReadLines(Path.Combine(output_folder, @"aggregate_OpenSearch.mytsv")).Count());
        }

        [Test]
        public static void TestProteinSplitAcrossFiles()
        {
            List<MassDiffAcceptor> massDiffAcceptors = new List<MassDiffAcceptor> { new OpenSearchMode() };
            SearchTask st = new SearchTask()
            {
                CommonParameters = new CommonParameters
                {
                    ScoreCutoff = 1,
                    Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"],
                    InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                    MinPeptideLength = 5,
                    ConserveMemory = false,
                    MaxMissedCleavages = 0,
                    MaxPeptideLength = null,
                    MaxModificationIsoforms = 4096,
                    BIons = true,
                    YIons = true,
                    TotalPartitions = 1,
                    LocalizeAll = true,
                    MaxDegreeOfParallelism = 1,
                    ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") },
                    ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") },
                    ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList(),
                    ProductMassTolerance = new AbsoluteTolerance(0.01),
                    ZdotIons = false,
                    CIons = false,

                    Max_mods_for_peptide = 3,

                    // Deconvolution stuff
                    DoPrecursorDeconvolution = true,
                    UseProvidedPrecursorInfo = true,
                    DeconvolutionIntensityRatio = 4,
                    DeconvolutionMaxAssumedChargeState = 10,
                    DeconvolutionMassTolerance = new PpmTolerance(5),
                },
                SearchParameters = new SearchParameters
                {
                    DoHistogramAnalysis = true,
                    MassDiffAcceptors = massDiffAcceptors,
                    SearchTarget = true,
                    SearchDecoy = false,
                    DoParsimony = true,
                    DoQuantification = true,
                    DisposeOfFileWhenDone = true,
                    AddCompIons = false,
                    NoOneHitWonders = false,
                    ModPeptidesAreUnique = true,
                    QuantifyPpmTol = 5,
                    DoLocalizationAnalysis = true,
                    WritePrunedDatabase = false,
                    KeepAllUniprotMods = true,

                },
            };

            string proteinDbFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles.xml");
            string mzmlFilePath1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles1.mzml");
            string mzmlFilePath2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles2.mzml");

            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            ModificationWithMass mod = new ModificationWithMass("mod1", "mt", motif, TerminusLocalization.Any, 10);

            IDictionary<int, List<Modification>> oneBasedModification = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification>{ mod } }
            };

            Protein prot1 = new Protein("MEDEEK", "prot1", oneBasedModifications: oneBasedModification);

            var possMod1 = prot1.Digest(st.CommonParameters.Protease, st.CommonParameters.MaxMissedCleavages, st.CommonParameters.MinPeptideLength, st.CommonParameters.MaxPeptideLength, st.CommonParameters.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep1 = possMod1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            var pep2 = possMod1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).Last();

            List<PeptideWithSetModifications> listForFile1 = new List<PeptideWithSetModifications> { pep1, pep2 };
            List<PeptideWithSetModifications> listForFile2 = new List<PeptideWithSetModifications> { pep2 };
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile1 = new TestDataFile(listForFile1);
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile2 = new TestDataFile(listForFile2);

            List<Protein> proteinList = new List<Protein> { prot1 };

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlFilePath1, false);
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile2, mzmlFilePath2, false);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList, proteinDbFilePath);

            string output_folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles");
            Directory.CreateDirectory(output_folder);

            st.RunTask(
                output_folder,
                new List<DbForTask> { new DbForTask(proteinDbFilePath, false) },
                new List<string> { mzmlFilePath1, mzmlFilePath2, },
                null);
        }

        #endregion Public Methods
    }
}