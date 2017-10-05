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

namespace Test
{
    [TestFixture]
    public static class BinGenerationTest
    {
        #region Public Methods

        [Test]
        public static void TestBinGeneration()
        {
            SearchTask st = new SearchTask
            {
                CommonParameters = new CommonParameters
                {
                    ScoreCutoff = 1,
                    DigestionParams = new DigestionParams
                    {
                        InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                    },
                    ConserveMemory = false,
                },
                SearchParameters = new SearchParameters
                {
                    DoHistogramAnalysis = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Open,
                    DecoyType = DecoyType.None,
                    DoParsimony = true,
                    DoQuantification = true
                },
            };

            string proteinDbFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.xml");
            string mzmlFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.mzml");

            Protein prot1 = new Protein("MEDEEK", "prot1");
            Protein prot2 = new Protein("MENEEK", "prot2");

            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            ModificationWithMass mod = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, 10);

            var possMod1 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>()).First();
            var pep1_0 = possMod1.GetPeptidesWithSetModifications(st.CommonParameters.DigestionParams, new List<ModificationWithMass>()).First();
            var pep1_10 = possMod1.GetPeptidesWithSetModifications(st.CommonParameters.DigestionParams, new List<ModificationWithMass> { mod }).Last();

            Protein prot3 = new Protein("MAAADAAAAAAAAAAAAAAA", "prot3");

            var possMod3 = prot3.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>()).First();
            var pep2_0 = possMod3.GetPeptidesWithSetModifications(st.CommonParameters.DigestionParams, new List<ModificationWithMass>()).First();
            var pep2_10 = possMod3.GetPeptidesWithSetModifications(st.CommonParameters.DigestionParams, new List<ModificationWithMass> { mod }).Last();

            Protein prot4 = new Protein("MNNDNNNN", "prot4");
            var possMod4 = prot4.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>()).First();
            var pep3_10 = possMod4.GetPeptidesWithSetModifications(st.CommonParameters.DigestionParams, new List<ModificationWithMass> { mod }).Last();

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
            SearchTask st = new SearchTask()
            {
                CommonParameters = new CommonParameters
                {
                    ScoreCutoff = 1,
                    DigestionParams = new DigestionParams
                    {
                        InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                        MaxMissedCleavages = 0
                    },
                    ConserveMemory = false,
                },
                SearchParameters = new SearchParameters
                {
                    DoHistogramAnalysis = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Open,
                    MatchBetweenRuns = true,
                    DoQuantification = true
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

            var possMod1 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>()).First();
            var pep1 = possMod1.GetPeptidesWithSetModifications(st.CommonParameters.DigestionParams, new List<ModificationWithMass>()).First();
            var pep2 = possMod1.GetPeptidesWithSetModifications(st.CommonParameters.DigestionParams, new List<ModificationWithMass>()).Last();

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