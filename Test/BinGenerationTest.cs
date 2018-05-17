﻿using EngineLayer;
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
                    DigestionParams = new DigestionParams(MinPeptideLength: 5, InitiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
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
            string mzmlFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.mzML");

            Protein prot1 = new Protein("MEDEEK", "prot1");
            Protein prot2 = new Protein("MENEEK", "prot2");

            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            ModificationWithMass mod = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, 10);

            var pep1_0 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            var pep1_10 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).Last();

            Protein prot3 = new Protein("MAAADAAAAAAAAAAAAAAA", "prot3");

            var pep2_0 = prot3.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            var pep2_10 = prot3.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass> { mod }).Last();

            Protein prot4 = new Protein("MNNDNNNN", "prot4");
            var pep3_10 = prot4.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass> { mod }).Last();

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

            Assert.AreEqual(3, File.ReadLines(Path.Combine(output_folder, @"aggregate.mytsv")).Count());
        }

        [Test]
        public static void TestProteinSplitAcrossFiles()
        {
            SearchTask st = new SearchTask()
            {
                CommonParameters = new CommonParameters
                {
                    ScoreCutoff = 1,
                    DigestionParams = new DigestionParams(MaxMissedCleavages: 0, MinPeptideLength: 5, InitiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
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
            string mzmlFilePath1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles1.mzML");
            string mzmlFilePath2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles2.mzML");

            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            ModificationWithMass mod = new ModificationWithMass("mod1", "mt", motif, TerminusLocalization.Any, 10);

            IDictionary<int, List<Modification>> oneBasedModification = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification>{ mod } }
            };

            Protein prot1 = new Protein("MEDEEK", "prot1", oneBasedModifications: oneBasedModification);

            var pep1 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            var pep2 = prot1.Digest(st.CommonParameters.DigestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).Last();

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