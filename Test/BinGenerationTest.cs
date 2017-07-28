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
    public class BinGenerationTest
    {

        #region Public Methods

        [Test]
        public static void TestBinGeneration()
        {
            List<MassDiffAcceptor> massDiffAcceptors = new List<MassDiffAcceptor> { new OpenSearchMode() };
            SearchTask st = new SearchTask()
            {
                DoHistogramAnalysis = true,
                MassDiffAcceptors = massDiffAcceptors,
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                SearchDecoy=false
            };

            string proteinDbFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.xml");
            string mzmlFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.mzml");

            Protein prot1 = new Protein("MDDDDK", "prot1");
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("D", out motif);
            Modification mod = new ModificationWithMass(null, null, motif, ModificationSites.Any, 10, null, new List<double> { 0 }, null, null);
            Protein prot2 = new Protein("MDDDDK", "prot2", oneBasedModifications: new Dictionary<int, List<Modification>> { { 3, new List<Modification> { mod } } });
            Protein prot3 = new Protein("MDNDDK", "prot3");

            var possMod1 = prot1.Digest(st.Protease, st.MaxMissedCleavages, st.MinPeptideLength, st.MaxPeptideLength, st.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep1 = possMod1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).Last();

            var possMod2 = prot2.Digest(st.Protease, st.MaxMissedCleavages, st.MinPeptideLength, st.MaxPeptideLength, st.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep2 = possMod2.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).Last();

            Protein prot4 = new Protein("MAAAAAAAAAAAAAAAAAAA", "prot4", oneBasedModifications: new Dictionary<int, List<Modification>> { { 3, new List<Modification> { mod } } });

            var possMod4 = prot4.Digest(st.Protease, st.MaxMissedCleavages, st.MinPeptideLength, st.MaxPeptideLength, st.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep4 = possMod4.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            var pep4mod = possMod4.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).Last();

            Protein prot5 = new Protein("MAAAAAAAAAAAAAAAAAAA", "prot5");



            List<PeptideWithSetModifications> pepsWithSetMods = new List<PeptideWithSetModifications> { pep1, pep1, pep2, pep2,pep4, pep4mod };
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepsWithSetMods);

            List<Protein> proteinList = new List<Protein> { prot1, prot3 , prot5 };

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlFilePath, false);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList, proteinDbFilePath);

            st.RunTask(
                TestContext.CurrentContext.TestDirectory,
                new List<DbForTask> { new DbForTask(proteinDbFilePath, false) },
                new List<string> { mzmlFilePath },
                null);
        }

        #endregion Public Methods

    }
}