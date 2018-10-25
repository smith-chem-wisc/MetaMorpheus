using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class SearchWithPeptidesAddedInParsimony
    {
        [Test]
        public static void SearchWithPeptidesAddedInParsimonyTest()
        {
            // Make sure can run the complete search task when multiple compact peptides may correspond to a single PWSM
            SearchTask st = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DecoyType = DecoyType.None,
                    ModPeptidesAreDifferent = false
                },
                CommonParameters = new CommonParameters(scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 2)),
            };

            string xmlName = "andguiaheow.xml";

            CommonParameters CommonParameters = new CommonParameters(
                scoreCutoff: 1,
                digestionParams: new DigestionParams(
                    maxMissedCleavages: 0,
                    minPeptideLength: 1,
                    maxModificationIsoforms: 2,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                    maxModsForPeptides: 1));

            ModificationMotif.TryGetMotif("A", out ModificationMotif motifA);
            Modification alanineMod = new Modification(_originalId: "111", _modificationType: "mt", _target: motifA, _locationRestriction: "Anywhere.", _monoisotopicMass: 111);

            var variableModifications = new List<Modification>();
            IDictionary<int, List<Modification>> oneBasedModifications1 = new Dictionary<int, List<Modification>>
                {
                    {2, new List<Modification>{ alanineMod } }
                };
            Protein protein1 = new Protein("MA", "protein1", oneBasedModifications: oneBasedModifications1);
            // Alanine = Glycine + CH2

            ModificationMotif.TryGetMotif("G", out ModificationMotif motif1);

            Modification glycineMod = new Modification(_originalId: "CH2 on Glycine", _modificationType: "mt", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: Chemistry.ChemicalFormula.ParseFormula("CH2").MonoisotopicMass);

            IDictionary<int, List<Modification>> oneBasedModifications2 = new Dictionary<int, List<Modification>>
                {
                    {2, new List<Modification>{glycineMod} }
                };
            Protein protein2 = new Protein("MG", "protein3", oneBasedModifications: oneBasedModifications2);

            PeptideWithSetModifications pepMA = protein1.Digest(CommonParameters.DigestionParams, new List<Modification>(), variableModifications).First();
            PeptideWithSetModifications pepMA111 = protein1.Digest(CommonParameters.DigestionParams, new List<Modification>(), variableModifications).Last();

            var pepMG = protein2.Digest(CommonParameters.DigestionParams, new List<Modification>(), variableModifications).First();

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein1, protein2 }, xmlName);

            string mzmlName = @"ajgdiu.mzML";

            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepMA, pepMG, pepMA111 }, true);

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSearchWithPeptidesAddedInParsimony");
            Directory.CreateDirectory(outputFolder);

            st.RunTask(outputFolder,
                new List<DbForTask> { new DbForTask(xmlName, false) },
                new List<string> { mzmlName }, "");
            Directory.Delete(outputFolder, true);
            File.Delete(mzmlName);
            File.Delete(xmlName);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }
    }
}