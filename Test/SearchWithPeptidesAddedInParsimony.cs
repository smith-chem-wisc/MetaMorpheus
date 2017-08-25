using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class SearchWithPeptidesAddedInParsimony
    {
        #region Public Methods

        [Test]
        public static void SearchWithPeptidesAddedInParsimonyTest()
        {
            // Make sure can run the complete search task when multiple compact peptides may correspond to a single PWSM
            SearchTask st = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchDecoy = false,
                    ModPeptidesAreUnique = false
                },
                CommonParameters = new CommonParameters
                {
                    ScoreCutoff = 1,
                    MinPeptideLength = 1
                }
            };

            string xmlName = "andguiaheow.xml";

            #region Generate protein and write to file

            CommonParameters CommonParameters = new CommonParameters
            {
                MaxMissedCleavages = 0,
                MinPeptideLength = null,
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                ScoreCutoff = 1
            };
            ModificationMotif.TryGetMotif("A", out ModificationMotif motifA);
            ModificationWithMass alanineMod = new ModificationWithMass("111", "mt", motifA, TerminusLocalization.Any, 111);

            var variableModifications = new List<ModificationWithMass>();
            IDictionary<int, List<Modification>> oneBasedModifications1 = new Dictionary<int, List<Modification>>
                {
                    {2, new List<Modification>{ alanineMod } }
                };
            Protein protein1 = new Protein("MA", "protein1", oneBasedModifications: oneBasedModifications1);
            // Alanine = Glycine + CH2

            ModificationMotif.TryGetMotif("G", out ModificationMotif motif1);

            ModificationWithMass glycineMod = new ModificationWithMass("CH2 on Glycine", "mt", motif1, TerminusLocalization.Any, Chemistry.ChemicalFormula.ParseFormula("CH2").MonoisotopicMass);

            IDictionary<int, List<Modification>> oneBasedModifications2 = new Dictionary<int, List<Modification>>
                {
                    {2, new List<Modification>{glycineMod} }
                };
            Protein protein2 = new Protein("MG", "protein3", oneBasedModifications: oneBasedModifications2);

            var prot1List = protein1.Digest(CommonParameters.Protease, CommonParameters.MaxMissedCleavages, CommonParameters.MinPeptideLength, CommonParameters.MaxPeptideLength, CommonParameters.InitiatorMethionineBehavior, new List<ModificationWithMass>());
            PeptideWithPossibleModifications pepWithPossibleModifications = prot1List.First();
            int maxModsForPeptide = 1;
            int maximumVariableModificationIsoforms = 2;
            var pep1list = pepWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maximumVariableModificationIsoforms, maxModsForPeptide);
            PeptideWithSetModifications pepMA = pep1list.First();
            PeptideWithSetModifications pepMA111 = pep1list.Last();

            var prot2List = protein2.Digest(CommonParameters.Protease, CommonParameters.MaxMissedCleavages, CommonParameters.MinPeptideLength, CommonParameters.MaxPeptideLength, CommonParameters.InitiatorMethionineBehavior, new List<ModificationWithMass>());
            pepWithPossibleModifications = prot2List.First();
            var pep2list = pepWithPossibleModifications.GetPeptidesWithSetModifications(variableModifications, maximumVariableModificationIsoforms, maxModsForPeptide);
            PeptideWithSetModifications pepMG = pep2list.First();

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein1, protein2 }, xmlName);

            #endregion Generate protein and write to file

            string mzmlName = @"ajgdiu.mzML";

            #region Generate and write the mzml

            {
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepMA, pepMG, pepMA111 }, true);

                IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            }

            #endregion Generate and write the mzml

            st.RunTask("",
                new List<DbForTask> { new DbForTask(xmlName, false) },
                new List<string> { mzmlName }, "");
        }

        #endregion Public Methods
    }
}