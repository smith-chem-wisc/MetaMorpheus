using EngineLayer;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using Proteomics.ProteolyticDigestion;
using MassSpectrometry;
using Chemistry;
using TaskLayer;
using ProteinGroup = EngineLayer.ProteinGroup;
using System.IO;
using Omics.Digestion;
using Omics.Modifications;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class ProteinGroupTest
    {
        [Test]
        public static void TestProteinGroupEquals()
        {
            Protein prot1 = new Protein("MEDEEK", "prot1");
            List<Protein> proteinList1 = new List<Protein> { prot1 };
            ProteinGroup proteinGroup1 = new ProteinGroup(new HashSet<Protein>(proteinList1),
                new HashSet<PeptideWithSetModifications>(), new HashSet<PeptideWithSetModifications>());

            List<Protein> proteinList2 = new List<Protein> { prot1 };
            ProteinGroup proteinGroup2 = new ProteinGroup(new HashSet<Protein>(proteinList2),
                new HashSet<PeptideWithSetModifications>(), new HashSet<PeptideWithSetModifications>());

            //two protein groups with the same protein should be equal
            Assert.IsTrue(proteinGroup1.Equals(proteinGroup2));

            Protein prot3 = new Protein("EDEEK", "prot3");
            List<Protein> proteinList3 = new List<Protein> { prot3 };
            ProteinGroup proteinGroup3 = new ProteinGroup(new HashSet<Protein>(proteinList3),
                new HashSet<PeptideWithSetModifications>(), new HashSet<PeptideWithSetModifications>());

            //two protein groups with different proteins should not be equal
            Assert.IsFalse(proteinGroup1.Equals(proteinGroup3));

            List<Protein> proteinList4 = new List<Protein> { prot1, prot3 };
            List<Protein> proteinList5 = new List<Protein> { prot3, prot1 };
            ProteinGroup proteinGroup4 = new ProteinGroup(new HashSet<Protein>(proteinList4),
                               new HashSet<PeptideWithSetModifications>(), new HashSet<PeptideWithSetModifications>());
            ProteinGroup proteinGroup5 = new ProteinGroup(new HashSet<Protein>(proteinList5),
                new HashSet<PeptideWithSetModifications>(), new HashSet<PeptideWithSetModifications>());

            //protein groups with the same proteins but in different order should be equal
            Assert.IsTrue(proteinGroup4.Equals(proteinGroup5));

            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(prot1,new DigestionParams(),1,3,CleavageSpecificity.Full,"",0,new Dictionary<int, Modification>(),0);
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(prot1, new DigestionParams(), 4, 6, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            ProteinGroup proteinGroup6 = new ProteinGroup(new HashSet<Protein>(proteinList1), new HashSet<PeptideWithSetModifications>(){pwsm1},
                new HashSet<PeptideWithSetModifications>(){pwsm1});
            ProteinGroup proteinGroup7 = new ProteinGroup(new HashSet<Protein>(proteinList1), new HashSet<PeptideWithSetModifications>() { pwsm2 },
                new HashSet<PeptideWithSetModifications>() { pwsm2 });

            //protein groups with the same proteins but different peptides should be equal
            Assert.IsTrue(proteinGroup6.Equals(proteinGroup7));

            //a protein group that is null should not be equal to a protein group that is not null
            ProteinGroup nullProteinGroup = null;
            Assert.IsFalse(proteinGroup1.Equals(nullProteinGroup));
        }

        [Test]
        public static void ProteinGroupToStringTest()
        {
            Protein prot1 = new Protein("MEDEEK", "prot1");
            Protein prot2 = new Protein("MENEEK", "prot2");

            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(prot1, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(prot2, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            List<Protein> proteinList1 = new List<Protein> { prot1, prot2 };

            ProteinGroup proteinGroup1 = new ProteinGroup(new HashSet<Protein>(proteinList1),
                new HashSet<PeptideWithSetModifications>() { pwsm1, pwsm2 }, new HashSet<PeptideWithSetModifications>() { pwsm1, pwsm2 });

            //string exectedProteinGroupToString = proteinGroup1.ToString();
            string exectedProteinGroupToString = "prot1|prot2\t|\t\t\t779.30073507823|778.3167194953201\t2\t\t\t2\t2\t\t\t\t\t\t0\tT\t0\t0\t0\t0\t0\t";
            Assert.AreEqual(exectedProteinGroupToString, proteinGroup1.ToString());


            Protein prot3 = new Protein("MAAADAAAAAAAAAAAAAAA", "prot3", isDecoy:true);
            List<Protein> proteinList3 = new List<Protein> { prot3 };
            ProteinGroup proteinGroup3 = new ProteinGroup(new HashSet<Protein>(proteinList3),
                               new HashSet<PeptideWithSetModifications>(), new HashSet<PeptideWithSetModifications>());
            string exectedProteinGroupWithDecoyToString = "prot1|prot2\t|\t\t\t779.30073507823|778.3167194953201\t2\t\t\t2\t2\t\t\t\t\t\t0\tT\t0\t0\t0\t0\t0\t";
            Assert.AreEqual(exectedProteinGroupWithDecoyToString, proteinGroup1.ToString());
        }

        [Test]
        public static void ProteinGroupMergeTest()
        {
            Protein prot1 = new Protein("MEDEEK", "prot1");
            Protein prot2 = new Protein("MENEEK", "prot2");
            Protein prot3 = new Protein("MAAADAAAAAAAAAAAAAAA", "prot3");
            Protein prot4 = new Protein("MNNDNNNN", "prot4");

            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(prot1, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(prot2, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pwsm3 = new PeptideWithSetModifications(prot3, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pwsm4 = new PeptideWithSetModifications(prot4, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            List<Protein> proteinList1 = new List<Protein> { prot1, prot2 };
            List<Protein> proteinList2 = new List<Protein> { prot3, prot4 };

            ProteinGroup proteinGroup1 = new ProteinGroup(new HashSet<Protein>(proteinList1),
                               new HashSet<PeptideWithSetModifications>(){pwsm1,pwsm2}, new HashSet<PeptideWithSetModifications>() { pwsm1,pwsm2 });
            ProteinGroup proteinGroup2 = new ProteinGroup(new HashSet<Protein>(proteinList2),
                               new HashSet<PeptideWithSetModifications>() { pwsm3,pwsm4 }, new HashSet<PeptideWithSetModifications>() { pwsm3,pwsm4 });

            proteinGroup1.MergeProteinGroupWith(proteinGroup2);

            //protein group 3 should have all proteins from protein group 1 and 2
            Assert.IsTrue(proteinGroup1.Proteins.Contains(prot1));
            Assert.IsTrue(proteinGroup1.Proteins.Contains(prot2));
            Assert.IsTrue(proteinGroup1.Proteins.Contains(prot3));
            Assert.IsTrue(proteinGroup1.Proteins.Contains(prot4));

            //protein group 3 should have no peptides
            Assert.AreEqual(4, proteinGroup1.AllPeptides.Count);
            Assert.AreEqual(4, proteinGroup1.UniquePeptides.Count);
        }

        [Test]
        public static void ProteinGroupDisplayModsTestWithGetIdentifiedPeptidesMethod()
        {
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif1);

            Dictionary<DissociationType, List<double>> NeutralLosses = new Dictionary<DissociationType, List<double>>();
            NeutralLosses.Add(DissociationType.HCD, new List<double> { 0 });

            Modification modFormula_C1 = new Modification(_originalId: "modC", _accession: "", _modificationType: "mt", _featureType: "", _target: motif1, _locationRestriction: "Anywhere.", _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("C1")), null, null, null, null, _neutralLosses: NeutralLosses, null, null);
            Modification modFormula_H1 = new Modification(_originalId: "modH", _accession: "", _modificationType: "mt", _featureType: "", _target: motif1, _locationRestriction: "Anywhere.", _chemicalFormula: new ChemicalFormula(ChemicalFormula.ParseFormula("H1")), null, null, null, null, _neutralLosses: NeutralLosses, null, null);

            IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification>{ modFormula_C1, modFormula_H1 }},
            };
            Protein protein1 = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: oneBasedModifications);

            Dictionary<int, Modification> allModsOneIsNterminus1 = new Dictionary<int, Modification>
            {
                {2, modFormula_C1},
            };

            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus1, 0);

            Dictionary<int, Modification> allModsOneIsNterminus2 = new Dictionary<int, Modification>
            {
                {2,modFormula_H1 },
            };

            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus2, 0);

            List<Protein> proteinList1 = new List<Protein> { protein1 };

            EngineLayer.ProteinGroup proteinGroup1 = new EngineLayer.ProteinGroup(new HashSet<Protein>(proteinList1),
                new HashSet<PeptideWithSetModifications>() { pwsm1, pwsm2 }, new HashSet<PeptideWithSetModifications>() { pwsm1, pwsm2 });

            proteinGroup1.DisplayModsOnPeptides = true;

            //This test just gets some lines in ProteinGroup covered. There is no accessible way to get the output of this method.
            Assert.DoesNotThrow(()=>proteinGroup1.GetIdentifiedPeptidesOutput(new List<SilacLabel>()));
        }

        [Test]
        public static void TestModificationInfoListInProteinGroupsOutput()
        {
            //Create GPTMD Task
            //Create Search Task
            GptmdTask task1 = new GptmdTask
            {
                CommonParameters = new CommonParameters(),
                GptmdParameters = new GptmdParameters
                {
                    ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                        b.ModificationType.Equals("Common Artifact")
                        || b.ModificationType.Equals("Common Biological")
                        || b.ModificationType.Equals("Metal")
                        || b.ModificationType.Equals("Less Common")
                        ).Select(b => (b.ModificationType, b.IdWithMotif)).ToList()
                }
            };

            SearchTask task2 = new SearchTask
            {
                CommonParameters = new CommonParameters(),

                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchTarget = true,
                    WritePrunedDatabase = true,
                    SearchType = SearchType.Classic
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("task1", task1), ("task2", task2) };
            string mzmlName = @"TestData\PrunedDbSpectra.mzml";
            string fastaName = @"TestData\DbForPrunedDb.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPrunedGeneration");

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();
            string final = Path.Combine(MySetUpClass.outputFolder, "task2", "DbForPrunedDbGPTMDproteinPruned.xml");
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(final, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out var ok);
            // ensures that protein out put contains the correct number of proteins to match the following conditions.
            // all proteins in DB have baseSequence!=null (not ambiguous)
            // all proteins that belong to a protein group are written to DB
            Assert.AreEqual(18, proteins.Count);
            int totalNumberOfMods = proteins.Sum(p => p.OneBasedPossibleLocalizedModifications.Count + p.SequenceVariations.Sum(sv => sv.OneBasedModifications.Count));

            //tests that modifications are being done correctly
            Assert.AreEqual(8, totalNumberOfMods);

            List<string> proteinGroupsOutput = File.ReadAllLines(Path.Combine(outputFolder, "task2", "AllQuantifiedProteinGroups.tsv")).ToList();
            string firstDataLine = proteinGroupsOutput[2];
            string modInfoListProteinTwo = firstDataLine.Split('\t')[14];

            Assert.AreEqual("#aa71[Oxidation on S,info:occupancy=0.33(1/3)];#aa72[Didehydro on Y,info:occupancy=0.33(1/3)]", modInfoListProteinTwo);
            Directory.Delete(outputFolder, true);
        }
    }
}

