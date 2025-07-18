﻿using EngineLayer;
using NUnit.Framework;
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
using System.Text.RegularExpressions;
using Omics;

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
            ProteinGroup proteinGroup1 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            List<Protein> proteinList2 = new List<Protein> { prot1 };
            ProteinGroup proteinGroup2 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList2),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            //two protein groups with the same protein should be equal
            Assert.That(proteinGroup1.Equals(proteinGroup2));

            Protein prot3 = new Protein("EDEEK", "prot3");
            List<Protein> proteinList3 = new List<Protein> { prot3 };
            ProteinGroup proteinGroup3 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList3),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            //two protein groups with different proteins should not be equal
            Assert.That(!proteinGroup1.Equals(proteinGroup3));

            List<Protein> proteinList4 = new List<Protein> { prot1, prot3 };
            List<Protein> proteinList5 = new List<Protein> { prot3, prot1 };
            ProteinGroup proteinGroup4 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList4),
                               new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            ProteinGroup proteinGroup5 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList5),
                new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());

            //protein groups with the same proteins but in different order should be equal
            Assert.That(proteinGroup4.Equals(proteinGroup5));

            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(prot1,new DigestionParams(),1,3,CleavageSpecificity.Full,"",0,new Dictionary<int, Modification>(),0);
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(prot1, new DigestionParams(), 4, 6, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            ProteinGroup proteinGroup6 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList1), new HashSet<IBioPolymerWithSetMods>(){pwsm1},
                new HashSet<IBioPolymerWithSetMods>(){pwsm1});
            ProteinGroup proteinGroup7 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList1), new HashSet<IBioPolymerWithSetMods>() { pwsm2 },
                new HashSet<IBioPolymerWithSetMods>() { pwsm2 });

            //protein groups with the same proteins but different peptides should be equal
            Assert.That(proteinGroup6.Equals(proteinGroup7));

            //a protein group that is null should not be equal to a protein group that is not null
            ProteinGroup nullProteinGroup = null;
            Assert.That(!proteinGroup1.Equals(nullProteinGroup));
        }

        [Test]
        public static void ProteinGroupToStringTest()
        {
            Protein prot1 = new Protein("MEDEEK", "prot1");
            Protein prot2 = new Protein("MENEEK", "prot2");

            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(prot1, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(prot2, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            List<Protein> proteinList1 = new List<Protein> { prot1, prot2 };

            ProteinGroup proteinGroup1 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods>() { pwsm1, pwsm2 }, new HashSet<IBioPolymerWithSetMods>() { pwsm1, pwsm2 });

            //string exectedProteinGroupToString = proteinGroup1.ToString();
            string exectedProteinGroupToString = "prot1|prot2\t|\t\t\t779.30073507823|778.3167194953201\t2\t\t\t2\t2\t\t\t\t\t\t0\tT\t0\t0\t0\t0\t0";
            Assert.That(proteinGroup1.ToString(), Is.EqualTo(exectedProteinGroupToString));


            Protein prot3 = new Protein("MAAADAAAAAAAAAAAAAAA", "prot3", isDecoy: true);
            List<Protein> proteinList3 = new List<Protein> { prot3 };
            ProteinGroup proteinGroup3 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList3),
                               new HashSet<IBioPolymerWithSetMods>(), new HashSet<IBioPolymerWithSetMods>());
            string exectedProteinGroupWithDecoyToString = "prot1|prot2\t|\t\t\t779.30073507823|778.3167194953201\t2\t\t\t2\t2\t\t\t\t\t\t0\tT\t0\t0\t0\t0\t0";
            Assert.That(proteinGroup1.ToString(), Is.EqualTo(exectedProteinGroupWithDecoyToString));
        }

        [Test]
        public static void TestProteinGroupStringAndHeaderHaveSameNumberOfTabs()
        {
            Protein prot1 = new Protein("MEDEEK", "prot1");
            Protein prot2 = new Protein("MENEEK", "prot2");

            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(prot1, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(prot2, new DigestionParams(), 1, 3, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            List<Protein> proteinList1 = new List<Protein> { prot1, prot2 };

            ProteinGroup proteinGroup1 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods>() { pwsm1, pwsm2 }, new HashSet<IBioPolymerWithSetMods>() { pwsm1, pwsm2 });

            string pgHeader = proteinGroup1.GetTabSeparatedHeader();
            string pgRow = proteinGroup1.ToString();
            string[] headerFields = pgHeader.Split('\t');
            string[] rowEntries = pgRow.Split("\t");
            Assert.That(headerFields.Length, Is.EqualTo(rowEntries.Length));
            Assert.That(Regex.Matches(pgHeader, @"\t").Count, Is.EqualTo(Regex.Matches(pgRow, @"\t").Count));
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

            ProteinGroup proteinGroup1 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList1),
                               new HashSet<IBioPolymerWithSetMods>(){pwsm1,pwsm2}, new HashSet<IBioPolymerWithSetMods>() { pwsm1,pwsm2 });
            ProteinGroup proteinGroup2 = new ProteinGroup(new HashSet<IBioPolymer>(proteinList2),
                               new HashSet<IBioPolymerWithSetMods>() { pwsm3,pwsm4 }, new HashSet<IBioPolymerWithSetMods>() { pwsm3,pwsm4 });

            proteinGroup1.MergeProteinGroupWith(proteinGroup2);

            //protein group 3 should have all proteins from protein group 1 and 2
            Assert.That(proteinGroup1.Proteins.Contains(prot1));
            Assert.That(proteinGroup1.Proteins.Contains(prot2));
            Assert.That(proteinGroup1.Proteins.Contains(prot3));
            Assert.That(proteinGroup1.Proteins.Contains(prot4));

            //protein group 3 should have no peptides
            Assert.That(proteinGroup1.AllPeptides.Count, Is.EqualTo(4));
            Assert.That(proteinGroup1.UniquePeptides.Count, Is.EqualTo(4));
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

            EngineLayer.ProteinGroup proteinGroup1 = new EngineLayer.ProteinGroup(new HashSet<IBioPolymer>(proteinList1),
                new HashSet<IBioPolymerWithSetMods>() { pwsm1, pwsm2 }, new HashSet<IBioPolymerWithSetMods>() { pwsm1, pwsm2 });

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
            Assert.That(proteins.Count, Is.EqualTo(18));
            int totalNumberOfMods = proteins.Sum(p => p.OneBasedPossibleLocalizedModifications.Count + p.SequenceVariations.Sum(sv => sv.OneBasedModifications.Count));

            //tests that modifications are being done correctly
            Assert.That(totalNumberOfMods, Is.EqualTo(8));

            List<string> proteinGroupsOutput = File.ReadAllLines(Path.Combine(outputFolder, "task2", "AllQuantifiedProteinGroups.tsv")).ToList();
            string firstDataLine = proteinGroupsOutput[2];
            string modInfoListProteinTwo = firstDataLine.Split('\t')[14];
            Assert.That(modInfoListProteinTwo, Is.EqualTo("#aa71[Oxidation on S,info:occupancy=0.33(1/3)];#aa72[Didehydro on Y,info:occupancy=0.33(1/3)]"));

            Directory.Delete(outputFolder, true);
        }
    }
}

