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
    public static class GptmdPrunedDbTests
    {
        // want a psm whose base sequence is not ambigous but full sequence is (ptm is not localized): make sure this does not make it in DB

        [Test]
        public static void TestPrunedGeneration()
        {
            //Create GPTMD Task
            //Create Search Task
            GptmdTask task1 = new GptmdTask
            {
                CommonParameters = new CommonParameters(),
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

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>
            {
                ("task1", task1),
                ("task2", task2)
            };

            string mzmlName = @"TestData\PrunedDbSpectra.mzml";
            string fastaName = @"TestData\DbForPrunedDb.fasta";
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, Environment.CurrentDirectory);
            engine.Run();
            string final = Path.Combine(MySetUpClass.outputFolder, "task2", "DbForPrunedDbGPTMDproteinPruned.xml");
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(final, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out var ok);

            // ensures that protein out put contins the correct number of proteins to match the folowing conditions.
            // all proteins in DB have baseSequence!=null (not ambiguous)
            // all proteins that belong to a protein group are written to DB
            Assert.AreEqual(18, proteins.Count); //used to be 20 before I changed the list of mods
            int totalNumberOfMods = 0;
            foreach (Protein p in proteins)
            {
                int numberOfMods = p.OneBasedPossibleLocalizedModifications.Count();
                totalNumberOfMods = totalNumberOfMods + numberOfMods;
            }

            //tests that modifications are being done correctly
            Assert.AreEqual(totalNumberOfMods, 0);
        }

        //test if prunedDatabase matches expected output
        [Test]
        public static void TestPrunedDatabase()
        {
            //#region setup

            //Create Search Task
            SearchTask task1 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    SearchTarget = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                    ModsToWriteSelection = new Dictionary<string, int>
                    {
                        {"ConnorModType", 1}
                    }
                },
                CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5))
            };

            //add task to task list
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>
            {
               ("task1", task1)
            };

            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);

            var connorMod = new ModificationWithMass("ConnorMod", "ConnorModType", motif, TerminusLocalization.Any, 10);

            GlobalVariables.AddMods(new List<ModificationWithLocation>
            {
                connorMod
            });

            //#endregion setup

            //#region Protein and Mod Creation

            //create modification lists
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where
                (b => task1.CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();

            //add modification to Protein object
            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationWithMass modToAdd = connorMod;
            ModificationWithMass modToAdd2 = connorMod;
            dictHere.Add(1, new List<Modification> { modToAdd });
            dictHere.Add(3, new List<Modification> { modToAdd2 });

            //protein Creation (One with mod and one without)
            Protein TestProteinWithMod = new Protein("PEPTID", "accession1", "organism", new List<Tuple<string, string>>(), dictHere);

            //#endregion Protein and Mod Creation

            //#region XML File

            //First Write XML Database

            string xmlName = "okkk.xml";

            //Add Mod to list and write XML input database
            Dictionary<string, HashSet<Tuple<int, Modification>>> modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var Hash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(3, modToAdd)
            };
            modList.Add("test", Hash);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProteinWithMod }, xmlName);

            //#endregion XML File

            //#region MZML File

            //now write MZML file
            var protein = ProteinDbLoader.LoadProteinXML(xmlName, true,
                DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);
            var digestedList = protein[0].Digest(task1.CommonParameters.DigestionParams, new List<ModificationWithMass> { },
                variableModifications).ToList();
            Assert.AreEqual(4, digestedList.Count);

            //Set Peptide with 1 mod at position 3
            PeptideWithSetModifications pepWithSetMods1 = digestedList[1];

            //Finally Write MZML file
            Assert.AreEqual("PEP[ConnorModType:ConnorMod]TID", pepWithSetMods1.Sequence);
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1 });
            string mzmlName = @"hello.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            //#endregion MZML File

            //run!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName },
                new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();

            string final = Path.Combine(MySetUpClass.outputFolder, "task1", "okkkpruned.xml");

            var proteins = ProteinDbLoader.LoadProteinXML(final, true,
                DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out ok);

            //check length
            Assert.AreEqual(1, proteins[0].OneBasedPossibleLocalizedModifications.Count);
            //check location (key)
            Assert.AreEqual(true, proteins[0].OneBasedPossibleLocalizedModifications.ContainsKey(3));
            List<Modification> listOfMods = proteins[0].OneBasedPossibleLocalizedModifications[3];
            //check Type, count, ID
            Assert.AreEqual("ConnorModType", listOfMods[0].modificationType);
            Assert.AreEqual("ConnorMod", listOfMods[0].id);
            Assert.AreEqual(1, listOfMods.Count);
        }

        [Test]
        public static void TestUserModSelectionInPrunedDB()
        {
            //#region setup

            List<(string, string)> listOfModsFixed = new List<(string, string)> { ("Common Fixed", "Carbamidomethyl of C"), ("Common Fixed", "Carbamidomethyl of U") };
            //Create Search Task
            SearchTask task5 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    SearchTarget = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                },
                CommonParameters = new CommonParameters(listOfModsFixed: listOfModsFixed)
            };

            task5.SearchParameters.ModsToWriteSelection["Mod"] = 0;
            task5.SearchParameters.ModsToWriteSelection["Common Fixed"] = 1;
            task5.SearchParameters.ModsToWriteSelection["Glycan"] = 2;
            task5.SearchParameters.ModsToWriteSelection["missing"] = 3;

            //add task 1 to task list
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("task5", task5) };
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif2);

            var connorMod = new ModificationWithMass("ModToNotAppear", "Mod", motif, TerminusLocalization.Any, 10);
            var connorMod2 = new ModificationWithMass("Default(Mod in DB and Observed)", "Common Fixed", motif, TerminusLocalization.Any, 10);
            var connorMod3 = new ModificationWithMass("ModToAlwaysAppear", "Glycan", motif, TerminusLocalization.Any, 10);
            var connorMod4 = new ModificationWithMass("ModObservedNotinDB", "missing", motif2, TerminusLocalization.Any, 5);

            GlobalVariables.AddMods(new List<ModificationWithLocation>
            {
                connorMod,
                connorMod2,
                connorMod3,
                connorMod4
            });

            //#endregion setup

            //#region Protein and Mod Creation

            //create modification lists
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => task5.CommonParameters.ListOfModsVariable.Contains
            ((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => task5.CommonParameters.ListOfModsFixed.Contains
            ((b.modificationType, b.id))).ToList();

            //add modification to Protein object
            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationWithMass modToAdd = connorMod;
            ModificationWithMass modToAdd2 = connorMod2;
            ModificationWithMass modToAdd3 = connorMod3;
            ModificationWithMass modToAdd4 = connorMod4;

            //add Fixed modifcation so can test if mod that is observed and not in DB
            fixedModifications.Add(connorMod4);
            listOfModsFixed.Add((connorMod4.modificationType, connorMod4.id));

            dictHere.Add(1, new List<Modification> { modToAdd });
            dictHere.Add(2, new List<Modification> { modToAdd2 }); //default
            dictHere.Add(3, new List<Modification> { modToAdd3 }); //Alway Appear

            var dictHere2 = new Dictionary<int, List<Modification>>();
            dictHere2.Add(1, new List<Modification> { modToAdd });
            dictHere2.Add(2, new List<Modification> { modToAdd2 }); //default
            dictHere2.Add(3, new List<Modification> { modToAdd3 }); //Alway Appear
            dictHere2.Add(4, new List<Modification> { modToAdd4 });//observed
            //protein Creation (One with mod and one without)
            Protein TestProteinWithModForDB = new Protein("PPPPPPPPPPE", "accession1", "organism", new List<Tuple<string, string>>(), dictHere);
            Protein TestProteinWithModObsevred = new Protein("PPPPPPPPPPE", "accession1", "organism", new List<Tuple<string, string>>(), dictHere2);

            //#endregion Protein and Mod Creation

            //#region XML File

            //First Write XML Database

            string xmlName = "selectedMods.xml";
            string xmlName2 = "selectedModsObvs.xml";

            //Add Mod to list and write XML input database
            Dictionary<string, HashSet<Tuple<int, Modification>>> modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var Hash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(1, modToAdd),
                new Tuple<int, Modification>(2, modToAdd2),
                new Tuple<int, Modification>(3, modToAdd3),
                new Tuple<int, Modification>(4, modToAdd4), //Observed Only
            };

            modList.Add("test", Hash);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProteinWithModForDB }, xmlName);

            //Add Observed Only
            modList.Add("test2", Hash);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProteinWithModObsevred }, xmlName2);

            //#endregion XML File

            //#region MZML File

            //now create MZML data
            var protein = ProteinDbLoader.LoadProteinXML(xmlName2, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);

            var digestedList = protein[0].Digest(task5.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            //Set Peptide with 1 mod at position 3
            PeptideWithSetModifications pepWithSetMods1 = digestedList[0];
            PeptideWithSetModifications pepWithSetMods2 = digestedList[1];
            PeptideWithSetModifications pepWithSetMods3 = digestedList[2];
            PeptideWithSetModifications pepWithSetMods4 = digestedList[3];
            PeptideWithSetModifications pepWithSetMods5 = digestedList[4];

            //CUSTOM PEP
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications>
            { pepWithSetMods1, pepWithSetMods2, pepWithSetMods3, pepWithSetMods4, pepWithSetMods5 });
            string mzmlName = @"newMzml.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            //#endregion MZML File

            //make sure this runs correctly
            //run!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();
            string final = Path.Combine(MySetUpClass.outputFolder, "task5", "selectedModspruned.xml");
            var proteins = ProteinDbLoader.LoadProteinXML(final, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out ok);
            var Dlist = proteins[0].Digest(task5.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
            Assert.AreEqual(Dlist[0].NumFixedMods, 1);

            //check length
            Assert.AreEqual(proteins[0].OneBasedPossibleLocalizedModifications.Count, 3);
            List<Modification> listOfLocalMods = new List<Modification>();
            listOfLocalMods.AddRange(proteins[0].OneBasedPossibleLocalizedModifications[2]);
            listOfLocalMods.AddRange(proteins[0].OneBasedPossibleLocalizedModifications[3]);
            listOfLocalMods.AddRange(proteins[0].OneBasedPossibleLocalizedModifications[11]);

            //check Type, count, ID
            Assert.AreEqual(listOfLocalMods[0].modificationType, "Common Fixed");
            Assert.AreEqual(listOfLocalMods[2].modificationType, "missing");
            Assert.IsFalse(listOfLocalMods.Contains(connorMod)); //make sure that mod set not to show up is not in mod list

            Assert.AreEqual(listOfLocalMods[0].id, "Default(Mod in DB and Observed)");
            Assert.AreEqual(listOfLocalMods[1].id, "ModToAlwaysAppear");
            //Makes sure Mod that was not in the DB but was observed is in pruned DB
            Assert.AreEqual(listOfLocalMods[2].id, "ModObservedNotinDB");
            Assert.AreEqual(listOfLocalMods.Count, 3);
        }
    }
}