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
    public class MyTaskTest
    {

        #region Public Fields

        public static string outputFolder = "";
        public static bool hasPrunedRun = false;

        #endregion Public Fields

        #region Public Methods

        [Test]
        public static void TestEverythingRunner()
        {
            #region Setup tasks

            foreach (var modFile in Directory.GetFiles(@"Mods"))
                GlobalTaskLevelSettings.AddMods(PtmListLoader.ReadModsFromFile(modFile));

            CalibrationTask task1 = new CalibrationTask();
            GptmdTask task2 = new GptmdTask();

            SearchTask task3 = new SearchTask()
            {
                DoParsimony = true
            };
            SearchTask task4 = new SearchTask()
            {
                ClassicSearch = false
            };
            List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>> {
                new Tuple<string, MetaMorpheusTask>("task1", task1),
                new Tuple<string, MetaMorpheusTask>("task2", task2),
                new Tuple<string, MetaMorpheusTask>("task3", task3),
                new Tuple<string, MetaMorpheusTask>("task4", task4),};

            #endregion Setup tasks

            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            Console.WriteLine("Size of variable Modificaitaons: " + variableModifications.Capacity);
            Console.WriteLine("Size of fixed Modificaitaons: " + fixedModifications.Capacity);
            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>(), new List<SequenceVariation>());

            var digestedList = ParentProtein.Digest(task1.Protease, 0, null, null, InitiatorMethionineBehavior.Retain, fixedModifications).ToList();

            Assert.AreEqual(2, digestedList.Count);

            PeptideWithPossibleModifications modPep1 = digestedList[0];
            var setList1 = modPep1.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();

            Assert.AreEqual(2, setList1.Count);

            PeptideWithSetModifications pepWithSetMods1 = setList1[0];

            PeptideWithPossibleModifications modPep2 = digestedList[1];
            var setList2 = modPep2.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();

            Assert.AreEqual(1, setList2.Count);

            PeptideWithSetModifications pepWithSetMods2 = setList2[0];

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("E", out motif);
            dictHere.Add(3, new List<Modification> { new ModificationWithMass("21", null, motif, ModificationSites.Any, 21.981943, null, new List<double> { 0 }, new List<double> { 21.981943 }, "") });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", new List<Tuple<string, string>>(), dictHere, new int?[0], new int?[0], new string[0], null, null, false, false, null, new List<SequenceVariation>());
            digestedList = ParentProteinToNotInclude.Digest(task1.Protease, 0, null, null, InitiatorMethionineBehavior.Retain, fixedModifications).ToList();
            var modPep3 = digestedList[0];
            Assert.AreEqual(1, digestedList.Count);
            var setList3 = modPep3.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();
            Assert.AreEqual(4, setList3.Count);

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, setList3[1] });

            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[] { 4 }, new int?[] { 8 }, new string[] { "chain" }, "name2", "fullname2", false, false, new List<DatabaseReference>(), new List<SequenceVariation>());

            #region Write the files

            string mzmlName = @"ok.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            string xmlName = "okk.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), new List<Protein> { ParentProtein, proteinWithChain }, xmlName);

            #endregion Write the files

            // RUN!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) });
            engine.Run();
        }

        [Test]
        public static void TestMultipleFilesRunner()
        {
            #region Setup tasks

            foreach (var modFile in Directory.GetFiles(@"Mods"))
                GlobalTaskLevelSettings.AddMods(PtmListLoader.ReadModsFromFile(modFile));

            CalibrationTask task1 = new CalibrationTask();
            GptmdTask task2 = new GptmdTask();

            SearchTask task3 = new SearchTask();
            SearchTask task4 = new SearchTask()
            {
                DoParsimony = true
            };
            List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>> {
                new Tuple<string, MetaMorpheusTask>("task1", task1),
                new Tuple<string, MetaMorpheusTask>("task2", task2),
                new Tuple<string, MetaMorpheusTask>("task3", task3),
                new Tuple<string, MetaMorpheusTask>("task4", task4),};

            #endregion Setup tasks

            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>(), new List<SequenceVariation>());

            var digestedList = ParentProtein.Digest(task1.Protease, 0, null, null, InitiatorMethionineBehavior.Retain, fixedModifications).ToList();

            Assert.AreEqual(2, digestedList.Count);

            PeptideWithPossibleModifications modPep1 = digestedList[0];
            var setList1 = modPep1.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();

            Assert.AreEqual(2, setList1.Count);

            PeptideWithSetModifications pepWithSetMods1 = setList1[0];

            PeptideWithPossibleModifications modPep2 = digestedList[1];
            var setList2 = modPep2.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();

            Assert.AreEqual(1, setList2.Count);

            PeptideWithSetModifications pepWithSetMods2 = setList2[0];

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("E", out motif);
            dictHere.Add(3, new List<Modification> { new ModificationWithMass("21", null, motif, ModificationSites.Any, 21.981943, null, new List<double> { 0 }, new List<double> { 21.981943 }, "") });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", new List<Tuple<string, string>>(), dictHere, new int?[0], new int?[0], new string[0], null, null, false, false, null, new List<SequenceVariation>());
            digestedList = ParentProteinToNotInclude.Digest(task1.Protease, 0, null, null, InitiatorMethionineBehavior.Retain, fixedModifications).ToList();
            var modPep3 = digestedList[0];
            Assert.AreEqual(1, digestedList.Count);
            var setList3 = modPep3.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();
            Assert.AreEqual(4, setList3.Count);

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile1 = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, setList3[1] });

            string mzmlName1 = @"ok1.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName1, false);

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile2 = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, setList3[1] });

            string mzmlName2 = @"ok2.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile2, mzmlName2, false);

            Protein proteinWithChain1 = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[] { 4 }, new int?[] { 8 }, new string[] { "chain" }, "name2", "fullname2", false, false, new List<DatabaseReference>(), new List<SequenceVariation>());
            Protein proteinWithChain2 = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[] { 4 }, new int?[] { 8 }, new string[] { "chain" }, "name2", "fullname2", false, false, new List<DatabaseReference>(), new List<SequenceVariation>());

            string xmlName = "okk.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), new List<Protein> { ParentProtein, proteinWithChain1, proteinWithChain2 }, xmlName);

            // RUN!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName1, mzmlName2 }, new List<DbForTask> { new DbForTask(xmlName, false) });
            engine.Run();
        }

        [Test]
        public static void MakeSureGptmdTaskMatchesExactMatches()
        {
            MetaMorpheusTask task1;

            #region Setup tasks

            {
                ModificationMotif motif;
                ModificationMotif.TryGetMotif("T", out motif);
                GlobalTaskLevelSettings.AddMods(new List<ModificationWithMass> { new ModificationWithMass("ok", null, motif, ModificationSites.Any, 229, null, null, null, "okType") });
                task1 = new GptmdTask()
                {
                    ListOfModsGptmd = new List<Tuple<string, string>> { new Tuple<string, string>("okType", "ok") },
                    ListOfModsVariable = new List<Tuple<string, string>>(),
                    ListOfModsFixed = new List<Tuple<string, string>>(),
                    PrecursorMassTolerance = new MzLibUtil.Tolerance(MzLibUtil.ToleranceUnit.Absolute, 1)
                };
            }

            #endregion Setup tasks

            string xmlName = "sweetness.xml";

            #region Generate protein and write to file

            {
                Protein theProtein = new Protein("MPEPTIDEKANTHE", "accession1", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>(), new List<SequenceVariation>());
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), new List<Protein> { theProtein }, xmlName);
            }

            #endregion Generate protein and write to file

            string mzmlName = @"ok.mzML";

            #region Generate and write the mzml

            {
                Dictionary<string, Modification> ok;
                var theProteins = ProteinDbLoader.LoadProteinXML(xmlName, true, new List<Modification>(), false, new List<string>(), out ok);

                List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();

                var targetDigested = theProteins[0].Digest(GlobalTaskLevelSettings.ProteaseDictionary["trypsin"], 1, null, null, InitiatorMethionineBehavior.Retain, fixedModifications).ToList();

                ModificationMotif motif;
                ModificationMotif.TryGetMotif("T", out motif);
                var okjhjf = targetDigested[0].GetPeptidesWithSetModifications(GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList(), 1, 0).ToList();
                PeptideWithSetModifications targetGood = okjhjf.First();

                var okjhj = targetDigested[1].GetPeptidesWithSetModifications(GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList(), 2, 1).ToList();
                PeptideWithSetModifications targetWithUnknownMod = okjhj.Last();
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { targetGood, targetWithUnknownMod }, true);

                IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            }

            #endregion Generate and write the mzml

            // RUN!
            var theStringResult = task1.RunTask(TestContext.CurrentContext.TestDirectory, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();
            Assert.IsTrue(theStringResult.Contains("Modifications added: 1"));
        }

        //test if prunedDatabase matches expected output
        [Test]
        public static void TestPrunedDatabase()
        {
            hasPrunedRun = true;

            #region setup

            //Create Search Task
            SearchTask task1 = new SearchTask()
            {
                WritePrunedDatabase = true,
                ListOfModsLocalize = new List<Tuple<string, string>> { new Tuple<string, string>("ConnorModType", "ConnorMod") }
            };

            //add task 1 to task list
            List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>> {
               new Tuple<string, MetaMorpheusTask>("task1", task1)};

            ModificationMotif motif;
            ModificationMotif.TryGetMotif("P", out motif);

            var connorMod = new ModificationWithMass("ConnorMod", null, motif, ModificationSites.Any, 10, null, null, null, "ConnorModType");

            GlobalTaskLevelSettings.AddMods(new List<ModificationWithLocation>
            {
                connorMod
            });

            #endregion setup

            #region Protein and Mod Creation

            //create modification lists
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            //add modification to Protein object
            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationWithMass modToAdd = connorMod;
            ModificationWithMass modToAdd2 = connorMod;
            dictHere.Add(1, new List<Modification> { modToAdd });
            dictHere.Add(3, new List<Modification> { modToAdd2 });

            //protein Creation (One with mod and one without)
            Protein TestProtein = new Protein("PEPTID", "accession1", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>(), new List<SequenceVariation>());
            Protein TestProteinWithMod = new Protein("PEPTID", "accession1", new List<Tuple<string, string>>(), dictHere, new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>(), new List<SequenceVariation>());

            #endregion Protein and Mod Creation

            #region XML File

            Console.WriteLine("hi");
            //First Write XML Database

            string xmlName = "okkk.xml";

            //Add Mod to list and write XML input database
            Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>> modList = new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>();
            var Hash = new HashSet<Tuple<int, ModificationWithMass>>
            {
                Tuple.Create(3, modToAdd)
            };
            modList.Add("test", Hash);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProteinWithMod }, xmlName);

            #endregion XML File

            #region MZML File

            //now write MZML file
            Dictionary<string, Modification> ok;
            var protein = ProteinDbLoader.LoadProteinXML(xmlName, true, new List<Modification>(), false, new List<string>(), out ok);
            var digestedList = protein[0].Digest(task1.Protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass> { }).ToList();
            Assert.AreEqual(1, digestedList.Count);

            PeptideWithPossibleModifications modPep1 = digestedList[0];

            var setList1 = modPep1.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();
            Assert.AreEqual(4, setList1.Count);

            //Set Peptide with 1 mod at position 3
            PeptideWithSetModifications pepWithSetMods1 = setList1[1];

            //Finally Write MZML file
            Assert.AreEqual("PEP[ConnorModType:ConnorMod]TID", pepWithSetMods1.Sequence);
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1 });
            string mzmlName = @"hello.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            #endregion MZML File

            //run!
            EverythingRunnerEngine.finishedAllTasksEngineHandler += SuccessfullyFinishedAllTasks;
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) });
            engine.Run();

            string outputFolderInThisTest = outputFolder;
            string final = Path.Combine(outputFolder, "task1", "okkkpruned.xml");
            //string[] files = Directory.GetFiles(fileAtPath);
            //string file = fileAtPath;
            var proteins = ProteinDbLoader.LoadProteinXML(final, true, new List<Modification>(), false, new List<string>(), out ok);
            //check length
            Assert.AreEqual(proteins[0].OneBasedPossibleLocalizedModifications.Count, 1);
            //check location (key)
            Assert.AreEqual(proteins[0].OneBasedPossibleLocalizedModifications.ContainsKey(3), true);
            List<Modification> listOfMods = new List<Modification>();
            listOfMods = proteins[0].OneBasedPossibleLocalizedModifications[3];
            //check Type, count, ID
            Assert.AreEqual(listOfMods[0].modificationType, "ConnorModType");
            Assert.AreEqual(listOfMods[0].id, "ConnorMod");
            Assert.AreEqual(listOfMods.Count, 1);
        }

        [Test]
        public static void TestUniquePeptideCount()
        {
            #region setup

            SearchTask testUnique = new SearchTask()
            {
                ListOfModsLocalize = new List<Tuple<string, string>> { new Tuple<string, string>("testUniqueModType", "testUniqueMod") }
            };

            List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>> {
               new Tuple<string, MetaMorpheusTask>("TestUnique", testUnique)};

            ModificationMotif motif;
            ModificationMotif.TryGetMotif("P", out motif);

            var testUniqeMod = new ModificationWithMass("testUniqeMod", null, motif, ModificationSites.Any, 10, null, null, null, "testUniqueModType");

            GlobalTaskLevelSettings.AddMods(new List<ModificationWithLocation>
            {
                testUniqeMod
            });

            #endregion setup

            #region mod setup and protein creation

            //create modification lists
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => testUnique.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => testUnique.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            //add modification to Protein object
            var modDictionary = new Dictionary<int, List<Modification>>();
            ModificationWithMass modToAdd = testUniqeMod;
            modDictionary.Add(1, new List<Modification> { modToAdd });
            modDictionary.Add(3, new List<Modification> { modToAdd });

            //protein Creation (One with mod and one without)
            Protein TestProtein = new Protein("PEPTID", "accession1", new List<Tuple<string, string>>(), modDictionary, new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>(), new List<SequenceVariation>());

            #endregion mod setup and protein creation

            #region XML setup

            //First Write XML Database

            string xmlName = "singleProteinWithTwoMods.xml";

            //Add Mod to list and write XML input database
            Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>> modList = new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>();
            var Hash = new HashSet<Tuple<int, ModificationWithMass>>
            {
                Tuple.Create(3, modToAdd)
            };
            modList.Add("test", Hash);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProtein }, xmlName);

            #endregion XML setup

            #region MZML setup

            //now write MZML file
            Dictionary<string, Modification> ok;
            var protein = ProteinDbLoader.LoadProteinXML(xmlName, true, new List<Modification>(), false, new List<string>(), out ok);
            var digestedList = protein[0].Digest(testUnique.Protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass> { }).ToList();
            Assert.AreEqual(1, digestedList.Count);

            PeptideWithPossibleModifications modPep1 = digestedList[0];

            var setList1 = modPep1.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();
            Assert.AreEqual(4, setList1.Count);

            //Finally Write MZML file
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { setList1[0], setList1[1], setList1[2], setList1[3], setList1[0], setList1[1] });
            string mzmlName = @"singleProteinWithRepeatedMods.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            #endregion MZML setup

            #region run

            EverythingRunnerEngine.finishedAllTasksEngineHandler += SuccessfullyFinishedAllTasks;
            string outputFolderInThisTest = outputFolder;
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) });
            engine.Run();

            List<string> found = new List<string>();
            string line;

            using (StreamReader file = new StreamReader(Path.Combine(outputFolder, "TestUnique", "results.txt")))
            {
                while ((line = file.ReadLine()) != null)
                {
                    if (line.Contains("Unique PSMS within 1% FDR"))
                    {
                        Console.WriteLine("adsf");
                        found.Add(line);
                        Assert.AreEqual(found[0], "Unique PSMS within 1% FDR: 4");
                    }
                }
            }

            #endregion run
        }

        #endregion Public Methods

        #region Private Methods

        private static void SuccessfullyFinishedAllTasks(object sender, string rootOutputFolderPath)
        {
            outputFolder = rootOutputFolderPath;
        }

        #endregion Private Methods

    }
}