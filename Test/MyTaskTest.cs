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

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>(), null);

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
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", new List<Tuple<string, string>>(), dictHere, new int?[0], new int?[0], new string[0], null, null, false, false, null, null);
            digestedList = ParentProteinToNotInclude.Digest(task1.Protease, 0, null, null, InitiatorMethionineBehavior.Retain, fixedModifications).ToList();
            var modPep3 = digestedList[0];
            Assert.AreEqual(1, digestedList.Count);
            var setList3 = modPep3.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();
            Assert.AreEqual(4, setList3.Count);

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, setList3[1] });

            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[] { 4 }, new int?[] { 8 }, new string[] { "chain" }, "name2", "fullname2", false, false, new List<DatabaseReference>(), null);

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
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>(), null);

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
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", new List<Tuple<string, string>>(), dictHere, new int?[0], new int?[0], new string[0], null, null, false, false, null, null);
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

            Protein proteinWithChain1 = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[] { 4 }, new int?[] { 8 }, new string[] { "chain" }, "name2", "fullname2", false, false, new List<DatabaseReference>(), null);
            Protein proteinWithChain2 = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[] { 4 }, new int?[] { 8 }, new string[] { "chain" }, "name2", "fullname2", false, false, new List<DatabaseReference>(), null);

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
                Protein theProtein = new Protein("MPEPTIDEKANTHE", "accession1", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>(), null);
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

        #endregion Public Methods

    }
}