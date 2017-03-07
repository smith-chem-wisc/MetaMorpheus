using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

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
            {
                var readMods = UsefulProteomicsDatabases.PtmListLoader.ReadMods(modFile).ToList();
                MetaMorpheusTask.AddModList(new ModList(modFile, readMods));
            }

            CalibrationTask task1 = new CalibrationTask();
            GptmdTask task2 = new GptmdTask();

            SearchTask task3 = new SearchTask();
            task3.ListOfModListsLocalize.Add(MetaMorpheusTask.AllModLists.First(b => b.FileName.EndsWith("m.txt")));
            task3.ListOfModListsLocalize.Add(MetaMorpheusTask.AllModLists.First(b => b.FileName.EndsWith("glyco.txt")));
            task3.DoParsimony = true;

            SearchTask task4 = new SearchTask();
            task4.ListOfModListsLocalize.Add(MetaMorpheusTask.AllModLists.First(b => b.FileName.EndsWith("m.txt")));
            task4.ListOfModListsLocalize.Add(MetaMorpheusTask.AllModLists.First(b => b.FileName.EndsWith("glyco.txt")));
            task4.ClassicSearch = false;
            List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>> {
                new Tuple<string, MetaMorpheusTask>("task1", task1),
                new Tuple<string, MetaMorpheusTask>("task2", task2),
                new Tuple<string, MetaMorpheusTask>("task3", task3),
                new Tuple<string, MetaMorpheusTask>("task4", task4),};

            #endregion Setup tasks

            List<ModificationWithMass> variableModifications = task1.ListOfModListsVariable.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> fixedModifications = task1.ListOfModListsFixed.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            //List<MorpheusModification> localizeableModifications = task1.ListOfModListsForCalibration.Where(b => b.Localize).SelectMany(b => b.Mods).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1", new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, false, false, null);

            var digestedList = ParentProtein.Digest(task1.Protease, 0, InitiatorMethionineBehavior.Retain, fixedModifications).ToList();

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
            dictHere.Add(3, new List<Modification> { new ModificationWithMass("21", null, motif, ModificationSites.Any, 21.981943, null, 0, new List<double> { 21.981943 }, null, "") });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", dictHere, new int?[0], new int?[0], new string[0], null, null, false, false, null);
            digestedList = ParentProteinToNotInclude.Digest(task1.Protease, 0, InitiatorMethionineBehavior.Retain, fixedModifications).ToList();
            var modPep3 = digestedList[0];
            Assert.AreEqual(1, digestedList.Count);
            var setList3 = modPep3.GetPeptidesWithSetModifications(variableModifications, 4096, 3).ToList();
            Assert.AreEqual(4, setList3.Count);
            Console.WriteLine(string.Join(",", setList3.Select(b => b.Sequence)));

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, setList3[1] });

            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", new Dictionary<int, List<Modification>>(), new int?[] { 4 }, new int?[] { 8 }, new string[] { "chain" }, null, null, false, false, null);

            #region Write the files

            string mzmlName = @"ok.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName);
            string xmlName = "okk.xml";
            GptmdTask.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), new List<Protein> { ParentProtein, proteinWithChain }, xmlName);

            #endregion Write the files

            // RUN!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) });
            var results = (EverythingRunnerResults)engine.Run();

            Assert.NotNull(results);
        }

        #endregion Public Methods

    }
}