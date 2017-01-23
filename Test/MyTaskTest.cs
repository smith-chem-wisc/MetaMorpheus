using InternalLogicEngineLayer;
using InternalLogicTaskLayer;
using MassSpectrometry;
using NUnit.Framework;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Chemistry;

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

			ModList modlist1 = new ModList("f.txt");
			ModList modlist2 = new ModList("v.txt");
			ModList modlist3 = new ModList("ptmlist.txt");
			ModList modlist4 = new ModList("m.txt");
			ObservableCollection<ModList> modList = new ObservableCollection<ModList> { modlist1, modlist2, modlist3 };
			CalibrationTask task1 = new CalibrationTask(modList);
			GptmdTask task2 = new GptmdTask(new ObservableCollection<ModList> { modlist1, modlist2, modlist3, modlist4 });

			IEnumerable<SearchMode> allSms = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("ye", 5) };

			SearchTask task3 = new SearchTask(new ObservableCollection<ModList> { modlist1, modlist2, modlist3, modlist4 }, allSms);
			task3.ListOfModListsForSearch[3].Localize = true;

			SearchTask task4 = new SearchTask(new ObservableCollection<ModList> { modlist1, modlist2, modlist3, modlist4 }, allSms);
			task4.ListOfModListsForSearch[3].Localize = true;
			task4.ClassicSearch = false;
			List<MyTaskEngine> taskList = new List<MyTaskEngine> { task1, task2, task3, task4 };

			#endregion Setup tasks

			List<MorpheusModification> variableModifications = task1.ListOfModListsForCalibration.Where(b => b.Variable).SelectMany(b => b.Mods).ToList();
			//List<MorpheusModification> fixedModifications = task1.ListOfModListsForCalibration.Where(b => b.Fixed).SelectMany(b => b.Mods).ToList();
			//List<MorpheusModification> localizeableModifications = task1.ListOfModListsForCalibration.Where(b => b.Localize).SelectMany(b => b.Mods).ToList();

			// Generate data for files
			Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1", new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);

			var digestedList = ParentProtein.Digest(task1.Protease, 0, InitiatorMethionineBehavior.Retain).ToList();

			Assert.AreEqual(2, digestedList.Count);

			PeptideWithPossibleModifications modPep1 = digestedList[0];
			var setList1 = modPep1.GetPeptideWithSetModifications(variableModifications, 4096, 3).ToList();

			Assert.AreEqual(2, setList1.Count);

			PeptideWithSetModifications pepWithSetMods1 = setList1[0];

			PeptideWithPossibleModifications modPep2 = digestedList[1];
			var setList2 = modPep2.GetPeptideWithSetModifications(variableModifications, 4096, 3).ToList();

			Assert.AreEqual(1, setList2.Count);

			PeptideWithSetModifications pepWithSetMods2 = setList2[0];


			var dictHere = new Dictionary<int, List<MorpheusModification>>();
			dictHere.Add(3, new List<MorpheusModification> { new MorpheusModification(null, ModificationType.AminoAcidResidue, 'E', 21.981943, null, null, '\0', double.NaN, false, null)});
			Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", dictHere, new int[0], new int[0], new string[0], null, null, 0, false, false);
			 digestedList = ParentProteinToNotInclude.Digest(task1.Protease, 0, InitiatorMethionineBehavior.Retain).ToList();
			var modPep3 = digestedList[0];
			Assert.AreEqual(1, digestedList.Count);
			var setList3 = modPep3.GetPeptideWithSetModifications(variableModifications, 4096, 3).ToList();
			Assert.AreEqual(4, setList3.Count);
			Console.WriteLine(string.Join(",", setList3.Select(b => b.Sequence)));

			IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, setList3[1] } );

			Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", new Dictionary<int, List<MorpheusModification>>(), new int[] { 4 }, new int[] { 8 }, new string[] { "chain" }, null, null, 0, false, false);

            #region Write the files

            string mzmlName = @"ok.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(myMsDataFile, mzmlName);
            string xmlName = "ok.xml";
            GptmdTask.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, string>>>(), new List<Protein> { ParentProtein, proteinWithChain }, xmlName);

            #endregion Write the files

            // RUN!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<XmlForTask> { new XmlForTask(xmlName, false) });
            var results = (EverythingRunnerResults)engine.Run();

            Assert.NotNull(results);
        }

        #endregion Public Methods

    }
}