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
            task1.InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain;
            GptmdTask task2 = new GptmdTask(new ObservableCollection<ModList> { modlist1, modlist2, modlist3, modlist4 });

            IEnumerable<SearchMode> allSms = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("ye", 5) };

            SearchTask task3 = new SearchTask(new ObservableCollection<ModList> { modlist1, modlist2, modlist3, modlist4 }, allSms);
            task3.ListOfModListsForSearch[3].Localize = true;

            SearchTask task4 = new SearchTask(new ObservableCollection<ModList> { modlist1, modlist2, modlist3, modlist4 }, allSms);
            task4.ListOfModListsForSearch[3].Localize = true;
            task4.ClassicSearch = false;
            List<MyTaskEngine> taskList = new List<MyTaskEngine> { task1, task2, task3, task4 };

            #endregion Setup tasks

            // Generate data for files
            Protein ParentProtein = new Protein("MAAAAAYYYYY", "accession", new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            PeptideWithPossibleModifications modPep = ParentProtein.Digest(task1.Protease, task1.MaxMissedCleavages, task1.InitiatorMethionineBehavior).First();
            Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            PeptideWithSetModifications pepWithSetMods = new PeptideWithSetModifications(modPep, twoBasedVariableAndLocalizeableModificationss);
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = new TestDataFile(pepWithSetMods);

            #region Write the files

            string mzmlName = @"ok.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(myMsDataFile, mzmlName);
            string xmlName = "ok.xml";
            GptmdTask.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, string>>>(), new List<Protein> { ParentProtein }, xmlName);

            #endregion Write the files

            // RUN!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<XmlForTask> { new XmlForTask(xmlName, false) });
            var results = (EverythingRunnerResults)engine.Run();

            Assert.NotNull(results);
        }

        #endregion Public Methods

    }
}