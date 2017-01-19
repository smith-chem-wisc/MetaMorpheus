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
            Console.WriteLine("Environment.CurrentDirectory is " + Environment.CurrentDirectory);

            MyEngine.OutLabelStatusHandler += MyEngine_outLabelStatusHandler;
            MyEngine.FinishedSingleEngineHandler += MyEngine_FinishedSingleEngineHandler;

            ModList modlist1 = new ModList("f.txt");
            ModList modlist2 = new ModList("v.txt");
            ModList modlist3 = new ModList("ptmlist.txt");
            ObservableCollection<ModList> modList = new ObservableCollection<ModList> { modlist1, modlist2, modlist3 };
            CalibrationTask task1 = new CalibrationTask(modList);
            task1.InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain;

            ModList modlist4 = new ModList("m.txt");
            GPTMDTask task2 = new GPTMDTask(new ObservableCollection<ModList> { modlist1, modlist2, modlist3, modlist4 });

            IEnumerable<SearchMode> allSms = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("ye", 5) };

            SearchTask task3 = new SearchTask(new ObservableCollection<ModList> { modlist1, modlist2, modlist3, modlist4 }, allSms);
            task3.ListOfModListsForSearch[3].Localize = true;

            SearchTask task4 = new SearchTask(new ObservableCollection<ModList> { modlist1, modlist2, modlist3, modlist4 }, allSms);
            task4.ListOfModListsForSearch[3].Localize = true;
            task4.ClassicSearch = false;


            string mzmlName = @"ok.mzML";
            Dictionary<int, List<MorpheusModification>> oneBasedPossibleLocalizedModifications = new Dictionary<int, List<MorpheusModification>>();
            Protein ParentProtein = new Protein("MAAAAAYYYYY", "accession", oneBasedPossibleLocalizedModifications, new int[0], new int[0], new string[0], null, null, 0, false, false);
            PeptideWithPossibleModifications modPep = ParentProtein.Digest(task1.Protease, task1.MaxMissedCleavages, task1.InitiatorMethionineBehavior).First();
            Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            PeptideWithSetModifications pepWithSetMods = new PeptideWithSetModifications(modPep, twoBasedVariableAndLocalizeableModificationss);

            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = new TestDataFile(pepWithSetMods);
            IO.MzML.MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(myMsDataFile, mzmlName);

            string xmlName = "ok.xml";
            var ye = new Dictionary<string, HashSet<Tuple<int, string>>>();
            GPTMDTask.WriteGPTMDdatabse(ye, new List<Protein> { ParentProtein }, xmlName);

            List<MyTaskEngine> taskList = new List<MyTaskEngine> { task1, task2, task3, task4 };
            List<string> startingRawFilenameList = new List<string> { mzmlName };
            List<XmlForTask> startingXmlDbFilenameList = new List<XmlForTask> { new XmlForTask(xmlName, false) };
            var engine = new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList);

            var results = (EverythingRunnerResults)engine.Run();

            Assert.NotNull(results);
        }

        #endregion Public Methods

        #region Private Methods

        private static void MyEngine_FinishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
        {
            Console.WriteLine(e.ToString());
        }

        private static void MyEngine_outLabelStatusHandler(object sender, StringEventArgs e)
        {
            Console.WriteLine(e.s);
        }

        #endregion Private Methods

    }
}