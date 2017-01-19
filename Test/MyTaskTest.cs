using InternalLogicEngineLayer;
using NUnit.Framework;
using System.Text;
using InternalLogicTaskLayer;
using System.Collections.Generic;
using System;
using MassSpectrometry;
using Spectra;
using System.Collections.ObjectModel;
using System.Collections;
using Proteomics;
using OldInternalLogic;

namespace Test
{
    [TestFixture]
    public class MyTaskTest
    {

        [Test]
        public static void TestEverythingRunner()
        {

            Console.WriteLine("Environment.CurrentDirectory is " + Environment.CurrentDirectory);
			MyEngine.OutLabelStatusHandler += MyEngine_outLabelStatusHandler;

			string mzmlName = "ok.mzML";
			Dictionary<int, List<MorpheusModification>> oneBasedPossibleLocalizedModifications = new Dictionary<int, List<MorpheusModification>>();
			Protein ParentProtein = new Protein("MQQQQQQQ", "accession", null, oneBasedPossibleLocalizedModifications,new int[0], new int[0], new string[0], null, null, 0, false);
			PeptideWithPossibleModifications modPep = new PeptideWithPossibleModifications(1, 8, ParentProtein, 0, "kk");
			Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
			PeptideWithSetModifications pepWithSetMods = new PeptideWithSetModifications(modPep, twoBasedVariableAndLocalizeableModificationss);

			IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = new TestDataFile(pepWithSetMods);
			IO.MzML.MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(myMsDataFile, mzmlName);

			string xmlName = "ok.xml";
			var ye = new Dictionary<string, HashSet<Tuple<int, string>>>();
			GPTMDTask.WriteGPTMDdatabse(ye, new List<Protein> { ParentProtein }, xmlName);

			ModList modlist1 = new ModList("f.txt");
			ModList modlist2 = new ModList("v.txt");
			ModList modlist3 = new ModList("ptmlist.txt");
			ObservableCollection<ModList> modList = new ObservableCollection<ModList> { modlist1, modlist2, modlist3 };
			CalibrationTask task1 = new CalibrationTask(modList);

			List<MyTaskEngine> taskList = new List<MyTaskEngine> { task1 };
			List<string> startingRawFilenameList = new List<string> { mzmlName };
			List<string> startingXmlDbFilenameList = new List<string> { xmlName };
			var engine =  new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList);

			var results = (EverythingRunnerResults)engine.Run();

			Assert.NotNull(results);
        }

		static void MyEngine_outLabelStatusHandler(object sender, StringEventArgs e)
		{
			Console.WriteLine(e.s);
		}
}
}