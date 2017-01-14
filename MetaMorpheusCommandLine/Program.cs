using InternalLogicEngineLayer;
using InternalLogicTaskLayer;
using IO.MzML;
using MassSpectrometry;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;

namespace MetaMorpheusCommandLine
{
	class Program
	{
		#region Private Fields

		const string elementsLocation = @"elements.dat";
		const string unimodLocation = @"unimod_tables.xml";
		const string uniprotLocation = @"ptmlist.txt";
		static bool inProgress;

		#endregion Private Fields

		#region Private Methods

		static void Main(string[] args)
		{
			var version = Assembly.GetExecutingAssembly().GetName().Version.ToString();
			if (version.Equals("1.0.0.0"))
				Console.WriteLine("Not a release version");
			else
				Console.WriteLine(version);

			if (args==null || args.Length==0)
			{
				Console.WriteLine("Usage:");
				return;
			}

			UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
			MyEngine.unimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation);
			MyEngine.uniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation);

			MyEngine.finishedSingleEngineHandler += MyEngine_finishedSingleEngineHandler;
			MyEngine.outLabelStatusHandler += MyEngine_outLabelStatusHandler;
			MyEngine.outProgressHandler += MyEngine_outProgressHandler;
			MyEngine.startingSingleEngineHander += MyEngine_startingSingleEngineHander;

			MyTaskEngine.finishedSingleTaskHandler += MyTaskEngine_finishedSingleTaskHandler;
			MyTaskEngine.finishedWritingFileHandler += MyTaskEngine_finishedWritingFileHandler;
			MyTaskEngine.startingSingleTaskHander += MyTaskEngine_startingSingleTaskHander;
		}

		static void RunModernSearchEngine()
		{
			List<CompactPeptide> peptideIndex;
			using (var file = File.OpenRead(Path.Combine(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\2017-01-13-09-57-20", "peptideIndex.ind")))
				peptideIndex = (List<CompactPeptide>)new NetSerializer.Serializer(SearchTask.GetSubclassesAndItself(typeof(List<CompactPeptide>))).Deserialize(file);
			Dictionary<float, List<int>> fragmentIndexDict;
			using (var file = File.OpenRead(Path.Combine(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\2017-01-13-09-57-20", "fragmentIndex.ind")))
				fragmentIndexDict = (Dictionary<float, List<int>>)new NetSerializer.Serializer(SearchTask.GetSubclassesAndItself(typeof(Dictionary<float, List<int>>))).Deserialize(file);
			var keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
			var fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();

			IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = new Mzml(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML", 400);
			myMsDataFile.Open();
			int spectraFileIndex = 0;
			double fragmentToleranceInDaltons = 0.01;
			var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("", 5) };

			var s = new ModernSearchEngine(myMsDataFile, spectraFileIndex, peptideIndex, keys, fragmentIndex, fragmentToleranceInDaltons, searchModes);
			s.Run();
		}

		static void RunSearchTask()
		{
			IEnumerable<ModList> modList = new List<ModList> { new ModList("f.txt"), new ModList("v.txt"), new ModList("p.txt") };
			IEnumerable<SearchMode> ism = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("", 5) };
			var s = new SearchTask(modList, ism);

			s.classicSearch = false;
			s.rawDataFilenameList = new List<string> { @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML" };
			s.xmlDbFilenameList = new List<string> { @"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-mouse-reviewed-12-23-2016.xml" };
			s.output_folder = Path.GetTempPath();

			s.Run();
		}

		static void MyTaskEngine_startingSingleTaskHander(object sender, SingleTaskEventArgs e)
		{
			if (inProgress)
				Console.WriteLine();
			inProgress = false;
			Console.WriteLine("Starting task:");
			Console.WriteLine(e.theTask);
		}

		static void MyTaskEngine_finishedWritingFileHandler(object sender, SingleFileEventArgs e)
		{
			if (inProgress)
				Console.WriteLine();
			inProgress = false;
			Console.WriteLine("Finished writing file: " + e.writtenFile);
		}

		static void MyTaskEngine_finishedSingleTaskHandler(object sender, SingleTaskEventArgs e)
		{
			if (inProgress)
				Console.WriteLine();
			inProgress = false;
			Console.WriteLine("Finished task: " + e.theTask.GetType().Name);
		}

		static void MyEngine_startingSingleEngineHander(object sender, SingleEngineEventArgs e)
		{
			if (inProgress)
				Console.WriteLine();
			inProgress = false;
			Console.WriteLine("Starting engine:" + e.myEngine.GetType().Name);
		}

		static void MyEngine_outProgressHandler(object sender, ProgressEventArgs e)
		{
			Console.Write(e.new_progress + " ");
			inProgress = true;
		}

		static void MyEngine_outLabelStatusHandler(object sender, string e)
		{
			if (inProgress)
				Console.WriteLine();
			inProgress = false;
			Console.WriteLine("Status: " + e);
		}

		static void MyEngine_finishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
		{
			if (inProgress)
				Console.WriteLine();
			inProgress = false;
			Console.WriteLine("Finished engine: ");
			Console.WriteLine(e);
		}

		#endregion Private Methods
	}
}
