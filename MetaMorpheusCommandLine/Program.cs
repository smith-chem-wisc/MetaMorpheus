using InternalLogicEngineLayer;
using InternalLogicTaskLayer;
using IO.MzML;
using MassSpectrometry;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace MetaMorpheusCommandLine
{
    internal class Program
    {

        #region Private Fields

        private const string elementsLocation = @"elements.dat";
        private const string unimodLocation = @"unimod_tables.xml";
        private const string uniprotLocation = @"ptmlist.txt";
        private static bool inProgress;

        #endregion Private Fields

        #region Private Methods

        private static void Main(string[] args)
        {
            if (MyEngine.MetaMorpheusVersion.Equals("1.0.0.0"))
                Console.WriteLine("Not a release version");
            else
                Console.WriteLine(MyEngine.MetaMorpheusVersion);

            if (args == null || args.Length == 0)
            {
                Console.WriteLine("Usage:");
                Console.WriteLine("\tmodern: runs the modern search engine");
                Console.WriteLine("\tsearch: runs the search task");
                return;
            }

            MyEngine.FinishedSingleEngineHandler += MyEngine_finishedSingleEngineHandler;
            MyEngine.OutLabelStatusHandler += MyEngine_outLabelStatusHandler;
            MyEngine.OutProgressHandler += MyEngine_outProgressHandler;
            MyEngine.StartingSingleEngineHander += MyEngine_startingSingleEngineHander;

            MyTaskEngine.FinishedSingleTaskHandler += MyTaskEngine_finishedSingleTaskHandler;
            MyTaskEngine.FinishedWritingFileHandler += MyTaskEngine_finishedWritingFileHandler;
            MyTaskEngine.StartingSingleTaskHander += MyTaskEngine_startingSingleTaskHander;

            switch (args[0])
            {
                case "modern":
                    RunModernSearchEngine();
                    break;

                case "search":
                    RunSearchTask();
                    break;
            }
        }

        private static void RunModernSearchEngine()
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
            double fragmentToleranceInDaltons = 0.01;
            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode(5) };

            var s = new ModernSearchEngine(myMsDataFile, peptideIndex, keys, fragmentIndex, fragmentToleranceInDaltons, searchModes);
            s.Run();
        }

        private static void RunSearchTask()
        {
            IEnumerable<ModList> modList = new List<ModList> { new ModList("f.txt"), new ModList("v.txt"), new ModList("p.txt") };
            IEnumerable<SearchMode> ism = new List<SearchMode> { new SinglePpmAroundZeroSearchMode(5) };
            var s = new SearchTask(modList, ism);

            s.ClassicSearch = false;
            s.rawDataFilenameList = new List<string> { @"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML" };
            s.dbFilenameList = new List<DbForTask> { new DbForTask(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-mouse-reviewed-12-23-2016.xml", false) };
            s.OutputFolder = Path.GetTempPath();

            s.Run();
        }

        private static void MyTaskEngine_startingSingleTaskHander(object sender, SingleTaskEventArgs e)
        {
            if (inProgress)
                Console.WriteLine();
            inProgress = false;
            Console.WriteLine("Starting task:");
            Console.WriteLine(e.TheTask);
        }

        private static void MyTaskEngine_finishedWritingFileHandler(object sender, SingleFileEventArgs e)
        {
            if (inProgress)
                Console.WriteLine();
            inProgress = false;
            Console.WriteLine("Finished writing file: " + e.writtenFile);
        }

        private static void MyTaskEngine_finishedSingleTaskHandler(object sender, SingleTaskEventArgs e)
        {
            if (inProgress)
                Console.WriteLine();
            inProgress = false;
            Console.WriteLine("Finished task: " + e.TheTask.GetType().Name);
        }

        private static void MyEngine_startingSingleEngineHander(object sender, SingleEngineEventArgs e)
        {
            if (inProgress)
                Console.WriteLine();
            inProgress = false;
            Console.WriteLine("Starting engine:" + e.myEngine.GetType().Name);
        }

        private static void MyEngine_outProgressHandler(object sender, ProgressEventArgs e)
        {
            Console.Write(e.new_progress + " ");
            inProgress = true;
        }

        private static void MyEngine_outLabelStatusHandler(object sender, StringEventArgs e)
        {
            if (inProgress)
                Console.WriteLine();
            inProgress = false;
            Console.WriteLine("Status: " + e.s);
        }

        private static void MyEngine_finishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
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