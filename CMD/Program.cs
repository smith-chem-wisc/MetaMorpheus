using EngineLayer;
using Fclp;
using Nett;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace MetaMorpheusCommandLine
{
    internal static class Program
    {

        #region Private Fields

        private static bool inProgress;

        #endregion Private Fields

        #region Private Methods

        private static void Main(string[] args)
        {
            if (MyEngine.MetaMorpheusVersion.Equals("1.0.0.0"))
                Console.WriteLine("Not a release version");
            else
                Console.WriteLine(MyEngine.MetaMorpheusVersion);

            var p = new FluentCommandLineParser<ApplicationArguments>();

            p.Setup(arg => arg.Tasks)
             .As('t', "tasks")
             .Required();

            p.Setup(arg => arg.Spectra)
             .As('s', "spectra")
             .Required();

            p.Setup(arg => arg.Databases)
             .As('d', "databases")
             .Required();

            var result = p.Parse(args);

            if (result.HasErrors == false)
            {
                MyEngine.FinishedSingleEngineHandler += MyEngine_finishedSingleEngineHandler;
                MyEngine.OutLabelStatusHandler += MyEngine_outLabelStatusHandler;
                MyEngine.OutProgressHandler += MyEngine_outProgressHandler;
                MyEngine.StartingSingleEngineHander += MyEngine_startingSingleEngineHander;

                MetaMorpheusTask.FinishedSingleTaskHandler += MyTaskEngine_finishedSingleTaskHandler;
                MetaMorpheusTask.FinishedWritingFileHandler += MyTaskEngine_finishedWritingFileHandler;
                MetaMorpheusTask.StartingSingleTaskHander += MyTaskEngine_startingSingleTaskHander;

                foreach (var modFile in Directory.GetFiles(@"Mods"))
                {
                    MetaMorpheusTask.AddModList(modFile);
                }

                List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>>();

                for (int i = 0; i < p.Object.Tasks.Count; i++)
                {
                    var draggedFilePath = p.Object.Tasks[i];

                    var uhum = Toml.ReadFile(draggedFilePath, MetaMorpheusTask.tomlConfig);

                    switch (uhum.Get<string>("TaskType"))
                    {
                        case "Search":
                            var ye1 = Toml.ReadFile<SearchTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                            taskList.Add(new Tuple<string, MetaMorpheusTask>("Task" + (i + 1) + "SearchTask", ye1));
                            break;

                        case "Calibrate":
                            var ye2 = Toml.ReadFile<CalibrationTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                            taskList.Add(new Tuple<string, MetaMorpheusTask>("Task" + (i + 1) + "CalibrationTask", ye2));
                            break;

                        case "Gptmd":
                            var ye3 = Toml.ReadFile<GptmdTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                            taskList.Add(new Tuple<string, MetaMorpheusTask>("Task" + (i + 1) + "GptmdTask", ye3));
                            break;
                    }
                }
                List<string> startingRawFilenameList = p.Object.Spectra;
                List<DbForTask> startingXmlDbFilenameList = p.Object.Databases.Select(b => new DbForTask(b, IsContaminant(b))).ToList();
                EverythingRunnerEngine a = new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList);
                a.Run();
            }
            Console.WriteLine("Usage:");
            Console.WriteLine("\t-t --tasks     List of task poml files");
            Console.WriteLine("\t-s --spectra   List of spectra files");
            Console.WriteLine("\t-e --databases List of database files");
            return;
        }

        private static bool IsContaminant(string b)
        {
            if (b.ToUpper().Contains("contaminant".ToUpper())
                || b.Contains("cRAP"))
                return true;
            return false;
        }

        private static void MyTaskEngine_startingSingleTaskHander(object sender, SingleTaskEventArgs e)
        {
            if (inProgress)
                Console.WriteLine();
            inProgress = false;
            Console.WriteLine("Starting task:");
            Console.WriteLine(e.TaskId);
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
            Console.WriteLine("Finished task: " + e.TaskId.GetType().Name);
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

        #region Public Classes

        public class ApplicationArguments
        {

            #region Public Properties

            public List<string> Tasks { get; set; }
            public List<string> Databases { get; set; }
            public List<string> Spectra { get; set; }

            #endregion Public Properties

        }

        #endregion Public Classes

    }
}