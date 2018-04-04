﻿using EngineLayer;
using Fclp;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using TaskLayer;

namespace MetaMorpheusCommandLine
{
    internal static class Program
    {
        #region Private Fields

        private static bool inProgress;

        private static System.CodeDom.Compiler.IndentedTextWriter myWriter = new System.CodeDom.Compiler.IndentedTextWriter(Console.Out, "\t");

        #endregion Private Fields

        #region Private Methods

        private static void Main(string[] args)
        {
            Console.WriteLine("Welcome to MetaMorpheus");
            Console.WriteLine(GlobalVariables.MetaMorpheusVersion);

            var p = new FluentCommandLineParser<ApplicationArguments>();

            p.Setup(arg => arg.Tasks)
             .As('t', "tasks")
             .SetDefault(new List<string>());

            p.Setup(arg => arg.OutputFolder)
             .As('o', "outputFolder")
             .SetDefault(null);

            p.Setup(arg => arg.MetaTasks)
             .As('m', "meta-task")
             .SetDefault(new List<string>());

            p.Setup(arg => arg.Spectra)
             .As('s', "spectra")
             .Required();

            p.Setup(arg => arg.Databases)
             .As('d', "databases")
             .Required();

            var result = p.Parse(args);

            if (p.Object.MetaTasks.Count != 0 || p.Object.Tasks.Count != 0)
            {
                if (!result.HasErrors)
                {
                    MetaMorpheusEngine.WarnHandler += WarnHandler;
                    MetaMorpheusEngine.OutProgressHandler += MyEngine_outProgressHandler;
                    MetaMorpheusEngine.StartingSingleEngineHander += MyEngine_startingSingleEngineHander;
                    MetaMorpheusEngine.FinishedSingleEngineHandler += MyEngine_finishedSingleEngineHandler;

                    MetaMorpheusTask.WarnHandler += WarnHandler;
                    MetaMorpheusTask.LogHandler += LogHandler;
                    MetaMorpheusTask.StartingSingleTaskHander += MyTaskEngine_startingSingleTaskHander;
                    MetaMorpheusTask.FinishedSingleTaskHandler += MyTaskEngine_finishedSingleTaskHandler;
                    MetaMorpheusTask.FinishedWritingFileHandler += MyTaskEngine_finishedWritingFileHandler;

                    foreach (var db in p.Object.Databases)
                        if (!Path.GetExtension(db).Equals(".fasta"))
                            GlobalVariables.AddMods(UsefulProteomicsDatabases.ProteinDbLoader.GetPtmListFromProteinXml(db).OfType<ModificationWithLocation>());

                    List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>();

                    for (int i = 0; i < p.Object.Tasks.Count; i++)
                    {
                        var filePath = p.Object.Tasks[i];

                        var uhum = Toml.ReadFile(filePath, MetaMorpheusTask.tomlConfig);

                        switch (uhum.Get<string>("TaskType"))
                        {
                            case "Search":
                                var ye1 = Toml.ReadFile<SearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                                taskList.Add(("Task" + (i + 1) + "SearchTask", ye1));
                                break;

                            case "Calibrate":
                                var ye2 = Toml.ReadFile<CalibrationTask>(filePath, MetaMorpheusTask.tomlConfig);
                                taskList.Add(("Task" + (i + 1) + "CalibrationTask", ye2));
                                break;

                            case "Gptmd":
                                var ye3 = Toml.ReadFile<GptmdTask>(filePath, MetaMorpheusTask.tomlConfig);
                                taskList.Add(("Task" + (i + 1) + "GptmdTask", ye3));
                                break;

                            case "XLSearch":
                                var ye4 = Toml.ReadFile<XLSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                                taskList.Add(("Task" + (i + 1) + "XLSearchTask", ye4));
                                break;

                            case "Neo":
                                Console.WriteLine("Neo tasks are meta-tasks that rely on several other tasks. Please use -m for meta instead of -t. Skipping.");
                                break;

                            default:
                                Console.WriteLine(uhum.Get<string>("TaskType") + " is not a known task type! Skipping.");
                                break;
                        }
                    }

                    for (int i = 0; i < p.Object.MetaTasks.Count; i++)
                    {
                        var filePath = p.Object.MetaTasks[i];
                        var uhum = Toml.ReadFile(filePath, MetaMorpheusTask.tomlConfig);
                        switch (uhum.Get<string>("TaskType"))
                        {
                            case "Search":
                                Console.WriteLine("Search tasks are individual tasks. Please use -t for task instead of -m. Skipping.");
                                break;

                            case "Calibrate":
                                Console.WriteLine("Calibrate tasks are individual tasks. Please use -t for task instead of -m. Skipping.");
                                break;

                            case "Gptmd":
                                Console.WriteLine("Gptmd tasks are individual tasks. Please use -t for task instead of -m. Skipping.");
                                break;

                            case "XLSearch":
                                Console.WriteLine("XLSearch tasks are individual tasks. Please use -t for task instead of -m. Skipping.");
                                break;

                            case "Neo":
                                var ye5 = Toml.ReadFile<NeoSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                                foreach (MetaMorpheusTask task in NeoLoadTomls.LoadTomls(ye5))
                                    taskList.Add(("Task" + (taskList.Count + 1) + ye5.TaskType, ye5));
                                break;

                            default:
                                Console.WriteLine(uhum.Get<string>("TaskType") + " is not a known task type! Skipping.");
                                break;
                        }
                    }

                    List<string> startingRawFilenameList = p.Object.Spectra.Select(b => Path.GetFullPath(b)).ToList();
                    List<DbForTask> startingXmlDbFilenameList = p.Object.Databases.Select(b => new DbForTask(Path.GetFullPath(b), IsContaminant(b))).ToList();
                    
                    string outputFolder = p.Object.OutputFolder;
                    if(outputFolder == null)
                    {
                        var pathOfFirstSpectraFile = Path.GetDirectoryName(startingRawFilenameList.First());
                        outputFolder = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
                    }
                    
                    EverythingRunnerEngine a = new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList, outputFolder);

                    try
                    {
                        a.Run();
                    }
                    catch (Exception e)
                    {
                        while (e.InnerException != null) e = e.InnerException;
                        var message = "Run failed, Exception: " + e.Message;
                        Console.WriteLine(message);
                    }
                }
                else
                {
                    Console.WriteLine("Error Text:" + result.ErrorText);
                }
            }
            else
            {
                Console.WriteLine("Error Text: No toml file was specified. Use -t for tasks or -m for meta-tasks.");
            }
        }

        private static void WriteMultiLineIndented(string toWrite)
        {
            string[] tokens = Regex.Split(toWrite, @"\r?\n|\r");
            foreach (var str in tokens)
            {
                myWriter.WriteLine(str);
            }
        }

        private static bool IsContaminant(string b)
        {
            if (b.ToUpper().Contains("contaminant".ToUpper())
                || b.ToUpper().Contains("CRAP"))
                return true;
            return false;
        }

        private static void MyTaskEngine_startingSingleTaskHander(object sender, SingleTaskEventArgs e)
        {
            if (inProgress)
                myWriter.WriteLine();
            inProgress = false;
            WriteMultiLineIndented("Starting task: " + e.DisplayName);
            myWriter.Indent++;
        }

        private static void MyTaskEngine_finishedWritingFileHandler(object sender, SingleFileEventArgs e)
        {
            if (inProgress)
                myWriter.WriteLine();
            inProgress = false;
            WriteMultiLineIndented("Finished writing file: " + e.writtenFile);
        }

        private static void MyTaskEngine_finishedSingleTaskHandler(object sender, SingleTaskEventArgs e)
        {
            if (inProgress)
                myWriter.WriteLine();
            inProgress = false;
            myWriter.Indent--;
            WriteMultiLineIndented("Finished task: " + e.DisplayName);
        }

        private static void MyEngine_startingSingleEngineHander(object sender, SingleEngineEventArgs e)
        {
            if (inProgress)
                myWriter.WriteLine();
            inProgress = false;
            WriteMultiLineIndented("Starting engine: " + e.myEngine.GetType().Name + " " + e.myEngine.GetId());
            myWriter.Indent++;
        }

        private static void MyEngine_finishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
        {
            if (inProgress)
                myWriter.WriteLine();
            inProgress = false;
            WriteMultiLineIndented("Engine results: " + e);
            myWriter.Indent--;
            WriteMultiLineIndented("Finished engine: " + e.myResults.MyEngine.GetType().Name + " " + e.myResults.MyEngine.GetId());
        }

        private static void MyEngine_outProgressHandler(object sender, ProgressEventArgs e)
        {
            myWriter.Write(e.new_progress + " ");
            inProgress = true;
        }

        private static void WarnHandler(object sender, StringEventArgs e)
        {
            if (inProgress)
                myWriter.WriteLine();
            inProgress = false;
            WriteMultiLineIndented("WARN: " + e.S);
        }

        private static void LogHandler(object sender, StringEventArgs e)
        {
            if (inProgress)
                myWriter.WriteLine();
            inProgress = false;
            WriteMultiLineIndented("Log: " + e.S);
        }

        #endregion Private Methods

        #region Public Classes

        public class ApplicationArguments
        {
            #region Public Properties

            public List<string> Tasks { get; set; }
            public List<string> Databases { get; set; }
            public List<string> Spectra { get; set; }
            public List<string> MetaTasks { get; set; }
            public string OutputFolder { get; set; }

            #endregion Public Properties
        }

        #endregion Public Classes
    }
}