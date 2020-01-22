using CommandLine;
using CommandLine.Text;
using EngineLayer;
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
    public static class Program
    {
        private static bool InProgress;

        private static System.CodeDom.Compiler.IndentedTextWriter MyWriter = new System.CodeDom.Compiler.IndentedTextWriter(Console.Out, "\t");

        public static void Main(string[] args)
        {
            Console.WriteLine("Welcome to MetaMorpheus");
            Console.WriteLine(GlobalVariables.MetaMorpheusVersion);

            var parser = new Parser(with => with.HelpWriter = null);
            var parserResult = parser.ParseArguments<CommandLineSettings>(args);

            parserResult
              .WithParsed<CommandLineSettings>(options => Run(options))
              .WithNotParsed(errs => DisplayHelp(parserResult, errs));
        }

        public static void DisplayHelp<T>(ParserResult<T> result, IEnumerable<Error> errs)
        {
            var helpText = HelpText.AutoBuild(result, h =>
            {
                h.AdditionalNewLineAfterOption = false;
                h.Copyright = "";
                return HelpText.DefaultParsingErrorsHandler(result, h);
            }, e => e);

            helpText.MaximumDisplayWidth = 300;

            helpText.AddPostOptionsLine("Example usage (Windows): ");
            helpText.AddPostOptionsLine("CMD.exe -d C:\\ExampleDatabase.fasta -s C:\\ExampleSpectra.mzML -t C:\\ExampleTask.toml");
            helpText.AddPostOptionsLine(Environment.NewLine);

            helpText.AddPostOptionsLine("Example usage (Linux): ");
            helpText.AddPostOptionsLine("dotnet CMD.dll -d home/mydata/ExampleDatabase.fasta -s home/mydata/ExampleSpectra.mzML -t home/mydata/ExampleTask.toml");
            helpText.AddPostOptionsLine(Environment.NewLine);

            Console.WriteLine(helpText);
        }

        private static void Run(CommandLineSettings settings)
        {
            try
            {
                settings.ValidateCommandLineSettings();
            }
            catch (Exception e)
            {
                Console.WriteLine("MetaMorpheus encountered the following error:" + Environment.NewLine + e.Message);
                return;
            }

            if (settings.GenerateDefaultTomls)
            {
                Console.WriteLine("Generating default tomls at location: " + settings.OutputFolder);
                CommandLineSettings.GenerateDefaultTaskTomls(settings.OutputFolder);
                return;
            }

            // set up microvignette
            if (settings.RunMicroVignette)
            {
                // set up the spectra file
                settings.Spectra.Clear();
                settings.Spectra.Add(Path.Combine(GlobalVariables.DataDir, @"Data", "SmallCalibratible_Yeast.mzML"));

                // set up the database
                settings.Databases.Clear();
                settings.Databases.Add(Path.Combine(GlobalVariables.DataDir, @"Data", "SmallYeast.fasta"));

                // set up the tasks (calibration, GPTMD, search)
                settings.Tasks.Clear();
                CommandLineSettings.GenerateDefaultTaskTomls(settings.OutputFolder);
                settings.Tasks.Add(Path.Combine(settings.OutputFolder, "CalibrationTask.toml"));
                settings.Tasks.Add(Path.Combine(settings.OutputFolder, "GptmdTask.toml"));
                settings.Tasks.Add(Path.Combine(settings.OutputFolder, "SearchTask.toml"));
            }

            MetaMorpheusEngine.WarnHandler += WarnHandler;
            MetaMorpheusEngine.OutProgressHandler += MyEngine_outProgressHandler;
            MetaMorpheusEngine.StartingSingleEngineHander += MyEngine_startingSingleEngineHander;
            MetaMorpheusEngine.FinishedSingleEngineHandler += MyEngine_finishedSingleEngineHandler;

            MetaMorpheusTask.WarnHandler += WarnHandler;
            MetaMorpheusTask.LogHandler += LogHandler;
            MetaMorpheusTask.StartingSingleTaskHander += MyTaskEngine_startingSingleTaskHander;
            MetaMorpheusTask.FinishedSingleTaskHandler += MyTaskEngine_finishedSingleTaskHandler;
            MetaMorpheusTask.FinishedWritingFileHandler += MyTaskEngine_finishedWritingFileHandler;

            bool containsRawFiles = settings.Spectra.Select(v => Path.GetExtension(v).ToLowerInvariant()).Any(v => v == ".raw");
            if (containsRawFiles && !GlobalVariables.GlobalSettings.UserHasAgreedToThermoRawFileReaderLicence)
            {
                // write the Thermo RawFileReader licence agreement
                Console.WriteLine(ThermoRawFileReader.ThermoRawFileReaderLicence.ThermoLicenceText);
                Console.WriteLine("\nIn order to search Thermo .raw files, you must agree to the above terms. Do you agree to the above terms? y/n\n");
                string res = Console.ReadLine().ToLowerInvariant();
                if (res == "y")
                {
                    var newGlobalSettings = new GlobalSettings
                    {
                        UserHasAgreedToThermoRawFileReaderLicence = true,
                        WriteExcelCompatibleTSVs = GlobalVariables.GlobalSettings.WriteExcelCompatibleTSVs
                    };

                    Toml.WriteFile<GlobalSettings>(newGlobalSettings, Path.Combine(GlobalVariables.DataDir, @"settings.toml"));
                    GlobalVariables.GlobalSettings = newGlobalSettings;
                }
                else
                {
                    Console.WriteLine("Thermo licence has been declined. Exiting MetaMorpheus. You can still search .mzML and .mgf files without agreeing to the Thermo licence.");
                    return;
                }
            }

            foreach (var db in settings.Databases)
            {
                if (!Path.GetExtension(db).Equals(".fasta"))
                {
                    GlobalVariables.AddMods(UsefulProteomicsDatabases.ProteinDbLoader.GetPtmListFromProteinXml(db).OfType<Modification>(), true);

                    // print any error messages reading the mods to the console
                    foreach (var error in GlobalVariables.ErrorsReadingMods)
                    {
                        Console.WriteLine(error);
                    }

                    GlobalVariables.ErrorsReadingMods.Clear();
                }
            }

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>();

            var tasks = settings.Tasks.ToList();
            for (int i = 0; i < tasks.Count; i++)
            {
                var filePath = tasks[i];

                var toml = Toml.ReadFile(filePath, MetaMorpheusTask.tomlConfig);

                switch (toml.Get<string>("TaskType"))
                {
                    case "Search":
                        var searchTask = Toml.ReadFile<SearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "SearchTask", searchTask));
                        break;

                    case "Calibrate":
                        var calibrationTask = Toml.ReadFile<CalibrationTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "CalibrationTask", calibrationTask));
                        break;

                    case "Gptmd":
                        var GptmdTask = Toml.ReadFile<GptmdTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "GptmdTask", GptmdTask));
                        break;

                    case "XLSearch":
                        var XlTask = Toml.ReadFile<XLSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "XLSearchTask", XlTask));
                        break;

                    default:
                        Console.WriteLine(toml.Get<string>("TaskType") + " is not a known task type! Skipping.");
                        break;
                }
            }

            List<string> startingRawFilenameList = settings.Spectra.Select(b => Path.GetFullPath(b)).ToList();
            List<DbForTask> startingXmlDbFilenameList = settings.Databases.Select(b => new DbForTask(Path.GetFullPath(b), IsContaminant(b))).ToList();

            EverythingRunnerEngine a = new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList, settings.OutputFolder);

            try
            {
                a.Run();
            }
            catch (Exception e)
            {
                while (e.InnerException != null)
                {
                    e = e.InnerException;
                }

                var message = "Run failed, Exception: " + e.Message;
                Console.WriteLine(message);
            }
        }

        private static void WriteMultiLineIndented(string toWrite)
        {
            string[] tokens = Regex.Split(toWrite, @"\r?\n|\r");
            foreach (var str in tokens)
            {
                MyWriter.WriteLine(str);
            }
        }

        private static bool IsContaminant(string b)
        {
            if (b.ToUpper().Contains("contaminant".ToUpper())
                || b.ToUpper().Contains("CRAP"))
            {
                return true;
            }

            return false;
        }

        private static void MyTaskEngine_startingSingleTaskHander(object sender, SingleTaskEventArgs e)
        {
            if (InProgress)
            {
                MyWriter.WriteLine();
            }

            InProgress = false;
            WriteMultiLineIndented("Starting task: " + e.DisplayName);
            MyWriter.Indent++;
        }

        private static void MyTaskEngine_finishedWritingFileHandler(object sender, SingleFileEventArgs e)
        {
            if (InProgress)
            {
                MyWriter.WriteLine();
            }

            InProgress = false;
            WriteMultiLineIndented("Finished writing file: " + e.WrittenFile);
        }

        private static void MyTaskEngine_finishedSingleTaskHandler(object sender, SingleTaskEventArgs e)
        {
            if (InProgress)
            {
                MyWriter.WriteLine();
            }

            InProgress = false;
            MyWriter.Indent--;
            WriteMultiLineIndented("Finished task: " + e.DisplayName);
        }

        private static void MyEngine_startingSingleEngineHander(object sender, SingleEngineEventArgs e)
        {
            if (InProgress)
            {
                MyWriter.WriteLine();
            }

            InProgress = false;
            WriteMultiLineIndented("Starting engine: " + e.MyEngine.GetType().Name + " " + e.MyEngine.GetId());
            MyWriter.Indent++;
        }

        private static void MyEngine_finishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
        {
            if (InProgress)
            {
                MyWriter.WriteLine();
            }

            InProgress = false;
            WriteMultiLineIndented("Engine results: " + e);
            MyWriter.Indent--;
            WriteMultiLineIndented("Finished engine: " + e.MyResults.MyEngine.GetType().Name + " " + e.MyResults.MyEngine.GetId());
        }

        private static void MyEngine_outProgressHandler(object sender, ProgressEventArgs e)
        {
            MyWriter.Write(e.NewProgress + " ");
            InProgress = true;
        }

        private static void WarnHandler(object sender, StringEventArgs e)
        {
            if (InProgress)
            {
                MyWriter.WriteLine();
            }

            InProgress = false;
            WriteMultiLineIndented("WARN: " + e.S);
        }

        private static void LogHandler(object sender, StringEventArgs e)
        {
            if (InProgress)
            {
                MyWriter.WriteLine();
            }

            InProgress = false;
            WriteMultiLineIndented("Log: " + e.S);
        }
    }
}