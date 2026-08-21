using CommandLine;
using CommandLine.Text;
using EngineLayer;
using IO.ThermoRawFileReader;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer.DatabaseLoading;
using Omics.Modifications;
using TaskLayer;

namespace MetaMorpheusCommandLine
{
    public static class Program
    {
        private static bool InProgress;
        private static CommandLineSettings CommandLineSettings;

        private static System.CodeDom.Compiler.IndentedTextWriter MyWriter = new System.CodeDom.Compiler.IndentedTextWriter(Console.Out, "\t");

        public static int Main(string[] args)
        {
            // an error code of 0 is returned if the program ran successfully.
            // otherwise, an error code of >0 is returned.
            // this makes it easier to determine via scripts when the program fails.
            int errorCode = 0;

            var parser = new Parser(with => with.HelpWriter = null);
            var parserResult = parser.ParseArguments<CommandLineSettings>(args);

            parserResult
              .WithParsed<CommandLineSettings>(options => errorCode = Run(options))
              .WithNotParsed(errs => errorCode = DisplayHelp(parserResult, errs));

            return errorCode;
        }

        public static int DisplayHelp<T>(ParserResult<T> result, IEnumerable<Error> errs)
        {
            Console.WriteLine("Welcome to MetaMorpheus");
            Console.WriteLine(GlobalVariables.MetaMorpheusVersion);

            int errorCode = 0;

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

            if (errs.Any(x => x.Tag != ErrorType.HelpRequestedError))
            {
                errorCode = 1;
            }

            return errorCode;
        }

        private static int Run(CommandLineSettings settings)
        {
            int errorCode = 0;

            if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                Console.WriteLine("Welcome to MetaMorpheus");
            }

            if (settings.CustomDataDirectory != null)
            {
                GlobalVariables.UserSpecifiedDataDir = settings.CustomDataDirectory;
            }

            GlobalVariables.SetUpGlobalVariables();

            if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                Console.WriteLine(GlobalVariables.MetaMorpheusVersion);
            }

            try
            {
                settings.ValidateCommandLineSettings();
                CommandLineSettings = settings;
            }
            catch (Exception e)
            {
                if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    Console.WriteLine("MetaMorpheus encountered the following error:" + Environment.NewLine + e.Message);
                }
                errorCode = 2;

                return errorCode;
            }

            if (settings.GenerateDefaultTomls)
            {
                if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    Console.WriteLine("Generating default tomls at location: " + settings.OutputFolder);
                }
                CommandLineSettings.GenerateDefaultTaskTomls(settings.OutputFolder);

                return errorCode;
            }

            // --acceptThermoLicence records the agreement without asking for it, for the situations that
            // cannot answer the prompt further down: a container, a scheduled cluster job, a CI runner,
            // anything with stdin closed. Those are also the situations most likely to need it, because
            // outside a Windows installation SetUpDataDirectory() leaves DataDir at the application
            // folder rather than a per-user one, so the agreement is recorded per extracted copy and
            // every new download, conda environment and container layer starts unagreed.
            if (AgreeToThermoLicence(settings, GlobalVariables.DataDir, Console.Out))
            {
                return errorCode;
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
                Console.WriteLine(ThermoRawFileReaderLicence.ThermoLicenceText);
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
                    errorCode = 3;
                    return errorCode;
                }
            }

            foreach (var db in settings.Databases)
            {
                if (Path.GetExtension(db).Equals(".xml"))
                {
                    GlobalVariables.AddMods(UsefulProteomicsDatabases.ProteinDbLoader.GetPtmListFromProteinXml(db).OfType<Modification>(), true);

                    // print any error messages reading the mods to the console
                    foreach (var error in GlobalVariables.ErrorsReadingMods)
                    {
                        if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                        {
                            Console.WriteLine(error);
                        }
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

                    case "GlycoSearch":
                        var GlycoTask = MetaMorpheusTask.ReadTaskTomlWithLowResFallback<GlycoSearchTask>(filePath);
                        taskList.Add(("Task" + (i + 1) + "GlycoSearchTask", GlycoTask));
                        break;

                    case "Average":
                        var AveragingTask = Toml.ReadFile<SpectralAveragingTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "AveragingTask", AveragingTask));
                        break;

                    default:
                        if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                        {
                            Console.WriteLine(toml.Get<string>("TaskType") + " is not a known task type! Skipping.");
                        }
                        break;
                }
            }

            List<string> startingRawFilenameList = settings.Spectra.Select(b => Path.GetFullPath(b)).ToList();
            List<DbForTask> startingXmlDbFilenameList = settings.Databases.Select(b => new DbForTask(Path.GetFullPath(b), IsContaminant(b))).ToList();

            // check that an experimental design is defined if normalization is enabled
            var searchTasks = taskList
                .Where(p => p.Item2.TaskType == MyTask.Search)
                .Select(p => (SearchTask)p.Item2);

            int designExitCode = ResolveExperimentalDesign(
                Directory.GetParent(startingRawFilenameList.First()).FullName,
                startingRawFilenameList,
                searchTasks.Any(p => p.SearchParameters.Normalize),
                settings.Verbosity == CommandLineSettings.VerbosityType.minimal
                    || settings.Verbosity == CommandLineSettings.VerbosityType.normal);

            if (designExitCode != 0)
            {
                return designExitCode;
            }

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

                if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    Console.WriteLine(message);
                }
                errorCode = 4;
            }

            return errorCode;
        }

        /// <summary>
        /// Works out which experimental design this run will use, reports whatever is wrong with it,
        /// and returns the process exit code: 0 to carry on, 5 to stop.
        ///
        /// Static, and parameterised over its console I/O, so the decision can be tested without
        /// driving a whole command line run. These branches decide whether a search starts at all,
        /// and before this they were reachable only by launching the CLI against real spectra.
        /// </summary>
        /// <param name="designDirectory">The folder beside the spectra files, where a design lives.</param>
        /// <param name="startingRawFilenameList">Full paths of the spectra files this run will search.</param>
        /// <param name="normalizationRequested">True when any search task asks for normalization, which is what makes a design mandatory rather than optional.</param>
        /// <param name="reportToConsole">True at minimal or normal verbosity. At "none" there is nobody to ask, so a recoverable problem is recovered from silently.</param>
        /// <param name="write">Console writer; defaults to <see cref="Console.WriteLine(string)"/>.</param>
        /// <param name="readLine">Console reader; defaults to <see cref="Console.ReadLine"/>.</param>
        public static int ResolveExperimentalDesign(
            string designDirectory,
            List<string> startingRawFilenameList,
            bool normalizationRequested,
            bool reportToConsole,
            Action<string> write = null,
            Func<string> readLine = null)
        {
            write ??= Console.WriteLine;
            readLine ??= Console.ReadLine;

            string pathToExperDesign = Path.Combine(designDirectory, GlobalVariables.ExperimentalDesignFileName);
            string pathToTmtDesign = Path.Combine(designDirectory, GlobalVariables.TmtExperimentalDesignFileName);

            bool hasClassicDesign = File.Exists(pathToExperDesign);
            bool hasTmtDesign = File.Exists(pathToTmtDesign);

            // Only one may be present. Two files carrying replicate structure can drift apart, and
            // nothing here could say which of them the user meant.
            if (hasClassicDesign && hasTmtDesign)
            {
                if (reportToConsole)
                {
                    write("Both ExperimentalDesign and TmtDesign.txt were found. Only one design file type may be present.");
                }

                return 5;
            }

            if (!hasClassicDesign && !hasTmtDesign)
            {
                // A design is optional until something needs it. Normalization is that something.
                if (normalizationRequested)
                {
                    if (reportToConsole)
                    {
                        write("No experimental design file present. Normalization requires a design (ExperimentalDesign.tsv or TmtDesign.txt).");
                    }

                    return 5;
                }

                return 0;
            }

            string designPath = hasClassicDesign ? pathToExperDesign : pathToTmtDesign;
            List<string> errors;

            if (hasClassicDesign)
            {
                ExperimentalDesign.ReadExperimentalDesign(designPath, startingRawFilenameList, out errors);
            }
            else
            {
                TmtExperimentalDesign.Read(designPath, startingRawFilenameList, out errors);
            }

            if (!errors.Any())
            {
                if (reportToConsole)
                {
                    write("Read " + Path.GetFileName(designPath) + " successfully");
                }

                return 0;
            }

            // Normalizing against a design that could not be read would silently produce different
            // numbers, so this is the one case that cannot be recovered from.
            if (normalizationRequested)
            {
                if (reportToConsole)
                {
                    foreach (string error in errors)
                    {
                        write(error);
                    }
                }

                return 5;
            }

            // Otherwise the run can continue without a design -- but only by discarding the one that
            // is there, which is the user's call whenever there is a user to ask.
            if (!reportToConsole)
            {
                File.Delete(designPath);
                return 0;
            }

            write((hasClassicDesign ? "An experimental design file" : "A TMT design file")
                + " was found, but an error occurred reading it. Delete and continue empty? y/n");
            write("First error: " + errors.First());

            string answer = readLine();
            if (answer?.ToLowerInvariant() == "y" || answer?.ToLowerInvariant() == "yes")
            {
                File.Delete(designPath);
                return 0;
            }

            return 5;
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
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Starting task: " + e.DisplayName);
                MyWriter.Indent++;
            }
        }

        private static void MyTaskEngine_finishedWritingFileHandler(object sender, SingleFileEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Finished writing file: " + e.WrittenFile);
            }
        }

        /// <summary>
        /// Carries out --acceptThermoLicence: prints the licence, records the agreement, and reports
        /// whether this was the flag on its own, in which case the caller has no run to continue into.
        /// Does nothing and reports false when the flag was not given.
        /// The data directory and the output writer are parameters rather than reached for through
        /// Console and GlobalVariables so that this can be exercised directly, Run() being private.
        /// </summary>
        public static bool AgreeToThermoLicence(CommandLineSettings settings, string dataDirectory, TextWriter output)
        {
            if (!settings.AcceptThermoLicence)
            {
                return false;
            }

            if (!GlobalVariables.GlobalSettings.UserHasAgreedToThermoRawFileReaderLicence)
            {
                // Print the licence whatever the verbosity. The flag is an affirmative act by
                // whoever typed it, and this is what puts the terms they agreed to in the log.
                output.WriteLine(ThermoRawFileReaderLicence.ThermoLicenceText);
                output.WriteLine("\nThe --acceptThermoLicence flag was given, which agrees to the above terms.");

                RecordThermoLicenceAgreement(dataDirectory);
            }

            if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                output.WriteLine("Agreed to the Thermo RawFileReader licence. Recorded in "
                    + Path.Combine(dataDirectory, @"settings.toml"));
            }

            // On its own the flag is a setup step, so there is no run to continue into.
            return settings.Tasks.Count + settings.Databases.Count + settings.Spectra.Count == 0;
        }

        /// <summary>
        /// Records agreement to the Thermo RawFileReader licence, in memory and in settings.toml,
        /// carrying the one other setting that file holds so recording the agreement does not reset it.
        /// </summary>
        private static void RecordThermoLicenceAgreement(string dataDirectory)
        {
            var newGlobalSettings = new GlobalSettings
            {
                UserHasAgreedToThermoRawFileReaderLicence = true,
                WriteExcelCompatibleTSVs = GlobalVariables.GlobalSettings.WriteExcelCompatibleTSVs
            };

            Toml.WriteFile<GlobalSettings>(newGlobalSettings, Path.Combine(dataDirectory, @"settings.toml"));
            GlobalVariables.GlobalSettings = newGlobalSettings;
        }

        private static void MyTaskEngine_finishedSingleTaskHandler(object sender, SingleTaskEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal || CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.minimal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;
            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal || CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.minimal)
            {
                MyWriter.Indent--;
                WriteMultiLineIndented("Finished task: " + e.DisplayName);
            }
        }

        private static void MyEngine_startingSingleEngineHander(object sender, SingleEngineEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Starting engine: " + e.MyEngine.GetType().Name + " " + e.MyEngine.GetId());
                MyWriter.Indent++;
            }
        }

        private static void MyEngine_finishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Engine results: " + e);
                MyWriter.Indent--;
                WriteMultiLineIndented("Finished engine: " + e.MyResults.MyEngine.GetType().Name + " " + e.MyResults.MyEngine.GetId());
            }
        }

        private static void MyEngine_outProgressHandler(object sender, ProgressEventArgs e)
        {
            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                MyWriter.Write(e.NewProgress + " ");
            }
            InProgress = true;
        }

        private static void WarnHandler(object sender, StringEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.minimal || CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.minimal || CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("WARN: " + e.S);
            }
        }

        private static void LogHandler(object sender, StringEventArgs e)
        {
            if (InProgress)
            {
                if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    MyWriter.WriteLine();
                }
            }

            InProgress = false;

            if (CommandLineSettings.Verbosity == CommandLineSettings.VerbosityType.normal)
            {
                WriteMultiLineIndented("Log: " + e.S);
            }
        }
    }
}
