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
using TaskLayer.ParallelSearch;
// Aliased rather than imported: `using Readers` collides with IO.ThermoRawFileReader over
// ThermoRawFileReaderLicence, which this file uses, and CS0104 is not worth a whole namespace.
using SdrfQuantAudit = Readers.SdrfQuantAudit;
using SdrfQuantAuditor = Readers.SdrfQuantAuditor;
using SdrfSearchDesign = Readers.SdrfSearchDesign;
using SdrfSearchDesignOptions = Readers.SdrfSearchDesignOptions;
using MzLibIsobaricMassTag = Omics.Modifications.IsobaricMassTag;

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

            // subscribed before startup so that what startup notices, such as a custom protease that
            // collided with a built-in, reaches the console; WarnHandler reads the verbosity from here.
            // Removed first so that calling Run more than once in a process does not print it twice.
            CommandLineSettings = settings;
            GlobalVariables.WarnHandler -= WarnHandler;
            GlobalVariables.WarnHandler += WarnHandler;
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

            if (settings.AuditSdrf != null)
            {
                return AuditSdrf(settings, Console.Out);
            }

            if (settings.SdrfDesign != null)
            {
                return WriteDesignFromSdrf(settings, Console.Out);
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

                    case "ParallelSearch":
                        var parallelSearchTask = Toml.ReadFile<ParallelSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "ParallelSearchTask", parallelSearchTask));
                        break;

                    case "Truncation":
                        var truncationTask = Toml.ReadFile<TruncationSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "TruncationSearchTask", truncationTask));
                        break;

                    case "DiaLibrarySearch":
                        var DiaTask = Toml.ReadFile<DiaLibrarySearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                        taskList.Add(("Task" + (i + 1) + "DiaLibrarySearchTask", DiaTask));
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

            bool reportToConsole = settings.Verbosity == CommandLineSettings.VerbosityType.minimal
                || settings.Verbosity == CommandLineSettings.VerbosityType.normal;

            // An SDRF beside the spectra files makes the design, before the design check reads it; the runner hands the
            // SDRF to the search tasks for this run only.
            int sdrfExitCode = ResolveSdrfDesign(
                Directory.GetParent(startingRawFilenameList.First()).FullName,
                startingRawFilenameList,
                taskList.Where(p => p.Item2.TaskType == MyTask.Search).Select(p => (p.Item1, (SearchTask)p.Item2)).ToList(),
                settings,
                searchTasks.Any(p => p.SearchParameters.Normalize),
                reportToConsole,
                out SdrfForRun runSdrf);

            if (sdrfExitCode != 0)
            {
                return sdrfExitCode;
            }

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

            // a task whose digestion gives seed peptides it cannot use would finish with far fewer, wrong identifications
            int digestionExitCode = RefuseSeedDigestion(taskList, Console.Out,
                settings.Verbosity == CommandLineSettings.VerbosityType.minimal
                    || settings.Verbosity == CommandLineSettings.VerbosityType.normal);
            if (digestionExitCode != 0)
            {
                return digestionExitCode;
            }

            EverythingRunnerEngine a = new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList, settings.OutputFolder,
                runSdrf);

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
        /// Refuses the run when a task's digestion settings ask for seed peptides the task cannot use (see
        /// <see cref="MetaMorpheusTask.GetSeedDigestionRefusal"/>). Checked for every task before any runs, as the GUI's Run
        /// button does, and returns the process exit code: 0 to carry on, 6 to stop.
        /// </summary>
        /// <param name="print">Whether to write the reason to <paramref name="output"/>; the run is refused either way.</param>
        public static int RefuseSeedDigestion(List<(string, MetaMorpheusTask)> taskList, TextWriter output, bool print)
        {
            foreach (var (taskName, task) in taskList)
            {
                string refusal = task.GetSeedDigestionRefusal(taskName);
                if (refusal != null)
                {
                    if (print)
                    {
                        output.WriteLine(refusal);
                    }
                    return 6;
                }
            }
            return 0;
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

            // The rule is shared with the GUI's Run button (EngineLayer.ExperimentalDesignCheck); this only shows it.
            var check = ExperimentalDesignCheck.Check(designDirectory, startingRawFilenameList, normalizationRequested);

            if (check.Verdict != ExperimentalDesignVerdict.AskToContinueWithoutDesign)
            {
                if (reportToConsole)
                {
                    foreach (string message in check.Messages)
                    {
                        write(message);
                    }
                }

                return check.Verdict == ExperimentalDesignVerdict.Proceed ? 0 : 5;
            }

            // An unreadable design nothing needs: the run may continue without it, and the file is left alone. It used to
            // be deleted -- silently at verbosity none, and on a y/n whose question was about continuing rather than about
            // deletion (#2256). At verbosity none there is nobody to ask, so the run continues.
            if (!reportToConsole)
            {
                return 0;
            }

            write(check.Question + " y/n");
            write("First error: " + check.Errors.First());

            string answer = readLine();
            if (answer?.ToLowerInvariant() == "y" || answer?.ToLowerInvariant() == "yes")
            {
                return 0;
            }

            return 5;
        }

        /// <summary>
        /// The automatic design for a command-line run (<see cref="SdrfAutomaticDesign.Resolve"/>, shared with the GUI),
        /// shown on the console. A design that cannot be made from the SDRF is handled as
        /// <see cref="ResolveExperimentalDesign"/> handles a design file it cannot read (QuantProject's oracle verdict,
        /// 2026-11-07): with normalization asked for, every reason and exit 5; at verbosity none, the run continues; otherwise
        /// every reason and a y/n prompt to continue without a design. A design MetaMorpheus wrote before that no longer fits
        /// and could not be removed (<see cref="SdrfDesignOutcome.StopsTheRun"/>) stops the run with exit 5 whatever its
        /// normalization and verbosity, with no prompt. Runs before <see cref="ResolveExperimentalDesign"/>,
        /// which then reads the design written here exactly as it reads one the user wrote. The search tasks are not
        /// changed: the SDRF comes back as <paramref name="runSdrf"/>, for the runner to hand them for this run only.
        /// </summary>
        /// <param name="designDirectory">The folder of the first spectra file, where a run reads its design.</param>
        /// <param name="startingRawFilenameList">Full paths of the spectra files this run will search.</param>
        /// <param name="searchTasks">The run's search tasks, by name; the design serves them all.</param>
        /// <param name="settings">The command line, for the --sdrf... switches.</param>
        /// <param name="normalizationRequested">True when any search task normalizes, which makes a design mandatory.</param>
        /// <param name="reportToConsole">True at minimal or normal verbosity.</param>
        /// <param name="runSdrf">The SDRF beside the spectra files, for <see cref="EverythingRunnerEngine"/>; null with none.</param>
        /// <param name="write">Console writer; defaults to <see cref="Console.WriteLine(string)"/>.</param>
        /// <param name="readLine">Console reader; defaults to <see cref="Console.ReadLine"/>.</param>
        /// <returns>0 to carry on, 5 to stop.</returns>
        public static int ResolveSdrfDesign(
            string designDirectory,
            List<string> startingRawFilenameList,
            IReadOnlyList<(string Name, SearchTask Task)> searchTasks,
            CommandLineSettings settings,
            bool normalizationRequested,
            bool reportToConsole,
            out SdrfForRun runSdrf,
            Action<string> write = null,
            Func<string> readLine = null)
        {
            write ??= Console.WriteLine;
            readLine ??= Console.ReadLine;

            var outcome = SdrfAutomaticDesign.Resolve(designDirectory, startingRawFilenameList, searchTasks,
                new SdrfDesignStatements
                {
                    ConditionColumns = settings.SdrfConditionColumns?.ToList() ?? new List<string>(),
                    SingleCondition = settings.SdrfSingleCondition,
                    UnknownFactorIsOwnCondition = settings.SdrfUnknownFactorIsCondition,
                    PlexColumn = settings.SdrfPlexColumn,
                    PlexPattern = settings.SdrfPlexPattern,
                    SinglePlex = settings.SdrfSinglePlex,
                },
                CommandLineWording);
            runSdrf = outcome.RunSdrf;

            if (reportToConsole)
            {
                foreach (string message in outcome.Messages)
                {
                    write(message);
                }
            }

            if (outcome.Refusal == null)
            {
                return 0;
            }

            // A design MetaMorpheus wrote before that no longer fits, still there: going on "without a design" would read it
            // and quantify with it, so there is nothing to ask (Alexander-Sol's review of #2980).
            if (outcome.StopsTheRun)
            {
                if (reportToConsole)
                {
                    write(outcome.Refusal + " The design MetaMorpheus wrote from the SDRF before could not be removed, and a run "
                        + "that went on would quantify with it, so the run stops.");
                }
                return 5;
            }

            if (normalizationRequested)
            {
                if (reportToConsole)
                {
                    write(outcome.Refusal + " Normalization needs an experimental design, so the run stops. No design was written.");
                }
                return 5;
            }

            if (!reportToConsole)
            {
                return 0;
            }

            write(outcome.Refusal + " Continue without an experimental design? y/n");
            string answer = readLine();
            return answer?.ToLowerInvariant() == "y" || answer?.ToLowerInvariant() == "yes" ? 0 : 5;
        }

        /// <summary>How the command line names the ways to state a grouping: its switches.</summary>
        private static readonly SdrfDesignWording CommandLineWording = new()
        {
            ReportHint = SwitchesTheReportNames,
            HowToStateAnotherGrouping = "Give --sdrfCondition, --sdrfSingleCondition or --sdrfUnknownFactorIsCondition to state another, "
                + "or delete the design file to make it again from every factor value column.",
        };

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
        /// Carries out --sdrfDesign: reads the SDRF for exactly the -s files with mzLib's SdrfSearchDesign, and writes
        /// ExperimentalDesign.tsv (label-free) or TmtDesign.txt (with --sdrfTag) beside the first spectra file, which is
        /// where a run given the same -s looks for it (ResolveExperimentalDesign).
        ///
        /// SdrfSearchDesign is the one reading every SDRF route shares, so a pipeline that checks an SDRF before a search
        /// (PXReprise) and this switch cannot disagree about it: it restricts the SDRF to the searched files with
        /// SdrfSearchScope (a row naming X.raw stands for a search of X.mzML, X-calib.mzML or X-averaged.mzML), then reads
        /// the kept rows with the label-free or isobaric reader, whose rules decide the design. Its report says what the
        /// join set aside and what the reader added, noted or refused.
        ///
        /// Refuses rather than writes anything MetaMorpheus would reject, because an invalid design is worse than none: a
        /// run finding one skips quantification with only a warning. The report is printed whatever the verbosity, as for
        /// --auditSdrf: it is what was asked for, and it is the only record of a biological replicate number added where
        /// the SDRF gives none.
        /// </summary>
        /// <returns>0 when the design was written; 4 when the SDRF could not be read; 5 when the design was refused or a
        /// design file is already there.</returns>
        public static int WriteDesignFromSdrf(CommandLineSettings settings, TextWriter output)
        {
            // -o is not honoured: a design anywhere else is one the run never reads. Said rather than
            // dropped silently, because "Wrote <path>" alone would read as -o having been obeyed.
            if (settings.OutputFolder != null)
            {
                output.WriteLine("-o is ignored with --sdrfDesign: the design is written beside the first spectra file, where a run looks for it.");
            }

            string designDirectory = Directory.GetParent(Path.GetFullPath(settings.Spectra.First())).FullName;
            bool isobaric = settings.SdrfTag != null;
            string designPath = Path.Combine(designDirectory,
                isobaric ? GlobalVariables.TmtExperimentalDesignFileName : GlobalVariables.ExperimentalDesignFileName);

            // Never overwrite a design: it may be the user's own, and a run with both kinds present exits 5.
            foreach (string existing in new[]
                     {
                         Path.Combine(designDirectory, GlobalVariables.ExperimentalDesignFileName),
                         Path.Combine(designDirectory, GlobalVariables.TmtExperimentalDesignFileName)
                     })
            {
                if (File.Exists(existing))
                {
                    output.WriteLine("A design file is already there, and it is not overwritten: " + existing);
                    return 5;
                }
            }

            // validated already, so a kit given is one MetaMorpheus knows
            MzLibIsobaricMassTag tag = null;
            if (isobaric)
            {
                MzLibIsobaricMassTag.TryGetIsobaricMassTag(settings.SdrfTag, out tag);
            }

            var design = SdrfSearchDesign.Read(settings.SdrfDesign, settings.Spectra.Select(Path.GetFullPath).ToList(),
                SdrfDesignOptions(settings, tag));

            output.WriteLine(design.Report());
            output.Write(SdrfAutomaticDesign.DefaultPlexNote(design, tag != null,
                SdrfAutomaticDesign.PlexIsStated(settings.SdrfPlexColumn, settings.SdrfPlexPattern, settings.SdrfSinglePlex)));

            if (!design.IsValid)
            {
                output.WriteLine("No design file was written.");

                output.Write(SwitchesTheReportNames(design));

                // Neither the join nor a reader ran only when the file itself could not be read: the settings were
                // checked before, so the other such refusals (no searched file, a plex source without a kit) cannot
                // arise here. Told apart by what ran, never by the report's words.
                return design.LabelFree is null && design.Isobaric is null && design.Scoping is null ? 4 : 5;
            }

            // MetaMorpheus has no TMT16 multiplex label: TMTpro 16-plex is searched with TMT18, the same
            // reagent. TmtDesign.txt names no kit, and TMT16's channels are TMT18's first sixteen, so the
            // file is one a TMT18 run reads as is. Said, so the search is not set to a label that is not there.
            if (tag?.TagType == Omics.Modifications.IsobaricMassTagType.TMT16)
            {
                output.WriteLine(Tmt16Note);
            }

            design.WriteDesign(designPath);
            output.WriteLine("Wrote " + designPath);
            return 0;
        }

        /// <summary>
        /// MetaMorpheus has no TMT16 multiplex label: TMTpro 16-plex is searched with TMT18, the same reagent. TmtDesign.txt
        /// names no kit, and TMT16's channels are TMT18's first sixteen, so the file is one a TMT18 run reads as is. Said, so
        /// the search is not set to a label that is not there.
        /// </summary>
        private const string Tmt16Note = SdrfAutomaticDesign.Tmt16Note;

        /// <summary>
        /// Where a refusal can be answered by stating a grouping, mzLib's report names its option; this says which switch
        /// that is, one line each, or nothing. Only whether a line is printed hangs on the report's words.
        /// </summary>
        private static string SwitchesTheReportNames(SdrfSearchDesign design)
        {
            var lines = new System.Text.StringBuilder();
            foreach (var (option, sw) in new[]
                     {
                         (nameof(SdrfSearchDesignOptions.SingleCondition), "--sdrfSingleCondition"),
                         (nameof(SdrfSearchDesignOptions.UnknownFactorIsOwnCondition), "--sdrfUnknownFactorIsCondition"),
                     })
            {
                if (design.Refusals.Any(r => r.Contains("(" + option + ")")))
                {
                    lines.AppendLine($"On the command line, {option} is {sw}.");
                }
            }
            return lines.ToString();
        }

        /// <summary>
        /// What the command line tells SdrfSearchDesign that the SDRF alone cannot say: the kit and where each file's plex
        /// is written (isobaric only), and how the condition is built. The condition columns are passed as given, and
        /// mzLib names them: 'genotype' is 'factor value[genotype]', the 'factor value[' prefix in any case is the
        /// specification's, and any other name with a '[' (such as 'characteristics[cell type]') is that column. With
        /// none, every factor value column makes the condition. A single condition and files whose factor nobody stated
        /// are groupings only the user can state (ruled 2026-11-11): never defaults.
        /// </summary>
        private static SdrfSearchDesignOptions SdrfDesignOptions(CommandLineSettings settings, MzLibIsobaricMassTag tag)
        {
            // An isobaric SDRF's plex: the switch given, else the SDRF's comment[plex] (ruled 2026-11-13).
            var plex = tag != null
                ? SdrfAutomaticDesign.PlexSource(settings.SdrfDesign, settings.SdrfPlexColumn, settings.SdrfPlexPattern, settings.SdrfSinglePlex)
                : (settings.SdrfPlexColumn, settings.SdrfPlexPattern, settings.SdrfSinglePlex);
            return new SdrfSearchDesignOptions
            {
                Tag = tag,
                PlexColumn = plex.Item1,
                PlexFileNamePattern = plex.Item2,
                SinglePlex = plex.Item3,
                ConditionColumns = settings.SdrfConditionColumns?.ToList(),
                SingleCondition = settings.SdrfSingleCondition,
                UnknownFactorIsOwnCondition = settings.SdrfUnknownFactorIsCondition,
            };
        }

        /// <summary>
        /// Prints what one SDRF document says about quantification, and runs nothing else.
        ///
        /// The failure this exists to catch is silent. A file whose comment[label] names the reagent kit
        /// ("TMT10plex" on every row) rather than the channel parses cleanly and validates cleanly, and
        /// then yields one channel for a ten-channel experiment. So does a file with two rows per
        /// (data file, channel), or one that declares no isobaric modification at all. None of that is a
        /// structural defect, so the SDRF validator is right not to report it -- it is a defect in what
        /// the annotation MEANS, and the only way a user finds out today is a wrong number downstream.
        ///
        /// All the reading is mzLib's SdrfQuantAuditor; this prints its report. Scoped to every SDRF
        /// rather than to TMT, because "this is not an isobaric experiment" is a useful answer to get.
        /// </summary>
        /// <returns>0 when the report was written, 4 when the document could not be read.</returns>
        public static int AuditSdrf(CommandLineSettings settings, TextWriter output)
        {
            SdrfQuantAudit audit;

            try
            {
                audit = SdrfQuantAuditor.Audit(settings.AuditSdrf, settings.AuditDataDirectory);
            }
            catch (Exception e)
            {
                // The auditor does not throw on a document it can read; getting here means it could not
                // be read at all -- not an SDRF, truncated, locked. Reported with the path, because the
                // usual cause is pointing at the wrong file.
                if (settings.Verbosity == CommandLineSettings.VerbosityType.minimal || settings.Verbosity == CommandLineSettings.VerbosityType.normal)
                {
                    output.WriteLine("The SDRF file could not be read: " + settings.AuditSdrf
                                     + Environment.NewLine + e.Message);
                }

                return 4;
            }

            // Written whatever the verbosity, including 'none'. -v suppresses MetaMorpheus talking about
            // itself; the report is the thing that was asked for, and a run that produced it silently
            // would have produced nothing at all.
            output.WriteLine(audit.ToReport());

            return 0;
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
