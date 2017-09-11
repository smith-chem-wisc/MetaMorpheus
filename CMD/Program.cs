using Chemistry;
using EngineLayer;
using Fclp;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace MetaMorpheusCommandLine
{
    internal static class Program
    {
        #region Private Methods

        private static void Main(string[] args)
        {
            Console.WriteLine(GlobalEngineLevelSettings.MetaMorpheusVersion);

            var goodAAsub = new List<string> {
            "KN",
"KT",
"KR",
"KI",
"KQ",
"KE",
"NK",
"NT",
"NL",
"NI",
"NH",
"ND",
"NY",
"KM",
"NS",
"TK",
"TR",
"TI",
"TP",
"TA",
"TS",
"TN",
"TL",
"TM",
"RK",
"RT",
"RL",
"RS",
"RI",
"RG",
"SN",
"ST",
"SR",
"SK",
"SI",
"SG",
"SC",
"RM",
"RW",
"SL",
"IK",
"IT",
"IR",
"IM",
"IL",
"IV",
"IN",
"IF",
"MK",
"MT",
"MI",
"ML",
"MV",
"IS",
"QK",
"QH",
"QP",
"QR",
"QL",
"QE",
"HN",
"HQ",
"HP",
"HR",
"HL",
"HD",
"HY",
"PT",
"PQ",
"PR",
"PL",
"PA",
"PS",
"PH",
"RQ",
"RP",
"RH",
"RC",
"LI",
"LQ",
"LP",
"LR",
"LV",
"LH",
"LF",
"LM",
"EK",
"EQ",
"ED",
"EA",
"EG",
"EV",
"DN",
"DH",
"DE",
"DA",
"DG",
"DV",
"DY",
"AT",
"AP",
"AE",
"AG",
"AV",
"AS",
"AD",
"GR",
"GE",
"GA",
"GV",
"GL",
"GD",
"GC",
"GK",
"GW",
"GS",
"VI",
"VL",
"VE",
"VA",
"VG",
"VD",
"VF",
"VM",
"YN",
"YH",
"YD",
"YS",
"YC",
"YF",
"SP",
"SA",
"SY",
"SF",
"SW",
"CL",
"CR",
"CG",
"CY",
"CS",
"CW",
"CF",
"WK",
"WR",
"WG",
"WS",
"WC",
"WL",
"LS",
"FI",
"FL",
"FV",
"FY",
"FS",
"FC",
"LW"};

            List<ModificationWithMass> subs = new List<ModificationWithMass>();
            string abc = "ABCDEFGHIJKLMNOPQRSTVWXYZ";
            foreach (var l1 in abc)
            {
                if (Residue.TryGetResidue(l1, out Residue r1))
                {
                    foreach (var l2 in abc)
                    {
                        if (!l1.Equals(l2) && Residue.TryGetResidue(l2, out Residue r2))
                        {
                            ChemicalFormula cf = new ChemicalFormula();
                            cf.Add(r2);
                            cf.Remove(r1);

                            ModificationMotif.TryGetMotif(r1.Letter.ToString(), out ModificationMotif motif);
                            if (goodAAsub.Contains(r1.Letter.ToString() + r2.Letter.ToString()))
                            {
                                subs.Add(new ModificationWithMassAndCf(r1.Letter.ToString() + "->" + r2.Letter.ToString(),
                                    "1 nucleotide substitution",
                                    motif,
                                    TerminusLocalization.Any,
                                    cf));
                            }
                            else
                            {
                                subs.Add(new ModificationWithMassAndCf(r1.Letter.ToString() + "->" + r2.Letter.ToString(),
                                    "2+ nucleotide substitution",
                                    motif,
                                    TerminusLocalization.Any,
                                    cf));
                            }
                        }
                    }
                }
            }

            foreach (var ok in subs.OrderBy(b => b.monoisotopicMass))
            {
                Console.WriteLine(ok);
                Console.WriteLine("//");
            }

            var p = new FluentCommandLineParser<ApplicationArguments>();

            Console.WriteLine(string.Join(" , ", args));

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
                MetaMorpheusEngine.WarnHandler += MyEngine_outLabelStatusHandler;
                MetaMorpheusEngine.OutProgressHandler += MyEngine_outProgressHandler;

                MetaMorpheusTask.WarnHandler += MyEngine_outLabelStatusHandler;

                foreach (var modFile in Directory.GetFiles(Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Mods")))
                    GlobalEngineLevelSettings.AddMods(UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(modFile));

                GlobalEngineLevelSettings.AddMods(GlobalEngineLevelSettings.UnimodDeserialized.OfType<ModificationWithLocation>());
                GlobalEngineLevelSettings.AddMods(GlobalEngineLevelSettings.UniprotDeseralized.OfType<ModificationWithLocation>());

                foreach (var db in p.Object.Databases)
                    if (!Path.GetExtension(db).Equals(".fasta"))
                        GlobalEngineLevelSettings.AddMods(UsefulProteomicsDatabases.ProteinDbLoader.GetPtmListFromProteinXml(db).OfType<ModificationWithLocation>());

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

                        case "XLSearch":
                            var ye4 = Toml.ReadFile<XLSearchTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                            taskList.Add(new Tuple<string, MetaMorpheusTask>("Task" + (i + 1) + "XLSearchTask", ye4));
                            break;
                    }
                }
                List<string> startingRawFilenameList = p.Object.Spectra.Select(b => Path.GetFullPath(b)).ToList();
                List<DbForTask> startingXmlDbFilenameList = p.Object.Databases.Select(b => new DbForTask(Path.GetFullPath(b), IsContaminant(b))).ToList();
                EverythingRunnerEngine a = new EverythingRunnerEngine(taskList, startingRawFilenameList, startingXmlDbFilenameList, null);
                a.Run();
            }
            Console.WriteLine("Error Text:" + result.ErrorText);
            Console.WriteLine("Version: {0}", Environment.Version.ToString());
            Console.WriteLine("OSVersion: {0}", Environment.OSVersion.ToString());
            Console.WriteLine("EmptyArgs:" + result.EmptyArgs);
            Console.WriteLine("EmptyArgs:" + string.Join(" , ", result.Errors.Select(b => b.Option.Description)));
            Console.WriteLine("Usage:");
            Console.WriteLine("\t-t --tasks     List of task poml files");
            Console.WriteLine("\t-s --spectra   List of spectra files");
            Console.WriteLine("\t-d --databases List of database files");
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
            //if (inProgress)
            //    Console.WriteLine();
            //inProgress = false;
            //Console.WriteLine("Starting task:");
            //Console.WriteLine(e.TaskId);
        }

        private static void MyTaskEngine_finishedWritingFileHandler(object sender, SingleFileEventArgs e)
        {
            //if (inProgress)
            //    Console.WriteLine();
            //inProgress = false;
            //Console.WriteLine("Finished writing file: " + e.writtenFile);
        }

        private static void MyTaskEngine_finishedSingleTaskHandler(object sender, SingleTaskEventArgs e)
        {
            //if (inProgress)
            //    Console.WriteLine();
            //inProgress = false;
            //Console.WriteLine("Finished task: " + e.TaskId.GetType().Name);
        }

        private static void MyEngine_startingSingleEngineHander(object sender, SingleEngineEventArgs e)
        {
            //if (inProgress)
            //    Console.WriteLine();
            //inProgress = false;
            //Console.WriteLine("Starting engine:" + e.myEngine.GetType().Name);
        }

        private static void MyEngine_outProgressHandler(object sender, ProgressEventArgs e)
        {
            //Console.Write(e.new_progress + " ");
            //inProgress = true;
        }

        private static void MyEngine_outLabelStatusHandler(object sender, StringEventArgs e)
        {
            //if (inProgress)
            //    Console.WriteLine();
            //inProgress = false;
            //Console.WriteLine("Status: " + e.s);
            Console.WriteLine(e.s);
        }

        private static void MyEngine_finishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
        {
            //if (inProgress)
            //    Console.WriteLine();
            //inProgress = false;
            //Console.WriteLine("Finished engine: ");
            //Console.WriteLine(e);
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