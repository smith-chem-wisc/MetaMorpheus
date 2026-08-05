using CommandLine;
using EngineLayer;
using Nett;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace MetaMorpheusCommandLine
{
    public class CommandLineSettings
    {
        public List<string> Spectra { get; private set; }
        public List<string> Tasks { get; private set; }
        public List<string> Databases { get; private set; }

        [Option('t', HelpText = "Single-task TOMLs (.toml file format); space-delimited")]
        public IEnumerable<string> _tasks { get; set; }

        [Option('d', HelpText = "Protein sequence databases (.fasta, .xml, .fasta.gz, .xml.gz file formats); space-delimited")]
        public IEnumerable<string> _databases { get; set; }

        [Option('s', HelpText = "Spectra to analyze (.raw, .mzML, .mgf, ms2.msalign file formats) or folder(s) containing spectra; space-delimited")]
        public IEnumerable<string> _spectra { get; set; }

        [Option('o', HelpText = "[Optional] Output folder")]
        public string OutputFolder { get; set; }

        [Option('g', HelpText = "[Optional] Generate default task tomls")]
        public bool GenerateDefaultTomls { get; set; }

        [Option('v', Default = VerbosityType.normal, HelpText = "[Optional] Determines how much text is written. Options are no output ('none'), minimal output and errors  ('minimal'), or normal ('normal')")]
        public VerbosityType Verbosity { get; set; }

        [Option("test", HelpText = "[Optional] Runs a small test search using a database and yeast data file included with this MetaMorpheus installation")]
        public bool RunMicroVignette { get; set; }

        [Option("mmsettings", HelpText = "[Optional] Path to MetaMorpheus settings")]
        public string CustomDataDirectory { get; set; }

        public enum VerbosityType { none, minimal, normal };

        public void ValidateCommandLineSettings()
        {
            Spectra = _spectra == null ? new List<string>() : _spectra.ToList();
            Tasks = _tasks == null ? new List<string>() : _tasks.ToList();
            Databases = _databases == null ? new List<string>() : _databases.ToList();

            if ((GenerateDefaultTomls || RunMicroVignette) && OutputFolder == null)
            {
                throw new MetaMorpheusException("An output path must be specified with the -o parameter.");
            }

            if (GenerateDefaultTomls || RunMicroVignette)
            {
                return;
            }

            if (OutputFolder == null && Spectra.Count < 1 && Tasks.Count < 1 && Databases.Count < 1 && Spectra.Count < 1)
            {
                throw new MetaMorpheusException("Use the --help parameter to view all parameters (e.g., \"CMD.exe --help\")");
            }

            if (OutputFolder == null && Spectra.Count > 0)
            {
                var pathOfFirstSpectraFile = Path.GetDirectoryName(Spectra.First());
                OutputFolder = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
            }

            if (Tasks.Count < 1)
            {
                throw new MetaMorpheusException("At least one task must be specified.");
            }

            if (Databases.Count < 1)
            {
                throw new MetaMorpheusException("At least one protein database must be specified.");
            }

            if (Spectra.Count < 1)
            {
                throw new MetaMorpheusException("At least one spectra file must be specified.");
            }

            // add spectra files from specified directories
            List<string> spectraFromDirectories = new List<string>();
            foreach (string item in Spectra)
            {
                if (Directory.Exists(item))
                {
                    FindSpectraFilesRecursive(item, spectraFromDirectories, Verbosity);
                }
                else if (!File.Exists(item))
                {
                    throw new MetaMorpheusException("The following is not a known file or directory: " + item);
                }
            }

            Spectra.AddRange(spectraFromDirectories);

            // Correct Bruker inner-file paths (.tdf/.tdf_bin for timsTOF, .baf for regular TOF) to their parent .d directory
            // This handles the case where a user provides a path to one of these files directly
            for (int i = 0; i < Spectra.Count; i++)
            {
                string ext = Path.GetExtension(Spectra[i]).ToLowerInvariant();
                if (ext == ".tdf" || ext == ".tdf_bin" || ext == ".baf")
                {
                    string parentDir = Path.GetDirectoryName(Spectra[i]);
                    if (parentDir != null &&
                        parentDir.EndsWith(".d", StringComparison.OrdinalIgnoreCase) &&
                        IsValidBrukerDirectory(parentDir))
                    {
                        Spectra[i] = parentDir;
                    }
                }
            }

            // Remove duplicate spectra entries (can occur if both .tdf and .tdf_bin were specified)
            Spectra = Spectra.Distinct().ToList();

            // remove spectra directories, after their spectra files have been added
            // but keep .d directories as they are valid Bruker spectra folders
            Spectra.RemoveAll(p => Directory.Exists(p) && !p.EndsWith(".d", StringComparison.OrdinalIgnoreCase));

            IEnumerable<string> fileNames = Tasks.Concat(Databases).Concat(Spectra);

            foreach (string filename in fileNames)
            {
                // .d folders are directories, so we need to check for both files and .d directories
                bool isDotDDirectory = filename.EndsWith(".d", StringComparison.OrdinalIgnoreCase) && Directory.Exists(filename);
                if (!File.Exists(filename) && !isDotDDirectory)
                {
                    throw new MetaMorpheusException("The following file does not exist: " + filename);
                }
            }

            foreach (string filename in Tasks)
            {
                string ext = Path.GetExtension(filename).ToLowerInvariant();

                if (ext != ".toml")
                {
                    throw new MetaMorpheusException("Tasks must be in .toml file format. Unrecognized file format: " + ext);
                }
            }

            foreach (string filename in Databases)
            {
                string ext = Path.GetExtension(filename).ToLowerInvariant();
                bool compressed = ext.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
                ext = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(filename)).ToLowerInvariant() : ext;

                if (!GlobalVariables.AcceptedDatabaseFormats.Contains(ext))
                {
                    throw new MetaMorpheusException("Unrecognized protein database file format: " + ext);
                }
            }

            foreach (string filename in Spectra)
            {
                string ext = Path.GetExtension(filename).ToLowerInvariant();

                if (!GlobalVariables.AcceptedSpectraFormats.Contains(ext))
                {
                    throw new MetaMorpheusException("Unrecognized spectra file format: " + ext);
                }
            }
        }

        public static void GenerateDefaultTaskTomls(string folderLocation)
        {
            try
            {
                if (!Directory.Exists(folderLocation))
                {
                    Directory.CreateDirectory(folderLocation);
                }

                CalibrationTask c = new CalibrationTask();
                Toml.WriteFile(c, Path.Combine(folderLocation, @"CalibrationTask.toml"), MetaMorpheusTask.tomlConfig);

                GptmdTask g = new GptmdTask();
                Toml.WriteFile(g, Path.Combine(folderLocation, @"GptmdTask.toml"), MetaMorpheusTask.tomlConfig);

                SearchTask s = new SearchTask();
                Toml.WriteFile(s, Path.Combine(folderLocation, @"SearchTask.toml"), MetaMorpheusTask.tomlConfig);

                XLSearchTask xl = new XLSearchTask();
                Toml.WriteFile(xl, Path.Combine(folderLocation, @"XLSearchTask.toml"), MetaMorpheusTask.tomlConfig);

                GlycoSearchTask glyco = new GlycoSearchTask();
                Toml.WriteFile(glyco, Path.Combine(folderLocation, @"GlycoSearchTask.toml"), MetaMorpheusTask.tomlConfig);
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException("Default tomls could not be written: " + e.Message);
            }
        }

        /// <summary>
        /// Recursively finds spectra files in a directory.
        /// For .d folders (timsTOF data), only adds them if they are valid (contain required tdf/tsf files).
        /// If a .d folder is invalid, recurses into it to find nested valid .d folders.
        /// </summary>
        private static void FindSpectraFilesRecursive(string path, List<string> spectraFiles, VerbosityType verbosity)
        {
            if (File.Exists(path))
            {
                if (GlobalVariables.AcceptedSpectraFormats.Contains(GlobalVariables.GetFileExtension(path).ToLowerInvariant()))
                {
                    // If a Bruker inner file is found (.tdf/.tdf_bin for timsTOF, .baf for regular TOF),
                    // check if the parent directory is a valid .d folder and add that instead
                    if (new[] { ".tdf", ".tdf_bin", ".baf" }.Contains(GlobalVariables.GetFileExtension(path).ToLowerInvariant()))
                    {
                        var parentDirectory = Directory.GetParent(path);
                        if (GlobalVariables.GetFileExtension(parentDirectory.FullName).ToLowerInvariant() == ".d" && IsValidBrukerDirectory(parentDirectory.FullName))
                            path = parentDirectory.FullName; // add the .d folder instead of the individual inner file
                    }

                    spectraFiles.Add(path);

                    if (verbosity == VerbosityType.normal)
                    {
                        Console.WriteLine("Found spectra file: " + path);
                    }
                }
                return;
            }

            if (Directory.Exists(path))
            {
                // Check if this is a .d folder
                if (path.EndsWith(".d", StringComparison.OrdinalIgnoreCase))
                {
                    // Check if it's a valid Bruker data folder (timsTOF or regular TOF)
                    if (IsValidBrukerDirectory(path))
                    {
                        spectraFiles.Add(path);

                        if (verbosity == VerbosityType.normal)
                        {
                            Console.WriteLine("Found spectra file: " + path);
                        }
                        return; // Valid .d folder - don't recurse into it
                    }
                    // Invalid .d folder - fall through to recurse and find nested valid .d folders
                }

                // Recurse into directory
                foreach (var entry in Directory.EnumerateFileSystemEntries(path))
                {
                    FindSpectraFilesRecursive(entry, spectraFiles, verbosity);
                }
            }
        }

        /// <summary>
        /// Checks if a .d directory is a valid Bruker data folder that MetaMorpheus can read,
        /// i.e. a timsTOF folder (analysis.tdf + analysis.tdf_bin) or a regular TOF folder (analysis.baf).
        /// </summary>
        private static bool IsValidBrukerDirectory(string directoryPath)
        {
            // timsTOF data: needs both the .tdf and .tdf_bin files
            if (File.Exists(Path.Combine(directoryPath, "analysis.tdf")) &&
                File.Exists(Path.Combine(directoryPath, "analysis.tdf_bin")))
            {
                return true;
            }

            // regular TOF data: has a .baf file
            if (File.Exists(Path.Combine(directoryPath, "analysis.baf")))
            {
                return true;
            }

            return false;
        }
    }
}
