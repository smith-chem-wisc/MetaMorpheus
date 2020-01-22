﻿using CommandLine;
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

        [Option('t', HelpText = "Single-task TOMLs (.toml file format)")]
        public IEnumerable<string> _tasks { get; set; }

        [Option('d', HelpText = "Protein sequence databases (.fasta, .xml, .fasta.gz, .xml.gz file formats)")]
        public IEnumerable<string> _databases { get; set; }

        [Option('s', HelpText = "Spectra to analyze (.raw, .mzML, .mgf file formats)")]
        public IEnumerable<string> _spectra { get; set; }

        [Option('o', HelpText = "Output folder")]
        public string OutputFolder { get; set; }

        [Option('g', HelpText = "Generate default task tomls")]
        public bool GenerateDefaultTomls { get; set; }

        [Option('v', HelpText = "Runs a small test search using a database and yeast data file included with this MetaMorpheus installation")]
        public bool RunMicroVignette { get; set; }
        
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

            IEnumerable<string> fileNames = Tasks.Concat(Databases).Concat(Spectra);

            foreach (string filename in fileNames)
            {
                if (!File.Exists(filename))
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
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException("Default tomls could not be written: " + e.Message);
            }
        }
    }
}
