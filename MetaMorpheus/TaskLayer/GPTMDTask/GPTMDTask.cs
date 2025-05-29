﻿using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Gptmd;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using Proteomics.ProteolyticDigestion;
using System.Globalization;
using Omics.Digestion;
using Omics.Modifications;
using System.Threading.Tasks;

namespace TaskLayer
{
    public class GptmdTask : MetaMorpheusTask
    {
        private const double tolForComboLoading = 1e-3;

        public GptmdTask() : base(MyTask.Gptmd)
        {
            CommonParameters = new CommonParameters();
            GptmdParameters = new GptmdParameters();
        }

        public GptmdParameters GptmdParameters { get; set; }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyFileManager myFileManager = new MyFileManager(true);

            // start loading first spectra file in the background
            Task<MsDataFile> nextFileLoadingTask = new(() => myFileManager.LoadFile(currentRawFileList[0], SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[0])));
            nextFileLoadingTask.Start();
            LoadModifications(taskId, out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);

            // start loading proteins in the background
            List<Protein> proteinList = null;
            Task<List<Protein>> proteinLoadingTask = new(() =>
            {
                var proteins = LoadProteins(taskId, dbFilenameList, true, DecoyType.Reverse,
                    localizeableModificationTypes,
                    CommonParameters);
                SanitizeProteinDatabase(proteins, TargetContaminantAmbiguity.RemoveContaminant);
                return proteins;
            });
            proteinLoadingTask.Start();

            // TODO: print error messages loading GPTMD mods
            List<Modification> gptmdModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => GptmdParameters.ListOfModsGptmd.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            IEnumerable<Tuple<double, double>> combos = LoadCombos(gptmdModifications).ToList();

            // write prose settings
            ProseCreatedWhileRunning.Append("The following G-PTM-D settings were used: "); 
            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.DigestionAgent + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxLength + "; ");
            if (CommonParameters.DigestionParams is DigestionParams digestionParams)
                ProseCreatedWhileRunning.Append("initiator methionine behavior = " + digestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("G-PTM-D modifications count = " + gptmdModifications.Count + "; ");

            // temporary search type for writing prose
            // the actual search type is technically file-specific but we don't allow file-specific notches, so it's safe to do this
            MassDiffAcceptor tempSearchMode = new DotMassDiffAcceptor("", GetAcceptableMassShifts(fixedModifications, variableModifications, gptmdModifications, combos), CommonParameters.PrecursorMassTolerance);
            ProseCreatedWhileRunning.Append("precursor mass tolerance(s) = {" + tempSearchMode.ToProseString() + "}; ");

            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + ". ");

            // start the G-PTM-D task
            Status("Running G-PTM-D...", new List<string> { taskId });
            MyTaskResults = new MyTaskResults(this)
            {
                NewDatabases = new List<DbForTask>()
            };

            var filePathToPrecursorMassTolerance = new Dictionary<string, Tolerance>();
            for (int i = 0; i < currentRawFileList.Count; i++)
            {
                string filePath = currentRawFileList[i];
                Tolerance fileTolerance = CommonParameters.PrecursorMassTolerance;
                if (fileSettingsList[i] != null && fileSettingsList[i].PrecursorMassTolerance != null)
                {
                    fileTolerance = fileSettingsList[i].PrecursorMassTolerance;
                }
                filePathToPrecursorMassTolerance.Add(filePath, fileTolerance);
            }

            // Store the psms from each file in one combined list
            List<SpectralMatch> allPsms = new();
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                var origDataFile = currentRawFileList[spectraFileIndex];

                // mark the file as in-progress
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);
                MassDiffAcceptor searchMode = new DotMassDiffAcceptor("", GetAcceptableMassShifts(fixedModifications, variableModifications, gptmdModifications, combos), combinedParams.PrecursorMassTolerance);

                NewCollection(Path.GetFileName(origDataFile), new List<string> { taskId, "Individual Spectra Files", origDataFile });

                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files", origDataFile });

                // ensure that the next file has finished loading from the async method
                nextFileLoadingTask.Wait();
                var myMsDataFile = nextFileLoadingTask.Result;
                // if another file exists, then begin loading it in while the previous is being searched
                if (origDataFile != currentRawFileList.Last())
                {
                    int nextFileIndex = spectraFileIndex + 1;
                    nextFileLoadingTask = new Task<MsDataFile>(() => myFileManager.LoadFile(currentRawFileList[nextFileIndex], SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[nextFileIndex])));
                    nextFileLoadingTask.Start();
                }
                Status("Getting ms2 scans...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy(b => b.PrecursorMass).ToArray();
                myFileManager.DoneWithFile(origDataFile);
                SpectralMatch[] psmArray = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];

                //spectral Library search and library generation have't applied to GPTMD yet
                bool writeSpctralLibrary = false;

                // ensure proteins are loaded in before proceeding with search
                switch (proteinLoadingTask.IsCompleted)
                {
                    case true when proteinList is null: // has finished loading but not been set
                        proteinList = proteinLoadingTask.Result;
                        break;
                    case true when proteinList.Any(): // has finished loading and already been set
                        break;
                    case false: // has not finished loading
                        proteinLoadingTask.Wait();
                        proteinList = proteinLoadingTask.Result;
                        break;
                }

                new ClassicSearchEngine(psmArray, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, null, null, null, 
                    proteinList, searchMode, combinedParams, this.FileSpecificParameters, null, new List<string> { taskId, "Individual Spectra Files", origDataFile }, writeSpctralLibrary).Run();


                var psmList = psmArray.Where(p => p != null).ToList();
                foreach(var psm in psmList)
                {
                    psm.SetMs2Scan(arrayOfMs2ScansSortedByMass[psm.ScanIndex].TheScan);
                }
                allPsms.AddRange(psmList);

                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
            }
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            // Run FDR analysis on allPsms
            // GPTMD doesn't work as well if you do FDR on a file-by-file basis. Presumably this is because it takes multiple files to get enough PSMs for all the different notches
            new FdrAnalysisEngine(allPsms, tempSearchMode.NumNotches, CommonParameters, this.FileSpecificParameters, new List<string> { taskId }, doPEP: false).Run();
            Dictionary<string, HashSet<Tuple<int, Modification>>> allModDictionary = new();
            new GptmdEngine(allPsms, gptmdModifications, combos, filePathToPrecursorMassTolerance, CommonParameters, this.FileSpecificParameters, new List<string> { taskId }, allModDictionary).Run();

            //Move this text after search because proteins don't get loaded until search begins.
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count(p => !p.IsDecoy) + " non-decoy protein entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            // run GPTMD engine
            Status("Creating the GPTMD Database", new List<string> { taskId });

            // write GPTMD databases
            if (dbFilenameList.Any(b => !b.IsContaminant))
            {
                List<string> databaseNames = new List<string>();
                foreach (var nonContaminantDb in dbFilenameList.Where(p => !p.IsContaminant))
                {
                    var dbName = Path.GetFileNameWithoutExtension(nonContaminantDb.FilePath);
                    var theExtension = Path.GetExtension(nonContaminantDb.FilePath).ToLowerInvariant();
                    bool compressed = theExtension.EndsWith("gz");
                    databaseNames.Add(compressed ? Path.GetFileNameWithoutExtension(dbName) : dbName);
                }
                string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", databaseNames) + "GPTMD.xml");

                var newModsActuallyWritten = ProteinDbWriter.WriteXmlDatabase(allModDictionary, proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                FinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });

                MyTaskResults.NewDatabases.Add(new DbForTask(outputXMLdbFullName, false));
                if(dbFilenameList.Any(p=>p.IsSpectralLibrary))
                {
                    MyTaskResults.NewDatabases.Add(dbFilenameList.Where(p => p.IsSpectralLibrary).First());
                }
                MyTaskResults.AddTaskSummaryText("Modifications added: " + newModsActuallyWritten.Select(b => b.Value).Sum());
                MyTaskResults.AddTaskSummaryText("Mods types and counts:");
                MyTaskResults.AddTaskSummaryText(string.Join(Environment.NewLine, newModsActuallyWritten.OrderByDescending(b => b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            }
            if (dbFilenameList.Any(b => b.IsContaminant))
            {
                // do NOT use this code (Path.GetFilenameWithoutExtension) because GPTMD on .xml.gz will result in .xml.xml file type being written
                //string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "GPTMD.xml");
                List<string> databaseNames = new List<string>();
                foreach (var contaminantDb in dbFilenameList.Where(p => p.IsContaminant))
                {
                    var dbName = Path.GetFileName(contaminantDb.FilePath);
                    int indexOfFirstDot = dbName.IndexOf(".");
                    databaseNames.Add(dbName.Substring(0, indexOfFirstDot));
                }
                string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", databaseNames) + "GPTMD.xml");

                var newModsActuallyWritten = ProteinDbWriter.WriteXmlDatabase(allModDictionary, proteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                FinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });
                MyTaskResults.NewDatabases.Add(new DbForTask(outputXMLdbFullNameContaminants, true));
                MyTaskResults.AddTaskSummaryText("Contaminant modifications added: " + newModsActuallyWritten.Select(b => b.Value).Sum());
                MyTaskResults.AddTaskSummaryText("Mods types and counts:");
                MyTaskResults.AddTaskSummaryText(string.Join(Environment.NewLine, newModsActuallyWritten.OrderByDescending(b => b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            }
            Status("Done", new List<string> { taskId });
            return MyTaskResults;
        }

        private static IEnumerable<Tuple<double, double>> LoadCombos(List<Modification> modificationsThatCanBeCombined)
        {
            using (StreamReader r = new StreamReader(Path.Combine(GlobalVariables.DataDir, "Data", @"combos.txt")))
            {
                while (r.Peek() >= 0)
                {
                    var line = r.ReadLine().Split(' ');
                    var mass1 = double.Parse(line[0], CultureInfo.InvariantCulture);
                    var mass2 = double.Parse(line[1], CultureInfo.InvariantCulture);
                    if (modificationsThatCanBeCombined.Any(b => b.ValidModification == true && Math.Abs((double)b.MonoisotopicMass - mass1) < tolForComboLoading) &&
                        modificationsThatCanBeCombined.Any(b => b.ValidModification == true && Math.Abs((double)b.MonoisotopicMass - mass2) < tolForComboLoading))
                        yield return new Tuple<double, double>(mass1, mass2);
                }
            }
        }

        private static IEnumerable<double> GetAcceptableMassShifts(List<Modification> fixedMods, List<Modification> variableMods, List<Modification> gptmdMods, IEnumerable<Tuple<double, double>> combos)
        {
            IEnumerable<double> gptmdNotches = gptmdMods.Where(b => b.ValidModification == true).Select(b => (double)b.MonoisotopicMass);
            IEnumerable<double> gptmdMinusOtherModsNotches = GetObservedMasses(variableMods.Concat(fixedMods), gptmdMods);
            IEnumerable<double> multipleGptmdNotches = combos.Select(b => b.Item1 + b.Item2);
            IEnumerable<double> zeroNotch = new List<double> { 0 };

            IEnumerable<double> allNotches = gptmdNotches.Concat(gptmdMinusOtherModsNotches).Concat(multipleGptmdNotches).Concat(zeroNotch);
            return allNotches.GroupBy(b => Math.Round(b, 5)).Select(b => b.FirstOrDefault()).OrderBy(b => b);
        }

        private static IEnumerable<double> GetObservedMasses(IEnumerable<Modification> enumerable, List<Modification> gptmdModifications)
        {
            foreach (var modOnPeptide in enumerable.Where(b => b.ValidModification == true))
            {
                foreach (var modToLocalize in gptmdModifications.Where(b => b.ValidModification == true))
                {
                    if (modOnPeptide.Target.Equals(modToLocalize.Target))
                    {
                        yield return (double)modToLocalize.MonoisotopicMass - (double)modOnPeptide.MonoisotopicMass;
                    }
                }
            }
        }
    }
}