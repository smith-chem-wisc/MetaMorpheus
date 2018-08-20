using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Gptmd;

#if NETFRAMEWORK

using IO.Thermo;

#else
#endif

using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using Proteomics.ProteolyticDigestion;

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
            // load modifications
            Status("Loading modifications...", taskId);
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
            List<string> localizeableModificationTypes = GlobalVariables.AllModTypesKnown.ToList();
            List<ModificationWithMass> gptmdModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => GptmdParameters.ListOfModsGptmd.Contains((b.modificationType, b.id))).ToList();
            IEnumerable<Tuple<double, double>> combos = LoadCombos(gptmdModifications).ToList();

            // what types of fragment ions to search for
            List<ProductType> ionTypes = new List<ProductType>();
            if (CommonParameters.BIons)
                ionTypes.Add(ProductType.BnoB1ions);
            if (CommonParameters.YIons)
                ionTypes.Add(ProductType.Y);
            if (CommonParameters.ZdotIons)
                ionTypes.Add(ProductType.Zdot);
            if (CommonParameters.CIons)
                ionTypes.Add(ProductType.C);

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, true, DecoyType.Reverse, localizeableModificationTypes, CommonParameters);

            List<PeptideSpectralMatch> allPsms = new List<PeptideSpectralMatch>();

            var numRawFiles = currentRawFileList.Count;

            // write prose settings
            ProseCreatedWhileRunning.Append("The following G-PTM-D settings were used: "); ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");
            ProseCreatedWhileRunning.Append("G-PTM-D modifications count = " + gptmdModifications.Count + "; ");

            // temporary search type for writing prose
            // the actual search type is technically file-specific but we don't allow file-specific notches, so it's safe to do this
            MassDiffAcceptor tempSearchMode = new DotMassDiffAcceptor("", GetAcceptableMassShifts(fixedModifications, variableModifications, gptmdModifications, combos), CommonParameters.PrecursorMassTolerance);
            ProseCreatedWhileRunning.Append("precursor mass tolerance(s) = {" + tempSearchMode.ToProseString() + "}; ");

            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + ". ");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count(p => !p.IsDecoy) + " non-decoy protein entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            // start the G-PTM-D task
            Status("Running G-PTM-D...", new List<string> { taskId });
            MyTaskResults = new MyTaskResults(this)
            {
                NewDatabases = new List<DbForTask>()
            };
            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));
            HashSet<DigestionParams> ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams));

            MyFileManager myFileManager = new MyFileManager(true);

            object lock1 = new object();
            object lock2 = new object();

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
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks, combinedParams);
                Status("Getting ms2 scans...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                myFileManager.DoneWithFile(origDataFile);
                PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];
                PeptideSpectralMatch[][] allDecoyPsmsArray = new PeptideSpectralMatch[CommonParameters.NumDecoyDatabases][];
                new ClassicSearchEngine(allPsmsArray, allDecoyPsmsArray, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ionTypes, searchMode, combinedParams, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();
                lock (lock2)
                {
                    allPsms.AddRange(allPsmsArray.Where(p => p != null));
                }
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
            }
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            // Group and order psms
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, ionTypes, ListOfDigestionParams, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string> { taskId });

            var resTest = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatchingTest = resTest.CompactPeptideToProteinPeptideMatching;

            foreach (var huh in allPsms)
            {
                huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatchingTest);
            }

            allPsms = allPsms.OrderByDescending(b => b.Score)
                .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                .GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass))
                .Select(b => b.First()).ToList();

            new FdrAnalysisEngine(allPsms, tempSearchMode.NumNotches, CommonParameters, new List<string> { taskId }).Run();

            var writtenFile = Path.Combine(OutputFolder, "GPTMD_Candidates.psmtsv");
            WritePsmsToTsv(allPsms, writtenFile, new Dictionary<string, int>(), double.PositiveInfinity);
            FinishedWritingFile(writtenFile, new List<string> { taskId });

            // get file-specific precursor mass tolerances for the GPTMD engine
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

            // run GPTMD engine
            var gptmdResults = (GptmdResults)new GptmdEngine(allPsms, gptmdModifications, combos, filePathToPrecursorMassTolerance, CommonParameters, new List<string> { taskId }).Run();

            // Stop if canceled
            if (GlobalVariables.StopLoops) { return MyTaskResults; }

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

                var newModsActuallyWritten = ProteinDbWriter.WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                FinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });

                MyTaskResults.NewDatabases.Add(new DbForTask(outputXMLdbFullName, false));
                MyTaskResults.AddNiceText("Modifications added: " + newModsActuallyWritten.Select(b => b.Value).Sum());
                MyTaskResults.AddNiceText("Mods types and counts:");
                MyTaskResults.AddNiceText(string.Join(Environment.NewLine, newModsActuallyWritten.OrderByDescending(b => b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
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

                var newModsActuallyWritten = ProteinDbWriter.WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                FinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });

                MyTaskResults.NewDatabases.Add(new DbForTask(outputXMLdbFullNameContaminants, true));
                MyTaskResults.AddNiceText("Contaminant modifications added: " + newModsActuallyWritten.Select(b => b.Value).Sum());
                MyTaskResults.AddNiceText("Mods types and counts:");
                MyTaskResults.AddNiceText(string.Join(Environment.NewLine, newModsActuallyWritten.OrderByDescending(b => b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            }
            return MyTaskResults;
        }

        private static IEnumerable<Tuple<double, double>> LoadCombos(List<ModificationWithMass> modificationsThatCanBeCombined)
        {
            using (StreamReader r = new StreamReader(Path.Combine(GlobalVariables.DataDir, "Data", @"combos.txt")))
            {
                while (r.Peek() >= 0)
                {
                    var line = r.ReadLine().Split(' ');
                    var mass1 = double.Parse(line[0]);
                    var mass2 = double.Parse(line[1]);
                    if (modificationsThatCanBeCombined.Any(b => Math.Abs(b.monoisotopicMass - mass1) < tolForComboLoading) &&
                        modificationsThatCanBeCombined.Any(b => Math.Abs(b.monoisotopicMass - mass2) < tolForComboLoading))
                        yield return new Tuple<double, double>(mass1, mass2);
                }
            }
        }

        private static IEnumerable<double> GetAcceptableMassShifts(List<ModificationWithMass> fixedMods, List<ModificationWithMass> variableMods, List<ModificationWithMass> gptmdMods, IEnumerable<Tuple<double, double>> combos)
        {
            IEnumerable<double> gptmdNotches = gptmdMods.Select(b => b.monoisotopicMass);
            IEnumerable<double> gptmdMinusOtherModsNotches = GetObservedMasses(variableMods.Concat(fixedMods), gptmdMods);
            IEnumerable<double> multipleGptmdNotches = combos.Select(b => b.Item1 + b.Item2);
            IEnumerable<double> zeroNotch = new List<double> { 0 };

            IEnumerable<double> allNotches = gptmdNotches.Concat(gptmdMinusOtherModsNotches).Concat(multipleGptmdNotches).Concat(zeroNotch);
            return allNotches.GroupBy(b => Math.Round(b, 5)).Select(b => b.FirstOrDefault()).OrderBy(b => b);
        }

        private static IEnumerable<double> GetObservedMasses(IEnumerable<ModificationWithMass> enumerable, List<ModificationWithMass> gptmdModifications)
        {
            foreach (var modOnPeptide in enumerable)
            {
                foreach (var modToLocalize in gptmdModifications)
                {
                    if (modOnPeptide.motif.Equals(modToLocalize.motif))
                    {
                        yield return modToLocalize.monoisotopicMass - modOnPeptide.monoisotopicMass;
                    }
                }
            }
        }
    }
}