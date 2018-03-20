using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Gptmd;

#if NETFRAMEWORK

using IO.Thermo;

#else
#endif

using MassSpectrometry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class GptmdTask : MetaMorpheusTask
    {
        #region Private Fields

        private const double tolForComboLoading = 1e-3;

        #endregion Private Fields

        #region Public Constructors

        public GptmdTask() : base(MyTask.Gptmd)
        {
            CommonParameters = new CommonParameters();
            GptmdParameters = new GptmdParameters();
        }

        #endregion Public Constructors

        #region Public Properties

        public GptmdParameters GptmdParameters { get; set; }

        #endregion Public Properties

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList)
        {
            myTaskResults = new MyTaskResults(this)
            {
                newDatabases = new List<DbForTask>()
            };
            Status("Loading modifications...", new List<string> { taskId });

            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
            List<string> localizeableModificationTypes = CommonParameters.LocalizeAll ? GlobalVariables.AllModTypesKnown.ToList() : CommonParameters.ListOfModTypesLocalize.ToList();

            List<ModificationWithMass> gptmdModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => GptmdParameters.ListOfModsGptmd.Contains((b.modificationType, b.id))).ToList();

            IEnumerable<Tuple<double, double>> combos = LoadCombos(gptmdModifications).ToList();

            List<Psm> allPsms = new List<Psm>();

            List<ProductType> ionTypes = new List<ProductType>();
            if (CommonParameters.BIons)
                ionTypes.Add(ProductType.B);
            if (CommonParameters.YIons)
                ionTypes.Add(ProductType.Y);
            if (CommonParameters.CIons)
                ionTypes.Add(ProductType.C);
            if (CommonParameters.ZdotIons)
                ionTypes.Add(ProductType.Zdot);

            Status("Loading proteins...", new List<string> { taskId });
            Dictionary<string, Modification> um = null;
            //Decoys are currently not being searched with DecoyType.None
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, DecoyType.Reverse, localizeableModificationTypes, b.IsContaminant, out um)).ToList();

            var numRawFiles = currentRawFileList.Count;

            proseCreatedWhileRunning.Append("The following G-PTM-D settings were used: "); proseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            proseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            proseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            if (CommonParameters.DigestionParams.MaxPeptideLength == null)
            {
                proseCreatedWhileRunning.Append("maximum peptide length = unspecified; ");
            }
            else
            {
                proseCreatedWhileRunning.Append("maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            }
            proseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            proseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");

            proseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("G-PTM-D modifications count = " + gptmdModifications.Count + "; ");
            //puppet searchmode for writing files. Actual searchmode is filespecific
            MassDiffAcceptor tempSearchMode = new DotMassDiffAcceptor("", GetAcceptableMassShifts(fixedModifications, variableModifications, gptmdModifications, combos), CommonParameters.PrecursorMassTolerance);
            proseCreatedWhileRunning.Append("parent mass tolerance(s) = {" + tempSearchMode.ToProseString() + "}; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            Status("Running G-PTM-D...", new List<string> { taskId });

            HashSet<IDigestionParams> ListOfDigestionParams = GetListOfDistinctDigestionParams(CommonParameters, fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b)));

            MyFileManager myFileManager = new MyFileManager(true);

            object lock1 = new object();
            object lock2 = new object();
            ParallelOptions parallelOptions = new ParallelOptions();
            if (CommonParameters.MaxParallelFilesToAnalyze.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxParallelFilesToAnalyze.Value;
            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);
                MassDiffAcceptor searchMode = new DotMassDiffAcceptor("", GetAcceptableMassShifts(fixedModifications, variableModifications, gptmdModifications, combos), combinedParams.PrecursorMassTolerance);

                NewCollection(Path.GetFileName(origDataFile), new List<string> { taskId, "Individual Spectra Files", origDataFile });
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                Status("Getting ms2 scans...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                myFileManager.DoneWithFile(origDataFile);
                Psm[] allPsmsArray = new Psm[arrayOfMs2ScansSortedByMass.Length];
                new ClassicSearchEngine(allPsmsArray, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ionTypes, searchMode, false, combinedParams, combinedParams.ProductMassTolerance, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();
                lock (lock2)
                {
                    allPsms.AddRange(allPsmsArray);
                }
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
            });
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            // Group and order psms

            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngineTest = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, ionTypes, ListOfDigestionParams, CommonParameters.ReportAllAmbiguity, new List<string> { taskId });

            var resTest = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngineTest.Run();
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatchingTest = resTest.CompactPeptideToProteinPeptideMatching;

            foreach (var huh in allPsms)
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatchingTest);

            allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

            new FdrAnalysisEngine(allPsms, tempSearchMode.NumNotches, false, new List<string> { taskId }).Run();

            var writtenFile = Path.Combine(OutputFolder, "GPTMD_Candidates.psmtsv");
            WritePsmsToTsv(allPsms, writtenFile, new Dictionary<string, int>());
            SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });

            var gptmdResults = (GptmdResults)new GptmdEngine(allPsms, gptmdModifications, combos, CommonParameters.PrecursorMassTolerance, new List<string> { taskId }).Run();

            if (dbFilenameList.Any(b => !b.IsContaminant))
            {
                // do NOT use this code (Path.GetFilenameWithoutExtension) because GPTMD on .xml.gz will result in .xml.xml file type being written
                //string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "GPTMD.xml");
                
                List<string> databaseNames = new List<string>();
                foreach (var nonContaminantDb in dbFilenameList.Where(p => !p.IsContaminant))
                {
                    var dbName = Path.GetFileName(nonContaminantDb.FilePath);
                    int indexOfFirstDot = dbName.IndexOf(".");
                    databaseNames.Add(dbName.Substring(0, indexOfFirstDot));
                }
                string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", databaseNames) + "GPTMD.xml");

                var newModsActuallyWritten = ProteinDbWriter.WriteXmlDatabase(gptmdResults.Mods, proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });

                myTaskResults.newDatabases.Add(new DbForTask(outputXMLdbFullName, false));
                myTaskResults.AddNiceText("Modifications added: " + newModsActuallyWritten.Select(b => b.Value).Sum());
                myTaskResults.AddNiceText("Mods types and counts:");
                myTaskResults.AddNiceText(string.Join(Environment.NewLine, newModsActuallyWritten.OrderByDescending(b => b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
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

                SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });

                myTaskResults.newDatabases.Add(new DbForTask(outputXMLdbFullNameContaminants, true));
                myTaskResults.AddNiceText("Contaminant modifications added: " + newModsActuallyWritten.Select(b => b.Value).Sum());
                myTaskResults.AddNiceText("Mods types and counts:");
                myTaskResults.AddNiceText(string.Join(Environment.NewLine, newModsActuallyWritten.OrderByDescending(b => b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            }
            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

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

        #endregion Private Methods
    }
}