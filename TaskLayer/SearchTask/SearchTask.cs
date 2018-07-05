using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public class SearchTask : MetaMorpheusTask
    {
        public SearchTask() : base(TaskType.Search)
        {
            CommonParameters = new CommonParameters();

            SearchParameters = new SearchParameters();
        }

        public SearchParameters SearchParameters { get; set; }

        private object IndexLock { get; } = new object();
        private object PsmLock { get; } = new object();
        private List<ModificationWithMass> VariableModifications { get; set; }
        private List<ModificationWithMass> FixedModifications { get; set; }
        private List<ProductType> IonTypes { get; set; }
        private HashSet<DigestionParams> ListOfDigestionParams { get; set; }
        private List<PeptideSpectralMatch> AllPsms { get; set; }

        public static MassDiffAcceptor GetMassDiffAcceptor(Tolerance precursorMassTolerance, MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact:
                    return precursorMassTolerance is PpmTolerance ?
                        new SinglePpmAroundZeroSearchMode(precursorMassTolerance.Value) as MassDiffAcceptor :
                        new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance.Value) as MassDiffAcceptor;

                case MassDiffAcceptorType.OneMM:
                    return new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, precursorMassTolerance);

                case MassDiffAcceptorType.TwoMM:
                    return new DotMassDiffAcceptor("2mm", new List<double> { 0, 1.0029, 2.0052 }, precursorMassTolerance);

                case MassDiffAcceptorType.ThreeMM:
                    return new DotMassDiffAcceptor("3mm", new List<double> { 0, 1.0029, 2.0052, 3.0077 }, precursorMassTolerance);

                case MassDiffAcceptorType.ModOpen:
                    return new IntervalMassDiffAcceptor("-187andUp", new List<DoubleRange> { new DoubleRange(-187, double.PositiveInfinity) });

                case MassDiffAcceptorType.Open:
                    return new OpenSearchMode();

                case MassDiffAcceptorType.Custom:
                    return ParseSearchMode(customMdac);

                default:
                    throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        protected override MyTaskResults RunSpecific(string outputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            // Load parameters
            MyFileManager myFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);
            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));
            ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams));

            // disable quantification if a .mgf is being used
            if (SearchParameters.DoQuantification && currentRawFileList.Any(x => Path.GetExtension(x).Equals(".mgf", StringComparison.OrdinalIgnoreCase)))
            {
                SearchParameters.DoQuantification = false;
            }

            // load modifications
            Status("Loading modifications...", taskId);
            VariableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            FixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
            List<string> localizeableModificationTypes = GlobalVariables.AllModTypesKnown.ToList();

            // what types of fragment ions to search for
            IonTypes = new List<ProductType>();
            if (CommonParameters.BIons && CommonParameters.AddCompIons) { IonTypes.Add(ProductType.B); }
            else if (CommonParameters.BIons) { IonTypes.Add(ProductType.BnoB1ions); }
            if (CommonParameters.YIons) { IonTypes.Add(ProductType.Y); }
            if (CommonParameters.ZdotIons) { IonTypes.Add(ProductType.Zdot); }
            if (CommonParameters.CIons) { IonTypes.Add(ProductType.C); }

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, SearchParameters.SearchTarget, SearchParameters.DecoyType, localizeableModificationTypes);

            // write prose settings
            ProseCreatedWhileRunning.Append("The following search settings were used: ");
            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", FixedModifications.Select(m => m.id)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", VariableModifications.Select(m => m.id)) + "; ");
            ProseCreatedWhileRunning.Append("max mods per peptide = " + CommonParameters.DigestionParams.MaxModsForPeptide + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            ProseCreatedWhileRunning.Append("precursor mass tolerance = " + CommonParameters.PrecursorMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("report PSM ambiguity = " + CommonParameters.ReportAllAmbiguity + ". ");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " non-decoy protein entries including " + proteinList.Count(p => p.IsContaminant) + " contaminant sequences. ");

            // start the search task
            MyTaskResults = new MyTaskResults(this);
            AllPsms = new List<PeptideSpectralMatch>();

            int completedFiles = 0;

            Status("Searching files...", taskId);
            Status("Searching files...", new List<string> { taskId, "Individual Spectra Files" });

            Dictionary<string, int[]> numMs2SpectraPerFile = new Dictionary<string, int[]>();
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                var origDataFile = currentRawFileList[spectraFileIndex];

                // mark the file as in-progress
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(combinedParams.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                numMs2SpectraPerFile.Add(Path.GetFileNameWithoutExtension(origDataFile), new int[] { myMsDataFile.GetAllScansList().Count(p => p.MsnOrder == 2), arrayOfMs2ScansSortedByMass.Length });
                myFileManager.DoneWithFile(origDataFile);

                // TODO: shuffle peptides and such
                var fileSpecificPsms = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];
                Search(taskId, proteinList, fileSpecificPsms, dbFilenameList, combinedParams, arrayOfMs2ScansSortedByMass, thisId, massDiffAcceptor);

                completedFiles++;
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            PostSearchAnalysisParameters parameters = new PostSearchAnalysisParameters();
            parameters.SearchTaskResults = MyTaskResults;
            parameters.SearchTaskId = taskId;
            parameters.SearchParameters = SearchParameters;
            parameters.CommonParameters = CommonParameters;
            parameters.ProteinList = proteinList;
            parameters.IonTypes = IonTypes;
            parameters.AllPsms = AllPsms;
            parameters.FixedModifications = FixedModifications;
            parameters.VariableModifications = VariableModifications;
            parameters.ListOfDigestionParams = ListOfDigestionParams;
            parameters.CurrentRawFileList = currentRawFileList;
            parameters.MyFileManager = myFileManager;
            parameters.NumNotches = GetNumNotches(SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            parameters.OutputFolder = outputFolder;
            parameters.FileSettingsList = fileSettingsList;
            parameters.NumMs2SpectraPerFile = numMs2SpectraPerFile;
            parameters.DatabaseFilenameList = dbFilenameList;
            PostSearchAnalysisTask postProcessing = new PostSearchAnalysisTask();
            postProcessing.Parameters = parameters;
            return postProcessing.Run();
        }

        private void Search(string taskId, List<Protein> proteinList, PeptideSpectralMatch[] fileSpecificPsms, List<DbForTask> dbFilenameList, CommonParameters combinedParams,
            Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass, List<string> thisId, MassDiffAcceptor massDiffAcceptor)
        {
            // modern search
            if (SearchParameters.SearchType == SearchType.Modern)
            {
                for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                {
                    List<CompactPeptide> peptideIndex = null;
                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count / combinedParams.TotalPartitions));

                    Status("Getting fragment dictionary...", new List<string> { taskId });
                    var indexEngine = new IndexingEngine(proteinListSubset, VariableModifications, FixedModifications, IonTypes, currentPartition, SearchParameters.DecoyType, ListOfDigestionParams, combinedParams, SearchParameters.MaxFragmentSize, new List<string> { taskId });
                    List<int>[] fragmentIndex = null;
                    lock (IndexLock)
                    {
                        GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, taskId);
                    }

                    Status("Searching files...", taskId);
                    new ModernSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, IonTypes, currentPartition, combinedParams, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();
                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                }
            }
            // nonspecific search
            else if (SearchParameters.SearchType == SearchType.NonSpecific)
            {
                List<List<ProductType>> terminusSeparatedIons = ProductTypeMethods.SeparateIonsByTerminus(IonTypes);
                foreach (List<ProductType> terminusSpecificIons in terminusSeparatedIons)
                {
                    for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                    {
                        List<CompactPeptide> peptideIndex = null;
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count / combinedParams.TotalPartitions));

                        List<int>[] fragmentIndex = new List<int>[1];

                        // Generate indices for nonspecifc search
                        Status("Getting fragment dictionary...", new List<string> { taskId });
                        var indexEngine = new IndexingEngine(proteinListSubset, VariableModifications, FixedModifications, terminusSpecificIons, currentPartition, SearchParameters.DecoyType, ListOfDigestionParams, combinedParams, SearchParameters.MaxFragmentSize, new List<string> { taskId });
                        lock (IndexLock)
                        {
                            GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, taskId);
                        }

                        // Generate indices for nonspecific search
                        Status("Getting precursor dictionary...", new List<string> { taskId });
                        List<CompactPeptide> peptideIndexPrecursor = null;
                        List<Protein> proteinListSubsetPrecursor = proteinList.GetRange(currentPartition * proteinList.Count / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count / combinedParams.TotalPartitions));
                        List<int>[] fragmentIndexPrecursor = new List<int>[1];
                        var indexEnginePrecursor = new PrecursorIndexingEngine(proteinListSubsetPrecursor, VariableModifications, FixedModifications, terminusSpecificIons, currentPartition, SearchParameters.DecoyType, ListOfDigestionParams, combinedParams, 0, new List<string> { taskId });
                        lock (IndexLock)
                        {
                            GenerateIndexes(indexEnginePrecursor, dbFilenameList, ref peptideIndexPrecursor, ref fragmentIndexPrecursor, taskId);
                        }

                        if (peptideIndex.Count != peptideIndexPrecursor.Count)
                        {
                            throw new MetaMorpheusException("peptideIndex not identical between indexing engines");
                        }

                        Status("Searching files...", taskId);
                        new NonSpecificEnzymeSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, fragmentIndexPrecursor, terminusSpecificIons, currentPartition, combinedParams, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();
                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                    }
                }
            }
            // classic search
            else
            {
                Status("Starting search...", thisId);
                new ClassicSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, VariableModifications, FixedModifications, proteinList, IonTypes, massDiffAcceptor, combinedParams, thisId).Run();

                ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
            }
            lock (PsmLock)
            {
                AllPsms.AddRange(fileSpecificPsms);
            }
        }

        private int GetNumNotches(MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact: return 1;
                case MassDiffAcceptorType.OneMM: return 2;
                case MassDiffAcceptorType.TwoMM: return 3;
                case MassDiffAcceptorType.ThreeMM: return 4;
                case MassDiffAcceptorType.ModOpen: return 1;
                case MassDiffAcceptorType.Open: return 1;
                case MassDiffAcceptorType.Custom: return ParseSearchMode(customMdac).NumNotches;

                default: throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        private static MassDiffAcceptor ParseSearchMode(string text)
        {
            MassDiffAcceptor ye = null;

            var split = text.Split(' ');

            switch (split[1])
            {
                case "dot":

                    var massShifts = Array.ConvertAll(split[4].Split(','), Double.Parse);
                    var newString = split[2].Replace("�", "");
                    var toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                    if (split[3].ToUpperInvariant().Equals("PPM"))
                    {
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new PpmTolerance(toleranceValue));
                    }
                    else if (split[3].ToUpperInvariant().Equals("DA"))
                    {
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new AbsoluteTolerance(toleranceValue));
                    }
                    break;

                case "interval":
                    IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[0], CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[1], CultureInfo.InvariantCulture)));
                    ye = new IntervalMassDiffAcceptor(split[0], doubleRanges);
                    break;

                case "OpenSearch":
                    ye = new OpenSearchMode();
                    break;

                case "daltonsAroundZero":
                    ye = new SingleAbsoluteAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                case "ppmAroundZero":
                    ye = new SinglePpmAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                default:
                    throw new MetaMorpheusException("Could not parse search mode string");
            }
            return ye;
        }

        private static IEnumerable<Type> GetSubclassesAndItself(Type type)
        {
            yield return type;
        }

        private static bool SameSettings(string pathToOldParamsFile, IndexingEngine indexEngine)
        {
            using (StreamReader reader = new StreamReader(pathToOldParamsFile))
            {
                if (reader.ReadToEnd().Equals(indexEngine.ToString()))
                {
                    return true;
                }
            }
            return false;
        }

        private static void WritePeptideIndex(List<CompactPeptide> peptideIndex, string peptideIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(peptideIndexFile))
            {
                ser.Serialize(file, peptideIndex);
            }
        }

        private static void WriteFragmentIndexNetSerializer(List<int>[] fragmentIndex, string fragmentIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<int>[]));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFile))
                ser.Serialize(file, fragmentIndex);
        }

        private static string GetExistingFolderWithIndices(IndexingEngine indexEngine, List<DbForTask> dbFilenameList)
        {
            // In every database location...
            foreach (var ok in dbFilenameList)
            {
                var baseDir = Path.GetDirectoryName(ok.FilePath);
                var directory = new DirectoryInfo(baseDir);
                DirectoryInfo[] directories = directory.GetDirectories();

                // Look at every subdirectory...
                foreach (DirectoryInfo possibleFolder in directories)
                {
                    if (File.Exists(Path.Combine(possibleFolder.FullName, "indexEngine.params")) &&
                        File.Exists(Path.Combine(possibleFolder.FullName, "peptideIndex.ind")) &&
                        File.Exists(Path.Combine(possibleFolder.FullName, "fragmentIndex.ind")) &&
                        SameSettings(Path.Combine(possibleFolder.FullName, "indexEngine.params"), indexEngine))
                    {
                        return possibleFolder.FullName;
                    }
                }
            }
            return null;
        }

        private static void WriteIndexEngineParams(IndexingEngine indexEngine, string fileName)
        {
            using (StreamWriter output = new StreamWriter(fileName))
            {
                output.Write(indexEngine);
            }
        }

        private static string GenerateOutputFolderForIndices(List<DbForTask> dbFilenameList)
        {
            var folder = Path.Combine(Path.GetDirectoryName(dbFilenameList.First().FilePath), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
            Directory.CreateDirectory(folder);
            return folder;
        }

        private void GenerateIndexes(IndexingEngine indexEngine, List<DbForTask> dbFilenameList, ref List<CompactPeptide> peptideIndex, ref List<int>[] fragmentIndex, string taskId)
        {
            string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);
            if (pathToFolderWithIndices == null)
            {
                var output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
                Status("Writing params...", new List<string> { taskId });
                var paramsFile = Path.Combine(output_folderForIndices, "indexEngine.params");
                WriteIndexEngineParams(indexEngine, paramsFile);
                SucessfullyFinishedWritingFile(paramsFile, new List<string> { taskId });

                Status("Running Index Engine...", new List<string> { taskId });
                var indexResults = (IndexingResults)indexEngine.Run();
                peptideIndex = indexResults.PeptideIndex;
                fragmentIndex = indexResults.FragmentIndex;

                Status("Writing peptide index...", new List<string> { taskId });
                var peptideIndexFile = Path.Combine(output_folderForIndices, "peptideIndex.ind");
                WritePeptideIndex(peptideIndex, peptideIndexFile);
                SucessfullyFinishedWritingFile(peptideIndexFile, new List<string> { taskId });

                Status("Writing fragment index...", new List<string> { taskId });
                var fragmentIndexFile = Path.Combine(output_folderForIndices, "fragmentIndex.ind");
                WriteFragmentIndexNetSerializer(fragmentIndex, fragmentIndexFile);
                SucessfullyFinishedWritingFile(fragmentIndexFile, new List<string> { taskId });
            }
            else
            {
                Status("Reading peptide index...", new List<string> { taskId });
                var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
                var ser = new NetSerializer.Serializer(messageTypes);
                using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "peptideIndex.ind")))
                {
                    peptideIndex = (List<CompactPeptide>)ser.Deserialize(file);
                }

                Status("Reading fragment index...", new List<string> { taskId });
                messageTypes = GetSubclassesAndItself(typeof(List<int>[]));
                ser = new NetSerializer.Serializer(messageTypes);
                using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "fragmentIndex.ind")))

                {
                    fragmentIndex = (List<int>[])ser.Deserialize(file);
                }
            }
        }
    }
}