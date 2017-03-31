using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;
using System.Threading;

namespace TaskLayer
{
    public class SearchTask : MetaMorpheusTask
    {

        #region Private Fields

        private const double binTolInDaltons = 0.003;

        #endregion Private Fields

        #region Public Constructors

        public SearchTask()
        {
            // Set default values here:
            ClassicSearch = true;
            DoParsimony = false;
            NoOneHitWonders = false;
            ModPeptidesAreUnique = false;
            Quantify = false;
            QuantifyRtTol = 3;
            QuantifyPpmTol = 5;
            SearchDecoy = true;
            DoHistogramAnalysis = false;
            MaxMissedCleavages = 2;
            Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            BIons = true;
            YIons = true;
            ZdotIons = false;
            CIons = false;

            ListOfModListsFixed = new List<string> { AllModLists.First(b => b.EndsWith("f.txt")) };
            ListOfModListsVariable = new List<string> { AllModLists.First(b => b.EndsWith("v.txt")) };
            ListOfModListsLocalize = new List<string> { AllModLists.First(b => b.EndsWith("ptmlist.txt")) };
            ListOfModListsToAlwaysKeep = new List<string> { AllModLists.First(b => b.EndsWith("ptmlist.txt")) };

            SearchModes = GlobalTaskLevelSettings.SearchModesKnown.Take(1).ToList();
            TaskType = MyTask.Search;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<string> ListOfModListsFixed { get; set; }
        public List<string> ListOfModListsVariable { get; set; }
        public List<string> ListOfModListsLocalize { get; set; }
        public List<string> ListOfModListsToAlwaysKeep { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public bool ClassicSearch { get; set; }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreUnique { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool Quantify { get; set; }
        public double QuantifyRtTol { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchDecoy { get; set; }
        public List<SearchMode> SearchModes { get; set; }
        public bool ConserveMemory { get; set; }

        public bool WritePrunedDatabase { get; set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string SpecificTaskInfo
        {
            get
            {
                var sb = new StringBuilder();
                sb.AppendLine("classicSearch: " + ClassicSearch);
                sb.AppendLine("doParsimony: " + DoParsimony);
                if (DoParsimony)
                {
                    sb.AppendLine("modifiedPeptidesAreUnique: " + ModPeptidesAreUnique);
                    sb.AppendLine("requireTwoPeptidesToIdProtein: " + NoOneHitWonders);
                }
                sb.AppendLine("quantify: " + Quantify);
                if (Quantify)
                {
                    sb.AppendLine("quantify rt tolerance: " + QuantifyRtTol);
                    sb.AppendLine("quantify ppm tolerance: " + QuantifyPpmTol);
                }
                sb.AppendLine("doHistogramAnalysis: " + DoHistogramAnalysis);
                sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModListsFixed));
                sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModListsVariable));
                sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModListsLocalize));
                sb.AppendLine("searchDecoy: " + SearchDecoy);
                sb.AppendLine("productMassTolerance: " + ProductMassTolerance);
                sb.AppendLine("searchModes: ");
                sb.Append(string.Join(Environment.NewLine, SearchModes.Select(b => "\t" + b.FileNameAddition)));
                return sb.ToString();
            }
        }

        #endregion Protected Properties

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId)
        {
            var mySearchTaskResults = new MyTaskResults(this);
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = ListOfModListsVariable.SelectMany(b => PtmListLoader.ReadModsFromFile(b)).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> fixedModifications = ListOfModListsFixed.SelectMany(b => PtmListLoader.ReadModsFromFile(b)).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> localizeableModifications = ListOfModListsLocalize.SelectMany(b => PtmListLoader.ReadModsFromFile(b)).OfType<ModificationWithMass>().ToList();

            List<PsmParent>[] allPsms = new List<PsmParent>[SearchModes.Count];
            for (int j = 0; j < SearchModes.Count; j++)
                allPsms[j] = new List<PsmParent>();

            Status("Loading proteins...", new List<string> { taskId });
            Dictionary<string, Modification> um;
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FileName, true, localizeableModifications, b.IsContaminant, null, out um)).ToList();

            List<CompactPeptide> peptideIndex = null;
            float[] keys = null;
            List<int>[] fragmentIndex = null;
            List<ProductType> lp = new List<ProductType>();
            if (BIons)
                lp.Add(ProductType.B);
            if (YIons)
                lp.Add(ProductType.Y);
            if (ZdotIons)
                lp.Add(ProductType.Zdot);
            if (CIons)
                lp.Add(ProductType.C);

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            if (!ClassicSearch)
            {
                Status("Getting fragment dictionary...", new List<string> { taskId });
                var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, InitiatorMethionineBehavior, MaxMissedCleavages, MaxModificationIsoforms, lp, new List<string> { taskId }, Thread.CurrentThread);
                string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);

                Dictionary<float, List<int>> fragmentIndexDict;
                if (pathToFolderWithIndices == null)
                {
                    Status("Generating indices...", new List<string> { taskId });
                    var output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
                    Status("Writing params...", new List<string> { taskId });
                    WriteIndexEngineParams(indexEngine, Path.Combine(output_folderForIndices, "indexEngine.params"), taskId);

                    var indexResults = (IndexingResults)indexEngine.Run();
                    mySearchTaskResults.AddResultText(indexResults);
                    peptideIndex = indexResults.PeptideIndex;
                    fragmentIndexDict = indexResults.FragmentIndexDict;

                    Status("Writing peptide index...", new List<string> { taskId });
                    WritePeptideIndex(peptideIndex, Path.Combine(output_folderForIndices, "peptideIndex.ind"), taskId);
                    Status("Writing fragment index...", new List<string> { taskId });
                    WriteFragmentIndexNetSerializer(fragmentIndexDict, Path.Combine(output_folderForIndices, "fragmentIndex.ind"), taskId);
                }
                else
                {
                    Status("Reading peptide index...", new List<string> { taskId });
                    var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
                    var ser = new NetSerializer.Serializer(messageTypes);
                    using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "peptideIndex.ind")))
                        peptideIndex = (List<CompactPeptide>)ser.Deserialize(file);

                    Status("Reading fragment index...", new List<string> { taskId });
                    messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
                    ser = new NetSerializer.Serializer(messageTypes);
                    using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "fragmentIndex.ind")))
                        fragmentIndexDict = (Dictionary<float, List<int>>)ser.Deserialize(file);
                }
                keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
                fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();
            }
            List<NewPsmWithFdr>[] allResultingIdentifications = null;

            Status("Searching files...", new List<string> { taskId });
            Thread taskThread = Thread.CurrentThread;
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                NewCollection(Path.GetFileName(origDataFile), new List<string> { taskId, "Individual Searches", origDataFile });
                Status("Loading spectra file...", new List<string> { taskId, "Individual Searches", origDataFile });
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                else
                    myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);
                Status("Opening spectra file...", new List<string> { taskId, "Individual Searches", origDataFile });

                if (ClassicSearch)
                {
                    var classicSearchResults = (ClassicSearchResults)new ClassicSearchEngine(MyEngine.GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray(), myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, SearchModes, MaxMissedCleavages, MaxModificationIsoforms, origDataFile, lp, new List<string> { taskId, "Individual Searches", origDataFile }, ConserveMemory, taskThread).Run();
                    mySearchTaskResults.AddResultText(classicSearchResults);
                    for (int i = 0; i < SearchModes.Count; i++)
                        allPsms[i].AddRange(classicSearchResults.OuterPsms[i]);
                    var analysisResults = new AnalysisEngine(classicSearchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups", new List<string> { taskId, "Individual Searches", origDataFile }), DoParsimony, NoOneHitWonders, ModPeptidesAreUnique, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFile }, Quantify, QuantifyRtTol, QuantifyPpmTol, taskThread).Run();
                    mySearchTaskResults.AddResultText(analysisResults);
                    allResultingIdentifications = ((AnalysisResults)analysisResults).AllResultingIdentifications;
                }
                else
                {
                    var modernSearchResults = (ModernSearchResults)new ModernSearchEngine(myMsDataFile, peptideIndex, keys, fragmentIndex, ProductMassTolerance, SearchModes, origDataFile, new List<string> { taskId, "Individual Searches", origDataFile }, taskThread).Run();
                    mySearchTaskResults.AddResultText(modernSearchResults);
                    for (int i = 0; i < SearchModes.Count; i++)
                        allPsms[i].AddRange(modernSearchResults.NewPsms[i]);
                    var analysisResults = new AnalysisEngine(modernSearchResults.NewPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups", new List<string> { taskId, "Individual Searches", origDataFile }), DoParsimony, NoOneHitWonders, ModPeptidesAreUnique, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFile }, Quantify, QuantifyRtTol, QuantifyPpmTol, null).Run();
                    mySearchTaskResults.AddResultText(analysisResults);
                    allResultingIdentifications = ((AnalysisResults)analysisResults).AllResultingIdentifications;
                }
                Status("Done!", new List<string> { taskId, "Individual Searches", origDataFile });
            }

            if (currentRawFileList.Count > 1)
            {
                var analysisResults = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + s, new List<string> { taskId }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, "aggregate" + s, new List<string> { taskId }), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, "aggregate_ProteinGroups" + s, new List<string> { taskId }), DoParsimony, NoOneHitWonders, ModPeptidesAreUnique, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId }, false, 0, 0, taskThread).Run();
                mySearchTaskResults.AddResultText(analysisResults);
                allResultingIdentifications = ((AnalysisResults)analysisResults).AllResultingIdentifications;
            }

            if (WritePrunedDatabase)
            {
                List<ModificationWithLocation> modificationsToAlwaysKeep = ListOfModListsToAlwaysKeep.SelectMany(b => PtmListLoader.ReadModsFromFile(b)).ToList();

                var goodPsmsForEachProtein = allResultingIdentifications.SelectMany(b => b).Where(b => b.QValueNotch < 0.01 && b.thisPSM.PeptidesWithSetModifications.Count == 1 && !b.IsDecoy).GroupBy(b => b.thisPSM.PeptidesWithSetModifications.First().Protein).ToDictionary(b => b.Key);

                foreach (var protein in proteinList)
                {
                    if (!protein.IsDecoy)
                    {
                        var modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>(goodPsmsForEachProtein[protein].SelectMany(b => b.thisPSM.PeptidesWithSetModifications.First().allModsOneIsNterminus.Select(c => new Tuple<int, ModificationWithMass>(c.Key - b.thisPSM.PeptidesWithSetModifications.First().OneBasedStartResidueInProtein + 1, c.Value))));

                        IDictionary<int, List<Modification>> modsToWrite = new Dictionary<int, List<Modification>>();
                        foreach (var modd in protein.OneBasedPossibleLocalizedModifications)
                            foreach (var mod in modd.Value)
                            {
                                if (modificationsToAlwaysKeep.Contains(mod as ModificationWithLocation)
                                    || modsObservedOnThisProtein.Contains(new Tuple<int, ModificationWithMass>(modd.Key, mod as ModificationWithMass)))
                                {
                                    if (!modsToWrite.ContainsKey(modd.Key))
                                        modsToWrite.Add(modd.Key, new List<Modification> { mod });
                                    else
                                        modsToWrite[modd.Key].Add(mod);
                                }
                            }
                        protein.OneBasedPossibleLocalizedModifications.Clear();
                        foreach (var kvp in modsToWrite)
                            protein.OneBasedPossibleLocalizedModifications.Add(kvp);
                    }
                }

                if (dbFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "pruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                    SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });
                }
                if (dbFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "pruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), proteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                    SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });
                }
            }

            return mySearchTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private static IEnumerable<Type> GetSubclassesAndItself(Type type)
        {
            yield return type;
        }

        private static bool SameSettings(string pathToOldParamsFile, IndexingEngine indexEngine)
        {
            using (StreamReader reader = new StreamReader(pathToOldParamsFile))
                if (reader.ReadToEnd().Equals(indexEngine.ToString()))
                    return true;
            return false;
        }

        private string GenerateOutputFolderForIndices(List<DbForTask> dbFilenameList)
        {
            var folder = Path.Combine(Path.GetDirectoryName(dbFilenameList.First().FileName), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
            if (!Directory.Exists(folder))
                Directory.CreateDirectory(folder);
            return folder;
        }

        private void WriteIndexEngineParams(IndexingEngine indexEngine, string fileName, string taskId)
        {
            using (StreamWriter output = new StreamWriter(fileName))
            {
                output.Write(indexEngine);
            }
            SucessfullyFinishedWritingFile(fileName, new List<string> { taskId });
        }

        private string GetExistingFolderWithIndices(IndexingEngine indexEngine, List<DbForTask> dbFilenameList)
        {
            // In every database location...
            foreach (var ok in dbFilenameList)
            {
                var baseDir = Path.GetDirectoryName(ok.FileName);
                var directory = new DirectoryInfo(baseDir);
                DirectoryInfo[] directories = directory.GetDirectories();

                // Look at every subdirectory...
                foreach (DirectoryInfo possibleFolder in directories)
                {
                    if (File.Exists(Path.Combine(possibleFolder.FullName, "indexEngine.params")) &&
                        File.Exists(Path.Combine(possibleFolder.FullName, "peptideIndex.ind")) &&
                        File.Exists(Path.Combine(possibleFolder.FullName, "fragmentIndex.ind")) &&
                        SameSettings(Path.Combine(possibleFolder.FullName, "indexEngine.params"), indexEngine))
                        return possibleFolder.FullName;
                }
            }
            return null;
        }

        private void WriteFragmentIndexNetSerializer(Dictionary<float, List<int>> fragmentIndex, string fragmentIndexFile, string taskId)
        {
            var messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFile))
                ser.Serialize(file, fragmentIndex);
            SucessfullyFinishedWritingFile(fragmentIndexFile, new List<string> { taskId });
        }

        private void WritePeptideIndex(List<CompactPeptide> peptideIndex, string peptideIndexFile, string taskId)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(peptideIndexFile))
            {
                ser.Serialize(file, peptideIndex);
            }

            SucessfullyFinishedWritingFile(peptideIndexFile, new List<string> { taskId });
        }

        #endregion Private Methods

    }
}