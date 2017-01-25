using InternalLogicEngineLayer;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace InternalLogicTaskLayer
{
    public class SearchTask : MyTaskEngine
    {

        #region Private Fields

        private const double binTolInDaltons = 0.003;

        #endregion Private Fields

        #region Public Constructors

        public SearchTask(IEnumerable<ModList> modList, IEnumerable<SearchMode> inputSearchModes)
        {
            // Set default values here:
            ClassicSearch = true;
            DoParsimony = false;
            SearchDecoy = true;
            DoHistogramAnalysis = false;
            MaxMissedCleavages = 2;
            Protease = ProteaseDictionary.Instance["trypsin (no proline rule)"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            BIons = true;
            YIons = true;
            ListOfModListsForSearch = new List<ModListForSearchTask>();
            foreach (var uu in modList)
                ListOfModListsForSearch.Add(new ModListForSearchTask(uu));
            ListOfModListsForSearch[0].Fixed = true;
            ListOfModListsForSearch[1].Variable = true;
            ListOfModListsForSearch[2].Localize = true;

            SearchModes = new List<SearchModeFoSearch>();
            foreach (var uu in inputSearchModes)
                SearchModes.Add(new SearchModeFoSearch(uu));
            SearchModes[0].Use = true;
            TaskType = MyTask.Search;
            MaxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public Tolerance ProductMassTolerance { get; set; }
        public bool ClassicSearch { get; set; }
        public bool DoParsimony { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public List<ModListForSearchTask> ListOfModListsForSearch { get; set; }
        public bool SearchDecoy { get; set; }
        public List<SearchModeFoSearch> SearchModes { get; set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string SpecificTaskInfo
        {
            get
            {
                var sb = new StringBuilder();
                sb.AppendLine("classicSearch: " + ClassicSearch);
                sb.AppendLine("doParsimony: " + DoParsimony);
                sb.AppendLine("doHistogramAnalysis: " + DoHistogramAnalysis);
                sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModListsForSearch.Where(b => b.Fixed).Select(b => b.FileName)));
                sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModListsForSearch.Where(b => b.Variable).Select(b => b.FileName)));
                sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModListsForSearch.Where(b => b.Localize).Select(b => b.FileName)));
                sb.AppendLine("searchDecoy: " + SearchDecoy);
                sb.AppendLine("productMassTolerance: " + ProductMassTolerance);
                sb.AppendLine("searchModes: ");
                sb.Append(string.Join(Environment.NewLine, SearchModes.Where(b => b.Use).Select(b => "\t" + b.SearchMode.FileNameAddition)));
                return sb.ToString();
            }
        }

        #endregion Protected Properties

        #region Public Methods

        public static IEnumerable<Type> GetSubclassesAndItself(Type type)
        {
            foreach (var ok in type.Assembly.GetTypes().Where(t => t.IsSubclassOf(type)))
                yield return ok;
            yield return type;
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...");
            List<MorpheusModification> variableModifications = ListOfModListsForSearch.Where(b => b.Variable).SelectMany(b => b.Mods).ToList();
            List<MorpheusModification> fixedModifications = ListOfModListsForSearch.Where(b => b.Fixed).SelectMany(b => b.Mods).ToList();
            List<MorpheusModification> localizeableModifications = ListOfModListsForSearch.Where(b => b.Localize).SelectMany(b => b.Mods).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            MatchXMLmodsToKnownMods(xmlDbFilenameList, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            List<SearchMode> searchModesS = SearchModes.Where(b => b.Use).Select(b => b.SearchMode).ToList();

            List<ParentSpectrumMatch>[] allPsms = new List<ParentSpectrumMatch>[searchModesS.Count];
            for (int j = 0; j < searchModesS.Count; j++)
                allPsms[j] = new List<ParentSpectrumMatch>();

            Status("Loading proteins...");
            var proteinList = xmlDbFilenameList.SelectMany(b => GetProteins(SearchDecoy, identifiedModsInXML, b)).ToList();

            List<CompactPeptide> peptideIndex = null;
            Dictionary<float, List<int>> fragmentIndexDict = null;
            float[] keys = null;
            List<int>[] fragmentIndex = null;
            List<ProductType> lp = new List<ProductType>();
            if (BIons)
                lp.Add(ProductType.B);
            if (YIons)
                lp.Add(ProductType.Y);

            if (!ClassicSearch)
            {
                Status("Getting fragment dictionary...");
                var indexEngine = new IndexEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, InitiatorMethionineBehavior, MaxMissedCleavages, MaxModificationIsoforms, lp);
                string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine);

                if (pathToFolderWithIndices == null)
                {
                    Status("Generating indices...");
                    var output_folderForIndices = GenerateOutputFolderForIndices();
                    Status("Writing params...");
                    writeIndexEngineParams(indexEngine, Path.Combine(output_folderForIndices, "indexEngine.params"));

                    var indexResults = (IndexResults)indexEngine.Run();
                    peptideIndex = indexResults.PeptideIndex;
                    fragmentIndexDict = indexResults.FragmentIndexDict;

                    Status("Writing peptide index...");
                    writePeptideIndex(peptideIndex, Path.Combine(output_folderForIndices, "peptideIndex.ind"));
                    Status("Writing fragment index...");
                    writeFragmentIndexNetSerializer(fragmentIndexDict, Path.Combine(output_folderForIndices, "fragmentIndex.ind"));
                }
                else
                {
                    Status("Reading peptide index...");
                    var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
                    var ser = new NetSerializer.Serializer(messageTypes);
                    using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "peptideIndex.ind")))
                        peptideIndex = (List<CompactPeptide>)ser.Deserialize(file);

                    Status("Reading fragment index...");
                    messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
                    ser = new NetSerializer.Serializer(messageTypes);
                    using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "fragmentIndex.ind")))
                        fragmentIndexDict = (Dictionary<float, List<int>>)ser.Deserialize(file);
                }
                keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
                fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();
            }

            var currentRawFileList = rawDataFilenameList;
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                Status("Loading spectra file...");
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile, MaxNumPeaksPerScan);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, MaxNumPeaksPerScan);
                Status("Opening spectra file...");
                myMsDataFile.Open();

                ClassicSearchEngine classicSearchEngine = null;
                ClassicSearchResults classicSearchResults = null;

                ModernSearchEngine modernSearchEngine = null;
                ModernSearchResults modernSearchResults = null;

                if (ClassicSearch)
                {
                    var listOfSortedms2Scans = GetMs2Scans(myMsDataFile).OrderBy(b => b.PrecursorMass).ToArray();
                    classicSearchEngine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModesS, MaxMissedCleavages, MaxModificationIsoforms, myMsDataFile.Name, lp);

                    classicSearchResults = (ClassicSearchResults)classicSearchEngine.Run();
                    for (int i = 0; i < searchModesS.Count; i++)
                        allPsms[i].AddRange(classicSearchResults.OuterPsms[i]);

                    AnalysisEngine analysisEngine = null;

                    analysisEngine = new AnalysisEngine(classicSearchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModesS, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups"), DoParsimony, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons);

                    analysisEngine.Run();
                }
                else
                {
                    modernSearchEngine = new ModernSearchEngine(myMsDataFile, peptideIndex, keys, fragmentIndex, ProductMassTolerance.Value, searchModesS, variableModifications, localizeableModifications);

                    modernSearchResults = (ModernSearchResults)modernSearchEngine.Run();
                    for (int i = 0; i < searchModesS.Count; i++)
                        allPsms[i].AddRange(modernSearchResults.NewPsms[i]);

                    var analysisEngine = new AnalysisEngine(modernSearchResults.NewPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModesS, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups"), DoParsimony, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons);

                    analysisEngine.Run();
                }
            }

            if (currentRawFileList.Count > 1)
            {
                var analysisEngine = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModesS, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, "aggregate" + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, "aggregate_ProteinGroups" + s), DoParsimony, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons);

                analysisEngine.Run();
            }
            return new MySearchTaskResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private static bool SameSettings(string pathToOldParamsFile, IndexEngine indexEngine)
        {
            using (StreamReader reader = new StreamReader(pathToOldParamsFile))
                if (reader.ReadToEnd().Equals(indexEngine.ToString()))
                    return true;
            return false;
        }

        private string GenerateOutputFolderForIndices()
        {
            var folder = Path.Combine(Path.GetDirectoryName(xmlDbFilenameList.First().FileName), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
            if (!Directory.Exists(folder))
                Directory.CreateDirectory(folder);
            return folder;
        }

        private void writeIndexEngineParams(IndexEngine indexEngine, string fileName)
        {
            using (StreamWriter output = new StreamWriter(fileName))
            {
                output.Write(indexEngine);
            }
            SucessfullyFinishedWritingFile(fileName);
        }

        private string GetExistingFolderWithIndices(IndexEngine indexEngine)
        {
            // In every database location...
            foreach (var ok in xmlDbFilenameList)
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

        private void writeFragmentIndexNetSerializer(Dictionary<float, List<int>> fragmentIndex, string fragmentIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFile))
                ser.Serialize(file, fragmentIndex);
            SucessfullyFinishedWritingFile(fragmentIndexFile);
        }

        private void writePeptideIndex(List<CompactPeptide> peptideIndex, string peptideIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(peptideIndexFile))
            {
                ser.Serialize(file, peptideIndex);
            }

            SucessfullyFinishedWritingFile(peptideIndexFile);
        }

        #endregion Private Methods

    }
}