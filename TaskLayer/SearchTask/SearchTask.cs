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
            SearchDecoy = true;
            DoHistogramAnalysis = false;
            MaxMissedCleavages = 2;
            Protease = ProteaseDictionary.Instance["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            BIons = true;
            YIons = true;
            ZdotIons = false;
            CIons = false;

            ListOfModListsFixed = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("f.txt")) };
            ListOfModListsVariable = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("v.txt")) };
            ListOfModListsLocalize = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("ptmlist.txt")) };

            SearchModes = SearchModesKnown.Take(1).ToList();
            TaskType = MyTask.Search;
            MaxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ModList> ListOfModListsFixed { get; set; }
        public List<ModList> ListOfModListsVariable { get; set; }
        public List<ModList> ListOfModListsLocalize { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public bool ClassicSearch { get; set; }
        public bool DoParsimony { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchDecoy { get; set; }
        public List<SearchMode> SearchModes { get; set; }

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
                sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModListsFixed.Select(b => b.FileName)));
                sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModListsVariable.Select(b => b.FileName)));
                sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModListsLocalize.Select(b => b.FileName)));
                sb.AppendLine("searchDecoy: " + SearchDecoy);
                sb.AppendLine("productMassTolerance: " + ProductMassTolerance);
                sb.AppendLine("searchModes: ");
                sb.Append(string.Join(Environment.NewLine, SearchModes.Select(b => "\t" + b.FileNameAddition)));
                return sb.ToString();
            }
        }

        #endregion Protected Properties

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            var mySearchTaskResults = new MySearchTaskResults(this);
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...");
            List<ModificationWithMass> variableModifications = ListOfModListsVariable.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> fixedModifications = ListOfModListsFixed.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> localizeableModifications = ListOfModListsLocalize.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();

            List<PsmParent>[] allPsms = new List<PsmParent>[SearchModes.Count];
            for (int j = 0; j < SearchModes.Count; j++)
                allPsms[j] = new List<PsmParent>();

            Status("Loading proteins...");
            Dictionary<string, Modification> um;
            var proteinList = dbFilenameList.SelectMany(b => ProteinDbLoader.LoadProteinDb(b.FileName, true, GetDict(localizeableModifications), b.IsContaminant, out um)).ToList();

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
                Status("Getting fragment dictionary...");
                var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, InitiatorMethionineBehavior, MaxMissedCleavages, MaxModificationIsoforms, lp);
                string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine);

                Dictionary<float, List<int>> fragmentIndexDict;
                if (pathToFolderWithIndices == null)
                {
                    Status("Generating indices...");
                    var output_folderForIndices = GenerateOutputFolderForIndices();
                    Status("Writing params...");
                    writeIndexEngineParams(indexEngine, Path.Combine(output_folderForIndices, "indexEngine.params"));

                    var indexResults = (IndexingResults)indexEngine.Run();
                    mySearchTaskResults.AddResultText(indexResults);
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
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile, MaxNumPeaksPerScan);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, MaxNumPeaksPerScan);
                Status("Opening spectra file...");
                myMsDataFile.Open();

                if (ClassicSearch)
                {
                    var classicSearchResults = (ClassicSearchResults)new ClassicSearchEngine(GetMs2Scans(myMsDataFile).OrderBy(b => b.PrecursorMass).ToArray(), myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, SearchModes, MaxMissedCleavages, MaxModificationIsoforms, myMsDataFile.Name, lp).Run();
                    mySearchTaskResults.AddResultText(classicSearchResults);
                    for (int i = 0; i < SearchModes.Count; i++)
                        allPsms[i].AddRange(classicSearchResults.OuterPsms[i]);
                    var analysisResults = new AnalysisEngine(classicSearchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups"), DoParsimony, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior).Run();
                    mySearchTaskResults.AddResultText(analysisResults);
                }
                else
                {
                    var modernSearchResults = (ModernSearchResults)new ModernSearchEngine(myMsDataFile, peptideIndex, keys, fragmentIndex, ProductMassTolerance, SearchModes).Run();
                    mySearchTaskResults.AddResultText(modernSearchResults);
                    for (int i = 0; i < SearchModes.Count; i++)
                        allPsms[i].AddRange(modernSearchResults.NewPsms[i]);
                    var analysisResults = new AnalysisEngine(modernSearchResults.NewPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups"), DoParsimony, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior).Run();
                    mySearchTaskResults.AddResultText(analysisResults);
                }
            }

            if (currentRawFileList.Count > 1)
            {
                var analysisResults = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, "aggregate" + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, "aggregate_ProteinGroups" + s), DoParsimony, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior).Run();
                mySearchTaskResults.AddResultText(analysisResults);
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

        private string GenerateOutputFolderForIndices()
        {
            var folder = Path.Combine(Path.GetDirectoryName(dbFilenameList.First().FileName), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
            if (!Directory.Exists(folder))
                Directory.CreateDirectory(folder);
            return folder;
        }

        private void writeIndexEngineParams(IndexingEngine indexEngine, string fileName)
        {
            using (StreamWriter output = new StreamWriter(fileName))
            {
                output.Write(indexEngine);
            }
            SucessfullyFinishedWritingFile(fileName);
        }

        private string GetExistingFolderWithIndices(IndexingEngine indexEngine)
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