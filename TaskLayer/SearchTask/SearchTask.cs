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
            Quantify = false;
            QuantifyRtTol = 0.5;
            QuantifyPpmTol = 10;
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

            SearchModes = MyEngine.SearchModesKnown.Take(1).ToList();
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
        public bool NoOneHitWonders { get; set; }
        public bool Quantify { get; set; }
        public double QuantifyRtTol { get; set; }
        public double QuantifyPpmTol { get; set; }
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
                sb.AppendLine("requireTwoPeptidesToIdProtein: " + NoOneHitWonders);
                sb.AppendLine("quantify: " + Quantify);
                if(Quantify)
                {
                    sb.AppendLine("quantify rt tolerance: " + QuantifyRtTol);
                    sb.AppendLine("quantify ppm tolerance: " + QuantifyPpmTol);
                }
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

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId)
        {
            var mySearchTaskResults = new MySearchTaskResults(this);
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = ListOfModListsVariable.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> fixedModifications = ListOfModListsFixed.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> localizeableModifications = ListOfModListsLocalize.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();

            List<PsmParent>[] allPsms = new List<PsmParent>[SearchModes.Count];
            for (int j = 0; j < SearchModes.Count; j++)
                allPsms[j] = new List<PsmParent>();

            Status("Loading proteins...", new List<string> { taskId });
            Dictionary<string, Modification> um;
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FileName, true, localizeableModifications, b.IsContaminant, out um)).ToList();

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
                var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, InitiatorMethionineBehavior, MaxMissedCleavages, MaxModificationIsoforms, lp, new List<string> { taskId });
                string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);

                Dictionary<float, List<int>> fragmentIndexDict;
                if (pathToFolderWithIndices == null)
                {
                    Status("Generating indices...", new List<string> { taskId });
                    var output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
                    Status("Writing params...", new List<string> { taskId });
                    writeIndexEngineParams(indexEngine, Path.Combine(output_folderForIndices, "indexEngine.params"), taskId);

                    var indexResults = (IndexingResults)indexEngine.Run();
                    mySearchTaskResults.AddResultText(indexResults);
                    peptideIndex = indexResults.PeptideIndex;
                    fragmentIndexDict = indexResults.FragmentIndexDict;

                    Status("Writing peptide index...", new List<string> { taskId });
                    writePeptideIndex(peptideIndex, Path.Combine(output_folderForIndices, "peptideIndex.ind"), taskId);
                    Status("Writing fragment index...", new List<string> { taskId });
                    writeFragmentIndexNetSerializer(fragmentIndexDict, Path.Combine(output_folderForIndices, "fragmentIndex.ind"), taskId);
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

            Status("Searching files...", new List<string> { taskId });
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
                    var classicSearchResults = (ClassicSearchResults)new ClassicSearchEngine(MyEngine.GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray(), myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, SearchModes, MaxMissedCleavages, MaxModificationIsoforms, origDataFile, lp, new List<string> { taskId, "Individual Searches", origDataFile }).Run();
                    mySearchTaskResults.AddResultText(classicSearchResults);
                    for (int i = 0; i < SearchModes.Count; i++)
                        allPsms[i].AddRange(classicSearchResults.OuterPsms[i]);
                    var analysisResults = new AnalysisEngine(classicSearchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups"), DoParsimony, NoOneHitWonders, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior, Quantify, QuantifyRtTol, QuantifyPpmTol).Run();
                    var analysisResults = new AnalysisEngine(classicSearchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups", new List<string> { taskId, "Individual Searches", origDataFile }), DoParsimony, NoOneHitWonders, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFile }).Run();
                    mySearchTaskResults.AddResultText(analysisResults);
                }
                else
                {
                    var modernSearchResults = (ModernSearchResults)new ModernSearchEngine(myMsDataFile, peptideIndex, keys, fragmentIndex, ProductMassTolerance, SearchModes, origDataFile, new List<string> { taskId, "Individual Searches", origDataFile }).Run();
                    mySearchTaskResults.AddResultText(modernSearchResults);
                    for (int i = 0; i < SearchModes.Count; i++)
                        allPsms[i].AddRange(modernSearchResults.NewPsms[i]);
                    var analysisResults = new AnalysisEngine(modernSearchResults.NewPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups", new List<string> { taskId, "Individual Searches", origDataFile }), DoParsimony, NoOneHitWonders, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFile }).Run();
                    var analysisResults = new AnalysisEngine(modernSearchResults.NewPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups"), DoParsimony, NoOneHitWonders, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior, Quantify, QuantifyRtTol, QuantifyPpmTol).Run();
                    mySearchTaskResults.AddResultText(analysisResults);
                }
                Status("Done!", new List<string> { taskId, "Individual Searches", origDataFile });
            }

            if (currentRawFileList.Count > 1)
            {
                var analysisResults = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, "aggregate" + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, "aggregate_ProteinGroups" + s), DoParsimony, NoOneHitWonders, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior, Quantify, QuantifyRtTol, QuantifyPpmTol).Run();
                var analysisResults = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, SearchModes, null, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate" + s, new List<string> { taskId }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, "aggregate" + s, new List<string> { taskId }), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTsv(h, OutputFolder, "aggregate_ProteinGroups" + s, new List<string> { taskId }), DoParsimony, NoOneHitWonders, MaxMissedCleavages, MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior, new List<string> { taskId }).Run();
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

        private string GenerateOutputFolderForIndices(List<DbForTask> dbFilenameList)
        {
            var folder = Path.Combine(Path.GetDirectoryName(dbFilenameList.First().FileName), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
            if (!Directory.Exists(folder))
                Directory.CreateDirectory(folder);
            return folder;
        }

        private void writeIndexEngineParams(IndexingEngine indexEngine, string fileName, string taskId)
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

        private void writeFragmentIndexNetSerializer(Dictionary<float, List<int>> fragmentIndex, string fragmentIndexFile, string taskId)
        {
            var messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFile))
                ser.Serialize(file, fragmentIndex);
            SucessfullyFinishedWritingFile(fragmentIndexFile, new List<string> { taskId });
        }

        private void writePeptideIndex(List<CompactPeptide> peptideIndex, string peptideIndexFile, string taskId)
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