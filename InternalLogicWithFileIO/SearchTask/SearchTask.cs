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

        #region Public Constructors

        public SearchTask(IEnumerable<ModList> modList, IEnumerable<SearchMode> inputSearchModes)
        {
            // Set default values here:
            classicSearch = true;
            doParsimony = false;
            searchDecoy = true;
            maxMissedCleavages = 2;
            protease = ProteaseDictionary.Instance["trypsin (no proline rule)"];
            maxModificationIsoforms = 4096;
            initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            bIons = true;
            yIons = true;
            listOfModListsForSearch = new List<ModListForSearchTask>();
            foreach (var uu in modList)
                listOfModListsForSearch.Add(new ModListForSearchTask(uu));
            listOfModListsForSearch[0].Fixed = true;
            listOfModListsForSearch[1].Variable = true;
            listOfModListsForSearch[2].Localize = true;

            searchModes = new List<SearchModeFoSearch>();
            foreach (var uu in inputSearchModes)
                searchModes.Add(new SearchModeFoSearch(uu));
            searchModes[0].Use = true;
            taskType = MyTaskEnum.Search;
            maxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public Tolerance productMassTolerance { get; set; }
        public bool classicSearch { get; set; }
        public bool doParsimony { get; set; }
        public List<ModListForSearchTask> listOfModListsForSearch { get; set; }
        public bool searchDecoy { get; set; }
        public List<SearchModeFoSearch> searchModes { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static IEnumerable<Type> GetSubclassesAndItself(Type type)
        {
            foreach (var ok in type.Assembly.GetTypes().Where(t => t.IsSubclassOf(type)))
                yield return ok;
            yield return type;
        }

        #endregion Public Methods

        #region Protected Methods

        protected override string GetSpecificTaskInfo()
        {
            var sb = new StringBuilder();
            sb.AppendLine("classicSearch: " + classicSearch);
            sb.AppendLine("doParsimony: " + doParsimony);
            sb.AppendLine("Fixed mod lists: " + string.Join(",", listOfModListsForSearch.Where(b => b.Fixed).Select(b => b.FileName)));
            sb.AppendLine("Variable mod lists: " + string.Join(",", listOfModListsForSearch.Where(b => b.Variable).Select(b => b.FileName)));
            sb.AppendLine("Localized mod lists: " + string.Join(",", listOfModListsForSearch.Where(b => b.Localize).Select(b => b.FileName)));
            sb.AppendLine("searchDecoy: " + searchDecoy);
            sb.AppendLine("productMassTolerance: " + productMassTolerance);
            sb.AppendLine("searchModes: ");
            sb.Append(string.Join(Environment.NewLine, searchModes.Where(b => b.Use).Select(b => "\t" + b.sm)));
            return sb.ToString();
        }

        protected override MyResults RunSpecific()
        {
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            status("Loading modifications...");
            List<MorpheusModification> variableModifications = listOfModListsForSearch.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> fixedModifications = listOfModListsForSearch.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> localizeableModifications = listOfModListsForSearch.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            MatchXMLmodsToKnownMods(xmlDbFilenameList, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            List<SearchMode> searchModesS = searchModes.Where(b => b.Use).Select(b => b.sm).ToList();

            List<ParentSpectrumMatch>[] allPsms = new List<ParentSpectrumMatch>[searchModesS.Count];
            for (int j = 0; j < searchModesS.Count; j++)
                allPsms[j] = new List<ParentSpectrumMatch>();

            status("Loading proteins...");
            var proteinList = xmlDbFilenameList.SelectMany(b => getProteins(searchDecoy, identifiedModsInXML, b)).ToList();

            List<CompactPeptide> peptideIndex = null;
            Dictionary<float, List<int>> fragmentIndexDict = null;
            float[] keys = null;
            List<int>[] fragmentIndex = null;

            if (!classicSearch)
            {
                status("Getting fragment dictionary...");
                var indexEngine = new IndexEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, protease);
                string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine);

                if (pathToFolderWithIndices == null)
                {
                    status("Generating indices...");
                    var output_folderForIndices = GenerateOutputFolderForIndices();
                    status("Writing params...");
                    writeIndexEngineParams(indexEngine, Path.Combine(output_folderForIndices, "indexEngine.params"));

                    var indexResults = (IndexResults)indexEngine.Run();
                    peptideIndex = indexResults.peptideIndex;
                    fragmentIndexDict = indexResults.fragmentIndexDict;

                    status("Writing peptide index...");
                    writePeptideIndex(peptideIndex, Path.Combine(output_folderForIndices, "peptideIndex.ind"));
                    status("Writing fragment index...");
                    writeFragmentIndexNetSerializer(fragmentIndexDict, Path.Combine(output_folderForIndices, "fragmentIndex.ind"));
                }
                else
                {
                    status("Reading peptide index...");
                    var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
                    var ser = new NetSerializer.Serializer(messageTypes);
                    using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "peptideIndex.ind")))
                        peptideIndex = (List<CompactPeptide>)ser.Deserialize(file);

                    status("Reading fragment index...");
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
                status("Loading spectra file...");
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile, maxNumPeaksPerScan);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, maxNumPeaksPerScan);
                status("Opening spectra file...");
                myMsDataFile.Open();

                ClassicSearchEngine classicSearchEngine = null;
                ClassicSearchResults classicSearchResults = null;

                ModernSearchEngine modernSearchEngine = null;
                ModernSearchResults modernSearchResults = null;

                if (classicSearch)
                {
                    var listOfSortedms2Scans = myMsDataFile.Where(b => b.MsnOrder == 2).Select(b => new LocalMs2Scan(b)).OrderBy(b => b.precursorMass).ToArray();

                    classicSearchEngine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, spectraFileIndex, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModesS);

                    classicSearchResults = (ClassicSearchResults)classicSearchEngine.Run();
                    for (int i = 0; i < searchModesS.Count; i++)
                        allPsms[i].AddRange(classicSearchResults.outerPsms[i]);

                    AnalysisEngine analysisEngine = null;

                    if (doParsimony)
                    {
                        analysisEngine = new AnalysisEngine(classicSearchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModesS, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFDR> h, string s) => WritePSMsToTSV(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTSV(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups"), doParsimony);
                    }
                    else
                        analysisEngine = new AnalysisEngine(classicSearchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModesS, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFDR> h, string s) => WritePSMsToTSV(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTSV(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups"), doParsimony);

                    analysisEngine.Run();
                }
                else
                {
                    modernSearchEngine = new ModernSearchEngine(myMsDataFile, spectraFileIndex, peptideIndex, keys, fragmentIndex, productMassTolerance.Value, searchModesS);

                    modernSearchResults = (ModernSearchResults)modernSearchEngine.Run();
                    for (int i = 0; i < searchModesS.Count; i++)
                        allPsms[i].AddRange(modernSearchResults.newPsms[i]);

                    var analysisEngine = new AnalysisEngine(modernSearchResults.newPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModesS, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFDR> h, string s) => WritePSMsToTSV(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTSV(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s + "_ProteinGroups"), doParsimony);

                    analysisEngine.Run();
                }
            }

            if (currentRawFileList.Count > 1)
            {
                var analysisEngine = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModesS, null, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, "aggregate"), (List<NewPsmWithFDR> h, string s) => WritePSMsToTSV(h, output_folder, "aggregate" + s), (List<ProteinGroup> h, string s) => WriteProteinGroupsToTSV(h, output_folder, "aggregate_ProteinGroups"), doParsimony);

                analysisEngine.Run();
            }
            return new MySearchTaskResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private string GenerateOutputFolderForIndices()
        {
            var folder = Path.Combine(Path.GetDirectoryName(xmlDbFilenameList.First()), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
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
                var baseDir = Path.GetDirectoryName(ok);
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

        private bool SameSettings(string pathToOldParamsFile, IndexEngine indexEngine)
        {
            using (StreamReader reader = new StreamReader(pathToOldParamsFile))
                if (reader.ReadToEnd().Equals(indexEngine.ToString()))
                    return true;
            return false;
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