using EngineLayer;
using EngineLayer.Indexing;
using EngineLayer.CrosslinkSearch;
using EngineLayer.CrosslinkAnalysis;
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
    public class XLSearchTask : MetaMorpheusTask
    {

        #region Private Fields

        private const double binTolInDaltons = 0.003;

        #endregion Private Fields

        #region Public Constructors

        public XLSearchTask() : base(MyTask.XLSearch)
        {
            // Set default values here:
            //Add crosslink para
            CrosslinkSearch = true;
            crosslinkerType = CrosslinkerType.DSSO;
            CrosslinkSearchTopNum = 50;
            CrosslinkSearchWithAllBeta = false;
            UdXLkerName = null;
            UdXLkerCleavable = false;
            UdXLkerShortMass = null;
            UdXLkerLongMass = null;
            UdXLkerTotalMass = null;
            UdXLkerResidue = null;
            XLprecusorMsTl = new Tolerance(ToleranceUnit.Absolute, 0.01);

            SearchDecoy = true;

            MaxMissedCleavages = 2;
            MinPeptideLength = null;
            MaxPeptideLength = null;
            Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            BIons = true;
            YIons = true;
            ZdotIons = false;
            CIons = false;

            LocalizeAll = true;

            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();

            KeepAllUniprotMods = true;

            SearchModes = GlobalTaskLevelSettings.SearchModesKnown.Take(1).ToList();

            UseProvidedPrecursorInfo = true;
            FindAllPrecursors = false;

        }

        #endregion Public Constructors

        #region Public Properties

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }

        public int MaxMissedCleavages { get; set; }

        public int? MinPeptideLength { get; set; }

        public int? MaxPeptideLength { get; set; }

        public int MaxModificationIsoforms { get; set; }

        public Protease Protease { get; set; }

        public bool BIons { get; set; }

        public bool YIons { get; set; }

        public bool ZdotIons { get; set; }

        public bool CIons { get; set; }
        public List<Tuple<string, string>> ListOfModsFixed { get; set; }
        public List<Tuple<string, string>> ListOfModsVariable { get; set; }
        public List<Tuple<string, string>> ListOfModsLocalize { get; set; }
        public Tolerance ProductMassTolerance { get; set; }

        //Crosslink search para
        public bool CrosslinkSearch { get; set; }
        public CrosslinkerType crosslinkerType { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public bool CrosslinkSearchWithAllBeta { get; set; }
        public string UdXLkerName { get; set; }
        public bool UdXLkerCleavable { get; set; }
        public double? UdXLkerTotalMass { get; set; }
        public double? UdXLkerShortMass { get; set; }
        public double? UdXLkerLongMass { get; set; }
        public string UdXLkerResidue { get; set; }
        public Tolerance XLprecusorMsTl { get; set; }

        public bool SearchDecoy { get; set; }
        public List<MassDiffAcceptor> SearchModes { get; set; }
        public bool LocalizeAll { get; set; }
        public bool KeepAllUniprotMods { get; set; }
        public bool UseProvidedPrecursorInfo { get; set; }
        public bool FindAllPrecursors { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(TaskType.ToString());
            sb.AppendLine("The initiator methionine behavior is set to "
                + InitiatorMethionineBehavior
                + " and the maximum number of allowed missed cleavages is "
                + MaxMissedCleavages);
            sb.AppendLine("MinPeptideLength: " + MinPeptideLength);
            sb.AppendLine("MaxPeptideLength: " + MaxPeptideLength);
            sb.AppendLine("maxModificationIsoforms: " + MaxModificationIsoforms);
            sb.AppendLine("protease: " + Protease);
            sb.AppendLine("bIons: " + BIons);
            sb.AppendLine("yIons: " + YIons);
            sb.AppendLine("cIons: " + CIons);
            sb.AppendLine("zdotIons: " + ZdotIons);
            if (CrosslinkSearch)
            {
                sb.AppendLine("CrosslinkSearch: " + CrosslinkSearch);
                sb.AppendLine("CrosslinkerType: " + crosslinkerType);
                sb.AppendLine("Crosslink Search Top Number: " + CrosslinkSearchTopNum);
                sb.AppendLine("CrosslinkSearchWithCrosslinkerMod: " + CrosslinkSearchWithAllBeta);
                sb.AppendLine("Crosslink peptides precusor tolerance: " + XLprecusorMsTl);
                if (crosslinkerType == CrosslinkerType.UserDefined)
                {
                    sb.AppendLine("User Defined Crosslinker Name: " + UdXLkerName
                        + "Cleavable: " + UdXLkerCleavable
                        + "Total Mass: " + UdXLkerTotalMass
                        + "Short Mod Mass: " + UdXLkerShortMass
                        + "Long Mod Mass: " + UdXLkerLongMass
                        + "Crosslink Residue: " + UdXLkerResidue);
                }
            }
            sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModsFixed));
            sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModsVariable));
            sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModsLocalize));
            sb.AppendLine("searchDecoy: " + SearchDecoy);
            sb.AppendLine("productMassTolerance: " + ProductMassTolerance);
            sb.AppendLine("searchModes: ");
            sb.Append(string.Join(Environment.NewLine, SearchModes.Select(b => "\t" + b.FileNameAddition)));
            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId)
        {
            myTaskResults = new MyTaskResults(this);
            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            List<ModificationWithMass> localizeableModifications;
            if (LocalizeAll)
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                if (!modsDictionary.ContainsKey(mod))
                    modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            List<PsmParent>[] allPsms = new List<PsmParent>[SearchModes.Count];
            for (int j = 0; j < SearchModes.Count; j++)
                allPsms[j] = new List<PsmParent>();
            List<Tuple<PsmCross, PsmCross>> allPsmsTopTuple = new List<Tuple<PsmCross, PsmCross>>();
            //List<PsmParent>[] allPsmsTop = new List<PsmParent>[arrayOfSortedMS2Scans.length];

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
            var crosslinker = new CrosslinkerTypeClass();
            if (CrosslinkSearch)
            {
                Status("Getting fragment dictionary...", new List<string> { taskId });
                var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, modsDictionary, Protease, InitiatorMethionineBehavior, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, new List<string> { taskId });
                string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);

                Dictionary<float, List<int>> fragmentIndexDict;
                if (pathToFolderWithIndices == null)
                {
                    var output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
                    Status("Writing params...", new List<string> { taskId });
                    WriteIndexEngineParams(indexEngine, Path.Combine(output_folderForIndices, "indexEngine.params"), taskId);

                    Status("Running Index Engine...", new List<string> { taskId });
                    var indexResults = (IndexingResults)indexEngine.Run();
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

                Status("Getting ms2 scans...", new List<string> { taskId, "Individual Searches", origDataFile });
                var intensityRatio = 4;
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, FindAllPrecursors, UseProvidedPrecursorInfo, intensityRatio, origDataFile).OrderBy(b => b.PrecursorMass).ToArray();

                if (CrosslinkSearch && !CrosslinkSearchWithAllBeta)
                {
                    crosslinker.SelectCrosslinker(crosslinkerType);

                    if (crosslinkerType == CrosslinkerType.UserDefined)
                    {
                        crosslinker.CrosslinkerName = UdXLkerName;
                        crosslinker.Cleavable = UdXLkerCleavable;
                        crosslinker.TotalMass = (double)UdXLkerTotalMass;
                        crosslinker.CleaveMassShort = (double)UdXLkerShortMass;
                        crosslinker.CleaveMassLong = (double)UdXLkerLongMass;
                        crosslinker.CrosslinkerModSite = UdXLkerResidue;
                    }

                    var crosslinkSearchResults = (CrosslinkSearchResults)new CrosslinkSearchEngine(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, SearchModes, crosslinker, CrosslinkSearchTopNum, CrosslinkSearchWithAllBeta, XLprecusorMsTl, lp, modsDictionary, nestedIds).Run();

                    allPsmsTopTuple.AddRange(crosslinkSearchResults.NewPsms);
                    MetaMorpheusEngineResults crosslinkAnalysisResults;

                    crosslinkAnalysisResults = new CrosslinkAnalysisEngine(crosslinkSearchResults.NewPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, Protease, SearchModes, arrayOfMs2ScansSortedByMass, ProductMassTolerance, (List<NewPsmWithFdr> h, string s, List<string> ss) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s, ss), (List<Tuple<NewPsmWithFdr, NewPsmWithFdr>> h, string s, List<string> ss) => WriteCrosslinkToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s, ss), MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFile }, modsDictionary, myMsDataFile, OutputFolder, crosslinker, CrosslinkSearchWithAllBeta).Run();

                }
                else if (CrosslinkSearch && CrosslinkSearchWithAllBeta)
                {
                    crosslinker.SelectCrosslinker(crosslinkerType);

                    if (crosslinkerType == CrosslinkerType.UserDefined)
                    {
                        crosslinker.CrosslinkerName = UdXLkerName;
                        crosslinker.Cleavable = UdXLkerCleavable;
                        crosslinker.TotalMass = (double)UdXLkerTotalMass;
                        crosslinker.CleaveMassShort = (double)UdXLkerShortMass;
                        crosslinker.CleaveMassLong = (double)UdXLkerLongMass;
                        crosslinker.CrosslinkerModSite = UdXLkerResidue;
                    }

                    var crosslinkSearchResults2 = (CrosslinkSearchResults)new CrosslinkSearchEngine2(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, SearchModes, new List<string> { taskId, "Individual Searches", origDataFile }, crosslinker, CrosslinkSearchTopNum, CrosslinkSearchWithAllBeta, XLprecusorMsTl, modsDictionary, lp, proteinList, Protease, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, variableModifications, fixedModifications, MaxModificationIsoforms).Run();
                    allPsmsTopTuple.AddRange(crosslinkSearchResults2.NewPsms);
                    MetaMorpheusEngineResults crosslinkAnalysisResults2;

                    crosslinkAnalysisResults2 = new CrosslinkAnalysisEngine(crosslinkSearchResults2.NewPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, Protease, SearchModes, arrayOfMs2ScansSortedByMass, ProductMassTolerance, (List<NewPsmWithFdr> h, string s, List<string> ss) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s, ss), (List<Tuple<NewPsmWithFdr, NewPsmWithFdr>> h, string s, List<string> ss) => WriteCrosslinkToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s, ss), MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFile }, modsDictionary, myMsDataFile, OutputFolder, crosslinker, CrosslinkSearchWithAllBeta).Run();

                }

                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Searches" }));
                if (currentRawFileList.Count > 1 && CrosslinkSearch)
                {
                    MetaMorpheusEngineResults allcrosslinkanalysisResults;
                    allcrosslinkanalysisResults = new CrosslinkAnalysisEngine(allPsmsTopTuple, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, Protease, SearchModes, null, ProductMassTolerance, (List<NewPsmWithFdr> h, string s, List<string> ss) => WritePsmsToTsv(h, OutputFolder, s, ss), (List<Tuple<NewPsmWithFdr, NewPsmWithFdr>> h, string s, List<string> ss) => WriteCrosslinkToTsv(h, OutputFolder, s, ss), MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, initiatorMethionineBehavior, new List<string> { taskId }, modsDictionary, null, OutputFolder, crosslinker, CrosslinkSearchWithAllBeta).Run();
                }
            }
            return myTaskResults;
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

        private int GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications peptideWithSetModifications)
        {
            if (oneIsNterminus == 1)
                return peptideWithSetModifications.OneBasedStartResidueInProtein;
            if (oneIsNterminus == peptideWithSetModifications.Length + 2)
                return peptideWithSetModifications.OneBasedEndResidueInProtein;
            return peptideWithSetModifications.OneBasedStartResidueInProtein + oneIsNterminus - 2;
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