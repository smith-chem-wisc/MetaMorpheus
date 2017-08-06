using EngineLayer;
using EngineLayer.CrosslinkAnalysis;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    public partial class XLSearchTask : MetaMorpheusTask
    {

        #region Private Fields

        private const double binTolInDaltons = 0.003;

        private FlashLFQEngine FlashLfqEngine;

        #endregion Private Fields

        #region Public Constructors

        public XLSearchTask() : base(MyTask.XLSearch)
        {
            // Set default values here:
            DoParsimony = false;
            NoOneHitWonders = false;
            ModPeptidesAreUnique = false;
            DoQuantification = false;
            QuantifyPpmTol = 5;
            SearchDecoy = true;
            DoHistogramAnalysis = false;
            MaxMissedCleavages = 2;
            MinPeptideLength = null;
            MaxPeptideLength = null;
            Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new PpmTolerance(10);
            BIons = true;
            YIons = true;
            ZdotIons = false;
            CIons = false;
            FlashLfqEngine = new FlashLFQEngine();

            LocalizeAll = true;
            DoLocalizationAnalysis = true;

            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();

            WritePrunedDatabase = false;
            KeepAllUniprotMods = true;

            ConserveMemory = false;
            MaxDegreeOfParallelism = null;
            CrosslinkerType = CrosslinkerType.DSS;
            CrosslinkSearchTopNum = 50;
            CrosslinkSearchWithAllBeta = false;
            UdXLkerName = null;
            UdXLkerCleavable = false;
            UdXLkerShortMass = null;
            UdXLkerLongMass = null;
            UdXLkerTotalMass = null;
            UdXLkerResidue = 'K';
            XLprecusorMsTl = new PpmTolerance(10);
            XLBetaPrecusorMsTl = new PpmTolerance(10);

            // Deconvolution stuff
            DoPrecursorDeconvolution = false;
            UseProvidedPrecursorInfo = true;
            DeconvolutionIntensityRatio = 4;
            DeconvolutionMaxAssumedChargeState = 10;
            DeconvolutionMassTolerance = new PpmTolerance(5);
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
        public Tolerance ProductMassTolerance { get; set; }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreUnique { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchDecoy { get; set; }

        public bool ConserveMemory { get; set; }

        public bool WritePrunedDatabase { get; set; }
        public bool KeepAllUniprotMods { get; set; }

        public bool DoLocalizationAnalysis { get; set; }
        public bool DoQuantification { get; set; }

        public CrosslinkerType CrosslinkerType { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public bool CrosslinkSearchWithAllBeta { get; set; }
        public string UdXLkerName { get; set; }
        public bool UdXLkerCleavable { get; set; }
        public double? UdXLkerTotalMass { get; set; }
        public double? UdXLkerShortMass { get; set; }
        public double? UdXLkerLongMass { get; set; }
        public char UdXLkerResidue { get; set; }
        public Tolerance XLprecusorMsTl { get; set; }
        public Tolerance XLBetaPrecusorMsTl { get; set; }

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

            sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModsFixed));
            sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModsVariable));
            sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModsLocalize));
            sb.AppendLine("searchDecoy: " + SearchDecoy);
            sb.AppendLine("productMassTolerance: " + ProductMassTolerance);

            sb.AppendLine("Crosslink Precusor mass tolerance: " + XLprecusorMsTl);
            sb.AppendLine("Beta Precusor mass tolerance: " + XLBetaPrecusorMsTl);

            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId)
        {
            myTaskResults = new MyTaskResults(this);
            List<Tuple<PsmCross, PsmCross>> allPsmsXLTuple = new List<Tuple<PsmCross, PsmCross>>();
            List<PsmCross> allPsmsXL = new List<PsmCross>();
            var compactPeptideToProteinPeptideMatch = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...", taskId);

            #region Load modifications

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

            #endregion Load modifications

            Status("Loading proteins...", new List<string> { taskId });
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, SearchDecoy, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> unknownModifications)).ToList();

            List<ProductType> ionTypes = new List<ProductType>();
            if (BIons)
                ionTypes.Add(ProductType.B);
            if (YIons)
                ionTypes.Add(ProductType.Y);
            if (ZdotIons)
                ionTypes.Add(ProductType.Zdot);
            if (CIons)
                ionTypes.Add(ProductType.C);

            List<CompactPeptide> peptideIndex = null;
            float[] keys = null;
            List<int>[] fragmentIndex = null;

            #region Generate indices for modern search

            Status("Getting fragment dictionary...", new List<string> { taskId });
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, Protease, InitiatorMethionineBehavior, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, ionTypes, new List<string> { taskId });
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

            #endregion Generate indices for modern search

            object lock2 = new object();
            MyFileManager myFileManager = new MyFileManager(true);
            Status("Searching files...", taskId);
            ParallelOptions parallelOptions = new ParallelOptions();
            if (MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = MaxDegreeOfParallelism.Value;
            double completedFiles = 0;

            var crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(CrosslinkerType);
            if (CrosslinkerType == CrosslinkerType.UserDefined)
            {
                crosslinker.CrosslinkerName = UdXLkerName;
                crosslinker.Cleavable = UdXLkerCleavable;
                crosslinker.TotalMass = (double)UdXLkerTotalMass;
                crosslinker.CleaveMassShort = (double)UdXLkerShortMass;
                crosslinker.CleaveMassLong = (double)UdXLkerLongMass;
                crosslinker.CrosslinkerModSite = UdXLkerResidue;
            }

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var origDataFile = currentRawFileList[spectraFileIndex];

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile);

                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

                Status("Starting search...", thisId);

                CrosslinkSearchResults xlsearchResults;
                if (CrosslinkSearchWithAllBeta)
                    xlsearchResults = ((CrosslinkSearchResults)new CrosslinkSearchEngine2(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, crosslinker, CrosslinkSearchTopNum, CrosslinkSearchWithAllBeta, XLprecusorMsTl, XLBetaPrecusorMsTl, modsDictionary, ionTypes, proteinList, Protease, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, variableModifications, fixedModifications, MaxModificationIsoforms, thisId).Run());
                else
                    xlsearchResults = ((CrosslinkSearchResults)new CrosslinkSearchEngine(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, crosslinker, CrosslinkSearchTopNum, CrosslinkSearchWithAllBeta, XLprecusorMsTl, XLBetaPrecusorMsTl, ionTypes, modsDictionary, thisId).Run());

                myFileManager.DoneWithFile(origDataFile);

                lock (lock2)
                {
                    allPsmsXLTuple.AddRange(xlsearchResults.NewPsms);
                    foreach (var item in xlsearchResults.NewPsms)
                    {
                        allPsmsXL.Add(item.Item1); allPsmsXL.Add(item.Item2);
                    }
                }
                ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                completedFiles++;
                ReportProgress(new ProgressEventArgs((int)completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }
            );
            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            Status("Crosslink analysis engine", taskId);
            MetaMorpheusEngineResults allcrosslinkanalysisResults;
            allcrosslinkanalysisResults = new CrosslinkAnalysisEngine(allPsmsXLTuple, compactPeptideToProteinPeptideMatch, proteinList, variableModifications, fixedModifications, Protease, null, ProductMassTolerance, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, ionTypes, InitiatorMethionineBehavior, modsDictionary, null, OutputFolder, crosslinker, new List<string> { taskId }).Run();
            allPsmsXLTuple = allPsmsXLTuple.OrderByDescending(p => p.Item1.XLTotalScore).ToList();

            //Write All Psms and pepxml
            var allPsmsXLTupleFDR = CrosslinkDoFalseDiscoveryRateAnalysis(allPsmsXLTuple, new OpenSearchMode());
            WriteCrosslinkToTsv(allPsmsXLTupleFDR, OutputFolder, "xl_fdr", new List<string> { taskId });
            if (allPsmsXLTupleFDR.Count != 0)
            {
                foreach (var fullFilePath in currentRawFileList)
                {
                    string fileNameNoExtension = Path.GetFileNameWithoutExtension(fullFilePath);
                    WritePepXML(allPsmsXLTupleFDR.Where(p => p.Item1.FullFilePath == fullFilePath).ToList(), dbFilenameList, variableModifications, fixedModifications, localizeableModifications, OutputFolder, fileNameNoExtension, new List<string> { taskId });
                }
            }

            //Write Inter Psms FDR
            var interPsmsXLTupleFDR = allPsmsXLTuple.Where(p => p.Item1.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() != p.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).ToList();
            interPsmsXLTupleFDR = CrosslinkDoFalseDiscoveryRateAnalysis(interPsmsXLTupleFDR, new OpenSearchMode()).Where(p => p.Item1.MostProbableProteinInfo.IsDecoy != true && p.Item2.MostProbableProteinInfo.IsDecoy != true && p.Item1.FdrInfo.QValue <= 0.01).ToList();
            WriteCrosslinkToTsv(interPsmsXLTupleFDR, OutputFolder, "xl_inter_fdr", new List<string> { taskId });
            WriteCrosslinkToTxtForCLMSVault(interPsmsXLTupleFDR, OutputFolder, "xl_inter_fdr_CLMSVault", crosslinker, new List<string> { taskId });

            var interPsmsXLTuple = allPsmsXLTuple.Where(p => p.Item1.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() != p.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).OrderBy(p => p.Item1.ScanNumber).ToList();
            WriteCrosslinkToTxtForPercolator(interPsmsXLTuple, OutputFolder, "xl_inter_perc", crosslinker, new List<string> { taskId });

            //Write Intra Psms FDR
            var intraPsmsXLTupleFDR = allPsmsXLTuple.Where(p => p.Item1.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() == p.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).ToList();
            intraPsmsXLTupleFDR = CrosslinkDoFalseDiscoveryRateAnalysis(intraPsmsXLTupleFDR, new OpenSearchMode()).Where(p => p.Item1.MostProbableProteinInfo.IsDecoy != true && p.Item2.MostProbableProteinInfo.IsDecoy != true && p.Item1.FdrInfo.QValue <= 0.01).ToList();
            WriteCrosslinkToTsv(intraPsmsXLTupleFDR, OutputFolder, "xl_intra_fdr", new List<string> { taskId });
            WriteCrosslinkToTxtForCLMSVault(intraPsmsXLTupleFDR, OutputFolder, "xl_intra_fdr_CLMSVault", crosslinker, new List<string> { taskId });

            var intraPsmsXLTuple = allPsmsXLTuple.Where(p => p.Item1.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() == p.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).OrderBy(p => p.Item1.ScanNumber).ToList();
            WriteCrosslinkToTxtForPercolator(intraPsmsXLTuple, OutputFolder, "xl_intra_perc", crosslinker, new List<string> { taskId });

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

        //Calculate the FDR of crosslinked peptide FP/(FP+TP)
        private static List<Tuple<PsmCross, PsmCross>> CrosslinkDoFalseDiscoveryRateAnalysis(List<Tuple<PsmCross, PsmCross>> items, MassDiffAcceptor sm)
        {
            var ids = new List<Tuple<PsmCross, PsmCross>>();
            foreach (var item in items)
            {
                ids.Add(new Tuple<PsmCross, PsmCross>(item.Item1, item.Item2));
            }

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            int[] cumulative_target_per_notch = new int[sm.NumNotches];
            int[] cumulative_decoy_per_notch = new int[sm.NumNotches];

            for (int i = 0; i < ids.Count; i++)
            {
                var item1 = ids[i].Item1; var item2 = ids[i].Item2;

                var isDecoy1 = item1.MostProbableProteinInfo.IsDecoy; var isDecoy2 = item1.MostProbableProteinInfo.IsDecoy;
                if (isDecoy1 || isDecoy2)
                    cumulative_decoy++;
                else
                    cumulative_target++;

                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                item1.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
                item2.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
            }

            double min_q_value = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                PsmCross id = ids[i].Item1;
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;
            }

            return ids;
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
            var folder = Path.Combine(Path.GetDirectoryName(dbFilenameList.First().FilePath), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
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