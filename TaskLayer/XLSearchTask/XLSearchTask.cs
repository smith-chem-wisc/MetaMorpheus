using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.CrosslinkAnalysis;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
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
    public class XLSearchTask : MetaMorpheusTask
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
            ProductMassTolerance = new AbsoluteTolerance(0.01);
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

            MassDiffAcceptors = GlobalTaskLevelSettings.SearchModesKnown.Take(1).ToList();

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
            XLprecusorMsTl = new AbsoluteTolerance(0.01);

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
        public List<MassDiffAcceptor> MassDiffAcceptors { get; set; }
        public bool ConserveMemory { get; set; }

        public bool WritePrunedDatabase { get; set; }
        public bool KeepAllUniprotMods { get; set; }

        public bool DoLocalizationAnalysis { get; set; }
        public bool DoQuantification { get; set; }

        public SearchType SearchType { get; set; }

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
            sb.AppendLine("SearchType: " + SearchType);
            sb.AppendLine("doParsimony: " + DoParsimony);
            if (DoParsimony)
            {
                sb.AppendLine("modifiedPeptidesAreUnique: " + ModPeptidesAreUnique);
                sb.AppendLine("requireTwoPeptidesToIdProtein: " + NoOneHitWonders);
            }
            sb.AppendLine("quantify: " + DoQuantification);
            if (DoQuantification)
                sb.AppendLine("quantify ppm tolerance: " + QuantifyPpmTol);
            sb.AppendLine("doHistogramAnalysis: " + DoHistogramAnalysis);
            sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModsFixed));
            sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModsVariable));
            sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModsLocalize));
            sb.AppendLine("searchDecoy: " + SearchDecoy);
            sb.AppendLine("productMassTolerance: " + ProductMassTolerance);
            sb.AppendLine("searchModes: ");
            sb.Append(string.Join(Environment.NewLine, MassDiffAcceptors.Select(b => "\t" + b.FileNameAddition)));
            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId)
        {
            myTaskResults = new MyTaskResults(this);

            List<Psm>[] allPsms = new List<Psm>[MassDiffAcceptors.Count];
            for (int searchModeIndex = 0; searchModeIndex < MassDiffAcceptors.Count; searchModeIndex++)
                allPsms[searchModeIndex] = new List<Psm>();
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
            Dictionary<string, Modification> unknownModifications;
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, SearchDecoy, localizeableModifications, b.IsContaminant, out unknownModifications)).ToList();

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
            MyFileManager myFileManager = new MyFileManager();
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
                SearchResults searchResults;
                if (SearchType == SearchType.Classic)
                    searchResults = ((SearchResults)new ClassicSearchEngine(arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, MassDiffAcceptors, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, ionTypes, thisId, ConserveMemory, InitiatorMethionineBehavior.Variable).Run());
                else
                    searchResults = ((SearchResults)(new ModernSearchEngine(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, MassDiffAcceptors, thisId).Run()));
                CrosslinkSearchResults xlsearchResults;
                if (CrosslinkSearchWithAllBeta)
                    xlsearchResults = ((CrosslinkSearchResults)new CrosslinkSearchEngine2(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, MassDiffAcceptors[0], crosslinker, CrosslinkSearchTopNum, CrosslinkSearchWithAllBeta, XLprecusorMsTl, modsDictionary, ionTypes, proteinList, Protease, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, variableModifications, fixedModifications, MaxModificationIsoforms, thisId).Run());
                else
                    xlsearchResults = ((CrosslinkSearchResults)new CrosslinkSearchEngine(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, MassDiffAcceptors[0], crosslinker, CrosslinkSearchTopNum, CrosslinkSearchWithAllBeta, XLprecusorMsTl, ionTypes, modsDictionary, thisId).Run());

                myFileManager.DoneWithFile(origDataFile);

                lock (lock2)
                {
                    for (int searchModeIndex = 0; searchModeIndex < MassDiffAcceptors.Count; searchModeIndex++)
                        allPsms[searchModeIndex].AddRange(searchResults.Psms[searchModeIndex]);
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
            allPsmsXLTuple.OrderByDescending(p => p.Item1.XLTotalScore);
            //WriteCrosslinkToTsv(allPsmsXLTuple, OutputFolder, "xl_all", new List<string> { taskId });

            var allPsmsXLTupleFDR = CrosslinkDoFalseDiscoveryRateAnalysis(allPsmsXLTuple, new OpenSearchMode());
            WriteCrosslinkToTsv(allPsmsXLTupleFDR, OutputFolder, "xl_fdr", new List<string> { taskId });

            var interPsmsXLTupleFDR = allPsmsXLTuple.Where(p => p.Item1.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() == p.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).ToList();
            interPsmsXLTupleFDR = CrosslinkDoFalseDiscoveryRateAnalysis(interPsmsXLTupleFDR, new OpenSearchMode()).Where(p => p.Item1.MostProbableProteinInfo.IsDecoy != true && p.Item2.MostProbableProteinInfo.IsDecoy != true && p.Item1.FdrInfo.QValue <= 0.01).ToList();
            WriteCrosslinkToTsv(interPsmsXLTupleFDR, OutputFolder, "xl_inter_fdr", new List<string> { taskId });

            var intraPsmsXLTupleFDR = allPsmsXLTuple.Where(p => p.Item1.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() != p.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).ToList();
            intraPsmsXLTupleFDR = CrosslinkDoFalseDiscoveryRateAnalysis(intraPsmsXLTupleFDR, new OpenSearchMode()).Where(p => p.Item1.MostProbableProteinInfo.IsDecoy != true && p.Item2.MostProbableProteinInfo.IsDecoy != true && p.Item1.FdrInfo.QValue <= 0.01).ToList();
            WriteCrosslinkToTsv(interPsmsXLTupleFDR, OutputFolder, "xl_inter_fdr", new List<string> { taskId });

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

        private void WriteCrosslinkToTsv(List<Tuple<PsmCross, PsmCross>> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".mytsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("File Name\tScan Numer\tPrecusor MZ\tPrecusor charge\tPrecusor mass" +
                    "\tPep1\tPep1 Protein Access\tPep1 Base sequence\tPep1 Full sequence\tPep1 mass\tPep1 XLBestScore\tPep1 topPos" +
                    "\tPep2\tPep2 Protein Access\tPep2 Base sequence\tPep2 Full sequence\tPep2 mass\tPep2 XLBestScore\tPep2 topPos\tTotalScore\tMass diff\tQValue");
                foreach (var item in items)
                {
                    var x = item.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession);
                    output.WriteLine(
                        item.Item1.FullFilePath
                        + "\t" + item.Item1.ScanNumber.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.ScanPrecursorMonoisotopicPeak.ToString() //CultureInfo.InvariantCulture
                        + "\t" + item.Item1.ScanPrecursorCharge.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.ScanPrecursorMass.ToString(CultureInfo.InvariantCulture)
                        + "\t"
                        + "\t" + item.Item1.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.MostProbableProteinInfo.BaseSequence
                        + "\t" + item.Item1.MostProbableProteinInfo.FullSequence
                        + "\t" + item.Item1.MostProbableProteinInfo.PeptideMonoisotopicMass.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.Score.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.topPosition[0].ToString(CultureInfo.InvariantCulture)
                        //+ "\t" + item.Item1.NScore.ToString(CultureInfo.InvariantCulture)
                        + "\t"
                        + "\t" + item.Item2.MostProbableProteinInfo.PeptidesWithSetModifications.Select(p => p.Protein.Accession).First().ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item2.MostProbableProteinInfo.BaseSequence
                        + "\t" + item.Item2.MostProbableProteinInfo.FullSequence
                        + "\t" + item.Item2.MostProbableProteinInfo.PeptideMonoisotopicMass.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item2.Score.ToString(CultureInfo.InvariantCulture)
                        + "\t" + item.Item1.topPosition[1].ToString(CultureInfo.InvariantCulture)
                        //+ "\t" + item.Item2.NScore.ToString(CultureInfo.InvariantCulture)

                        + "\t" + item.Item1.XLTotalScore.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (item.Item1.ScanPrecursorMass - item.Item2.MostProbableProteinInfo.PeptideMonoisotopicMass - item.Item1.MostProbableProteinInfo.PeptideMonoisotopicMass)
                        + "\t" + item.Item1.FdrInfo.QValue.ToString(CultureInfo.InvariantCulture));
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
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
                int notch1 = item1.MostProbableProteinInfo.Notch; int notch2 = item1.MostProbableProteinInfo.Notch;
                if (isDecoy1 || isDecoy2)
                    cumulative_decoy++;
                else
                    cumulative_target++;


                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                item1.SetValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
                item2.SetValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
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
        
        #endregion Private Methods

    }
}
