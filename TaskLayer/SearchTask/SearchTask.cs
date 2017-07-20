using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
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
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class SearchTask : MetaMorpheusTask
    {

        #region Private Fields

        private const double binTolInDaltons = 0.003;

        private FlashLFQEngine FlashLfqEngine;

        #endregion Private Fields

        #region Public Constructors

        public SearchTask() : base(MyTask.Search)
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
            MinPeptideLength = 5;
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
            MaxDegreeOfParallelism = 1;

            // Deconvolution stuff
            DoPrecursorDeconvolution = true;
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

        public bool addCompIons { get; set; }

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

            Status("Loading modifications...", taskId);

            #region Load modifications

            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> localizeableModifications;
            if (LocalizeAll)
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

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
            if (SearchType == SearchType.Modern)
            {
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
            }

            object lock2 = new object();
            MyFileManager myFileManager = new MyFileManager();
            Status("Searching files...", taskId);
            ParallelOptions parallelOptions = new ParallelOptions();
            if (MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = MaxDegreeOfParallelism.Value;
            double completedFiles = 0;
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
                        searchResults = ((SearchResults)new ClassicSearchEngine(arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, MassDiffAcceptors, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, ionTypes, thisId, ConserveMemory, InitiatorMethionineBehavior, this.addCompIons).Run());
                    else
                        searchResults = ((SearchResults)(new ModernSearchEngine(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, MassDiffAcceptors, thisId, this.addCompIons, ionTypes).Run()));

                    myFileManager.DoneWithFile(origDataFile);

                    lock (lock2)
                    {
                        for (int searchModeIndex = 0; searchModeIndex < MassDiffAcceptors.Count; searchModeIndex++)
                            allPsms[searchModeIndex].AddRange(searchResults.Psms[searchModeIndex]);
                    }
                    ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                    completedFiles++;
                    ReportProgress(new ProgressEventArgs((int)completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
                }
            );
            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            Status("Ordering and grouping psms...", taskId);
            for (int j = 0; j < allPsms.Length; j++)
                allPsms[j] = allPsms[j].Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => Math.Abs(b.ScanPrecursorMass - b.PeptideMonoisotopicMass)).GroupBy(b => new Tuple<string, int, double>(b.FullFilePath, b.ScanNumber, b.PeptideMonoisotopicMass)).Select(b => b.First()).ToList();

            // Group and order psms
            Status("Matching peptides to proteins...", taskId);
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, MassDiffAcceptors, Protease, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, InitiatorMethionineBehavior, fixedModifications, variableModifications, MaxModificationIsoforms, new List<string> { taskId });
            var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;

            if (DoLocalizationAnalysis)
            {
                Status("Running localization analysis...", taskId);
                Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    Status("Running localization analysis...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile);
                    var localizationEngine = new LocalizationEngine(allPsms.SelectMany(b => b).Where(b => b != null && b.FullFilePath.Equals(origDataFile)).ToList(), ionTypes, myMsDataFile, ProductMassTolerance, new List<string> { taskId, "Individual Spectra Files", origDataFile }, this.addCompIons);
                    localizationEngine.Run();
                    myFileManager.DoneWithFile(origDataFile);
                    ReportProgress(new ProgressEventArgs(100, "Done with localization analysis!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
                });
                Status("Grouping by matched ions...", taskId);
                for (int j = 0; j < allPsms.Length; j++)
                    allPsms[j] = allPsms[j].GroupBy(b => new Tuple<string, int, MatchedIonMassesListPositiveIsMatch>(b.FullFilePath, b.ScanNumber, b.LocalizationResults.MatchedIonMassesListPositiveIsMatch)).Select(b => b.First()).ToList();
            }

            Status("Running FDR analysis...", taskId);

            var fdrAnalysisResults = new FdrAnalysisEngine(allPsms, MassDiffAcceptors, new List<string> { taskId }).Run();

            new ModificationAnalysisEngine(allPsms, MassDiffAcceptors, new List<string> { taskId }).Run();

            ProteinAnalysisResults proteinAnalysisResults = null;
            if (DoParsimony)
                proteinAnalysisResults = (ProteinAnalysisResults)(new ProteinAnalysisEngine(allPsms, compactPeptideToProteinPeptideMatching, MassDiffAcceptors, NoOneHitWonders, ModPeptidesAreUnique, new List<string> { taskId }).Run());

            if (DoQuantification)
            {
                // use FlashLFQ to quantify peaks across all files
                Status("Quantifying...", taskId);
                FlashLfqEngine.PassFilePaths(currentRawFileList.ToArray());

                if (!FlashLfqEngine.ReadPeriodicTable())
                    throw new Exception("Quantification error - could not find periodic table file");

                if (!FlashLfqEngine.ParseArgs(new string[] {
                        "--ppm " + QuantifyPpmTol,
                        "--sil true",
                        "--pau false",
                        "--mbr " + MatchBetweenRuns }
                    ))
                    throw new Exception("Quantification error - Could not pass parameters to quantification engine");

                var psmsBelowOnePercentFdr = allPsms.SelectMany(v => v).Where(p => p.FdrInfo.QValue < 0.01 && !p.MostProbableProteinInfo.IsDecoy);
                foreach (var psm in psmsBelowOnePercentFdr)
                    FlashLfqEngine.AddIdentification(Path.GetFileNameWithoutExtension(psm.FullFilePath), psm.MostProbableProteinInfo.BaseSequence, psm.MostProbableProteinInfo.FullSequence, psm.MostProbableProteinInfo.PeptideMonoisotopicMass, psm.ScanRetentionTime, psm.ScanPrecursorCharge, string.Join("|", psm.MostProbableProteinInfo.PeptidesWithSetModifications.Select(v => v.Protein.Accession).Distinct().OrderBy(v => v)));

                FlashLfqEngine.ConstructBinsFromIdentifications();

                Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    Status("Quantifying...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile);
                    FlashLfqEngine.Quantify(myMsDataFile, origDataFile);
                    myFileManager.DoneWithFile(origDataFile);
                    GC.Collect();
                    ReportProgress(new ProgressEventArgs(100, "Done quantifying!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
                });

                if (FlashLfqEngine.mbr)
                    FlashLfqEngine.RetentionTimeCalibrationAndErrorCheckMatchedFeatures();

                // assign quantities to PSMs
                Dictionary<string, List<Psm>> baseseqToPsm = new Dictionary<string, List<Psm>>();
                List<Psm> list;
                foreach (var psm in psmsBelowOnePercentFdr)
                {
                    if (baseseqToPsm.TryGetValue(psm.MostProbableProteinInfo.BaseSequence, out list))
                        list.Add(psm);
                    else
                        baseseqToPsm.Add(psm.MostProbableProteinInfo.BaseSequence, new List<Psm>() { psm });
                }

                var summedPeaks = FlashLfqEngine.SumFeatures(FlashLfqEngine.allFeaturesByFile.SelectMany(p => p).ToList(), "BaseSequence");
                foreach (var summedPeak in summedPeaks)
                {
                    if (baseseqToPsm.TryGetValue(summedPeak.BaseSequence, out list))
                    {
                        var psmsForThisBaseSeqAndFile = list.GroupBy(p => p.FullFilePath);
                        foreach (var file in psmsForThisBaseSeqAndFile)
                        {
                            int j = Array.IndexOf(FlashLfqEngine.filePaths, file.Key);

                            foreach (var psm in file)
                                psm.QuantIntensity[0] = summedPeak.intensitiesByFile[j];
                        }
                    }
                }
            }

            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            if (DoHistogramAnalysis)
            {
                for (int j = 0; j < MassDiffAcceptors.Count; j++)
                {
                    var limitedpsms_with_fdr = allPsms[j].Where(b => (b.FdrInfo.QValue <= 0.01)).ToList();
                    if (limitedpsms_with_fdr.Any(b => !b.MostProbableProteinInfo.IsDecoy))
                    {
                        Status("Running histogram analysis...", new List<string> { taskId });
                        var myTreeStructure = new BinTreeStructure();
                        myTreeStructure.GenerateBins(limitedpsms_with_fdr, binTolInDaltons);
                        WriteTree(myTreeStructure, OutputFolder, "aggregate_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });
                    }
                }
            }

            // Now that we are done with fdr analysis and localization analysis, can write the results!
            Status("Writing results...", taskId);
            for (int j = 0; j < MassDiffAcceptors.Count; j++)
            {
                WritePsmsToTsv(allPsms[j], OutputFolder, "aggregatePSMs_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });

                myTaskResults.AddNiceText("All target PSMS within 1% FDR " + MassDiffAcceptors[j].FileNameAddition + ": " + allPsms[j].Count(a => a.FdrInfo.QValue <= .01 && a.MostProbableProteinInfo.IsDecoy == false));

                var uniquePeptides = allPsms[j].GroupBy(b => b.MostProbableProteinInfo.FullSequence).Select(b => b.FirstOrDefault()).ToList();

                WritePsmsToTsv(uniquePeptides, OutputFolder, "aggregateUniquePeptides_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });

                myTaskResults.AddNiceText("Unique peptides within 1% FDR " + MassDiffAcceptors[j].FileNameAddition + ": " + uniquePeptides.Count(a => a.FdrInfo.QValue <= .01 && a.MostProbableProteinInfo.IsDecoy == false));

                var psmsGroupedByFile = allPsms[j].GroupBy(p => p.FullFilePath);

                // individual psm files (with global psm fdr, global parsimony)
                foreach (var group in psmsGroupedByFile)
                {
                    var psmsForThisFile = group.ToList();

                    var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);

                    myTaskResults.AddNiceText("PSMs within 1% FDR " + MassDiffAcceptors[j].FileNameAddition + " " + strippedFileName + ": " + psmsForThisFile.Count(a => a.FdrInfo.QValue <= .01 && a.MostProbableProteinInfo.IsDecoy == false));

                    WritePsmsToTsv(psmsForThisFile, OutputFolder, strippedFileName + "_PSMs_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
                    var uniquePeptidesForFile = psmsForThisFile.GroupBy(b => b.MostProbableProteinInfo.FullSequence).Select(b => b.FirstOrDefault()).ToList();
                    WritePsmsToTsv(uniquePeptidesForFile, OutputFolder, strippedFileName + "_UniquePeptides_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
                    myTaskResults.AddNiceText("Unique peptides within 1% FDR " + MassDiffAcceptors[j].FileNameAddition + " " + strippedFileName + ": " + uniquePeptidesForFile.Count(a => a.FdrInfo.QValue <= .01 && a.MostProbableProteinInfo.IsDecoy == false));
                }

                if (DoParsimony)
                {
                    // aggregate protein group file
                    foreach (var pg in proteinAnalysisResults.ProteinGroups[j])
                    {
                        if (pg.ProteinGroupScore != 0)
                            pg.AggregateQuantifyHelper(currentRawFileList);
                    }

                    WriteProteinGroupsToTsv(proteinAnalysisResults.ProteinGroups[j], OutputFolder, "aggregateProteinGroups_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });

                    // individual protein group files (local protein fdr, global parsimony, global psm fdr)
                    foreach (var fullFilePath in currentRawFileList)
                    {
                        var strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                        var subsetProteinGroupsForThisFile = new List<ProteinGroup>();
                        foreach (var pg in proteinAnalysisResults.ProteinGroups[j])
                        {
                            var subsetPg = pg.ConstructSubsetProteinGroup(fullFilePath);
                            subsetPg.Score();

                            if (subsetPg.ProteinGroupScore != 0)
                            {
                                subsetPg.CalculateSequenceCoverage();
                                subsetPg.Quantify();
                                subsetProteinGroupsForThisFile.Add(subsetPg);
                            }
                        }
                        subsetProteinGroupsForThisFile = new ProteinAnalysisEngine(allPsms, compactPeptideToProteinPeptideMatching, MassDiffAcceptors, NoOneHitWonders, ModPeptidesAreUnique, new List<string> { taskId, "Individual Spectra Files", fullFilePath }).DoProteinFdr(subsetProteinGroupsForThisFile);
                        WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, OutputFolder, strippedFileName + "_" + MassDiffAcceptors[j].FileNameAddition + "_ProteinGroups", new List<string> { taskId, "Individual Spectra Files", fullFilePath });

                        Status("Writing mzid...", new List<string> { taskId, "Individual Spectra Files", fullFilePath });
                        WriteMzidentml(psmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList(), subsetProteinGroupsForThisFile, variableModifications, fixedModifications, new List<Protease> { Protease }, 0.01, MassDiffAcceptors[j], ProductMassTolerance, MaxMissedCleavages, OutputFolder, strippedFileName + "_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId, "Individual Spectra Files", fullFilePath });
                        ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", fullFilePath }));
                    }
                }

                if (DoQuantification)
                {
                    foreach (var fullFilePath in currentRawFileList)
                    {
                        var strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);
                        var peaksForThisFile = FlashLfqEngine.allFeaturesByFile[Array.IndexOf(FlashLfqEngine.filePaths, fullFilePath)];

                        WritePeakQuantificationResultsToTsv(peaksForThisFile, OutputFolder, strippedFileName + "_" + MassDiffAcceptors[j].FileNameAddition + "_QuantifiedPeaks", new List<string> { taskId, "Individual Spectra Files", fullFilePath });
                    }

                    var summedPeaksByPeptide = FlashLfqEngine.SumFeatures(FlashLfqEngine.allFeaturesByFile.SelectMany(p => p).ToList(), "BaseSequence");
                    WritePeptideQuantificationResultsToTsv(summedPeaksByPeptide.ToList(), OutputFolder, "aggregateQuantifiedPeptides_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });
                }
            }

            if (WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { taskId });

                List<Modification> modificationsToAlwaysKeep = new List<Modification>();
                if (KeepAllUniprotMods)
                    modificationsToAlwaysKeep.AddRange(GlobalTaskLevelSettings.AllModsKnown.Where(b => b.modificationType.Equals("Uniprot")));

                var goodPsmsForEachProtein = allPsms.SelectMany(b => b).Where(b => b.FdrInfo.QValueNotch < 0.01 && b.MostProbableProteinInfo.PeptidesWithSetModifications.Count == 1 && !b.MostProbableProteinInfo.IsDecoy).GroupBy(b => b.MostProbableProteinInfo.PeptidesWithSetModifications.First().Protein).ToDictionary(b => b.Key);

                foreach (var protein in proteinList)
                {
                    if (!protein.IsDecoy)
                    {
                        HashSet<Tuple<int, ModificationWithMass>> modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>();
                        if (goodPsmsForEachProtein.ContainsKey(protein))
                            modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>(goodPsmsForEachProtein[protein].SelectMany(b => b.MostProbableProteinInfo.PeptidesWithSetModifications.First().allModsOneIsNterminus.Select(c => new Tuple<int, ModificationWithMass>(GetOneBasedIndexInProtein(c.Key, b.MostProbableProteinInfo.PeptidesWithSetModifications.First()), c.Value))));

                        IDictionary<int, List<Modification>> modsToWrite = new Dictionary<int, List<Modification>>();
                        foreach (var modd in protein.OneBasedPossibleLocalizedModifications)
                            foreach (var mod in modd.Value)
                            {
                                if (modificationsToAlwaysKeep.Contains(mod as Modification)
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

                //writes all proteins
                if (dbFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                    SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });
                }
                if (dbFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                    SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });
                }

                //writes only detected proteins
                if (dbFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), goodPsmsForEachProtein.Keys.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                    SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });
                }
                if (dbFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), goodPsmsForEachProtein.Keys.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                    SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });
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