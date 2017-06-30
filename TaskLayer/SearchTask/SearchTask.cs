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
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class SearchTask : MetaMorpheusTask
    {

        #region Private Fields

        private const double binTolInDaltons = 0.003;
        private const int max_mods_for_peptide = 3;

        #endregion Private Fields

        #region Public Constructors

        public SearchTask() : base(MyTask.Search)
        {
            // Set default values here:
            DoParsimony = false;
            NoOneHitWonders = false;
            ModPeptidesAreUnique = false;
            Quantify = false;
            QuantifyPpmTol = 5;
            SearchDecoy = true;
            DoHistogramAnalysis = false;
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
            DoLocalizationAnalysis = true;

            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();

            WritePrunedDatabase = false;
            KeepAllUniprotMods = true;

            MassDiffAcceptors = GlobalTaskLevelSettings.SearchModesKnown.Take(1).ToList();

            FindAllPrecursors = true;
            UseProvidedPrecursorInfo = true;

            ConserveMemory = false;
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
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreUnique { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool Quantify { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchDecoy { get; set; }
        public List<MassDiffAcceptor> MassDiffAcceptors { get; set; }
        public bool ConserveMemory { get; set; }

        public bool WritePrunedDatabase { get; set; }
        public bool LocalizeAll { get; set; }
        public bool KeepAllUniprotMods { get; set; }

        public bool FindAllPrecursors { get; set; }
        public bool UseProvidedPrecursorInfo { get; set; }
        public bool DoLocalizationAnalysis { get; set; }

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
            sb.AppendLine("quantify: " + Quantify);
            if (Quantify)
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

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> localizeableModifications;
            if (LocalizeAll)
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            #region Populate modsDictionary

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            {
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
            }

            #endregion Populate modsDictionary

            List<PsmParent>[] allPsms = new List<PsmParent>[MassDiffAcceptors.Count()];
            for (int j = 0; j < MassDiffAcceptors.Count; j++)
                allPsms[j] = new List<PsmParent>();

            Status("Loading proteins...", new List<string> { taskId });
            Dictionary<string, Modification> um;
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FileName, SearchDecoy, localizeableModifications, b.IsContaminant, out um)).ToList();

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

            if (SearchType == SearchType.Modern)
            {
                #region Generate indices for modern search

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

                #endregion Generate indices for modern search
            }

            //FlashLFQEngine FlashLfqEngine = null;
            //if (Quantify)
            //{
            //    FlashLfqEngine = new FlashLFQEngine();
            //    FlashLfqEngine.PassFilePaths(currentRawFileList.ToArray());

            //    if (!FlashLfqEngine.ReadPeriodicTable())
            //        throw new Exception("Quantification error - could not find periodic table file");

            //    if (!FlashLfqEngine.ParseArgs(new string[] {
            //    "--ppm " + QuantifyPpmTol,
            //    "--sil true",
            //    "--pau false",
            //    "--mbr " + MatchBetweenRuns }
            //    ))
            //        throw new Exception("Quantification error - Could not pass parameters");
            //}

            // individual file analysis
            object lock1 = new object();

            Status("Searching files...", new List<string> { taskId });
            Parallel.For(0, currentRawFileList.Count, spectraFileIndex =>
                {
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    NewCollection(Path.GetFileName(origDataFile), new List<string> { taskId, "Individual Spectra Files", origDataFile });

                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                    lock (lock1) // Lock because reading is sequential
                    {
                        Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        if (Path.GetExtension(origDataFile).Equals(".mzML"))
                            myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                        else
                            myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);
                    }

                    Status("Getting ms2 scans...", new List<string> { taskId, "Individual Spectra Files", origDataFile });

                    Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, FindAllPrecursors, UseProvidedPrecursorInfo, 4, origDataFile).OrderBy(b => b.PrecursorMass).ToArray();

                    if (SearchType == SearchType.Classic)
                    {
                        var classicSearchResults = (SearchResults)new ClassicSearchEngine(arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, MassDiffAcceptors, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, new List<string> { taskId, "Individual Spectra Files", origDataFile }, ConserveMemory).Run();
                        for (int searchModeIndex = 0; searchModeIndex < MassDiffAcceptors.Count(); searchModeIndex++)
                            allPsms[searchModeIndex].AddRange(classicSearchResults.Psms[searchModeIndex]);
                    }
                    else
                    {
                        var modernSearchResults = (SearchResults)new ModernSearchEngine(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, MassDiffAcceptors, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();
                        for (int searchModeIndex = 0; searchModeIndex < MassDiffAcceptors.Count; searchModeIndex++)
                            allPsms[searchModeIndex].AddRange(modernSearchResults.Psms[searchModeIndex]);
                    }
                    ReportProgress(new ProgressEventArgs(100, "Done with search!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
                }
            );
            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            // Group and order psms
            for (int j = 0; j < MassDiffAcceptors.Count; j++)
                allPsms[j] = allPsms[j].Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => Math.Abs(b.ScanPrecursorMass - b.PeptideMonoisotopicMass)).GroupBy(b => new Tuple<string, int, double>(b.FullFilePath, b.ScanNumber, b.PeptideMonoisotopicMass)).Select(b => b.First()).ToList();
            List<ProteinGroup>[] ProteinGroups = null;
            bool matchSequencesToActualProteinPeptides = true;
            if (matchSequencesToActualProteinPeptides)
            {
                SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsms, modsDictionary, proteinList, MassDiffAcceptors, Protease, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, InitiatorMethionineBehavior, fixedModifications, variableModifications, MaxModificationIsoforms);
                var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;

                Status("Running FDR analysis...", new List<string> { taskId });

                var analysisResults = new FdrAnalysisEngine(allPsms,
                    compactPeptideToProteinPeptideMatching,
                    MassDiffAcceptors,
                    DoParsimony, NoOneHitWonders, ModPeptidesAreUnique,
                    new List<string> { taskId }).Run();
                ProteinGroups = ((FdrAnalysisResults)analysisResults).ProteinGroups;

                bool doQuantification = false;
                bool doHistogramAnalysis = true;
                bool needToBackToLookAtSpectraAgain = DoLocalizationAnalysis || doQuantification;
                if (needToBackToLookAtSpectraAgain)
                {
                    Status("Analyzing Spectra...", new List<string> { taskId });
                    Parallel.For(0, currentRawFileList.Count, spectraFileIndex =>
                    {
                        var origDataFile = currentRawFileList[spectraFileIndex];
                        IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                        lock (lock1) // Lock because reading is sequential
                        {
                            Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                            if (Path.GetExtension(origDataFile).Equals(".mzML"))
                                myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                            else
                                myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);
                        }
                        if (DoLocalizationAnalysis)
                        {
                            Status("Running localization analysis...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                            var localizationEngine = new LocalizationEngine(allPsms.SelectMany(b => b).Where(b => b.FullFilePath.Equals(origDataFile)), lp, myMsDataFile, ProductMassTolerance);
                            localizationEngine.Run();
                        }

                        if (doQuantification)
                        {
                            //if (FlashLfqEngine != null && ((FlashLfqEngine.filePaths.Count() == 1) || (FlashLfqEngine.mbr && myMsDataFile == null) || (!FlashLfqEngine.mbr && myMsDataFile != null)))
                            //{
                            //    Status("Quantifying peptides...", new List<string> { taskId });
                            //    // use FlashLFQ to quantify peaks
                            //    var psmsBelowOnePercentFdr = orderedPsmsWithFDR.Where(p => p.QValue < 0.01 && !p.IsDecoy);

                            //    foreach (var psm in psmsBelowOnePercentFdr)
                            //        FlashLfqEngine.AddIdentification(psm.FileName, psm.Pli.BaseSequence, psm.Pli.FullSequence, psm.Pli.PeptideMonoisotopicMass, psm.ScanRetentionTime, psm.ScanPrecursorCharge, string.Join("|", psm.Pli.PeptidesWithSetModifications.Select(v => v.Protein.Accession).Distinct().OrderBy(v => v)));

                            //    FlashLfqEngine.ConstructBinsFromIdentifications();

                            //    if (myMsDataFile != null)
                            //    {
                            //        string fileName = orderedPsmsWithFDR.First().FileName;
                            //        string path = FlashLfqEngine.filePaths.Where(p => System.IO.Path.GetFileNameWithoutExtension(p).Equals(fileName)).First();
                            //        FlashLfqEngine.Quantify(myMsDataFile, path);
                            //    }
                            //    else
                            //    {
                            //        for (int i = 0; i < FlashLfqEngine.filePaths.Length; i++)
                            //        {
                            //            FlashLfqEngine.Quantify(null, FlashLfqEngine.filePaths[i]);
                            //            GC.Collect();
                            //        }
                            //    }

                            //    if (FlashLfqEngine.mbr)
                            //        FlashLfqEngine.RetentionTimeCalibrationAndErrorCheckMatchedFeatures();

                            //    var allPeaks = FlashLfqEngine.allFeaturesByFile.SelectMany(p => p).ToList();
                            //    writeQuantifiedPeaksAction(allPeaks, "QuantifiedPeaks_" + searchModes[j].FileNameAddition, nestedIds);

                            //    var summedPeaks = FlashLfqEngine.SumFeatures(allPeaks);
                            //    writeQuantifiedPeptidesAction(summedPeaks.ToList(), "QuantifiedPeptides_" + searchModes[j].FileNameAddition, nestedIds);

                            //    // assign quantities to PSMs
                            //    string[] fileNames = new string[FlashLfqEngine.filePaths.Length];
                            //    for (int i = 0; i < FlashLfqEngine.filePaths.Length; i++)
                            //        fileNames[i] = System.IO.Path.GetFileNameWithoutExtension(FlashLfqEngine.filePaths[i]);

                            //    Dictionary<string, List<PsmParent>> baseseqToPsm = new Dictionary<string, List<PsmParent>>();
                            //    List<PsmParent> list;
                            //    foreach (var psm in orderedPsmsWithFDR)
                            //    {
                            //        if (baseseqToPsm.TryGetValue(psm.Pli.BaseSequence, out list))
                            //            list.Add(psm);
                            //        else
                            //            baseseqToPsm.Add(psm.Pli.BaseSequence, new List<PsmParent>() { psm });
                            //    }

                            //    foreach (var summedPeak in summedPeaks)
                            //    {
                            //        if (baseseqToPsm.TryGetValue(summedPeak.BaseSequence, out list))
                            //        {
                            //            var psmsForThisBaseSeqAndFile = list.GroupBy(p => p.FileName);
                            //            foreach (var file in psmsForThisBaseSeqAndFile)
                            //            {
                            //                int i = Array.IndexOf(fileNames, file.Key);

                            //                foreach (var psm in file)
                            //                    psm.QuantIntensity[0] = summedPeak.intensitiesByFile[i];
                            //            }
                            //        }
                            //    }
                            //}

                            //if (FlashLfqEngine != null)
                            //    foreach (var proteinGroup in proteinGroups)
                            //        proteinGroup.Quantify();

                            //if (FlashLfqEngine != null)
                            //{
                            //    // call multifile protein quantification helper function (need all the filenames to organize results properly)
                            //    var files = orderedPsmsWithFDR.Select(p => p.FileName).Distinct().ToList();
                            //    foreach (var pg in proteinGroups[j])
                            //        pg.AggregateQuantifyHelper(files);
                            //}
                        }
                        ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
                    });
                }
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

                if (doHistogramAnalysis)
                {
                    for (int j = 0; j < MassDiffAcceptors.Count; j++)
                    {
                        var limitedpsms_with_fdr = allPsms[j].Where(b => (b.FdrInfo.QValue <= 0.01)).ToList();
                        if (limitedpsms_with_fdr.Any(b => !b.Pli.IsDecoy))
                        {
                            Status("Running histogram analysis...", new List<string> { taskId });
                            var myTreeStructure = new BinTreeStructure();
                            myTreeStructure.GenerateBins(limitedpsms_with_fdr, binTolInDaltons);
                            WriteTree(myTreeStructure, OutputFolder, "aggregate_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });
                        }
                    }
                }
            }

            // Now that we are done with fdr analysis and localization analysis, can write the results!
            for (int j = 0; j < MassDiffAcceptors.Count; j++)
            {
                WritePsmsToTsv(allPsms[j], OutputFolder, "aggregatePSMs_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });

                //myAnalysisResults.AddText("Unique peptides within 1% FDR: " + peptidesWithFDR.Count(a => a.QValue <= .01 && a.IsDecoy == false));

                WritePsmsToTsv(allPsms[j].GroupBy(b => b.Pli.FullSequence).Select(b => b.FirstOrDefault()), OutputFolder, "aggregateUniquePeptides_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId });

                var psmsGroupedByFile = allPsms[j].GroupBy(p => p.FullFilePath);

                // individual psm files (with global psm fdr, global parsimony)
                foreach (var group in psmsGroupedByFile)
                {
                    var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);
                    WritePsmsToTsv(group.ToList(), OutputFolder, strippedFileName + "_allPSMs_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
                }

                // individual protein group files (local protein fdr, global parsimony, global psm fdr)
                var fullFilePaths = psmsGroupedByFile.Select(p => p.Key).Distinct();
                if (DoParsimony)
                    foreach (var fullFilePath in fullFilePaths)
                    {
                        var strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                        var subsetProteinGroupsForThisFile = new List<ProteinGroup>();
                        foreach (var pg in ProteinGroups[j])
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
                        WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, OutputFolder, strippedFileName + "_" + MassDiffAcceptors[j].FileNameAddition + "_ProteinGroups", new List<string> { taskId, "Individual Spectra Files", fullFilePath });

                        WriteMzidentml(psmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList(), subsetProteinGroupsForThisFile, variableModifications, fixedModifications, new List<Protease> { Protease }, 0.01, MassDiffAcceptors[j], ProductMassTolerance, MaxMissedCleavages, OutputFolder, strippedFileName + "_" + MassDiffAcceptors[j].FileNameAddition, new List<string> { taskId, "Individual Spectra Files", fullFilePath });
                    }
            }

            if (WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { taskId });

                List<Modification> modificationsToAlwaysKeep = new List<Modification>();
                if (KeepAllUniprotMods)
                    modificationsToAlwaysKeep.AddRange(GlobalTaskLevelSettings.AllModsKnown.Where(b => b.modificationType.Equals("Uniprot")));

                var goodPsmsForEachProtein = allPsms.SelectMany(b => b).Where(b => b.FdrInfo.QValueNotch < 0.01 && b.Pli.PeptidesWithSetModifications.Count == 1 && !b.Pli.IsDecoy).GroupBy(b => b.Pli.PeptidesWithSetModifications.First().Protein).ToDictionary(b => b.Key);

                foreach (var protein in proteinList)
                {
                    if (!protein.IsDecoy)
                    {
                        HashSet<Tuple<int, ModificationWithMass>> modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>();
                        if (goodPsmsForEachProtein.ContainsKey(protein))
                            modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>(goodPsmsForEachProtein[protein].SelectMany(b => b.Pli.PeptidesWithSetModifications.First().allModsOneIsNterminus.Select(c => new Tuple<int, ModificationWithMass>(GetOneBasedIndexInProtein(c.Key, b.Pli.PeptidesWithSetModifications.First()), c.Value))));

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

                if (dbFilenameList.Any(b => !b.IsContaminant))
                {
                    string outputXMLdbFullName = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "pruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                    SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });
                }
                if (dbFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FileName))) + "pruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

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