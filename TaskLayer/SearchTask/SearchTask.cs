using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using FlashLFQ;
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

        public SearchTask() : base(MyTask.Search)
        {
            // Set default values here:
            ClassicSearch = true;
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

            SearchModes = GlobalTaskLevelSettings.SearchModesKnown.Take(1).ToList();

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
        public bool ClassicSearch { get; set; }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreUnique { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool Quantify { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchDecoy { get; set; }
        public List<SearchMode> SearchModes { get; set; }
        public bool ConserveMemory { get; set; }

        public bool WritePrunedDatabase { get; set; }
        public bool LocalizeAll { get; set; }
        public bool KeepAllUniprotMods { get; set; }

        public bool FindAllPrecursors { get; set; }
        public bool UseProvidedPrecursorInfo { get; set; }
        public bool DoLocalizationAnalysis { get; set; }

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
            sb.AppendLine("classicSearch: " + ClassicSearch);
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

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            if (!ClassicSearch)
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
            List<NewPsmWithFdr>[] allResultingIdentifications = null;
            List<NewPsmWithFdr>[] allResultingPeptides = null;

            FlashLFQEngine FlashLfqEngine = null;
            if (Quantify)
            {
                FlashLfqEngine = new FlashLFQEngine();
                FlashLfqEngine.PassFilePaths(currentRawFileList.ToArray());

                if (!FlashLfqEngine.ReadPeriodicTable())
                    throw new Exception("Quantification error - could not find periodic table file");

                if (!FlashLfqEngine.ParseArgs(new string[] {
                "--ppm " + QuantifyPpmTol,
                "--sil true",
                "--pau false",
                "--mbr " + MatchBetweenRuns }
                ))
                    throw new Exception("Quantification error - Could not pass parameters");
            }

            // individual file analysis
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
                    myMsDataFile = ThermoDynamicData.InitiateDynamicConnection(origDataFile);

                Status("Getting ms2 scans...", new List<string> { taskId, "Individual Searches", origDataFile });
                var intensityRatio = 4;
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, FindAllPrecursors, UseProvidedPrecursorInfo, intensityRatio, origDataFile).OrderBy(b => b.PrecursorMass).ToArray();

                if (ClassicSearch)
                {
                    var classicSearchResults = (ClassicSearchResults)new ClassicSearchEngine(arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, SearchModes, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, new List<string> { taskId, "Individual Searches", origDataFile }, ConserveMemory).Run();
                    for (int searchModeIndex = 0; searchModeIndex < SearchModes.Count; searchModeIndex++)
                        allPsms[searchModeIndex].AddRange(classicSearchResults.OuterPsms[searchModeIndex]);
                }
                else
                {
                    var modernSearchResults = (ModernSearchResults)new ModernSearchEngine(arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, ProductMassTolerance, SearchModes, new List<string> { taskId, "Individual Searches", origDataFile }).Run();
                    for (int searchModeIndex = 0; searchModeIndex < SearchModes.Count; searchModeIndex++)
                        allPsms[searchModeIndex].AddRange(modernSearchResults.OuterPsms[searchModeIndex]);
                }
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Searches", origDataFile }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Searches" }));

            var analysisResults = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(),
                compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications,
                Protease, SearchModes, null, ProductMassTolerance,
                (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, "aggregate_" + s, new List<string> { taskId }),
                (List<NewPsmWithFdr> h, string s, List<string> ss) => WritePsmsToTsv(h, OutputFolder, s, ss),
                (List<ProteinGroup> h, string s, List<string> ss) => WriteProteinGroupsToTsv(h, OutputFolder, s, ss),
                (List<NewPsmWithFdr> h, List<ProteinGroup> g, SearchMode m, string s, List<string> ss) => WriteMzidentml(h, g, variableModifications, fixedModifications, new List<Protease> { Protease }, 0.01, m, ProductMassTolerance, MaxMissedCleavages, OutputFolder, s, ss),
                DoParsimony, NoOneHitWonders, ModPeptidesAreUnique, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength,
                MaxModificationIsoforms, DoHistogramAnalysis, lp, binTolInDaltons, initiatorMethionineBehavior,
                new List<string> { taskId }, modsDictionary, null, currentRawFileList).Run();

            allResultingIdentifications = ((AnalysisResults)analysisResults).AllResultingIdentifications;
            allResultingPeptides = ((AnalysisResults)analysisResults).allResultingPeptides;
            var ProteinGroups = ((AnalysisResults)analysisResults).ProteinGroups;

            bool doQuantification = false;
            bool doHistogramAnalysis = true;
            bool needToBackToLookAtSpectraAgain = DoLocalizationAnalysis || doQuantification;
            if (needToBackToLookAtSpectraAgain)
            {
                for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
                {
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                    if (Path.GetExtension(origDataFile).Equals(".mzML"))
                        myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                    else
                        myMsDataFile = ThermoDynamicData.InitiateDynamicConnection(origDataFile);
                    if (DoLocalizationAnalysis)
                    {
                        var localizationEngine = new LocalizationEngine(allResultingIdentifications.SelectMany(b => b).Select(b => b.thisPSM).Where(b => b.FileName.Equals(origDataFile)), lp, myMsDataFile, ProductMassTolerance);
                        var localizationResults = localizationEngine.Run();
                    }

                    if (doQuantification)
                    {
                        //if (FlashLfqEngine != null && ((FlashLfqEngine.filePaths.Count() == 1) || (FlashLfqEngine.mbr && myMsDataFile == null) || (!FlashLfqEngine.mbr && myMsDataFile != null)))
                        //{
                        //    Status("Quantifying peptides...", new List<string> { taskId });
                        //    // use FlashLFQ to quantify peaks
                        //    var psmsBelowOnePercentFdr = orderedPsmsWithFDR.Where(p => p.QValue < 0.01 && !p.IsDecoy);

                        //    foreach (var psm in psmsBelowOnePercentFdr)
                        //        FlashLfqEngine.AddIdentification(psm.thisPSM.FileName, psm.thisPSM.Pli.BaseSequence, psm.thisPSM.Pli.FullSequence, psm.thisPSM.Pli.PeptideMonoisotopicMass, psm.thisPSM.ScanRetentionTime, psm.thisPSM.ScanPrecursorCharge, string.Join("|", psm.thisPSM.Pli.PeptidesWithSetModifications.Select(v => v.Protein.Accession).Distinct().OrderBy(v => v)));

                        //    FlashLfqEngine.ConstructBinsFromIdentifications();

                        //    if (myMsDataFile != null)
                        //    {
                        //        string fileName = orderedPsmsWithFDR.First().thisPSM.FileName;
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

                        //    Dictionary<string, List<NewPsmWithFdr>> baseseqToPsm = new Dictionary<string, List<NewPsmWithFdr>>();
                        //    List<NewPsmWithFdr> list;
                        //    foreach (var psm in orderedPsmsWithFDR)
                        //    {
                        //        if (baseseqToPsm.TryGetValue(psm.thisPSM.Pli.BaseSequence, out list))
                        //            list.Add(psm);
                        //        else
                        //            baseseqToPsm.Add(psm.thisPSM.Pli.BaseSequence, new List<NewPsmWithFdr>() { psm });
                        //    }

                        //    foreach (var summedPeak in summedPeaks)
                        //    {
                        //        if (baseseqToPsm.TryGetValue(summedPeak.BaseSequence, out list))
                        //        {
                        //            var psmsForThisBaseSeqAndFile = list.GroupBy(p => p.thisPSM.FileName);
                        //            foreach (var file in psmsForThisBaseSeqAndFile)
                        //            {
                        //                int i = Array.IndexOf(fileNames, file.Key);

                        //                foreach (var psm in file)
                        //                    psm.thisPSM.QuantIntensity[0] = summedPeak.intensitiesByFile[i];
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
                        //    var files = orderedPsmsWithFDR.Select(p => p.thisPSM.FileName).Distinct().ToList();
                        //    foreach (var pg in proteinGroups[j])
                        //        pg.AggregateQuantifyHelper(files);
                        //}
                    }
                }
            }

            if (DoLocalizationAnalysis && doHistogramAnalysis)
            {
                for (int j = 0; j < SearchModes.Count; j++)
                {
                    var limitedpsms_with_fdr = allResultingIdentifications[j].Where(b => (b.QValue <= 0.01)).ToList();
                    if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                    {
                        Status("Running histogram analysis...", new List<string> { taskId });
                        var myTreeStructure = new BinTreeStructure();
                        myTreeStructure.GenerateBins(limitedpsms_with_fdr, binTolInDaltons);
                        WriteTree(myTreeStructure, OutputFolder, "aggregate_" + SearchModes[j].FileNameAddition, new List<string> { taskId });
                    }
                }
            }

            // Now that we are done with fdr analysis and localization analysis, can write the results!
            for (int j = 0; j < SearchModes.Count; j++)
            {
                WritePsmsToTsv(allResultingIdentifications[j], OutputFolder, "allPSMs_" + SearchModes[j].FileNameAddition, new List<string> { taskId });

                WritePsmsToTsv(allResultingPeptides[j], OutputFolder, "uniquePeptides_" + SearchModes[j].FileNameAddition, new List<string> { taskId });

                var psmsGroupedByFilename = allResultingIdentifications[j].GroupBy(p => p.thisPSM.FileName);

                // individual psm files (with global psm fdr, global parsimony)
                foreach (var group in psmsGroupedByFilename)
                {
                    var fileName = System.IO.Path.GetFileNameWithoutExtension(group.First().thisPSM.FileName);
                    WritePsmsToTsv(group.ToList(), OutputFolder, fileName + "_allPSMs_" + SearchModes[j].FileNameAddition, new List<string> { taskId, "Individual Searches", group.First().thisPSM.FileName });
                }

                // individual protein group files (local protein fdr, global parsimony, global psm fdr)
                var fileNames = psmsGroupedByFilename.Select(p => p.Key).Distinct();
                if (DoParsimony)
                    foreach (var file in fileNames)
                    {
                        var subsetProteinGroupsForThisFile = new List<ProteinGroup>();
                        foreach (var pg in ProteinGroups[j])
                        {
                            var subsetPg = pg.ConstructSubsetProteinGroup(file);
                            subsetPg.Score();

                            if (subsetPg.ProteinGroupScore != 0)
                            {
                                subsetPg.CalculateSequenceCoverage();
                                subsetPg.Quantify();
                                subsetProteinGroupsForThisFile.Add(subsetPg);
                            }
                        }
                        WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, OutputFolder, file + "_" + SearchModes[j].FileNameAddition + "_ProteinGroups", new List<string> { taskId, "Individual Searches", file });

                        WriteMzidentml(psmsGroupedByFilename.Where(p => p.Key == file).SelectMany(g => g).ToList(), subsetProteinGroupsForThisFile, variableModifications, fixedModifications, new List<Protease> { Protease }, 0.01, SearchModes[j], ProductMassTolerance, MaxMissedCleavages, OutputFolder, file + "_" + SearchModes[j].FileNameAddition, new List<string> { taskId, "Individual Searches", file });
                    }
            }

            if (WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { taskId });

                List<Modification> modificationsToAlwaysKeep = new List<Modification>();
                if (KeepAllUniprotMods)
                    modificationsToAlwaysKeep.AddRange(GlobalTaskLevelSettings.AllModsKnown.Where(b => b.modificationType.Equals("Uniprot")));

                var goodPsmsForEachProtein = allResultingIdentifications.SelectMany(b => b).Where(b => b.QValueNotch < 0.01 && b.thisPSM.Pli.PeptidesWithSetModifications.Count == 1 && !b.IsDecoy).GroupBy(b => b.thisPSM.Pli.PeptidesWithSetModifications.First().Protein).ToDictionary(b => b.Key);

                foreach (var protein in proteinList)
                {
                    if (!protein.IsDecoy)
                    {
                        HashSet<Tuple<int, ModificationWithMass>> modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>();
                        if (goodPsmsForEachProtein.ContainsKey(protein))
                            modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>(goodPsmsForEachProtein[protein].SelectMany(b => b.thisPSM.Pli.PeptidesWithSetModifications.First().allModsOneIsNterminus.Select(c => new Tuple<int, ModificationWithMass>(GetOneBasedIndexInProtein(c.Key, b.thisPSM.Pli.PeptidesWithSetModifications.First()), c.Value))));

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