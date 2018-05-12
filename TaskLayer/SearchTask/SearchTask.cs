using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using EngineLayer.Indexing;
using EngineLayer.Localization;
using EngineLayer.ModernSearch;
using EngineLayer.ModificationAnalysis;
using EngineLayer.NonSpecificEnzymeSearch;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Serialization;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class SearchTask : MetaMorpheusTask
    {
        #region Public Constructors

        public SearchTask() : base(MyTask.Search)
        {
            CommonParameters = new CommonParameters();

            SearchParameters = new SearchParameters();
        }

        #endregion Public Constructors

        #region Public Properties

        public SearchParameters SearchParameters { get; set; }

        #endregion Public Properties

        #region Public Methods
      
        public static MassDiffAcceptor GetMassDiffAcceptor(Tolerance precursorMassTolerance, MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact:
                    if (precursorMassTolerance is PpmTolerance)
                        return new SinglePpmAroundZeroSearchMode(precursorMassTolerance.Value);
                    else
                        return new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance.Value);

                case MassDiffAcceptorType.OneMM:
                    return new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, precursorMassTolerance);

                case MassDiffAcceptorType.TwoMM:
                    return new DotMassDiffAcceptor("2mm", new List<double> { 0, 1.0029, 2.0052 }, precursorMassTolerance);

                case MassDiffAcceptorType.ThreeMM:
                    return new DotMassDiffAcceptor("3mm", new List<double> { 0, 1.0029, 2.0052, 3.0077 }, precursorMassTolerance);

                case MassDiffAcceptorType.ModOpen:
                    return new IntervalMassDiffAcceptor("-187andUp", new List<DoubleRange> { new DoubleRange(-187, double.PositiveInfinity) });

                case MassDiffAcceptorType.Open:
                    return new OpenSearchMode();

                case MassDiffAcceptorType.Custom:
                    return ParseSearchMode(customMdac);

                default:
                    throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        #endregion Public Methods

        #region Protected Methods

        protected static void WriteTree(BinTreeStructure myTreeStructure, string writtenFile)
        {
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tUnimodDiffs\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tprotNtermLocFrac\tpepNtermLocFrac\tpepCtermLocFrac\tprotCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot");
                foreach (Bin bin in myTreeStructure.FinalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString("F4", CultureInfo.InvariantCulture)
                        + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.Count == 0 ? double.NaN : (double)bin.CountDecoy / bin.Count).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : (double)bin.LocalizeableTarget / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.Mine
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.UnimodDiffs
                        + "\t" + bin.AA
                        + "\t" + bin.Combos
                        + "\t" + string.Join(",", bin.modsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.PepClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ProtClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.FracWithSingle).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + ((double)bin.Overlapping / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.MedianLength).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.UniprotID);
                }
            }
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            // load modifications
            Status("Loading modifications...", taskId);
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
            List<string> localizeableModificationTypes = CommonParameters.LocalizeAll ? GlobalVariables.AllModTypesKnown.ToList() : CommonParameters.ListOfModTypesLocalize.ToList();

            // what types of fragment ions to search for
            List<ProductType> ionTypes = new List<ProductType>();
            if (CommonParameters.BIons && SearchParameters.AddCompIons)
                ionTypes.Add(ProductType.B);
            else if (CommonParameters.BIons)
                ionTypes.Add(ProductType.BnoB1ions);
            if (CommonParameters.YIons)
                ionTypes.Add(ProductType.Y);
            if (CommonParameters.ZdotIons)
                ionTypes.Add(ProductType.Zdot);
            if (CommonParameters.CIons)
                ionTypes.Add(ProductType.C);

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, SearchParameters.SearchTarget, SearchParameters.DecoyType, localizeableModificationTypes);

            // write prose settings
            proseCreatedWhileRunning.Append("The following search settings were used: ");
            proseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            proseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            proseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            proseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            proseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            proseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("max mods per peptide = " + CommonParameters.DigestionParams.MaxModsForPeptide + "; ");
            proseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            proseCreatedWhileRunning.Append("precursor mass tolerance = " + CommonParameters.PrecursorMassTolerance + "; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; ");
            proseCreatedWhileRunning.Append("report PSM ambiguity = " + CommonParameters.ReportAllAmbiguity + ". ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count(p => !p.IsDecoy) + " non-decoy protein entries including " + proteinList.Count(p => p.IsContaminant) + " contaminant sequences. ");
            
            // start the search task
            myTaskResults = new MyTaskResults(this);
            List<PeptideSpectralMatch> allPsms = new List<PeptideSpectralMatch>();
            FlashLFQResults flashLfqResults = null;

            ParallelOptions parallelOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = CommonParameters.MaxParallelFilesToAnalyze
            };

            MyFileManager myFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);

            HashSet<IDigestionParams> ListOfDigestionParams = GetListOfDistinctDigestionParams(CommonParameters, fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b)));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);
            Status("Searching files...", new List<string> { taskId, "Individual Spectra Files" });

            Dictionary<string, int[]> numMs2SpectraPerFile = new Dictionary<string, int[]>();
            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var origDataFile = currentRawFileList[spectraFileIndex];

                // mark the file as in-progress
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);
                
                MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(combinedParams.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                numMs2SpectraPerFile.Add(Path.GetFileNameWithoutExtension(origDataFile), new int[] { myMsDataFile.Count(p => p.MsnOrder == 2), arrayOfMs2ScansSortedByMass.Length });
                myFileManager.DoneWithFile(origDataFile);

                var fileSpecificPsms = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];

                // modern search
                if (SearchParameters.SearchType == SearchType.Modern)
                {
                    for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                    {
                        List<CompactPeptide> peptideIndex = null;
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                        #region Generate indices for modern search

                        Status("Getting fragment dictionary...", new List<string> { taskId });
                        var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, ionTypes, currentPartition, SearchParameters.DecoyType, ListOfDigestionParams, combinedParams, SearchParameters.MaxFragmentSize, new List<string> { taskId });
                        List<int>[] fragmentIndex = null;
                        lock (indexLock)
                            GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, taskId);

                        #endregion Generate indices for modern search

                        Status("Searching files...", taskId);

                        new ModernSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, ionTypes, currentPartition, combinedParams, SearchParameters.AddCompIons, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                    }
                }
                // nonspecific search
                else if (SearchParameters.SearchType == SearchType.NonSpecific)
                {
                    List<List<ProductType>> terminusSeparatedIons = ProductTypeMethod.SeparateIonsByTerminus(ionTypes);
                    foreach (List<ProductType> terminusSpecificIons in terminusSeparatedIons)
                    {
                        for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                        {
                            List<CompactPeptide> peptideIndex = null;
                            List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                            List<int>[] fragmentIndex = new List<int>[1];

                            #region Generate indices for nonspecifc search

                            Status("Getting fragment dictionary...", new List<string> { taskId });
                            var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, terminusSpecificIons, currentPartition, SearchParameters.DecoyType, ListOfDigestionParams, combinedParams, SearchParameters.MaxFragmentSize, new List<string> { taskId });
                            lock (indexLock)
                                GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, taskId);

                            #endregion Generate indices for nonspecifc search

                            #region Generate indices for nonspecific search

                            Status("Getting precursor dictionary...", new List<string> { taskId });
                            List<CompactPeptide> peptideIndexPrecursor = null;
                            List<Protein> proteinListSubsetPrecursor = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));
                            List<int>[] fragmentIndexPrecursor = new List<int>[1];
                            var indexEnginePrecursor = new PrecursorIndexingEngine(proteinListSubsetPrecursor, variableModifications, fixedModifications, terminusSpecificIons, currentPartition, SearchParameters.DecoyType, ListOfDigestionParams, combinedParams, 0, new List<string> { taskId });
                            lock (indexLock)
                                GenerateIndexes(indexEnginePrecursor, dbFilenameList, ref peptideIndexPrecursor, ref fragmentIndexPrecursor, taskId);

                            if (peptideIndex.Count != peptideIndexPrecursor.Count)
                                throw new MetaMorpheusException("peptideIndex not identical between indexing engines");

                            #endregion Generate indices for nonspecific search

                            Status("Searching files...", taskId);

                            new NonSpecificEnzymeSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, fragmentIndexPrecursor, terminusSpecificIons, currentPartition, combinedParams, SearchParameters.AddCompIons, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                            ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                        }
                    }
                }
                // classic search
                else
                {
                    Status("Starting search...", thisId);
                    new ClassicSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, ionTypes, massDiffAcceptor, SearchParameters.AddCompIons, combinedParams, thisId).Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                }
                lock (psmLock)
                {
                    allPsms.AddRange(fileSpecificPsms);
                }

                completedFiles++;
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            });
            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            // Group and order psms
            Status("Matching peptides to proteins...", taskId);
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            if (proteinList.Any())
            {
                if (SearchParameters.SearchType == SearchType.NonSpecific)
                {
                    List<List<ProductType>> terminusSeparatedIons = ProductTypeMethod.SeparateIonsByTerminus(ionTypes);
                    MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
                    foreach (List<ProductType> terminusSpecificIons in terminusSeparatedIons)
                        new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, allPsms, proteinList, fixedModifications, variableModifications, terminusSpecificIons, ListOfDigestionParams, massDiffAcceptor, CommonParameters.ReportAllAmbiguity, new List<string> { taskId }).Run();
                }
                else
                {
                    SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, ionTypes, ListOfDigestionParams, CommonParameters.ReportAllAmbiguity, new List<string> { taskId });
                    var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                    compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;
                }
            }

            ProteinParsimonyResults proteinAnalysisResults = null;
            if (SearchParameters.DoParsimony)
            {
                proteinAnalysisResults = (ProteinParsimonyResults)(new ProteinParsimonyEngine(compactPeptideToProteinPeptideMatching, SearchParameters.ModPeptidesAreDifferent, new List<string> { taskId }).Run());
            }

            Status("Resolving most probable peptide...", new List<string> { taskId });
            foreach (var huh in allPsms)
            {
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
            }
            
            Status("Running FDR analysis...", taskId);
            int massDiffAcceptorNumNotches = GetNumNotches(SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            var fdrAnalysisResults = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsms, massDiffAcceptorNumNotches, CommonParameters, new List<string> { taskId }).Run());

            //sort the list of psms by the score used for fdr calculation
            if (fdrAnalysisResults.DeltaScoreImprovement)
            {
                allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.DeltaScore).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
            }
            else
            {
                allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
            }

            List<EngineLayer.ProteinGroup> proteinGroups = null;

            if (SearchParameters.DoParsimony)
            {
                var proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(proteinAnalysisResults.ProteinGroups, allPsms, SearchParameters.NoOneHitWonders, SearchParameters.ModPeptidesAreDifferent, true, new List<string> { taskId }).Run();
                proteinGroups = proteinScoringAndFdrResults.sortedAndScoredProteinGroups;
            }

            if (SearchParameters.DoLocalizationAnalysis)
            {
                Status("Running localization analysis...", taskId);
                Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);
                    
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    Status("Running localization analysis...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                    new LocalizationEngine(allPsms.Where(b => b.FullFilePath.Equals(origDataFile)).ToList(), ionTypes, myMsDataFile, combinedParams.ProductMassTolerance, new List<string> { taskId, "Individual Spectra Files", origDataFile }, SearchParameters.AddCompIons).Run();
                    myFileManager.DoneWithFile(origDataFile);
                    ReportProgress(new ProgressEventArgs(100, "Done with localization analysis!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
                });
            }

            new ModificationAnalysisEngine(allPsms, new List<string> { taskId }).Run();

            if (SearchParameters.DoQuantification)
            {
                // pass quantification parameters to FlashLFQ
                Status("Quantifying...", taskId);

                // construct file info for FlashLFQ
                var rawfileinfos = new List<RawFileInfo>();
                foreach (var file in currentRawFileList)
                {
                    if (myFileManager.SeeIfOpen(file))
                        rawfileinfos.Add(new RawFileInfo(file, myFileManager.LoadFile(file, null, null, false, false)));
                    else
                        rawfileinfos.Add(new RawFileInfo(file));
                    myFileManager.DoneWithFile(file);
                }

                // get PSMs to pass to FlashLFQ
                var unambiguousPsmsBelowOnePercentFdr = allPsms.Where(p => p.FdrInfo.QValue < 0.01 && !p.IsDecoy && p.FullSequence != null);
                var mypsmsGroupedByFile = unambiguousPsmsBelowOnePercentFdr.GroupBy(p => p.FullFilePath);

                // pass protein group info for each PSM
                Dictionary<PeptideSpectralMatch, List<string>> psmToProteinGroupNames = new Dictionary<PeptideSpectralMatch, List<string>>();
                if (proteinGroups != null)
                {
                    EngineLayer.ProteinGroup.FilesForQuantification = rawfileinfos.Select(p => p.fullFilePathWithExtension).ToArray();

                    foreach (var proteinGroup in proteinGroups)
                    {
                        foreach (var psm in proteinGroup.AllPsmsBelowOnePercentFDR)
                        {
                            if (psmToProteinGroupNames.TryGetValue(psm, out List<string> proteinGroupNames))
                            {
                                proteinGroupNames.Add(proteinGroup.ProteinGroupName);
                            }
                            else
                            {
                                psmToProteinGroupNames.Add(psm, new List<string> { proteinGroup.ProteinGroupName });
                            }
                        }
                    }
                }
                else
                {
                    // if protein groups were not constructed, just use accession numbers
                    foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
                    {
                        var proteinsAccessionString = psm.CompactPeptides.SelectMany(b => b.Value.Item2).Select(b => b.Protein.Accession).Distinct();
                        psmToProteinGroupNames.Add(psm, proteinsAccessionString.ToList());
                    }
                }

                // some PSMs may not have protein groups (if 2 peptides are required to construct a protein group, some PSMs will be left over)
                // the peptides should still be quantified but not considered for protein quantification
                foreach (var psm in unambiguousPsmsBelowOnePercentFdr)
                {
                    if (!psmToProteinGroupNames.ContainsKey(psm))
                    {
                        psmToProteinGroupNames.Add(psm, new List<string>() { "" });
                    }
                }

                // pass PSM info to FlashLFQ
                var flashLFQIdentifications = new List<Identification>();
                foreach (var file in mypsmsGroupedByFile)
                {
                    var rawfileinfo = rawfileinfos.Where(p => p.fullFilePathWithExtension.Equals(file.Key)).First();
                    foreach (var psm in file)
                    {
                        flashLFQIdentifications.Add(new Identification(rawfileinfo, psm.BaseSequence, psm.FullSequence, (double)psm.PeptideMonisotopicMass, psm.ScanRetentionTime, psm.ScanPrecursorCharge, psmToProteinGroupNames[psm]));
                    }
                }

                // run FlashLFQ
                var FlashLfqEngine = new FlashLFQEngine(flashLFQIdentifications, SearchParameters.QuantifyPpmTol, 5.0, SearchParameters.MatchBetweenRuns, 5.0, false, 2, false, true, true, GlobalVariables.ElementsLocation);
                if (flashLFQIdentifications.Any())
                {
                    flashLfqResults = FlashLfqEngine.Run();
                }

                // get protein intensity back from FlashLFQ
                if (proteinGroups != null && flashLfqResults != null)
                {
                    Dictionary<string, EngineLayer.ProteinGroup> proteinGroupNameToProteinGroup = new Dictionary<string, EngineLayer.ProteinGroup>();
                    foreach (var proteinGroup in proteinGroups)
                    {
                        proteinGroup.IntensitiesByFile = new double[EngineLayer.ProteinGroup.FilesForQuantification.Length];
                        if (!proteinGroupNameToProteinGroup.ContainsKey(proteinGroup.ProteinGroupName))
                        {
                            proteinGroupNameToProteinGroup.Add(proteinGroup.ProteinGroupName, proteinGroup);
                        }
                    }

                    foreach (var flashLfqProteinGroup in flashLfqResults.proteinGroups)
                    {
                        if (proteinGroupNameToProteinGroup.TryGetValue(flashLfqProteinGroup.Key, out EngineLayer.ProteinGroup metamorpheusProteinGroup))
                        {
                            for (int i = 0; i < EngineLayer.ProteinGroup.FilesForQuantification.Length; i++)
                            {
                                metamorpheusProteinGroup.IntensitiesByFile[i] = flashLfqProteinGroup.Value.intensities[rawfileinfos[i]];
                            }
                        }
                    }
                }
            }

            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            if (SearchParameters.DoHistogramAnalysis)
            {
                var limitedpsms_with_fdr = allPsms.Where(b => (b.FdrInfo.QValue <= 0.01)).ToList();
                if (limitedpsms_with_fdr.Any(b => !b.IsDecoy))
                {
                    Status("Running histogram analysis...", new List<string> { taskId });
                    var myTreeStructure = new BinTreeStructure();
                    myTreeStructure.GenerateBins(limitedpsms_with_fdr, SearchParameters.HistogramBinTolInDaltons);
                    var writtenFile = Path.Combine(OutputFolder, "aggregate.mytsv");
                    WriteTree(myTreeStructure, writtenFile);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });
                }
            }

            Status("Writing results...", taskId);
            {
                if (currentRawFileList.Count > 1)
                {
                    var writtenFile = Path.Combine(OutputFolder, "aggregatePSMs.psmtsv");
                    WritePsmsToTsv(allPsms, writtenFile, SearchParameters.ModsToWriteSelection);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });

                    var writtenFileForPercolator = Path.Combine(OutputFolder, "forPercolator.tsv");
                    WritePsmsForPercolator(allPsms, writtenFileForPercolator);
                    SucessfullyFinishedWritingFile(writtenFileForPercolator, new List<string> { taskId });
                }
                myTaskResults.AddNiceText("All target PSMS within 1% FDR: " + allPsms.Count(a => a.FdrInfo.QValue < .01 && !a.IsDecoy));
            }

            var peptides = allPsms.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
            {
                if (currentRawFileList.Count > 1)
                {
                    var writtenFile = Path.Combine(OutputFolder, "aggregatePeptides.psmtsv");
                    WritePsmsToTsv(peptides, writtenFile, SearchParameters.ModsToWriteSelection);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId });
                }
                myTaskResults.AddNiceText("Target peptides within 1% FDR: " + peptides.Count(a => a.FdrInfo.QValue < 0.01 && !a.IsDecoy));
            }

            if (SearchParameters.DoParsimony)
            {
                myTaskResults.AddNiceText("Target protein groups within 1% FDR: " + proteinGroups.Count(b => b.QValue < 0.01 && !b.isDecoy) + Environment.NewLine);
            }

            var psmsGroupedByFile = allPsms.GroupBy(p => p.FullFilePath);

            // individual psm files (with global psm fdr, global parsimony)
            foreach (var group in psmsGroupedByFile)
            {
                var psmsForThisFile = group.ToList();

                var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);

                {
                    var writtenFile = Path.Combine(OutputFolder, strippedFileName + "_PSMs.psmtsv");
                    WritePsmsToTsv(psmsForThisFile, writtenFile, SearchParameters.ModsToWriteSelection);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
                    myTaskResults.AddNiceText("Target PSMs within 1% FDR in " + strippedFileName + ": " + psmsForThisFile.Count(a => a.FdrInfo.QValue < .01 && a.IsDecoy == false));
                    myTaskResults.AddNiceText("MS2 spectra in " + strippedFileName + ": " + numMs2SpectraPerFile[strippedFileName][0]);
                    myTaskResults.AddNiceText("Number of precursor species fragmented in " + strippedFileName + ": " + numMs2SpectraPerFile[strippedFileName][1]);
                }

                var writtenFileForPercolator = Path.Combine(OutputFolder, strippedFileName + "_forPercolator.tsv");
                WritePsmsForPercolator(psmsForThisFile, writtenFileForPercolator);
                SucessfullyFinishedWritingFile(writtenFileForPercolator, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
            }
            foreach (var group in psmsGroupedByFile)
            {
                var psmsForThisFile = group.ToList();

                var strippedFileName = Path.GetFileNameWithoutExtension(group.First().FullFilePath);

                {
                    var peptidesForFile = psmsForThisFile.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
                    var writtenFile = Path.Combine(OutputFolder, strippedFileName + "_Peptides.psmtsv");
                    WritePsmsToTsv(peptidesForFile, writtenFile, SearchParameters.ModsToWriteSelection);
                    SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId, "Individual Spectra Files", group.First().FullFilePath });
                    myTaskResults.AddNiceText("Target peptides within 1% FDR in " + strippedFileName + ": " + peptidesForFile.Count(a => a.FdrInfo.QValue < .01 && a.IsDecoy == false));
                }
            }

            if (SearchParameters.DoParsimony)
            {
                if (currentRawFileList.Count > 1)
                    WriteProteinGroupsToTsv(proteinGroups, OutputFolder, "aggregateProteinGroups", new List<string> { taskId }, psmsGroupedByFile.Select(b => b.Key).ToList());

                // individual protein group files (local protein fdr, global parsimony, global psm fdr)
                foreach (var fullFilePath in currentRawFileList)
                {
                    List<PeptideSpectralMatch> psmsForThisFile = psmsGroupedByFile.Where(p => p.Key == fullFilePath).SelectMany(g => g).ToList();

                    var strippedFileName = Path.GetFileNameWithoutExtension(fullFilePath);

                    var subsetProteinGroupsForThisFile = new List<EngineLayer.ProteinGroup>();
                    foreach (var pg in proteinGroups)
                        subsetProteinGroupsForThisFile.Add(pg.ConstructSubsetProteinGroup(fullFilePath));

                    var subsetProteinScoringAndFdrResults = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(subsetProteinGroupsForThisFile, psmsForThisFile, SearchParameters.NoOneHitWonders, SearchParameters.ModPeptidesAreDifferent, false, new List<string> { taskId, "Individual Spectra Files", fullFilePath }).Run();
                    subsetProteinGroupsForThisFile = subsetProteinScoringAndFdrResults.sortedAndScoredProteinGroups;

                    myTaskResults.AddNiceText("Target protein groups within 1 % FDR in " + strippedFileName + ": " + subsetProteinGroupsForThisFile.Count(b => b.QValue < 0.01 && !b.isDecoy));

                    WriteProteinGroupsToTsv(subsetProteinGroupsForThisFile, OutputFolder, strippedFileName + "_ProteinGroups", new List<string> { taskId, "Individual Spectra Files", fullFilePath }, new List<string> { fullFilePath });

                    Status("Writing mzid...", new List<string> { taskId, "Individual Spectra Files", fullFilePath });
                    var mzidFilePath = Path.Combine(OutputFolder, strippedFileName + ".mzid");
                    MzIdentMLWriter.WriteMzidentml(psmsForThisFile, subsetProteinGroupsForThisFile, variableModifications, fixedModifications, new List<Protease> { CommonParameters.DigestionParams.Protease }, 0.01, CommonParameters.ProductMassTolerance, CommonParameters.PrecursorMassTolerance, CommonParameters.DigestionParams.MaxMissedCleavages, mzidFilePath);
                    SucessfullyFinishedWritingFile(mzidFilePath, new List<string> { taskId, "Individual Spectra Files", fullFilePath });

                    ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", fullFilePath }));
                }
            }

            if (SearchParameters.DoQuantification && flashLfqResults != null)
            {
                foreach (var file in flashLfqResults.peaks)
                    WritePeakQuantificationResultsToTsv(file.Value, OutputFolder, file.Key.filenameWithoutExtension + "_QuantifiedPeaks", new List<string> { taskId, "Individual Spectra Files", file.Key.fullFilePathWithExtension });

                if (currentRawFileList.Count > 1)
                    WritePeakQuantificationResultsToTsv(flashLfqResults.peaks.SelectMany(p => p.Value).ToList(), OutputFolder, "aggregateQuantifiedPeaks", new List<string> { taskId });

                WritePeptideQuantificationResultsToTsv(flashLfqResults.peptideBaseSequences.Select(p => p.Value).OrderBy(p => p.Sequence).ToList(), OutputFolder, "aggregateQuantifiedPeptidesByBaseSeq", new List<string> { taskId });
                WritePeptideQuantificationResultsToTsv(flashLfqResults.peptideModifiedSequences.Select(p => p.Value).OrderBy(p => p.Sequence).ToList(), OutputFolder, "aggregateQuantifiedPeptidesByFullSeq", new List<string> { taskId });
            }

            // write pruned database
            if (SearchParameters.WritePrunedDatabase)
            {
                Status("Writing Pruned Database...", new List<string> { taskId });
                List<Modification> modificationsToWriteIfBoth = new List<Modification>();
                List<Modification> modificationsToWriteIfInDatabase = new List<Modification>();
                List<Modification> modificationsToWriteIfObserved = new List<Modification>();

                var confidentPsms = allPsms.Where(b => b.FdrInfo.QValueNotch < 0.01 && b.FdrInfo.QValue < 0.01 && !b.IsDecoy && b.BaseSequence != null).ToList();
                var proteinToConfidentBaseSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                // associate all confident PSMs with all possible proteins they could be digest products of (before or after parsimony)
                foreach (PeptideSpectralMatch psm in confidentPsms)
                {
                    var myPepsWithSetMods = psm.CompactPeptides.SelectMany(p => p.Value.Item2);

                    foreach (var peptide in myPepsWithSetMods)
                    {
                        if (proteinToConfidentBaseSequences.TryGetValue(peptide.Protein, out var myPepList))
                            myPepList.Add(peptide);
                        else
                            proteinToConfidentBaseSequences.Add(peptide.Protein, new List<PeptideWithSetModifications> { peptide });
                    }
                }
                
                // Add user mod selection behavours to Pruned DB
                foreach (var modType in SearchParameters.ModsToWriteSelection)
                {
                    if (modType.Value == 1) // Write if observed and in database
                        modificationsToWriteIfBoth.AddRange(GlobalVariables.AllModsKnown.Where(b => b.modificationType.Equals(modType.Key)));
                    if (modType.Value == 2) // Write if in database
                        modificationsToWriteIfInDatabase.AddRange(GlobalVariables.AllModsKnown.Where(b => b.modificationType.Equals(modType.Key)));
                    if (modType.Value == 3) // Write if observed
                        modificationsToWriteIfObserved.AddRange(GlobalVariables.AllModsKnown.Where(b => b.modificationType.Equals(modType.Key)));
                }
                //generates dictionary of proteins with only localized modifications
                var ModPsms = allPsms.Where(b => b.FdrInfo.QValueNotch < 0.01 && b.FdrInfo.QValue < 0.01 && !b.IsDecoy && b.BaseSequence != null && b.FullSequence!=null).ToList();
                var proteinToConfidentModifiedSequences = new Dictionary<Protein, List<PeptideWithSetModifications>>();

                foreach (PeptideSpectralMatch psm in ModPsms)
                {
                    var myPepsWithSetMods = psm.CompactPeptides.SelectMany(p => p.Value.Item2);

                    foreach (var peptide in myPepsWithSetMods)
                    {
                        if (proteinToConfidentModifiedSequences.TryGetValue(peptide.Protein, out var myPepList))
                            myPepList.Add(peptide);
                        else
                            proteinToConfidentModifiedSequences.Add(peptide.Protein, new List<PeptideWithSetModifications> { peptide });
                    }
                }

                // mods included in pruned database will only be confidently localized mods (peptide's FullSequence != null)
                foreach (var protein in proteinList)
                {
                    if (!protein.IsDecoy)
                    {
                        HashSet<Tuple<int, ModificationWithMass>> modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>();
                        if (proteinToConfidentModifiedSequences.ContainsKey(protein))
                        {
                                modsObservedOnThisProtein = new HashSet<Tuple<int, ModificationWithMass>>(proteinToConfidentModifiedSequences[protein].SelectMany(b => b.allModsOneIsNterminus.Select(c => new Tuple<int, ModificationWithMass>(GetOneBasedIndexInProtein(c.Key, b), c.Value))));
                        }
                                                
                        IDictionary<int, List<Modification>> modsToWrite = new Dictionary<int, List<Modification>>();

                        foreach (var observedMod in modsObservedOnThisProtein)
                        {
                            //Add if observed (regardless if in database)
                            var tempMod = observedMod.Item2;

                            if (modificationsToWriteIfObserved.Contains(tempMod as Modification))
                            {
                                if (!modsToWrite.ContainsKey(observedMod.Item1))
                                    modsToWrite.Add(observedMod.Item1, new List<Modification> { observedMod.Item2 as Modification });
                                else
                                    modsToWrite[observedMod.Item1].Add(observedMod.Item2 as Modification);
                                continue;
                            }
                        }

                        // Add if in database (two cases: always or if observed)
                        foreach (var modd in protein.OneBasedPossibleLocalizedModifications)
                            foreach (var mod in modd.Value)
                            {
                                //Add if always In Database or if was observed and in database and not set to not include
                                if (modificationsToWriteIfInDatabase.Contains(mod as Modification) ||
                                (modsObservedOnThisProtein.Contains(new Tuple<int, ModificationWithMass>(modd.Key, mod as ModificationWithMass)) && modificationsToWriteIfBoth.Contains(mod as Modification)))
                                {
                                    if (!modsToWrite.ContainsKey(modd.Key))
                                        modsToWrite.Add(modd.Key, new List<Modification> { mod });
                                    else
                                        modsToWrite[modd.Key].Add(mod);
                                }
                            }
                      

                        
                       if (proteinToConfidentBaseSequences.TryGetValue(protein, out var peptideSequences))
                       { 
                            // removes all annotated mods on proteins                           
                            if (protein.Accession == protein.Accession)
                            {
                                protein.OneBasedPossibleLocalizedModifications.Clear();
                                // adds confidently localized and identified mods
                                foreach (var kvp in modsToWrite)
                                    protein.OneBasedPossibleLocalizedModifications.Add(kvp);
                            }
                            
                            
                       }
                        
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

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);

                    SucessfullyFinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });
                }
                if (dbFilenameList.Any(b => b.IsContaminant))
                {
                    string outputXMLdbFullNameContaminants = Path.Combine(OutputFolder, string.Join("-", dbFilenameList.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");

                    ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);

                    SucessfullyFinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });
                }
            }

            return myTaskResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private static MassDiffAcceptor ParseSearchMode(string text)
        {
            MassDiffAcceptor ye = null;

            var split = text.Split(' ');

            switch (split[1])
            {
                case "dot":

                    var massShifts = Array.ConvertAll(split[4].Split(','), Double.Parse);
                    var newString = split[2].Replace("�", "");
                    var toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                    if (split[3].ToUpperInvariant().Equals("PPM"))
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new PpmTolerance(toleranceValue));
                    else if (split[3].ToUpperInvariant().Equals("DA"))
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new AbsoluteTolerance(toleranceValue));
                    break;

                case "interval":
                    IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[0], CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[1], CultureInfo.InvariantCulture)));
                    ye = new IntervalMassDiffAcceptor(split[0], doubleRanges);
                    break;

                case "OpenSearch":
                    ye = new OpenSearchMode();
                    break;

                case "daltonsAroundZero":
                    ye = new SingleAbsoluteAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                case "ppmAroundZero":
                    ye = new SinglePpmAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                default:
                    throw new MetaMorpheusException("Could not parse search mode string");
            }
            return ye;
        }

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

        private static void WritePeptideIndex(List<CompactPeptide> peptideIndex, string peptideIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(peptideIndexFile))
            {
                ser.Serialize(file, peptideIndex);
            }
        }

        private static void WriteFragmentIndexNetSerializer(List<int>[] fragmentIndex, string fragmentIndexFile)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<int>[]));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFile))
                ser.Serialize(file, fragmentIndex);
        }

        private static string GetExistingFolderWithIndices(IndexingEngine indexEngine, List<DbForTask> dbFilenameList)
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

        private static void WriteIndexEngineParams(IndexingEngine indexEngine, string fileName)
        {
            using (StreamWriter output = new StreamWriter(fileName))
            {
                output.Write(indexEngine);
            }
        }

        private static string GenerateOutputFolderForIndices(List<DbForTask> dbFilenameList)
        {
            var folder = Path.Combine(Path.GetDirectoryName(dbFilenameList.First().FilePath), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
            if (!Directory.Exists(folder))
                Directory.CreateDirectory(folder);
            return folder;
        }

        private static int GetOneBasedIndexInProtein(int oneIsNterminus, PeptideWithSetModifications peptideWithSetModifications)
        {
            if (oneIsNterminus == 1)
                return peptideWithSetModifications.OneBasedStartResidueInProtein;
            if (oneIsNterminus == peptideWithSetModifications.Length + 2)
                return peptideWithSetModifications.OneBasedEndResidueInProtein;
            return peptideWithSetModifications.OneBasedStartResidueInProtein + oneIsNterminus - 2;
        }

        private void WritePsmsForPercolator(List<PeptideSpectralMatch> psmList, string writtenFileForPercolator)
        {
            using (StreamWriter output = new StreamWriter(writtenFileForPercolator))
            {
                output.WriteLine("SpecId\tLabel\tScanNr\tF1\tF2\tPeptide\tProteins");
                output.WriteLine("DefaultDirection\t-\t-\t1\t1\t\t");
                for (int i = 0; i < psmList.Count; i++)
                {
                    var heh = psmList[i];

                    output.Write(i.ToString());
                    output.Write('\t' + (heh.IsDecoy ? -1 : 1).ToString());
                    output.Write('\t' + heh.ScanNumber.ToString());

                    // Features
                    {
                        output.Write('\t' + string.Join("\t", heh.Features));
                    }

                    // HACKY: Ignores all ambiguity
                    var pwsm = heh.CompactPeptides.First().Value.Item2.First();

                    output.Write('\t' + (pwsm.PreviousAminoAcid + "." + pwsm.Sequence + "." + pwsm.NextAminoAcid).ToString());
                    output.Write('\t' + (pwsm.Protein.Accession).ToString());
                    output.WriteLine();
                }
            }
        }

        private int GetNumNotches(MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact:
                    return 1;

                case MassDiffAcceptorType.OneMM:
                    return 2;

                case MassDiffAcceptorType.TwoMM:
                    return 3;

                case MassDiffAcceptorType.ThreeMM:
                    return 4;

                case MassDiffAcceptorType.ModOpen:
                    return 1;

                case MassDiffAcceptorType.Open:
                    return 1;

                case MassDiffAcceptorType.Custom:
                    return ParseSearchMode(customMdac).NumNotches;

                default:
                    throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        private void GenerateIndexes(IndexingEngine indexEngine, List<DbForTask> dbFilenameList, ref List<CompactPeptide> peptideIndex, ref List<int>[] fragmentIndex, string taskId)
        {
            string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);
            if (pathToFolderWithIndices == null)
            {
                var output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
                Status("Writing params...", new List<string> { taskId });
                var paramsFile = Path.Combine(output_folderForIndices, "indexEngine.params");
                WriteIndexEngineParams(indexEngine, paramsFile);
                SucessfullyFinishedWritingFile(paramsFile, new List<string> { taskId });

                Status("Running Index Engine...", new List<string> { taskId });
                var indexResults = (IndexingResults)indexEngine.Run();
                peptideIndex = indexResults.PeptideIndex;
                fragmentIndex = indexResults.FragmentIndex;

                Status("Writing peptide index...", new List<string> { taskId });
                var peptideIndexFile = Path.Combine(output_folderForIndices, "peptideIndex.ind");
                WritePeptideIndex(peptideIndex, peptideIndexFile);
                SucessfullyFinishedWritingFile(peptideIndexFile, new List<string> { taskId });

                Status("Writing fragment index...", new List<string> { taskId });
                var fragmentIndexFile = Path.Combine(output_folderForIndices, "fragmentIndex.ind");
                WriteFragmentIndexNetSerializer(fragmentIndex, fragmentIndexFile);
                SucessfullyFinishedWritingFile(fragmentIndexFile, new List<string> { taskId });
            }
            else
            {
                Status("Reading peptide index...", new List<string> { taskId });
                var messageTypes = GetSubclassesAndItself(typeof(List<CompactPeptide>));
                var ser = new NetSerializer.Serializer(messageTypes);
                using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "peptideIndex.ind")))
                {
                    peptideIndex = (List<CompactPeptide>)ser.Deserialize(file);
                }

                Status("Reading fragment index...", new List<string> { taskId });
                messageTypes = GetSubclassesAndItself(typeof(List<int>[]));
                ser = new NetSerializer.Serializer(messageTypes);
                using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "fragmentIndex.ind")))
                    fragmentIndex = (List<int>[])ser.Deserialize(file);
            }
        }

        private void WriteProteinGroupsToTsv(List<EngineLayer.ProteinGroup> items, string outputFolder, string strippedFileName, List<string> nestedIds, List<string> FileNames)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, strippedFileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(EngineLayer.ProteinGroup.GetTabSeparatedHeader(FileNames.Count == 1));
                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        private void WritePeptideQuantificationResultsToTsv(List<FlashLFQ.Peptide> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(FlashLFQ.Peptide.TabSeparatedHeader);

                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        private void WritePeakQuantificationResultsToTsv(List<FlashLFQ.ChromatographicPeak> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(FlashLFQ.ChromatographicPeak.TabSeparatedHeader);

                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }
                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        #endregion Private Methods
    }
}