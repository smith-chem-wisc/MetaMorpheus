using EngineLayer;
using EngineLayer.CrosslinkAnalysis;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
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

        #endregion Private Fields

        #region Public Constructors

        public XLSearchTask() : base(MyTask.XLSearch)
        {
            CommonParameters = new CommonParameters();
            CommonParameters.DoPrecursorDeconvolution = false;
            CommonParameters.ConserveMemory = false;
            XlSearchParameters = new XlSearchParameters();
        }

        #endregion Public Constructors

        #region Public Properties

        public XlSearchParameters XlSearchParameters { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(TaskType.ToString());
            sb.AppendLine("The initiator methionine behavior is set to "
                + CommonParameters.DigestionParams.InitiatorMethionineBehavior
                + " and the maximum number of allowed missed cleavages is "
                + CommonParameters.DigestionParams.MaxMissedCleavages);
            sb.AppendLine("MinPeptideLength: " + CommonParameters.DigestionParams.MinPeptideLength);
            sb.AppendLine("MaxPeptideLength: " + CommonParameters.DigestionParams.MaxPeptideLength);
            sb.AppendLine("maxModificationIsoforms: " + CommonParameters.DigestionParams.MaxModificationIsoforms);
            sb.AppendLine("protease: " + CommonParameters.DigestionParams.Protease);
            sb.AppendLine("bIons: " + CommonParameters.BIons);
            sb.AppendLine("yIons: " + CommonParameters.YIons);
            sb.AppendLine("cIons: " + CommonParameters.CIons);
            sb.AppendLine("zdotIons: " + CommonParameters.ZdotIons);

            sb.AppendLine("Fixed mod lists: " + string.Join(",", CommonParameters.ListOfModsFixed));
            sb.AppendLine("Variable mod lists: " + string.Join(",", CommonParameters.ListOfModsVariable));
            sb.AppendLine("Localized mod lists: " + string.Join(",", CommonParameters.ListOfModsLocalize));
            sb.AppendLine("searchDecoy: " + XlSearchParameters.SearchDecoy);
            sb.AppendLine("productMassTolerance: " + CommonParameters.ProductMassTolerance);

            sb.AppendLine("Crosslink Precusor mass tolerance: " + XlSearchParameters.XlPrecusorMsTl);
            sb.AppendLine("Beta Precusor mass tolerance: " + XlSearchParameters.XlBetaPrecusorMsTl);

            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList)
        {
            myTaskResults = new MyTaskResults(this);
            List<PsmCross> allPsms = new List<PsmCross>();
            var compactPeptideToProteinPeptideMatch = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...", taskId);

            #region Load modifications

            List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> localizeableModifications;
            if (CommonParameters.LocalizeAll)
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
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
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, XlSearchParameters.SearchDecoy, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> unknownModifications)).ToList();

            List<ProductType> ionTypes = new List<ProductType>();
            if (CommonParameters.BIons)
                ionTypes.Add(ProductType.B);
            if (CommonParameters.YIons)
                ionTypes.Add(ProductType.Y);
            if (CommonParameters.ZdotIons)
                ionTypes.Add(ProductType.Zdot);
            if (CommonParameters.CIons)
                ionTypes.Add(ProductType.C);
            TerminusType terminusType = ProductTypeToTerminusType.IdentifyTerminusType(ionTypes);

            var crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(XlSearchParameters.CrosslinkerType);
            if (XlSearchParameters.CrosslinkerType == CrosslinkerType.UserDefined)
            {
                crosslinker.CrosslinkerName = XlSearchParameters.UdXLkerName;
                crosslinker.Cleavable = XlSearchParameters.UdXLkerCleavable;
                crosslinker.TotalMass = (double)XlSearchParameters.UdXLkerTotalMass;
                crosslinker.CleaveMassShort = (double)XlSearchParameters.UdXLkerShortMass;
                crosslinker.CleaveMassLong = (double)XlSearchParameters.UdXLkerLongMass;
                crosslinker.CrosslinkerModSite = XlSearchParameters.UdXLkerResidue;
            }
            ParallelOptions parallelOptions = new ParallelOptions();
            if (CommonParameters.MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxDegreeOfParallelism.Value;
            MyFileManager myFileManager = new MyFileManager(XlSearchParameters.DisposeOfFileWhenDone);

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);

            proseCreatedWhileRunning.Append("The following crosslink discovery were used: ");
            proseCreatedWhileRunning.Append("crosslinker name = " + crosslinker.CrosslinkerName + "; ");
            proseCreatedWhileRunning.Append("crosslinker type = " + crosslinker.Cleavable + "; ");
            proseCreatedWhileRunning.Append("crosslinker mass = " + crosslinker.TotalMass + "; ");
            proseCreatedWhileRunning.Append("crosslinker modification site(s) = " + crosslinker.CrosslinkerModSite + "; ");

            proseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            proseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            proseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            if (CommonParameters.DigestionParams.MaxPeptideLength == null)
            {
                proseCreatedWhileRunning.Append("maximum peptide length = unspecified; ");
            }
            else
            {
                proseCreatedWhileRunning.Append("maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            }
            proseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            proseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");

            proseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");

            proseCreatedWhileRunning.Append("parent mass tolerance(s) = UNKNOWN; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                List<PsmCross> newPsms = new List<PsmCross>();
                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile, null, null, false, false);
                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, CommonParameters.DoPrecursorDeconvolution, CommonParameters.UseProvidedPrecursorInfo, CommonParameters.DeconvolutionIntensityRatio, CommonParameters.DeconvolutionMaxAssumedChargeState, CommonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<CompactPeptide> peptideIndex = null;
                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / CommonParameters.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / CommonParameters.TotalPartitions) - (currentPartition * proteinList.Count() / CommonParameters.TotalPartitions));

                    float[] keys = null;
                    List<int>[] fragmentIndex = null;

                    #region Generate indices for modern search

                    Status("Getting fragment dictionary...", new List<string> { taskId });
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, ionTypes, currentPartition, XlSearchParameters.SearchDecoy, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.TotalPartitions, new List<string> { taskId });

                    Dictionary<float, List<int>> fragmentIndexDict;
                    lock (indexLock)
                    {
                        string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);

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
                    }
                    keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
                    fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();

                    #endregion Generate indices for modern search

                    Status("Searching files...", taskId);

                    if (XlSearchParameters.CrosslinkSearchWithAllBeta)
                        new CrosslinkSearchEngine2(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, crosslinker, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.CrosslinkSearchWithAllBeta, XlSearchParameters.XlPrecusorMsTl, XlSearchParameters.XlBetaPrecusorMsTl, modsDictionary, ionTypes, proteinList, variableModifications, fixedModifications, CommonParameters, thisId).Run();
                    else
                        new CrosslinkSearchEngine(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, CommonParameters.ProductMassTolerance, crosslinker, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.CrosslinkSearchWithAllBeta, XlSearchParameters.XlPrecusorMsTl, XlSearchParameters.XlBetaPrecusorMsTl, ionTypes, modsDictionary, thisId).Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                }

                lock (psmLock)
                {
                    allPsms.AddRange(newPsms);
                }

                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            });

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            Status("Crosslink analysis engine", taskId);
            MetaMorpheusEngineResults allcrosslinkanalysisResults;
            allcrosslinkanalysisResults = new CrosslinkAnalysisEngine(allPsms, compactPeptideToProteinPeptideMatch, proteinList, variableModifications, fixedModifications, ionTypes, modsDictionary, OutputFolder, crosslinker, terminusType, CommonParameters, new List<string> { taskId }).Run();
            if (XlSearchParameters.XlOutAll)
            {
                WriteAllToTsv(allPsms, OutputFolder, "allPsms", new List<string> { taskId });

            }
            var allPsmsXL = allPsms.Where(p => p.CrossType == PsmCrossType.Cross).Where(p => p.XLBestScore >= 3 && p.BetaPsmCross.XLBestScore >= 3).OrderByDescending(p => p.XLTotalScore).ToList();
            var allPsmsXLFDR = CrosslinkFDRAnalysis(allPsmsXL).ToList();
            //Write Inter Psms FDR
            var interPsmsXL = allPsmsXL.Where(p => p.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() != p.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).OrderByDescending(p => p.XLTotalScore).ToList();
            //var interPsmsXLFDR = interPsmsXL.Where(p => p.XLBestScore >= 2 && p.BetaPsmCross.XLBestScore >= 2).ToList();
            //interPsmsXLFDR = CrosslinkDoFalseDiscoveryRateAnalysis(interPsmsXLFDR).ToList();
            var interPsmsXLFDR = allPsmsXLFDR.Where(p => p.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() != p.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).ToList();
            if (XlSearchParameters.XlOutCrosslink)
            {
                WriteCrosslinkToTsv(interPsmsXLFDR, OutputFolder, "xl_inter_fdr", new List<string> { taskId });
            }
            if (XlSearchParameters.XlOutCLMSVault)
            {
                var interPsmsXLFDR_CLMSVault = interPsmsXLFDR.Where(p => p.MostProbableProteinInfo.IsDecoy != true && p.BetaPsmCross.MostProbableProteinInfo.IsDecoy != true).ToList();
                WriteCrosslinkToTxtForCLMSVault(interPsmsXLFDR_CLMSVault, OutputFolder, "xl_inter_fdr_CLMSVault", crosslinker, new List<string> { taskId });
            }
            if (XlSearchParameters.XlOutPepXML && interPsmsXLFDR.Count != 0)
            {
                var interPsmsXLFDR_PepXML = interPsmsXLFDR.Where(p => p.MostProbableProteinInfo.IsDecoy != true && p.BetaPsmCross.MostProbableProteinInfo.IsDecoy != true).ToList();
                foreach (var fullFilePath in currentRawFileList)
                {
                    string fileNameNoExtension = Path.GetFileNameWithoutExtension(fullFilePath);
                    WritePepXML_xl(interPsmsXLFDR_PepXML.Where(p => p.FullFilePath == fullFilePath).ToList(), dbFilenameList, variableModifications, fixedModifications, localizeableModifications, OutputFolder, fileNameNoExtension, new List<string> { taskId });
                }
            }
            if (XlSearchParameters.XlOutPercolator)
            {
                var interPsmsXLPercolator = interPsmsXL.Where(p => p.XLBestScore >= 2 && p.BetaPsmCross.XLBestScore >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteCrosslinkToTxtForPercolator(interPsmsXLPercolator, OutputFolder, "xl_inter_perc", crosslinker, new List<string> { taskId });
            }

            //Write Intra Psms FDR
            var intraPsmsXL = allPsmsXL.Where(p => p.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() == p.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).OrderByDescending(p => p.XLTotalScore).ToList();
            //var intraPsmsXLFDR = intraPsmsXL.Where(p => p.XLBestScore >= 2 && p.BetaPsmCross.XLBestScore >= 2).ToList();
            //intraPsmsXLFDR = CrosslinkDoFalseDiscoveryRateAnalysis(intraPsmsXLFDR).ToList();
            var intraPsmsXLFDR = allPsmsXLFDR.Where(p => p.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() == p.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).ToList();
            if (XlSearchParameters.XlOutCrosslink)
            {
                WriteCrosslinkToTsv(intraPsmsXLFDR, OutputFolder, "xl_intra_fdr", new List<string> { taskId });
            }
            if (XlSearchParameters.XlOutCLMSVault)
            {
                var intraPsmsXLFDR_CLMSVault = intraPsmsXLFDR.Where(p => p.MostProbableProteinInfo.IsDecoy != true && p.BetaPsmCross.MostProbableProteinInfo.IsDecoy != true).ToList();
                WriteCrosslinkToTxtForCLMSVault(intraPsmsXLFDR_CLMSVault, OutputFolder, "xl_intra_fdr_CLMSVault", crosslinker, new List<string> { taskId });
            }
            if (XlSearchParameters.XlOutPepXML && intraPsmsXLFDR.Count != 0)
            {
                var intraPsmsXLFDR_PepXML = intraPsmsXLFDR.Where(p => p.MostProbableProteinInfo.IsDecoy != true && p.BetaPsmCross.MostProbableProteinInfo.IsDecoy != true).ToList();
                foreach (var fullFilePath in currentRawFileList)
                {
                    string fileNameNoExtension = Path.GetFileNameWithoutExtension(fullFilePath);
                    WritePepXML_xl(intraPsmsXLFDR_PepXML.Where(p => p.FullFilePath == fullFilePath).ToList(), dbFilenameList, variableModifications, fixedModifications, localizeableModifications, OutputFolder, fileNameNoExtension, new List<string> { taskId });
                }
            }
            if (XlSearchParameters.XlOutPercolator)
            {
                var intraPsmsXLPercolator = intraPsmsXL.Where(p => p.XLBestScore >= 2 && p.BetaPsmCross.XLBestScore >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteCrosslinkToTxtForPercolator(intraPsmsXLPercolator, OutputFolder, "xl_intra_perc", crosslinker, new List<string> { taskId });
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

        //Calculate the FDR of crosslinked peptide FP/(FP+TP)
        private static List<PsmCross> CrosslinkDoFalseDiscoveryRateAnalysis(List<PsmCross> items)
        {
            var ids = new List<PsmCross>();
            foreach (var item in items)
            {
                ids.Add(item);
            }

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            for (int i = 0; i < ids.Count; i++)
            {
                var item1 = ids[i]; var item2 = ids[i].BetaPsmCross;

                var isDecoy1 = item1.MostProbableProteinInfo.IsDecoy; var isDecoy2 = item2.MostProbableProteinInfo.IsDecoy;
                if (isDecoy1 || isDecoy2)
                    cumulative_decoy++;
                else
                    cumulative_target++;

                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                item1.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
                //item2.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
            }

            double min_q_value = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                PsmCross id = ids[i];
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;
            }

            return ids;
        }

        //Calculate the FDR of crosslinked peptide (D - 2DD) / T
        private static List<PsmCross> CrosslinkFDRAnalysis(List<PsmCross> items)
        {
            var ids = new List<PsmCross>();
            foreach (var item in items)
            {
                ids.Add(item);
            }
            int cumulative_target = 0;
            int cumulative_decoy = 0;
            int cumulative_decoy_decoy = 0;

            for (int i = 0; i < ids.Count; i++)
            {
                var item1 = ids[i]; var item2 = ids[i].BetaPsmCross;

                var isDecoy1 = item1.MostProbableProteinInfo.IsDecoy; var isDecoy2 = item2.MostProbableProteinInfo.IsDecoy;
                if (isDecoy1 || isDecoy2)
                    cumulative_decoy++;
                else
                    cumulative_target++;

                if (isDecoy1 && isDecoy2)
                {
                    cumulative_decoy_decoy++;
                }

                double temp_q_value = (double)(cumulative_decoy - 2*cumulative_decoy_decoy) / (cumulative_target + cumulative_decoy);
                item1.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
// item2.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
            }

            double min_q_value = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                PsmCross id = ids[i];
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;
            }
            return ids;
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