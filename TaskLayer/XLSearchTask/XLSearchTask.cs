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

        #endregion Private Fields

        #region Public Constructors

        public XLSearchTask() : base(MyTask.XLSearch)
        {
            CommonParameters = new CommonParameters
            {
                MaxMissedCleavages = 2,
                MinPeptideLength = 5,
                MaxPeptideLength = null,
                MaxModificationIsoforms = 4096,
                Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"],
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable,
                ProductMassTolerance = new AbsoluteTolerance(0.01),
                BIons = true,
                YIons = true,
                ZdotIons = false,
                CIons = false,

                TotalPartitions = 1,
                LocalizeAll = true,

                ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") },
                ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") },
                ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList(),

                Max_mods_for_peptide = 3,

                ConserveMemory = true,
                MaxDegreeOfParallelism = 1,
                ScoreCutoff = 5,

                // Deconvolution stuff
                DoPrecursorDeconvolution = true,
                UseProvidedPrecursorInfo = true,
                DeconvolutionIntensityRatio = 4,
                DeconvolutionMaxAssumedChargeState = 10,
                DeconvolutionMassTolerance = new PpmTolerance(5),

            };
            XlSearchParameters = new XlSearchParameters
            {
                SearchDecoy = true,
                CrosslinkerType = CrosslinkerType.DSS,
                CrosslinkSearchTopNum = 50,
                CrosslinkSearchWithAllBeta = false,
                UdXLkerName = null,
                UdXLkerCleavable = false,
                UdXLkerShortMass = null,
                UdXLkerLongMass = null,
                UdXLkerTotalMass = null,
                UdXLkerResidue = 'K',
                XlPrecusorMsTl = new PpmTolerance(10),
                XlBetaPrecusorMsTl = new PpmTolerance(10),
            };

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
                + CommonParameters.InitiatorMethionineBehavior
                + " and the maximum number of allowed missed cleavages is "
                + CommonParameters.MaxMissedCleavages);
            sb.AppendLine("MinPeptideLength: " + CommonParameters.MinPeptideLength);
            sb.AppendLine("MaxPeptideLength: " + CommonParameters.MaxPeptideLength);
            sb.AppendLine("maxModificationIsoforms: " + CommonParameters.MaxModificationIsoforms);
            sb.AppendLine("protease: " + CommonParameters.Protease);
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

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId)
        {
            myTaskResults = new MyTaskResults(this);
            List<PsmCross> allPsms = new List<PsmCross>();
            var compactPeptideToProteinPeptideMatch = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();


            Status("Loading modifications...", taskId);

            #region Load modifications

            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> localizeableModifications;
            if (CommonParameters.LocalizeAll)
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
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

            proseCreatedWhileRunning.Append("protease = " + CommonParameters.Protease + "; ");
            proseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.MaxMissedCleavages + "; ");
            proseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.MinPeptideLength + "; ");
            if (CommonParameters.MaxPeptideLength == null)
            {
                proseCreatedWhileRunning.Append("maximum peptide length = unspecified; ");
            }
            else
            {
                proseCreatedWhileRunning.Append("maximum peptide length = " + CommonParameters.MaxPeptideLength + "; ");
            }
            proseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.InitiatorMethionineBehavior + "; ");
            proseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");

            proseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.MaxModificationIsoforms + "; ");
            proseCreatedWhileRunning.Append("parent mass tolerance(s) = UNKNOWN; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                //Psm[][] fileSpecificPsms = new Psm[MassDiffAcceptors.Count()][];
                List<PsmCross> newPsms = new List<PsmCross>();
                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile);
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
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, ionTypes, currentPartition, XlSearchParameters.SearchDecoy, CommonParameters, new List<string> { taskId });

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
                        new CrosslinkSearchEngine2(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, keys, fragmentIndex, CommonParameters.ProductMassTolerance, crosslinker, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.CrosslinkSearchWithAllBeta, XlSearchParameters.XlPrecusorMsTl, XlSearchParameters.XlBetaPrecusorMsTl, modsDictionary, ionTypes, proteinList, CommonParameters.Protease, CommonParameters.MaxMissedCleavages, CommonParameters.MinPeptideLength, CommonParameters.MaxPeptideLength, variableModifications, fixedModifications, CommonParameters.MaxModificationIsoforms, thisId).Run();
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
            allcrosslinkanalysisResults = new CrosslinkAnalysisEngine(allPsms, compactPeptideToProteinPeptideMatch, proteinList, variableModifications, fixedModifications, CommonParameters.Protease, null, CommonParameters.ProductMassTolerance, CommonParameters.MaxMissedCleavages, CommonParameters.MinPeptideLength, CommonParameters.MaxPeptideLength, CommonParameters.MaxModificationIsoforms, ionTypes, CommonParameters.InitiatorMethionineBehavior, modsDictionary, null, OutputFolder, crosslinker, new List<string> { taskId }, terminusType).Run();

            var allPsmsXL = allPsms.Where(p=>p.CrossType == PsmCrossType.Cross).OrderByDescending(p => p.ScanNumber).ToList();
            //Write Inter Psms FDR
            var interPsmsXLFDR = allPsmsXL.Where(p => p.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() != p.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).OrderByDescending(p=>p.XLTotalScore).ToList();
            interPsmsXLFDR = CrosslinkDoFalseDiscoveryRateAnalysis(interPsmsXLFDR).Where(p => p.MostProbableProteinInfo.IsDecoy != true && p.BetaPsmCross.MostProbableProteinInfo.IsDecoy != true && p.FdrInfo.QValue <= 0.01).ToList();
            WriteCrosslinkToTsv(interPsmsXLFDR, OutputFolder, "xl_inter_fdr", new List<string> { taskId });
            WriteCrosslinkToTxtForCLMSVault(interPsmsXLFDR, OutputFolder, "xl_inter_fdr_CLMSVault", crosslinker, new List<string> { taskId });
            //
            if (interPsmsXLFDR.Count != 0)
            {
                foreach (var fullFilePath in currentRawFileList)
                {
                    string fileNameNoExtension = Path.GetFileNameWithoutExtension(fullFilePath);
                    WritePepXML(interPsmsXLFDR.Where(p => p.FullFilePath == fullFilePath).ToList(), dbFilenameList, variableModifications, fixedModifications, localizeableModifications, OutputFolder, fileNameNoExtension, new List<string> { taskId });
                }
            }
            //
            var interPsmsXLPercolator = allPsmsXL.Where(p => p.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() != p.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).OrderBy(p => p.ScanNumber).ToList();
            WriteCrosslinkToTxtForPercolator(interPsmsXLPercolator, OutputFolder, "xl_inter_perc", crosslinker, new List<string> { taskId });

            //Write Intra Psms FDR
            var intraPsmsXLFDR = allPsmsXL.Where(p => p.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() == p.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).OrderByDescending(p => p.XLTotalScore).ToList();
            intraPsmsXLFDR = CrosslinkDoFalseDiscoveryRateAnalysis(intraPsmsXLFDR).Where(p => p.MostProbableProteinInfo.IsDecoy != true && p.BetaPsmCross.MostProbableProteinInfo.IsDecoy != true).ToList();
            WriteCrosslinkToTsv(intraPsmsXLFDR, OutputFolder, "xl_intra_fdr", new List<string> { taskId });
            WriteCrosslinkToTxtForCLMSVault(intraPsmsXLFDR, OutputFolder, "xl_intra_fdr_CLMSVault", crosslinker, new List<string> { taskId });
            //
            if (intraPsmsXLFDR.Count != 0)
            {
                foreach (var fullFilePath in currentRawFileList)
                {
                    string fileNameNoExtension = Path.GetFileNameWithoutExtension(fullFilePath);
                    WritePepXML(intraPsmsXLFDR.Where(p => p.FullFilePath == fullFilePath).ToList(), dbFilenameList, variableModifications, fixedModifications, localizeableModifications, OutputFolder, fileNameNoExtension, new List<string> { taskId });
                }
            }
            //
            var intraPsmsXLPercolator = allPsmsXL.Where(p => p.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First() == p.BetaPsmCross.MostProbableProteinInfo.PeptidesWithSetModifications.Select(b => b.Protein.Accession).First()).OrderBy(p => p.ScanNumber).ToList();
            WriteCrosslinkToTxtForPercolator(intraPsmsXLPercolator, OutputFolder, "xl_intra_perc", crosslinker, new List<string> { taskId });

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
                item2.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0);
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