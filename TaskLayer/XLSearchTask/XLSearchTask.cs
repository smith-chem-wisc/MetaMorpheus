using EngineLayer;
using EngineLayer.CrosslinkAnalysis;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public partial class XLSearchTask : MetaMorpheusTask
    {
        private const double BinToleranceInDaltons = 0.003;

        public XLSearchTask() : base(TaskType.XLSearch)
        {
            CommonParameters = new CommonParameters();
            XlSearchParameters = new XlSearchParameters();
        }

        public XlSearchParameters XlSearchParameters { get; set; }

        private object IndexLock { get; set; } = new object();
        private object PsmLock { get; set; } = new object();

        protected override MyTaskResults RunSpecific(string outputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            List<PsmCross> allPsms = new List<PsmCross>();
            var compactPeptideToProteinPeptideMatch = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();

            Status("Loading modifications...", taskId);
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();
            List<string> localizeableModificationTypes = GlobalVariables.AllModTypesKnown.ToList();

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, true, XlSearchParameters.DecoyType, localizeableModificationTypes);

            List<ProductType> ionTypes = new List<ProductType>();
            if (CommonParameters.BIons) { ionTypes.Add(ProductType.BnoB1ions); }
            if (CommonParameters.YIons) { ionTypes.Add(ProductType.Y); }
            if (CommonParameters.ZdotIons) { ionTypes.Add(ProductType.Zdot); }
            if (CommonParameters.CIons) { ionTypes.Add(ProductType.C); }

            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(ionTypes);

            var crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(XlSearchParameters.CrosslinkerType);
            if (XlSearchParameters.CrosslinkerType == CrosslinkerType.UserDefined)
            {
                crosslinker = GenerateUserDefinedCrosslinker(XlSearchParameters);
            }

            MyFileManager myFileManager = new MyFileManager(XlSearchParameters.DisposeOfFileWhenDone);

            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));
            HashSet<DigestionParams> ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams));

            int completedFiles = 0;

            Status("Searching files...", taskId);

            ProseCreatedWhileRunning.Append("The following crosslink discovery were used: ");
            ProseCreatedWhileRunning.Append("crosslinker name = " + crosslinker.CrosslinkerName + "; ");
            ProseCreatedWhileRunning.Append("crosslinker type = " + crosslinker.Cleavable + "; ");
            ProseCreatedWhileRunning.Append("crosslinker mass = " + crosslinker.TotalMass + "; ");
            ProseCreatedWhileRunning.Append("crosslinker modification site(s) = " + crosslinker.CrosslinkerModSites + "; ");

            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");

            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");

            ProseCreatedWhileRunning.Append("parent mass tolerance(s) = " + XlSearchParameters.XlPrecusorMsTl + "; ");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);
                List<PsmCross> newPsms = new List<PsmCross>();

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<CompactPeptide> peptideIndex = null;
                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count / combinedParams.TotalPartitions));

                    // Generate indices for modern search

                    Status("Getting fragment dictionary...", new List<string> { taskId });
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, ionTypes, currentPartition, XlSearchParameters.DecoyType, ListOfDigestionParams, combinedParams, 30000.0, new List<string> { taskId });
                    List<int>[] fragmentIndex = null;
                    lock (IndexLock)
                        GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, taskId);

                    Status("Searching files...", taskId);

                    new TwoPassCrosslinkSearchEngine(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, ionTypes, currentPartition, combinedParams, false, XlSearchParameters.XlPrecusorMsTl, crosslinker, XlSearchParameters.CrosslinkSearchTop, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.XlQuench_H2O, XlSearchParameters.XlQuench_NH2, XlSearchParameters.XlQuench_Tris, XlSearchParameters.XlCharge_2_3, XlSearchParameters.XlCharge_2_3_PrimeFragment, thisId).Run();
                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                }

                lock (PsmLock)
                {
                    allPsms.AddRange(newPsms);
                }

                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            Status("Crosslink analysis engine", taskId);
            MetaMorpheusEngineResults allcrosslinkanalysisResults;
            allcrosslinkanalysisResults = new CrosslinkAnalysisEngine(allPsms, compactPeptideToProteinPeptideMatch, proteinList, variableModifications, fixedModifications, ionTypes, outputFolder, crosslinker, terminusType, CommonParameters, new List<string> { taskId }).Run();
            allPsms = allPsms.Where(p => p != null).ToList();
            if (XlSearchParameters.XlOutAll)
            {
                try
                {
                    WriteAllToTsv(allPsms, outputFolder, "allPsms", new List<string> { taskId });
                }
                catch (Exception)
                {
                    throw;
                }
            }
            var allPsmsXL = allPsms.Where(p => p.CrossType == PsmCrossType.Cross).Where(p => p.XLBestScore >= CommonParameters.ScoreCutoff && p.BetaPsmCross.XLBestScore >= CommonParameters.ScoreCutoff).ToList();
            foreach (var item in allPsmsXL)
            {
                if (item.OneBasedStartResidueInProtein.HasValue)
                {
                    item.XlProteinPos = item.OneBasedStartResidueInProtein.Value + item.XlPos - 1;
                }
                if (item.BetaPsmCross.OneBasedStartResidueInProtein.HasValue)
                {
                    item.BetaPsmCross.XlProteinPos = item.BetaPsmCross.OneBasedStartResidueInProtein.Value + item.BetaPsmCross.XlPos - 1;
                }
            }

            // Inter Crosslink

            //Write Inter Psms FDR
            var interPsmsXL = allPsmsXL.Where(p => !p.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First().Contains(p.BetaPsmCross.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First()) &&
            !p.BetaPsmCross.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First().Contains(p.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First())).OrderByDescending(p => p.XLQvalueTotalScore).ToList();
            foreach (var item in interPsmsXL)
            {
                item.CrossType = PsmCrossType.Inter;
            }
            var interPsmsXLFDR = CrosslinkDoFalseDiscoveryRateAnalysis(interPsmsXL).ToList();
            //var interPsmsXLFDR = CrosslinkFDRAnalysis(interPsmsXL).ToList();
            if (XlSearchParameters.XlOutCrosslink)
            {
                WriteCrosslinkToTsv(interPsmsXLFDR, outputFolder, "xl_inter_fdr", new List<string> { taskId });
            }

            if (XlSearchParameters.XlOutPercolator)
            {
                try
                {
                    var interPsmsXLPercolator = interPsmsXL.Where(p => p.XLBestScore >= 2 && p.BetaPsmCross.XLBestScore >= 2).OrderBy(p => p.ScanNumber).ToList();
                    WriteCrosslinkToTxtForPercolator(interPsmsXLPercolator, outputFolder, "xl_inter_perc", crosslinker, new List<string> { taskId });
                }
                catch (Exception)
                {
                    throw;
                }
            }

            // Intra Cross-link

            //Write Intra Psms FDR
            var intraPsmsXL = allPsmsXL.Where(p => p.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First() == p.BetaPsmCross.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First() ||
            p.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First().Contains(p.BetaPsmCross.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First()) ||
            p.BetaPsmCross.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First().Contains(p.CompactPeptides.First().Value.Item2.Select(b => b.Protein.Accession).First())).OrderByDescending(p => p.XLQvalueTotalScore).ToList();
            foreach (var item in intraPsmsXL)
            {
                item.CrossType = PsmCrossType.Intra;
            }
            var intraPsmsXLFDR = CrosslinkDoFalseDiscoveryRateAnalysis(intraPsmsXL).ToList();
            //var intraPsmsXLFDR = CrosslinkFDRAnalysis(intraPsmsXL).ToList();
            if (XlSearchParameters.XlOutCrosslink)
            {
                WriteCrosslinkToTsv(intraPsmsXLFDR, outputFolder, "xl_intra_fdr", new List<string> { taskId });
            }

            if (XlSearchParameters.XlOutPercolator)
            {
                try
                {
                    var intraPsmsXLPercolator = intraPsmsXL.Where(p => p.XLBestScore >= 2 && p.BetaPsmCross.XLBestScore >= 2).OrderBy(p => p.ScanNumber).ToList();
                    WriteCrosslinkToTxtForPercolator(intraPsmsXLPercolator, outputFolder, "xl_intra_perc", crosslinker, new List<string> { taskId });
                }
                catch (Exception)
                {
                    throw;
                }
            }

            // Single peptide

            var singlePsms = allPsms.Where(p => p.CrossType == PsmCrossType.Singe && p.FullSequence != null && !p.FullSequence.Contains("Crosslink")).OrderByDescending(p => p.Score).ToList();
            var singlePsmsFDR = SingleFDRAnalysis(singlePsms).ToList();
            if (XlSearchParameters.XlOutAll)
            {
                WriteSingleToTsv(singlePsmsFDR, outputFolder, "single_fdr", new List<string> { taskId });
            }

            // Loop peptide

            var loopPsms = allPsms.Where(p => p.CrossType == PsmCrossType.Loop).OrderByDescending(p => p.XLTotalScore).ToList();
            var loopPsmsFDR = SingleFDRAnalysis(loopPsms).ToList();
            if (XlSearchParameters.XlOutAll)
            {
                WriteSingleToTsv(loopPsmsFDR, outputFolder, "loop_fdr", new List<string> { taskId });
            }

            // deadend peptide

            var deadendPsms = allPsms.Where(p => p.CrossType == PsmCrossType.DeadEnd || p.CrossType == PsmCrossType.DeadEndH2O || p.CrossType == PsmCrossType.DeadEndNH2 || p.CrossType == PsmCrossType.DeadEndTris).OrderByDescending(p => p.XLTotalScore).ToList();

            deadendPsms.AddRange(allPsms.Where(p => p.CrossType == PsmCrossType.Singe && p.FullSequence != null && p.FullSequence.Contains("Crosslink")).ToList());
            var deadendPsmsFDR = SingleFDRAnalysis(deadendPsms).ToList();
            if (XlSearchParameters.XlOutAll)
            {
                WriteSingleToTsv(deadendPsmsFDR, outputFolder, "deadend_fdr", new List<string> { taskId });
            }

            if (XlSearchParameters.XlOutPepXML)
            {
                List<PsmCross> allPsmsFDR = new List<PsmCross>();
                allPsmsFDR.AddRange(intraPsmsXLFDR.Where(p => p.IsDecoy != true && p.BetaPsmCross.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR.AddRange(interPsmsXLFDR.Where(p => p.IsDecoy != true && p.BetaPsmCross.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR.AddRange(singlePsmsFDR.Where(p => p.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR.AddRange(loopPsmsFDR.Where(p => p.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR.AddRange(deadendPsmsFDR.Where(p => p.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR = allPsmsFDR.OrderBy(p => p.ScanNumber).ToList();
                foreach (var fullFilePath in currentRawFileList)
                {
                    string fileNameNoExtension = Path.GetFileNameWithoutExtension(fullFilePath);
                    WritePepXML_xl(allPsmsFDR.Where(p => p.FullFilePath == fullFilePath).ToList(), proteinList, dbFilenameList[0].FilePath, variableModifications, fixedModifications, localizeableModificationTypes, outputFolder, fileNameNoExtension, new List<string> { taskId });
                }
            }

            if (XlSearchParameters.XlOutAll)
            {
                List<PsmCross> allPsmsXLFDR = new List<PsmCross>();
                allPsmsXLFDR.AddRange(intraPsmsXLFDR.Where(p => p.IsDecoy != true && p.BetaPsmCross.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsXLFDR.AddRange(interPsmsXLFDR.Where(p => p.IsDecoy != true && p.BetaPsmCross.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                try
                {
                    allPsmsXLFDR = allPsmsXLFDR.OrderByDescending(p => p.XLQvalueTotalScore).ToList();
                    var allPsmsXLFDRGroup = FindCrosslinks(allPsmsXLFDR);
                    WriteCrosslinkToTsv(allPsmsXLFDRGroup, outputFolder, "allPsmsXLFDRGroup", new List<string> { taskId });
                }
                catch (Exception)
                {
                    throw;
                }
            }

            return MyTaskResults;
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

        //Calculate the FDR of single peptide FP/TP
        private static List<PsmCross> SingleFDRAnalysis(List<PsmCross> items)
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
                var item1 = ids[i];

                var isDecoy1 = item1.IsDecoy;
                if (isDecoy1)
                    cumulative_decoy++;
                else
                    cumulative_target++;

                double temp_q_value = (double)cumulative_decoy / cumulative_target;
                item1.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0, 0, 0, 0, false);
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

        //Calculate the FDR of crosslinked peptide FP/TP
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

                var isDecoy1 = item1.IsDecoy; var isDecoy2 = item2.IsDecoy;
                if (isDecoy1 || isDecoy2)
                    cumulative_decoy++;
                else
                    cumulative_target++;

                double temp_q_value = (double)cumulative_decoy / cumulative_target;
                item1.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, 0, 0, 0, 0, 0, 0, false);
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

        //Find crosslinks
        private static List<PsmCross> FindCrosslinks(List<PsmCross> items)
        {
            List<PsmCross> psmCrossCrosslinks = new List<PsmCross>();
            List<string> allString = new List<string>();
            foreach (var item in items)
            {
                if (item.ProteinAccesion != null && item.BetaPsmCross.ProteinAccesion != null)
                {
                    string st;
                    if (item.ProteinAccesion.CompareTo(item.BetaPsmCross.ProteinAccesion) > 0)
                    {
                        st = item.BetaPsmCross.ProteinAccesion + "(" + item.BetaPsmCross.XlProteinPos + ")" + item.ProteinAccesion + "(" + item.XlProteinPos + ")";
                    }
                    else if (item.ProteinAccesion.CompareTo(item.BetaPsmCross.ProteinAccesion) == 0 && item.XlProteinPos.CompareTo(item.BetaPsmCross.XlProteinPos) > 0)
                    {
                        st = item.BetaPsmCross.ProteinAccesion + "(" + item.BetaPsmCross.XlProteinPos + ")" + item.ProteinAccesion + "(" + item.XlProteinPos + ")";
                    }
                    else { st = item.ProteinAccesion + "(" + item.XlProteinPos + ")" + item.BetaPsmCross.ProteinAccesion + "(" + item.BetaPsmCross.XlProteinPos + ")"; }

                    if (!allString.Contains(st))
                    {
                        allString.Add(st);
                        psmCrossCrosslinks.Add(item);
                    }
                }
            }
            return psmCrossCrosslinks;
        }

        //Generate user defined crosslinker
        public static CrosslinkerTypeClass GenerateUserDefinedCrosslinker(XlSearchParameters xlSearchParameters)
        {
            var crosslinker = new CrosslinkerTypeClass(
            xlSearchParameters.UdXLkerResidues,
        xlSearchParameters.UdXLkerResidues2,
        xlSearchParameters.UdXLkerName,
        xlSearchParameters.UdXLkerCleavable,
        (xlSearchParameters.UdXLkerTotalMass.HasValue ? (double)xlSearchParameters.UdXLkerTotalMass : 9999),
        (xlSearchParameters.UdXLkerShortMass.HasValue ? (double)xlSearchParameters.UdXLkerShortMass : 9999),
        (xlSearchParameters.UdXLkerLongMass.HasValue ? (double)xlSearchParameters.UdXLkerLongMass : 9999),
        (xlSearchParameters.UdXLkerLoopMass.HasValue ? (double)xlSearchParameters.UdXLkerLoopMass : 9999),
        (xlSearchParameters.UdXLkerDeadendMassH2O.HasValue ? (double)xlSearchParameters.UdXLkerDeadendMassH2O : 9999),
        (xlSearchParameters.UdXLkerDeadendMassNH2.HasValue ? (double)xlSearchParameters.UdXLkerDeadendMassNH2 : 9999),
        (xlSearchParameters.UdXLkerDeadendMassTris.HasValue ? (double)xlSearchParameters.UdXLkerDeadendMassTris : 9999)
            );
            return crosslinker;
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

        private static void WriteIndexEngineParams(IndexingEngine indexEngine, string fileName)
        {
            using (StreamWriter output = new StreamWriter(fileName))
            {
                output.Write(indexEngine);
            }
        }

        private void WriteFragmentIndexNetSerializer(List<int>[] fragmentIndex, string fragmentIndexFile, string taskId)
        {
            var messageTypes = GetSubclassesAndItself(typeof(Dictionary<float, List<int>>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFile))
                ser.Serialize(file, fragmentIndex);
            SucessfullyFinishedWritingFile(fragmentIndexFile, new List<string> { taskId });
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
                    peptideIndex = (List<CompactPeptide>)ser.Deserialize(file);

                Status("Reading fragment index...", new List<string> { taskId });
                messageTypes = GetSubclassesAndItself(typeof(List<int>[]));
                ser = new NetSerializer.Serializer(messageTypes);
                using (var file = File.OpenRead(Path.Combine(pathToFolderWithIndices, "fragmentIndex.ind")))
                    fragmentIndex = (List<int>[])ser.Deserialize(file);
            }
        }
    }
}