using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;

namespace TaskLayer
{
    public partial class XLSearchTask : MetaMorpheusTask
    {
        public XLSearchTask() : base(MyTask.XLSearch)
        {
            CommonParameters = new CommonParameters(
                precursorMassTolerance: new PpmTolerance(10), 
                scoreCutoff: 3
            );

            XlSearchParameters = new XlSearchParameters();
        }

        public XlSearchParameters XlSearchParameters { get; set; }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            List<CrosslinkSpectralMatch> allPsms = new List<CrosslinkSpectralMatch>();

            LoadModifications(taskId, out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, true, XlSearchParameters.DecoyType, localizeableModificationTypes, CommonParameters);
            
            var crosslinker = new Crosslinker();
            crosslinker = crosslinker.SelectCrosslinker(XlSearchParameters.CrosslinkerType);
            if (XlSearchParameters.CrosslinkerType == CrosslinkerType.UserDefined)
            {
                crosslinker = GenerateUserDefinedCrosslinker(XlSearchParameters);
            }

            MyFileManager myFileManager = new MyFileManager(XlSearchParameters.DisposeOfFileWhenDone);

            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));
            HashSet<DigestionParams> ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

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

            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif) + "; "));
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");

            ProseCreatedWhileRunning.Append("parent mass tolerance(s) = " + CommonParameters.PrecursorMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);

                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks, combinedParams);

                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

                CrosslinkSpectralMatch[] newPsms = new CrosslinkSpectralMatch[arrayOfMs2ScansSortedByMass.Length];
                for (int currentPartition = 0; currentPartition < CommonParameters.TotalPartitions; currentPartition++)
                {
                    List<PeptideWithSetModifications> peptideIndex = null;
                    List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                    Status("Getting fragment dictionary...", new List<string> { taskId });
                    var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, currentPartition, UsefulProteomicsDatabases.DecoyType.Reverse, ListOfDigestionParams, combinedParams, 30000.0, new List<string> { taskId });
                    List<int>[] fragmentIndex = null;

                    GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, proteinList, GlobalVariables.AllModsKnown.ToList(), taskId);

                    Status("Searching files...", taskId);
                    new TwoPassCrosslinkSearchEngine(newPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition, combinedParams, crosslinker, XlSearchParameters.RestrictToTopNHits, XlSearchParameters.CrosslinkSearchTopNum, XlSearchParameters.XlQuench_H2O, XlSearchParameters.XlQuench_NH2, XlSearchParameters.XlQuench_Tris, XlSearchParameters.XlCharge_2_3, false, thisId).Run();
                    ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + CommonParameters.TotalPartitions + "!", thisId));
                }

                allPsms.AddRange(newPsms.Where(p => p != null));

                completedFiles++;
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            allPsms = allPsms.Where(p => p != null).ToList();
            if (XlSearchParameters.XlOutAll)
            {
                WriteAllToTsv(allPsms, OutputFolder, "allPsms", new List<string> { taskId });
            }

            var allPsmsXL = allPsms.Where(p => p.CrossType == PsmCrossType.Cross).ToList();
            foreach (var csm in allPsmsXL)
            {
                // alpha peptide crosslink residue in the protein
                csm.XlProteinPos = csm.OneBasedStartResidueInProtein.Value + csm.ModPositions.First() - 1;

                // beta crosslink residue in protein
                csm.BetaPeptide.XlProteinPos = csm.BetaPeptide.BestMatchingPeptideWithSetMods.First().Pwsm.OneBasedStartResidueInProtein + csm.BetaPeptide.ModPositions.First() - 1;
            }

            //Write Inter Psms FDR
            var interPsmsXL = allPsmsXL.Where(p => !p.BestMatchingPeptideWithSetMods.First().Pwsm.Protein.Accession.Equals(p.BetaPeptide.BestMatchingPeptideWithSetMods.First().Pwsm.Protein.Accession)).OrderByDescending(p => p.XLQvalueTotalScore).ToList();
            foreach (var item in interPsmsXL)
            {
                item.CrossType = PsmCrossType.Inter;
            }
            var interPsmsXLFDR = CrosslinkDoFalseDiscoveryRateAnalysis(interPsmsXL).ToList();

            if (XlSearchParameters.XlOutCrosslink)
            {
                string file = Path.Combine(OutputFolder, "xl_inter_fdr.tsv");
                WritePsmCrossToTsv(interPsmsXLFDR, file, 2);
                FinishedWritingFile(file, new List<string> { taskId });
            }

            if (XlSearchParameters.XlOutPercolator)
            {
                var interPsmsXLPercolator = interPsmsXL.Where(p => p.BestScore >= 2 && p.BetaPeptide.BestScore >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteCrosslinkToTxtForPercolator(interPsmsXLPercolator, OutputFolder, "xl_inter_perc", crosslinker, new List<string> { taskId });
            }

            //Write Intra Psms FDR
            var intraPsmsXL = allPsmsXL.Where(p => p.BestMatchingPeptideWithSetMods.First().Pwsm.Protein.Accession.Equals(p.BetaPeptide.BestMatchingPeptideWithSetMods.First().Pwsm.Protein.Accession)).OrderByDescending(p => p.XLQvalueTotalScore).ToList();
            foreach (var item in intraPsmsXL)
            {
                item.CrossType = PsmCrossType.Intra;
            }
            var intraPsmsXLFDR = CrosslinkDoFalseDiscoveryRateAnalysis(intraPsmsXL).ToList();

            if (XlSearchParameters.XlOutCrosslink)
            {
                string file = Path.Combine(OutputFolder, "xl_intra_fdr.tsv");
                WritePsmCrossToTsv(intraPsmsXLFDR, file, 2);
                FinishedWritingFile(file, new List<string> { taskId });
            }

            if (XlSearchParameters.XlOutPercolator)
            {
                var intraPsmsXLPercolator = intraPsmsXL.Where(p => p.BestScore >= 2 && p.BetaPeptide.BestScore >= 2).OrderBy(p => p.ScanNumber).ToList();
                WriteCrosslinkToTxtForPercolator(intraPsmsXLPercolator, OutputFolder, "xl_intra_perc", crosslinker, new List<string> { taskId });
            }

            // single peptide FDR
            var writtenFileSingle = Path.Combine(OutputFolder, "single_fdr" + ".tsv");
            var singlePsms = allPsms.Where(p => p.CrossType == PsmCrossType.Single && p.XLTotalScore >= 2 && (string.IsNullOrEmpty(p.FullSequence) ? true : !p.FullSequence.Contains("Crosslink"))).OrderByDescending(p => p.Score).ToList();
            var singlePsmsFDR = SingleFDRAnalysis(singlePsms).ToList();
            WritePsmCrossToTsv(singlePsmsFDR, writtenFileSingle, 1);
            FinishedWritingFile(writtenFileSingle, new List<string> { taskId });

            // loop FDR
            var writtenFileLoop = Path.Combine(OutputFolder, "loop_fdr" + ".tsv");
            var loopPsms = allPsms.Where(p => p.CrossType == PsmCrossType.Loop && p.XLTotalScore >= 2).OrderByDescending(p => p.XLTotalScore).ToList();
            var loopPsmsFDR = SingleFDRAnalysis(loopPsms).ToList();
            WritePsmCrossToTsv(loopPsmsFDR, writtenFileLoop, 1);
            FinishedWritingFile(writtenFileLoop, new List<string> { taskId });

            // deadend FDRd
            var writtenFileDeadend = Path.Combine(OutputFolder, "deadend_fdr" + ".tsv");
            var deadendPsms = allPsms.Where(p => p.BestScore > 2 && (p.CrossType == PsmCrossType.DeadEnd || p.CrossType == PsmCrossType.DeadEndH2O || p.CrossType == PsmCrossType.DeadEndNH2 || p.CrossType == PsmCrossType.DeadEndTris)).OrderByDescending(p => p.XLTotalScore).ToList();
            //If parameter.modification contains crosslinker.deadend as variable mod, then the deadend will be in the following form. 
            deadendPsms.AddRange(allPsms.Where(p => p.CrossType == PsmCrossType.Single && p.XLTotalScore >= 2 && (string.IsNullOrEmpty(p.FullSequence) ? true : p.FullSequence.Contains("Crosslink"))).ToList());
            var deadendPsmsFDR = SingleFDRAnalysis(deadendPsms).ToList();
            WritePsmCrossToTsv(deadendPsmsFDR, writtenFileDeadend, 1);
            FinishedWritingFile(writtenFileDeadend, new List<string> { taskId });

            if (XlSearchParameters.XlOutPepXML)
            {
                List<CrosslinkSpectralMatch> allPsmsFDR = new List<CrosslinkSpectralMatch>();
                allPsmsFDR.AddRange(intraPsmsXLFDR.Where(p => p.IsDecoy != true && p.BetaPeptide.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR.AddRange(interPsmsXLFDR.Where(p => p.IsDecoy != true && p.BetaPeptide.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR.AddRange(singlePsmsFDR.Where(p => p.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR.AddRange(loopPsmsFDR.Where(p => p.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR.AddRange(deadendPsmsFDR.Where(p => p.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsFDR = allPsmsFDR.OrderBy(p => p.ScanNumber).ToList();
                allPsms.ForEach(p => p.ResolveAllAmbiguities());
                foreach (var fullFilePath in currentRawFileList)
                {
                    string fileNameNoExtension = Path.GetFileNameWithoutExtension(fullFilePath);
                    WritePepXML_xl(allPsmsFDR.Where(p => p.FullFilePath == fullFilePath).ToList(), proteinList, dbFilenameList[0].FilePath, variableModifications, fixedModifications, localizeableModificationTypes, OutputFolder, fileNameNoExtension, new List<string> { taskId });
                }
            }

            if (XlSearchParameters.XlOutAll)
            {
                List<CrosslinkSpectralMatch> allPsmsXLFDR = new List<CrosslinkSpectralMatch>();
                allPsmsXLFDR.AddRange(intraPsmsXLFDR.Where(p => p.IsDecoy != true && p.BetaPeptide.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());
                allPsmsXLFDR.AddRange(interPsmsXLFDR.Where(p => p.IsDecoy != true && p.BetaPeptide.IsDecoy != true && p.FdrInfo.QValue <= 0.05).ToList());

                allPsmsXLFDR = allPsmsXLFDR.OrderByDescending(p => p.XLQvalueTotalScore).ToList();
                var allPsmsXLFDRGroup = FindCrosslinks(allPsmsXLFDR);
                WritePsmCrossToTsv(allPsmsXLFDRGroup, Path.Combine(OutputFolder, "allPsmsXLFDRGroup.tsv"), 2);
            }
            return MyTaskResults;
        }
        
        //Calculate the FDR of single peptide FP/TP
        private static List<CrosslinkSpectralMatch> SingleFDRAnalysis(List<CrosslinkSpectralMatch> items)
        {
            var ids = new List<CrosslinkSpectralMatch>();
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
                CrosslinkSpectralMatch id = ids[i];
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;
            }

            return ids;
        }

        //Calculate the FDR of crosslinked peptide FP/TP
        private static List<CrosslinkSpectralMatch> CrosslinkDoFalseDiscoveryRateAnalysis(List<CrosslinkSpectralMatch> items)
        {
            var ids = new List<CrosslinkSpectralMatch>();
            foreach (var item in items)
            {
                ids.Add(item);
            }

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            for (int i = 0; i < ids.Count; i++)
            {
                var item1 = ids[i];
                var item2 = ids[i].BetaPeptide;

                var isDecoy1 = item1.IsDecoy;
                var isDecoy2 = item2.IsDecoy;

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
                CrosslinkSpectralMatch id = ids[i];
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;
            }

            return ids;
        }

        //Find crosslinks
        private static List<CrosslinkSpectralMatch> FindCrosslinks(List<CrosslinkSpectralMatch> csms)
        {
            List<CrosslinkSpectralMatch> psmCrossCrosslinks = new List<CrosslinkSpectralMatch>();
            List<string> allString = new List<string>();
            foreach (var csm in csms)
            {
                if (csm.ProteinAccession != null && csm.BetaPeptide.ProteinAccession != null)
                {
                    string st;
                    if (csm.ProteinAccession.CompareTo(csm.BetaPeptide.ProteinAccession) > 0)
                    {
                        st = csm.BetaPeptide.ProteinAccession + "(" + csm.BetaPeptide.XlProteinPos + ")" + csm.ProteinAccession + "(" + csm.XlProteinPos + ")";
                    }
                    else if (csm.ProteinAccession.CompareTo(csm.BetaPeptide.ProteinAccession) == 0 && csm.XlProteinPos.CompareTo(csm.BetaPeptide.XlProteinPos) > 0)
                    {
                        st = csm.BetaPeptide.ProteinAccession + "(" + csm.BetaPeptide.XlProteinPos + ")" + csm.ProteinAccession + "(" + csm.XlProteinPos + ")";
                    }
                    else { st = csm.ProteinAccession + "(" + csm.XlProteinPos + ")" + csm.BetaPeptide.ProteinAccession + "(" + csm.BetaPeptide.XlProteinPos + ")"; }

                    if (!allString.Contains(st))
                    {
                        allString.Add(st);
                        psmCrossCrosslinks.Add(csm);
                    }
                }
            }
            return psmCrossCrosslinks;
        }

        //Generate user defined crosslinker
        public static Crosslinker GenerateUserDefinedCrosslinker(XlSearchParameters xlSearchParameters)
        {
            var crosslinker = new Crosslinker(
                xlSearchParameters.CrosslinkerResidues,
                xlSearchParameters.CrosslinkerResidues2,
                xlSearchParameters.CrosslinkerName,
                xlSearchParameters.IsCleavable,
                xlSearchParameters.CrosslinkerTotalMass ?? double.NaN,
                xlSearchParameters.CrosslinkerShortMass ?? double.NaN,
                xlSearchParameters.CrosslinkerLongMass ?? double.NaN,
                xlSearchParameters.CrosslinkerLoopMass ?? double.NaN,
                xlSearchParameters.CrosslinkerDeadEndMassH2O ?? double.NaN,
                xlSearchParameters.CrosslinkerDeadEndMassNH2 ?? double.NaN,
                xlSearchParameters.CrosslinkerDeadEndMassTris ?? double.NaN
            );

            return crosslinker;
        }
    }
}