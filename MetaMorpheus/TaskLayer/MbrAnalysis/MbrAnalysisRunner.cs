using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using FlashLFQ;
using MassSpectrometry;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Collections.Concurrent;
using EngineLayer.PsmTsv;
using Microsoft.ML;
using Microsoft.ML.Data;

namespace TaskLayer.MbrAnalysis
{
    
    public static class MbrAnalysisRunner
    {
        public const string mbrAnalysisFolder = "MbrAnalysis";

        /// <summary>
        /// Performs secondary analysis of MBR results by searching acceptor files for candidate spectra,
        /// and comparing those spectra to a spectral library. Results (PEP model and .psmtsv) are written to
        /// unique MbrAnalysis folder.
        /// </summary>
        /// <param name="parameters"></param>
        /// <param name="commonParameters"></param>
        /// <param name="fileSpecificParameters"></param>
        public static MbrAnalysisResults RunMbrAnalysis(
            PostSearchAnalysisParameters parameters,
            CommonParameters commonParameters,
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            Dictionary<Identification, PeptideWithSetModifications> idsToPwsms = null)
        {

            if (!parameters.SearchParameters.DoMbrAnalysis | !parameters.SearchParameters.MatchBetweenRuns)
            {
                return null;
            }

            List<SpectraFileInfo> spectraFiles = parameters.FlashLfqResults.Peaks.Select(p => p.Key).ToList();
            List<PeptideSpectralMatch> allPeptides = GetAllPeptides(parameters, commonParameters, fileSpecificParameters);

            int maxThreadsPerFile = commonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();

            ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> bestMbrMatches = new();

            //Dictionary<SpectraFileInfo, Dictionary<int, double>> ms1TicDictionary = new();
            foreach (SpectraFileInfo spectraFile in spectraFiles)
            {
                List<ChromatographicPeak> fileSpecificMbrPeaks =
                    parameters.FlashLfqResults.Peaks[spectraFile].Where(p => p.IsMbrPeak).ToList();
                if (fileSpecificMbrPeaks == null || (!fileSpecificMbrPeaks.Any())) break;

                MyFileManager myFileManager = new(true);
                MsDataFile myMsDataFile =
                    myFileManager.LoadFile(spectraFile.FullFilePathWithExtension, commonParameters);
                //ms1TicDictionary.TryAdd(spectraFile, RetrieveMs1TicInfo(myMsDataFile));

                MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(
                    commonParameters.PrecursorMassTolerance,
                    parameters.SearchParameters.MassDiffAcceptorType, parameters.SearchParameters.CustomMdac);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByRT = MetaMorpheusTask
                    .GetMs2Scans(myMsDataFile, spectraFile.FullFilePathWithExtension, commonParameters)
                    .OrderBy(b => b.RetentionTime).ToArray();

                var list = Directory.GetFiles(parameters.OutputFolder, "*.*", SearchOption.AllDirectories);
                string matchingvalue = list.Where(p => p.Contains("spectralLibrary")).First().ToString();
                string spectralLibraryPath = Path.Combine(parameters.OutputFolder, matchingvalue);
                //var updatedLib = new SpectralLibrary(new List<string> { Path.Combine(thisTaskOutputFolder, matchingvalue) });
                //string spectralLibraryPath = Path.Combine(parameters.OutputFolder, @"spectralLibrary.msp");
                SpectralLibrary library = new(new List<string>() { spectralLibraryPath });

                MiniClassicSearchEngine mcse = new(
                    arrayOfMs2ScansSortedByRT,
                    massDiffAcceptor,
                    commonParameters,
                    library,
                    fileSpecificParameters.
                        Where(t => t.FileName.Equals(spectraFile.FilenameWithoutExtension)).
                        Select(t => t.Parameters).FirstOrDefault());

                Parallel.ForEach(threads, (i) =>
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    for (; i < fileSpecificMbrPeaks.Count; i += maxThreadsPerFile)
                    {
                        ChromatographicPeak mbrPeak = fileSpecificMbrPeaks[i];
                        PeptideSpectralMatch bestDonorPsm = allPeptides.Where(p =>
                            p.FullSequence == mbrPeak.Identifications.First().ModifiedSequence).FirstOrDefault();
                        if (bestDonorPsm == null)
                        {
                            break;
                        }
                        PeptideWithSetModifications bestDonorPwsm = bestDonorPsm.BestMatchingPeptides.First().Peptide;

                        IEnumerable<PeptideSpectralMatch> peptideSpectralMatches =
                            mcse.SearchAroundPeak(bestDonorPwsm, mbrPeak.Apex.IndexedPeak.RetentionTime);

                        if (peptideSpectralMatches == null || !peptideSpectralMatches.Any())
                        {
                            bestMbrMatches.TryAdd(mbrPeak, new MbrSpectralMatch(null, mbrPeak));
                        }
                        else
                        {
                            bestMbrMatches.TryAdd(mbrPeak,
                                new MbrSpectralMatch(BestPsmForMbrPeak(peptideSpectralMatches), mbrPeak));
                        }
                    }
                });

                library.CloseConnections();
            }

            if (bestMbrMatches.Any())
            {
                List<PeptideSpectralMatch> allPsms = parameters.AllPsms.
                    OrderByDescending(p => p.Score).
                    ThenBy(p => p.FdrInfo.QValue).
                    ThenBy(p => p.FullFilePath).
                    ThenBy(x => x.ScanNumber).
                    ThenBy(p => p.FullSequence).
                    ThenBy(p => p.ProteinAccession).ToList();

                AssignEstimatedPsmQvalue(bestMbrMatches, allPsms);
                FDRAnalysisOfMbrPsms(bestMbrMatches, allPsms, parameters, fileSpecificParameters);
                AssignEstimatedPsmPepQValue(bestMbrMatches, allPsms);
                foreach (MbrSpectralMatch match in bestMbrMatches.Values) match.FindOriginalPsm(allPsms);
            }

            Directory.CreateDirectory(Path.Join(parameters.OutputFolder, mbrAnalysisFolder));
            WriteMbrPsmResults(bestMbrMatches, parameters);

            return new MbrAnalysisResults(bestMbrMatches, parameters.FlashLfqResults, idsToPwsms);
        }

        //private static Dictionary<int, double> RetrieveMs1TicInfo(MsDataFile msDataFile)
        //{
        //    Dictionary<int, double> zeroBasedScanToTicDictionary = new();
        //    //Dictionary<double, double> retentionTimeToTicDictionary = new();
        //    IEnumerable<MsDataScan> ms1Scans = msDataFile.GetAllScansList().Where(x => x.MsnOrder == 1);
        //    int zeroBasedScanIndex = 0;
        //    foreach (MsDataScan scan in ms1Scans)
        //    {
        //        //retentionTimeToTicDictionary.TryAdd(scan.RetentionTime, scan.TotalIonCurrent);
        //        zeroBasedScanToTicDictionary.TryAdd(zeroBasedScanIndex, scan.TotalIonCurrent);
        //        zeroBasedScanIndex++;
        //    }
        //    return zeroBasedScanToTicDictionary;
        //}

        /// <summary>
        /// Performs secondary analysis of MBR results by searching acceptor files for candidate spectra,
        /// and comparing those spectra to a spectral library. Results (PEP model and .psmtsv) are written to
        /// unique MbrAnalysis folder.
        /// </summary>
        /// <param name="parameters"></param>
        /// <param name="commonParameters"></param>
        /// <param name="fileSpecificParameters"></param>
        /// <param name ="calibratedFilePaths"> Optional param to link spectraFileInfo objects to the file path of a calibrated version of the same file</param>
        public static MbrAnalysisResults RunMbrAnalysisFromMaxQuant(
            List<SpectraFileInfo> spectraFiles,
            Dictionary<string, List<ChromatographicPeak>> mbrPeaksDict,
            Dictionary<string, PsmFromTsv> donorPsmDict,
            string spectralLibraryPath,
            string outputFolder,
            CommonParameters commonParameters,
            Dictionary<SpectraFileInfo, string> calibratedFilePaths = null)
        {

            int maxThreadsPerFile = commonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();

            ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> bestMbrMatches = new();

            foreach (SpectraFileInfo spectraFile in spectraFiles)
            {
                MyFileManager myFileManager = new(true);
                string spectraFilePath = calibratedFilePaths != null && calibratedFilePaths.TryGetValue(spectraFile, out string calibFilePath)
                    ? calibFilePath
                    : spectraFile.FullFilePathWithExtension;
                MsDataFile myMsDataFile =
                    myFileManager.LoadFile(spectraFilePath, commonParameters);
                // These parameters are hardcoded for now
                MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(
                    commonParameters.PrecursorMassTolerance,
                    MassDiffAcceptorType.OneMM, customMdac: "");
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByRT = MetaMorpheusTask
                    .GetMs2Scans(myMsDataFile, spectraFile.FullFilePathWithExtension, commonParameters)
                    .OrderBy(b => b.RetentionTime).ToArray();
                SpectralLibrary library = new(new List<string>() { spectralLibraryPath });

                MiniClassicSearchEngine mcse = new(
                    arrayOfMs2ScansSortedByRT,
                    massDiffAcceptor,
                    commonParameters,
                    library,
                    fileSpecificParameters: null,
                    maxQuantAnalysis: true);

                // This is a cursed LINQ statement
                var fileSpecificMbrPeaks =
                    mbrPeaksDict.ToDictionary(
                        kvp => kvp.Key, 
                        kvp => kvp.Value.Where(p => p.SpectraFileInfo == spectraFile).
                            FirstOrDefault()
                    ).
                        Where(kvp => kvp.Value != null).
                        ToList();

                Parallel.ForEach(threads, (i) =>
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    for (; i < fileSpecificMbrPeaks.Count; i += maxThreadsPerFile)
                    {
                        ChromatographicPeak mbrPeak = fileSpecificMbrPeaks[i].Value;
                        string maxQuantFullSeq = fileSpecificMbrPeaks[i].Key;

                        if (donorPsmDict.TryGetValue(maxQuantFullSeq, out PsmFromTsv donorPsm))
                        {
                            IEnumerable<PeptideSpectralMatch> peptideSpectralMatches =
                                mcse.SearchAroundPeak(
                                    donorPsm.PeptideWithSetModifications,
                                    mbrPeak.IsotopicEnvelopes.First().IndexedPeak.RetentionTime
                                    );

                            if (peptideSpectralMatches == null || !peptideSpectralMatches.Any())
                            {
                                bestMbrMatches.TryAdd(mbrPeak, new MbrSpectralMatch(null, mbrPeak));
                            }
                            else
                            {
                                PeptideSpectralMatch bestMatch = BestPsmForMbrPeak(peptideSpectralMatches);
                                if (bestMatch != null && mbrPeak is MaxQuantChromatographicPeak mqPeak)
                                {
                                    mqPeak.SpectralContrastAngle = bestMatch.SpectralAngle;
                                    mbrPeak = mqPeak;
                                }
                                bestMbrMatches.TryAdd(mbrPeak, 
                                    new MbrSpectralMatch(bestMatch, mbrPeak));
                            }
                        }
                    }
                });

                library.CloseConnections();
            }

            // Currently, this type of analysis doesn't make sense when analyzing MaxQuant results
            //if (bestMbrMatches.Any())
            //{
            //    List<PeptideSpectralMatch> allPsms = parameters.AllPsms.
            //        OrderByDescending(p => p.Score).
            //        ThenBy(p => p.FullFilePath).
            //        ThenBy(x => x.ScanNumber).
            //        ThenBy(p => p.FullSequence).
            //        ThenBy(p => p.ProteinAccession).ToList();

            //    AssignEstimatedPsmQvalue(bestMbrMatches, allPsms);
            //    FDRAnalysisOfMbrPsms(bestMbrMatches, allPsms, parameters, fileSpecificParameters);
            //    AssignEstimatedPsmPepQValue(bestMbrMatches, allPsms);
            //    foreach (MbrSpectralMatch match in bestMbrMatches.Values) match.FindOriginalPsm(allPsms);
            //}

            WriteMbrPsmResults(bestMbrMatches, null, outputFolder);

            return new MbrAnalysisResults(bestMbrMatches);
        }

        private static List<PeptideSpectralMatch> GetAllPeptides(
            PostSearchAnalysisParameters parameters,
            CommonParameters commonParameters,
            List<(string, CommonParameters)> fileSpecificParameters)
        {
            List<PeptideSpectralMatch> peptides = new();
            peptides = parameters.AllPsms.Where(b => b.FullSequence != null).GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();

            new FdrAnalysisEngine(peptides, parameters.NumNotches, commonParameters, fileSpecificParameters, new List<string> { parameters.SearchTaskId }, "Peptide").Run();

            if (!parameters.SearchParameters.WriteDecoys)
            {
                peptides.RemoveAll(b => b.IsDecoy);
            }
            if (!parameters.SearchParameters.WriteContaminants)
            {
                peptides.RemoveAll(b => b.IsContaminant);
            }

            double qValueCutoff = 0.01;
            if (parameters.AllPsms.Count > 100)//PEP is not computed when there are fewer than 100 psms
            {
                peptides.RemoveAll(p => p.FdrInfo.PEP_QValue > qValueCutoff);
            }
            else
            {
                peptides.RemoveAll(p => p.FdrInfo.QValue > qValueCutoff);
            }

            return peptides;
        }

        private static PeptideSpectralMatch BestPsmForMbrPeak(IEnumerable<PeptideSpectralMatch> peptideSpectralMatches)
        {
            List<PeptideSpectralMatch> nonNullPsms = peptideSpectralMatches.Where(p => p != null).ToList();

            if (nonNullPsms.Any())
            {
                // Setting Qvalue, QValueNotch, PEP, and PEP_Qvalue equal to 0 is necessary for MetaDraw to read the .psmtsv
                foreach (PeptideSpectralMatch psm in nonNullPsms)
                {
                    psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                }
                if (nonNullPsms.Select(p => p.SpectralAngle).Any(g => g != double.NaN))
                {
                    double maxSpectralAngle = nonNullPsms.Select(p => p.SpectralAngle).Max();
                    return nonNullPsms.Where(p => p.SpectralAngle == maxSpectralAngle).FirstOrDefault();
                }
                return nonNullPsms.FirstOrDefault();
            }
            return null;
        }
        private static void AssignEstimatedPsmQvalue(ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> bestMbrMatches, List<PeptideSpectralMatch> allPsms)
        {
            double[] allScores = allPsms.Select(s => s.Score).OrderByDescending(s => s).ToArray();
            double[] allQValues = allPsms.OrderByDescending(s => s.Score).Select(q => q.FdrInfo.QValue).ToArray();

            foreach (MbrSpectralMatch match in bestMbrMatches.Values)
            {
                if (match.spectralLibraryMatch != null)
                {
                    int myIndex = 0;
                    while (myIndex < (allScores.Length - 1) && allScores[myIndex] >= match.spectralLibraryMatch.Score)
                    {
                        myIndex++;
                    }
                    if (myIndex == allScores.Length - 1)
                    {
                        match.spectralLibraryMatch.FdrInfo.QValue = allQValues.Last();
                    }
                    else
                    {
                        double estimatedQ = (allQValues[myIndex] + allQValues[myIndex + 1]) / 2;
                        match.spectralLibraryMatch.FdrInfo.QValue = estimatedQ;
                    }
                }
            }
        }
        private static void FDRAnalysisOfMbrPsms(ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> bestMbrMatches, List<PeptideSpectralMatch> allPsms,
            PostSearchAnalysisParameters parameters, List<(string, CommonParameters)> fileSpecificParameters)
        {
            List<PeptideSpectralMatch> psms = bestMbrMatches.
                Select(p => p.Value.spectralLibraryMatch).
                Where(v => v != null).
                ToList();
            List<int>[] psmGroupIndices = PEP_Analysis_Cross_Validation.Get_PSM_Group_Indices(psms, 1);
            MLContext mlContext = new MLContext();
            IEnumerable<PsmData>[] PSMDataGroups = new IEnumerable<PsmData>[1];

            string searchType = "standard";
            if (psms[0].DigestionParams.Protease.Name == "top-down")
            {
                searchType = "top-down";
            }

            Dictionary<string, int> sequenceToPsmCount = PEP_Analysis_Cross_Validation.GetSequenceToPSMCount(allPsms);
            int chargeStateMode = PEP_Analysis_Cross_Validation.GetChargeStateMode(allPsms);

            Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = PEP_Analysis_Cross_Validation.ComputeHydrophobicityValues(allPsms, fileSpecificParameters, false);
            Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = PEP_Analysis_Cross_Validation.ComputeHydrophobicityValues(allPsms, fileSpecificParameters, true);
            PEP_Analysis_Cross_Validation.ComputeMobilityValues(allPsms, fileSpecificParameters);

            Dictionary<string, float> fileSpecificMedianFragmentMassErrors = PEP_Analysis_Cross_Validation.GetFileSpecificMedianFragmentMassError(allPsms);

            PSMDataGroups[0] = PEP_Analysis_Cross_Validation.CreatePsmData(searchType, fileSpecificParameters, psms, psmGroupIndices[0], sequenceToPsmCount, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode);

            string[] trainingVariables = PsmData.trainingInfos[searchType];

            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[1];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features", numberOfTrees: 400);
            var pipeline = mlContext.Transforms.Concatenate("Features", trainingVariables).Append(trainer);

            IDataView dataView = mlContext.Data.LoadFromEnumerable(PSMDataGroups[0]);

            string outputFolder = parameters.OutputFolder;

            trainedModels[0] = pipeline.Fit(dataView);

            PEP_Analysis_Cross_Validation.Compute_PSM_PEP(psms, psmGroupIndices[0], mlContext, trainedModels[0], searchType, fileSpecificParameters, sequenceToPsmCount, fileSpecificMedianFragmentMassErrors, chargeStateMode, outputFolder);
        }

        private static void AssignEstimatedPsmPepQValue(ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> bestMbrMatches, List<PeptideSpectralMatch> allPsms)
        {
            List<double> pepValues = bestMbrMatches. 
                Select(p => p.Value.spectralLibraryMatch).
                Where(p => p != null).
                OrderBy(p => p.FdrInfo.PEP).
                Select(p => p.FdrInfo.PEP).
                ToList();

            foreach (MbrSpectralMatch match in bestMbrMatches.Values)
            {
                if (match.spectralLibraryMatch == null) continue;

                int myIndex = 0;
                while (myIndex < (pepValues.Count - 1) && pepValues[myIndex] <= match.spectralLibraryMatch.FdrInfo.PEP)
                {
                    myIndex++;
                }
                if (myIndex == pepValues.Count - 1)
                {
                    match.spectralLibraryMatch.FdrInfo.PEP_QValue = pepValues.Last();
                }
                else
                {
                    double estimatedQ = (pepValues[myIndex - 1] + pepValues[myIndex]) / 2;
                    match.spectralLibraryMatch.FdrInfo.PEP_QValue = estimatedQ;
                }
            }
        }

        private static void WriteMbrPsmResults(ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> bestMbrMatches, PostSearchAnalysisParameters parameters, 
            string outputFolder = null)
        {
            string mbrOutputFolder = outputFolder ?? Path.Join(parameters.OutputFolder, mbrAnalysisFolder);
            using (var output = new StreamWriter(Path.Combine(mbrOutputFolder, @"MbrAnalysis.psmtsv")))
            {
                output.WriteLine(MbrSpectralMatch.TabSeparatedHeader);
                IEnumerable<MbrSpectralMatch> orderedMatches = bestMbrMatches.Select(p => p.Value).
                    Where(p => p.spectralLibraryMatch != null).OrderByDescending(p => p.spectralLibraryMatch.Score);
                foreach (MbrSpectralMatch match in orderedMatches)
                {
                    output.WriteLine(match.ToString());
                }
            }
        }

    }
}
