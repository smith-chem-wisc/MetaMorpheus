using EngineLayer;
using EngineLayer.SpectralRecovery;
using EngineLayer.FdrAnalysis;
using FlashLFQ;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Collections.Concurrent;
using Microsoft.ML;
using Microsoft.ML.Data;

namespace TaskLayer.SpectralRecovery
{
    public static class SpectralRecoveryRunner
    {
        public const string SpectralRecoveryFolder = "SpectralRecovery";

        /// <summary>
        /// Performs secondary analysis of MBR results by searching acceptor files for candidate spectra,
        /// and comparing those spectra to a spectral library. Results (PEP model and .psmtsv) are written to
        /// unique MbrAnalysis folder.
        /// </summary>
        /// <param name="parameters"></param>
        /// <param name="commonParameters"></param>
        /// <param name="fileSpecificParameters"></param>
        public static SpectralRecoveryResults RunSpectralRecoveryAlgorithm(
            PostSearchAnalysisParameters parameters,
            CommonParameters commonParameters,
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters)
        {

            if (!parameters.SearchParameters.DoSpectralRecovery | !parameters.SearchParameters.MatchBetweenRuns)
            {
                return null;
            }

            List<SpectraFileInfo> spectraFiles = parameters.FlashLfqResults.Peaks.Select(p => p.Key).ToList();
            List<PeptideSpectralMatch> allPeptides = GetAllPeptides(parameters, commonParameters, fileSpecificParameters);

            int maxThreadsPerFile = commonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();

            ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches = new();

            foreach (SpectraFileInfo spectraFile in spectraFiles)
            {
                List<ChromatographicPeak> fileSpecificMbrPeaks =
                    parameters.FlashLfqResults.Peaks[spectraFile].Where(p => p.IsMbrPeak).ToList();
                if (!fileSpecificMbrPeaks.Any()) break;

                MyFileManager myFileManager = new(true);

                string spectralLibraryFile = Directory
                    .GetFiles(parameters.OutputFolder, "*.*", SearchOption.AllDirectories)
                    .First(p => p.Contains("spectralLibrary"))
                    .ToString();
                string spectralLibraryPath = Path.Combine(parameters.OutputFolder, spectralLibraryFile);
                SpectralLibrary library = new(new List<string>() { spectralLibraryPath });

                MiniClassicSearchEngine mcse = new(
                    dataFile: myFileManager.LoadFile(spectraFile.FullFilePathWithExtension, commonParameters),
                    spectralLibrary: library,
                    commonParameters,
                    fileSpecificParameters.
                        Where(t => t.FileName.Equals(spectraFile.FilenameWithoutExtension)).
                        Select(t => t.Parameters).FirstOrDefault(),
                    spectraFile.FullFilePathWithExtension);

                Parallel.ForEach(threads, (i) =>
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    for (; i < fileSpecificMbrPeaks.Count; i += maxThreadsPerFile)
                    {
                        ChromatographicPeak mbrPeak = fileSpecificMbrPeaks[i];
                        PeptideSpectralMatch bestDonorPsm = allPeptides.FirstOrDefault(p =>
                            p.FullSequence == mbrPeak.Identifications.First().ModifiedSequence);
                        if (bestDonorPsm == null)
                        {
                            break;
                        }
                        PeptideWithSetModifications bestDonorPwsm = bestDonorPsm.BestMatchingPeptides.First().Peptide;

                        List<SpectralRecoveryPSM> peptideSpectralMatches = mcse.SearchAroundPeak(bestDonorPwsm, mbrPeak);

                        if (peptideSpectralMatches == null || !peptideSpectralMatches.Any())
                        {
                            bestMbrMatches.TryAdd(mbrPeak, null);
                        }
                        else
                        {
                            bestMbrMatches.TryAdd(mbrPeak, BestPsmForMbrPeak(peptideSpectralMatches));
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
                foreach (SpectralRecoveryPSM match in bestMbrMatches.Values) match.FindOriginalPsm(allPsms);
            }

            Directory.CreateDirectory(Path.Join(parameters.OutputFolder, SpectralRecoveryFolder));
            WriteSpectralRecoveryPsmResults(bestMbrMatches, parameters);

            return new SpectralRecoveryResults(bestMbrMatches, parameters.FlashLfqResults);
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

        private static SpectralRecoveryPSM BestPsmForMbrPeak(IEnumerable<SpectralRecoveryPSM> SpectralRecoveryPSMes)
        {
            List<SpectralRecoveryPSM> nonNullPsms = SpectralRecoveryPSMes.Where(p => p != null).ToList();

            if (nonNullPsms.Any())
            {
                // Setting Qvalue, QValueNotch, PEP, and PEP_Qvalue equal to 0 is necessary for MetaDraw to read the .psmtsv
                foreach (SpectralRecoveryPSM psm in nonNullPsms)
                {
                    psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                }

                return nonNullPsms.MaxBy(p => p.SpectralAngle);
            }
            return null;
        }

        private static void AssignEstimatedPsmQvalue(ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches, List<PeptideSpectralMatch> allPsms)
        {
            double[] allScores = allPsms.Select(s => s.Score).OrderByDescending(s => s).ToArray();
            double[] allQValues = allPsms.OrderByDescending(s => s.Score).Select(q => q.FdrInfo.QValue).ToArray();

            foreach (SpectralRecoveryPSM match in bestMbrMatches.Values)
            {
                if (match != null)
                {
                    int myIndex = 0;
                    while (myIndex < (allScores.Length - 1) && allScores[myIndex] >= match.Score)
                    {
                        myIndex++;
                    }
                    if (myIndex == allScores.Length - 1)
                    {
                        match.FdrInfo.QValue = allQValues.Last();
                    }
                    else
                    {
                        double estimatedQ = (allQValues[myIndex] + allQValues[myIndex + 1]) / 2;
                        match.FdrInfo.QValue = estimatedQ;
                    }
                }
            }
        }
        private static void FDRAnalysisOfMbrPsms(ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches, List<PeptideSpectralMatch> allPsms,
            PostSearchAnalysisParameters parameters, List<(string, CommonParameters)> fileSpecificParameters)
        {
            List<PeptideSpectralMatch> psms = bestMbrMatches.
                Select(p => p.Value as PeptideSpectralMatch).
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

        private static void AssignEstimatedPsmPepQValue(ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches, List<PeptideSpectralMatch> allPsms)
        {
            List<double> pepValues = bestMbrMatches. 
                Select(p => p.Value).
                Where(p => p != null).
                OrderBy(p => p.FdrInfo.PEP).
                Select(p => p.FdrInfo.PEP).
                ToList();

            foreach (SpectralRecoveryPSM match in bestMbrMatches.Values)
            {
                if (match == null) continue;

                int myIndex = 0;
                while (myIndex < (pepValues.Count - 1) && pepValues[myIndex] <= match.FdrInfo.PEP)
                {
                    myIndex++;
                }
                if (myIndex == pepValues.Count - 1)
                {
                    match.FdrInfo.PEP_QValue = pepValues.Last();
                }
                else
                {
                    double estimatedQ = (pepValues[myIndex - 1] + pepValues[myIndex]) / 2;
                    match.FdrInfo.PEP_QValue = estimatedQ;
                }
            }
        }

        private static void WriteSpectralRecoveryPsmResults(ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches, PostSearchAnalysisParameters parameters)
        {
            string mbrOutputPath = Path.Combine(Path.Join(parameters.OutputFolder, SpectralRecoveryFolder), @"RecoveredSpectra.psmtsv");
            using (var output = new StreamWriter(mbrOutputPath))
            {
                output.WriteLine(SpectralRecoveryPSM.TabSeparatedHeader);
                IEnumerable<SpectralRecoveryPSM> orderedMatches = bestMbrMatches.Select(p => p.Value).
                    Where(p => p != null).OrderByDescending(p => p.Score);
                foreach (SpectralRecoveryPSM match in orderedMatches)
                {
                    output.WriteLine(match.ToString());
                }
            }
        }

    }
}
