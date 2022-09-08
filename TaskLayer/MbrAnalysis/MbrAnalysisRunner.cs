using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using EngineLayer.Localization;
using EngineLayer.ModificationAnalysis;
using FlashLFQ;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics.Distributions;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;
using System.Collections.Concurrent;
using Microsoft.ML;
using Microsoft.ML.Data;
using TaskLayer.MbrAnalysis;

namespace TaskLayer.MbrAnalysis
{
    
    public static class MbrAnalysisRunner
    {
        public const string mbrAnalysisFolder = "MbrAnalysis";

        public static void RunMbrAnalysis(
            PostSearchAnalysisParameters parameters,
            CommonParameters commonParameters,
            List<(string, CommonParameters)> fileSpecificParameters)
        {

            if (!parameters.SearchParameters.DoMbrAnalysis | !parameters.SearchParameters.MatchBetweenRuns)
            {
                return;
            }

            List<SpectraFileInfo> spectraFiles = parameters.FlashLfqResults.Peaks.Select(p => p.Key).ToList();
            List<PeptideSpectralMatch> allPeptides = GetAllPeptides(parameters, commonParameters, fileSpecificParameters);

            int maxThreadsPerFile = commonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();

            ConcurrentBag<MbrSpectralMatch> bestMbrMatches = new();

            foreach (SpectraFileInfo spectraFile in spectraFiles)
            {
                List<ChromatographicPeak> fileSpecificMbrPeaks = parameters.FlashLfqResults.Peaks[spectraFile].Where(p => p.IsMbrPeak).ToList();
                if (fileSpecificMbrPeaks == null || (!fileSpecificMbrPeaks.Any())) break;

                MyFileManager myFileManager = new(true);
                MsDataFile myMsDataFile = myFileManager.LoadFile(spectraFile.FullFilePathWithExtension, commonParameters);
                MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(commonParameters.PrecursorMassTolerance,
                    parameters.SearchParameters.MassDiffAcceptorType, parameters.SearchParameters.CustomMdac);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByRT = MetaMorpheusTask.GetMs2Scans(myMsDataFile, spectraFile.FullFilePathWithExtension, commonParameters)
                    .OrderBy(b => b.RetentionTime).ToArray();
                double[] arrayOfRTs = arrayOfMs2ScansSortedByRT.Select(p => p.TheScan.RetentionTime).ToArray();
                string spectralLibraryPath = Path.Combine(parameters.OutputFolder, @"spectralLibrary.msp");
                SpectralLibrary library = new(new List<string>() { spectralLibraryPath });

                Parallel.ForEach(threads, (i) =>
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }
                    for (; i < fileSpecificMbrPeaks.Count; i += maxThreadsPerFile)
                    {
                        ChromatographicPeak mbrPeak = fileSpecificMbrPeaks[i];
                        PeptideSpectralMatch bestDonorPsm = allPeptides.
                            Where(p => p.FullSequence == mbrPeak.Identifications.First().ModifiedSequence).First();
                        if (bestDonorPsm == null) break;
                        PeptideWithSetModifications bestDonorPwsm = bestDonorPsm.BestMatchingPeptides.First().Peptide;

                        Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass =
                            GetScansInWindow(mbrPeak, arrayOfRTs, arrayOfMs2ScansSortedByRT);
                        PeptideSpectralMatch[] peptideSpectralMatches = new PeptideSpectralMatch[arrayOfMs2ScansSortedByRT.Count()];

                        MiniClassicSearchEngine mcse = new(bestDonorPwsm, peptideSpectralMatches, arrayOfMs2ScansSortedByMass,
                            parameters.VariableModifications, parameters.FixedModifications, massDiffAcceptor,
                            commonParameters, fileSpecificParameters, library, new List<string> { parameters.SearchTaskId });
                        mcse.Run();

                        if (peptideSpectralMatches.Any())
                        {
                            bestMbrMatches.Add(new MbrSpectralMatch(BestPsmForMbrPeak(peptideSpectralMatches), mbrPeak));
                        }
                        else
                        {
                            bestMbrMatches.Add(new MbrSpectralMatch(null, mbrPeak));
                        }
                    }
                });
            }
            List<PeptideSpectralMatch> allPsms = parameters.AllPsms.OrderByDescending(p => p.Score).ThenBy(p => p.FdrInfo.QValue).
                ThenBy(p => p.FullFilePath).ThenBy(x => x.ScanNumber).ThenBy(p => p.FullSequence).ThenBy(p => p.ProteinAccession).ToList();
            
            AssignEstimatedPsmQvalue(bestMbrMatches, allPsms);
            FDRAnalysisOfMbrPsms(bestMbrMatches, allPsms, parameters, fileSpecificParameters);
            AssignEstimatedPsmPepQValue(bestMbrMatches, allPsms);
            foreach (MbrSpectralMatch match in bestMbrMatches) match.FindOriginalPsm(allPsms);
            WriteMbrPsmResults(bestMbrMatches, parameters);
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
            peptides.RemoveAll(p => p.FdrInfo.QValue > commonParameters.QValueOutputFilter);
            return peptides;
        }

        private static PeptideSpectralMatch BestPsmForMbrPeak(PeptideSpectralMatch[] peptideSpectralMatches)
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

        /// <summary>
        /// Finds MS2 scans falling within the relevant time window
        /// </summary>
        /// <param name="mbrPeak"> Acceptor peak </param>
        /// <param name="arrayOfRTs"> </param>
        /// <param name="arrayOfMs2ScansSortedByRT"> </param>
        /// <returns> An array of MS2 scans falling within a 2 minute retention time window of the Acceptor peak apex.
        ///           This array is sorted by precursor mass. </returns>
        private static Ms2ScanWithSpecificMass[] GetScansInWindow(ChromatographicPeak mbrPeak, double[] arrayOfRTs, Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByRT)
        {
            double apexRT = mbrPeak.Apex.IndexedPeak.RetentionTime;
            double peakHalfWidth = 1.0; //Placeholder value to determine retention time window
            int startIndex = Array.BinarySearch(arrayOfRTs, apexRT - peakHalfWidth);
            if (startIndex < 0)
                startIndex = ~startIndex;
            int endIndex = Array.BinarySearch(arrayOfRTs, apexRT + peakHalfWidth);
            if (endIndex < 0)
                endIndex = ~endIndex;

            return arrayOfMs2ScansSortedByRT[startIndex..endIndex].OrderBy(b => b.PrecursorMass).ToArray();
        }

        private static void AssignEstimatedPsmQvalue(ConcurrentBag<MbrSpectralMatch> bestMbrMatches, List<PeptideSpectralMatch> allPsms)
        {
            double[] allScores = allPsms.Select(s => s.Score).OrderByDescending(s => s).ToArray();
            double[] allQValues = allPsms.OrderByDescending(s => s.Score).Select(q => q.FdrInfo.QValue).ToArray();

            foreach (MbrSpectralMatch match in bestMbrMatches)
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
        private static void FDRAnalysisOfMbrPsms(ConcurrentBag<MbrSpectralMatch> bestMbrMatches, List<PeptideSpectralMatch> allPsms,
            PostSearchAnalysisParameters parameters, List<(string, CommonParameters)> fileSpecificParameters)
        {
            List<PeptideSpectralMatch> psms = bestMbrMatches.Select(p => p.spectralLibraryMatch).Where(v => v != null).ToList();
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
            Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE = PEP_Analysis_Cross_Validation.ComputeMobilityValues(allPsms, fileSpecificParameters);

            Dictionary<string, float> fileSpecificMedianFragmentMassErrors = PEP_Analysis_Cross_Validation.GetFileSpecificMedianFragmentMassError(allPsms);

            PSMDataGroups[0] = PEP_Analysis_Cross_Validation.CreatePsmData(searchType, fileSpecificParameters, psms, psmGroupIndices[0], sequenceToPsmCount, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode);

            string[] trainingVariables = PsmData.trainingInfos[searchType];

            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[1];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features", numberOfTrees: 400);
            var pipeline = mlContext.Transforms.Concatenate("Features", trainingVariables).Append(trainer);

            IDataView dataView = mlContext.Data.LoadFromEnumerable(PSMDataGroups[0]);

            string outputFolder = Path.Join(parameters.OutputFolder, mbrAnalysisFolder);

            trainedModels[0] = pipeline.Fit(dataView);

            int ambiguousPeptidesResolved = PEP_Analysis_Cross_Validation.Compute_PSM_PEP(psms, psmGroupIndices[0], mlContext, trainedModels[0], searchType, fileSpecificParameters, sequenceToPsmCount, fileSpecificMedianFragmentMassErrors, chargeStateMode, outputFolder);

        }

        private static void AssignEstimatedPsmPepQValue(ConcurrentBag<MbrSpectralMatch> bestMbrMatches, List<PeptideSpectralMatch> allPsms)
        {
            List<double> pepValues = bestMbrMatches.Select(p => p.spectralLibraryMatch).Where(p => p != null).OrderBy(p => p.FdrInfo.PEP).Select(p => p.FdrInfo.PEP).ToList();
            foreach (MbrSpectralMatch match in bestMbrMatches)
            {
                if (match.spectralLibraryMatch != null)
                {
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
        }

        private static void WriteMbrPsmResults(ConcurrentBag<MbrSpectralMatch> bestMbrMatches, PostSearchAnalysisParameters parameters)
        {
            string mbrOutputPath = Path.Combine(Path.Join(parameters.OutputFolder, mbrAnalysisFolder), @"MbrAnalysis.psmtsv");
            using (var output = new StreamWriter(mbrOutputPath))
            {
                output.WriteLine(MbrSpectralMatch.TabSeparatedHeader);
                IEnumerable<MbrSpectralMatch> orderedMatches = bestMbrMatches.Where(p => p.spectralLibraryMatch != null).OrderByDescending(p => p.spectralLibraryMatch.Score);
                foreach (MbrSpectralMatch match in orderedMatches)
                {
                    output.WriteLine(match.ToString());
                }
            }
        }

    }
}
