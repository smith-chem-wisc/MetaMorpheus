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
using Microsoft.ML;
using Microsoft.ML.Data;
using Omics;

namespace TaskLayer.MbrAnalysis
{
    
    public static class SpectralRecoveryRunner
    {
        public const string spectralRecoveryFolder = "SpectralRecovery";

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
            List<SpectralMatch> allPeptides = GetAllPeptides(parameters, commonParameters, fileSpecificParameters);

            int maxThreadsPerFile = commonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();

            ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches = new();

            foreach (SpectraFileInfo spectraFile in spectraFiles)
            {
                List<ChromatographicPeak> fileSpecificMbrPeaks =
                    parameters.FlashLfqResults.Peaks[spectraFile].Where(p => p.IsMbrPeak).ToList();
                if (fileSpecificMbrPeaks == null || (!fileSpecificMbrPeaks.Any())) break;

                MyFileManager myFileManager = new(true);
                MsDataFile myMsDataFile =
                    myFileManager.LoadFile(spectraFile.FullFilePathWithExtension, commonParameters);
                MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(
                    commonParameters.PrecursorMassTolerance,
                    parameters.SearchParameters.MassDiffAcceptorType, parameters.SearchParameters.CustomMdac);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByRT = MetaMorpheusTask
                    .GetMs2Scans(myMsDataFile, spectraFile.FullFilePathWithExtension, commonParameters)
                    .OrderBy(b => b.RetentionTime).ToArray();

                var list = Directory.GetFiles(parameters.OutputFolder, "*.*", SearchOption.AllDirectories);
                string matchingvalue = list.Where(p => p.Contains("SpectralLibrary")).First().ToString();
                string spectralLibraryPath = Path.Combine(parameters.OutputFolder, matchingvalue);
                //var updatedLib = new SpectralLibrary(new List<string> { Path.Combine(thisTaskOutputFolder, matchingvalue) });
                //string spectralLibraryPath = Path.Combine(parameters.OutputFolder, @"SpectralLibrary.msp");
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
                        SpectralMatch bestDonorPsm = allPeptides.Where(p =>
                            p.FullSequence == mbrPeak.Identifications.First().ModifiedSequence).FirstOrDefault();
                        if (bestDonorPsm == null)
                        {
                            break;
                        }
                        IBioPolymerWithSetMods bestDonorPwsm = bestDonorPsm.BestMatchingBioPolymersWithSetMods.First().Peptide;

                        IEnumerable<SpectralMatch> peptideSpectralMatches =
                            mcse.SearchAroundPeak(bestDonorPwsm, mbrPeak.Apex.IndexedPeak.RetentionTime);

                        if (peptideSpectralMatches == null || !peptideSpectralMatches.Any())
                        {
                            bestMbrMatches.TryAdd(mbrPeak, new SpectralRecoveryPSM(null, mbrPeak));
                        }
                        else
                        {
                            bestMbrMatches.TryAdd(mbrPeak,
                                new SpectralRecoveryPSM(BestPsmForMbrPeak(peptideSpectralMatches), mbrPeak));
                        }
                    }
                });

                library.CloseConnections();
            }

            if (bestMbrMatches.Any())
            {
                List<SpectralMatch> allPsms = parameters.AllPsms.
                    OrderByDescending(p => p).ToList();

                FDRAnalysisOfMbrPsms(bestMbrMatches, allPsms, parameters, fileSpecificParameters);

                foreach (SpectralRecoveryPSM match in bestMbrMatches.Values) match.FindOriginalPsm(allPsms);
            }

            Directory.CreateDirectory(Path.Join(parameters.OutputFolder, spectralRecoveryFolder));
            WriteSpectralRecoveryPsmResults(bestMbrMatches, parameters);

            return new SpectralRecoveryResults(bestMbrMatches, parameters.FlashLfqResults);
        }

        private static List<SpectralMatch> GetAllPeptides(
            PostSearchAnalysisParameters parameters,
            CommonParameters commonParameters,
            List<(string, CommonParameters)> fileSpecificParameters)
        {
            var peptides = parameters.AllPsms;
            PostSearchAnalysisTask postProcessing = new PostSearchAnalysisTask
            {
                Parameters = parameters,
                FileSpecificParameters = fileSpecificParameters,
                CommonParameters = commonParameters
            };

            FilteredPsms.Filter(peptides,
                commonParameters,
                includeDecoys: false,
                includeContaminants: false,
                includeAmbiguous: false,
                includeHighQValuePsms: false);

            return peptides;
        }

        private static SpectralMatch BestPsmForMbrPeak(IEnumerable<SpectralMatch> peptideSpectralMatches)
        {
            List<SpectralMatch> nonNullPsms = peptideSpectralMatches.Where(p => p != null).ToList();

            if (nonNullPsms.Any())
            {
                // Setting Qvalue, QValueNotch, PEP, and PEP_Qvalue equal to 0 is necessary for MetaDraw to read the .psmtsv
                foreach (SpectralMatch psm in nonNullPsms)
                {
                    psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
                    psm.PeptideFdrInfo = psm.PsmFdrInfo;
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
        private static void AssignEstimatedPsmQvalue(ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches, List<SpectralMatch> allPsms)
        {
            double[] allScores = allPsms.Select(s => s.Score).OrderByDescending(s => s).ToArray();
            double[] allQValues = allPsms.OrderByDescending(s => s.Score).Select(q => q.PsmFdrInfo.QValue).ToArray();

            foreach (SpectralRecoveryPSM match in bestMbrMatches.Values)
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
        private static void FDRAnalysisOfMbrPsms(ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches, List<SpectralMatch> allPsms,
            PostSearchAnalysisParameters parameters, List<(string, CommonParameters)> fileSpecificParameters)
        {
            List<SpectralMatch> psms = bestMbrMatches.
                Select(p => p.Value.spectralLibraryMatch).
                Where(v => v != null).
                ToList();

            new FdrAnalysisEngine(psms, parameters.NumNotches, fileSpecificParameters.First().Item2, fileSpecificParameters,
                new List<string> { parameters.SearchTaskId }, analysisType: "PSM", doPEP: true, outputFolder: parameters.OutputFolder).Run();

        }

        private static void WriteSpectralRecoveryPsmResults(ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches, PostSearchAnalysisParameters parameters)
        {
            string mbrOutputPath = Path.Combine(Path.Join(parameters.OutputFolder, spectralRecoveryFolder), @"RecoveredSpectra.psmtsv");
            using (var output = new StreamWriter(mbrOutputPath))
            {
                output.WriteLine(SpectralRecoveryPSM.TabSeparatedHeader);
                IEnumerable<SpectralRecoveryPSM> orderedMatches = bestMbrMatches.Select(p => p.Value).
                    Where(p => p.spectralLibraryMatch != null).OrderByDescending(p => p.spectralLibraryMatch.Score);
                foreach (SpectralRecoveryPSM match in orderedMatches)
                {
                    output.WriteLine(match.ToString());
                }
            }
        }

    }
}
