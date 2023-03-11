using EngineLayer.PsmTsv;
using FlashLFQ;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;
using Chemistry;
using EngineLayer;
using MassSpectrometry.MzSpectra;
using Proteomics.ProteolyticDigestion;

namespace TaskLayer.MbrAnalysis
{
    public class MbrAnalysisResults
    {
        public readonly ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> BestMbrMatches;
        public readonly FlashLfqResults FlashLfqResults;
        private Dictionary<string, List<string>> PeptideScoreDict;
        public Dictionary<Identification, PeptideWithSetModifications> IdsToPwsms { get; }
        public Dictionary<SpectraFileInfo, Dictionary<int, double>> Ms1TicDictionary { get; }
        public bool MaxQuantAnalysis { get; }

        /// <summary>
        /// A Tab Separated Header that is similar to the QuantifiedPeptides header,
        /// but with a new column appended at the zero-indexed 16th position
        /// </summary>
        public string PeptidesTabSeparatedHeader
        {
            get
            {
                List<string> scoreFiles = PeptideScoreDict.TryGetValue("fileName", out var scoreFileList)
                    ? scoreFileList
                    : new List<string>();

                StringBuilder sb = new();
                sb.Append(FlashLFQ.Peptide.TabSeparatedHeader(FlashLfqResults.SpectraFiles).Trim());
                foreach (string file in scoreFiles)
                {
                    sb.Append('\t');
                    sb.Append("Spectral Angle " + file);
                }

                return sb.ToString();
            }
        }

        public MbrAnalysisResults(ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> bestMbrMatches,
            FlashLfqResults flashLfqResults,
            Dictionary<Identification, PeptideWithSetModifications> idsToPwsms = null,
            Dictionary<SpectraFileInfo, Dictionary<int, double>> ms1TicDictionary = null)
        {
            BestMbrMatches = bestMbrMatches;
            FlashLfqResults = flashLfqResults;
            IdsToPwsms = idsToPwsms;
            MaxQuantAnalysis = false;
            PopulatePeptideScoreDict();
        }

        public MbrAnalysisResults(ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> bestMbrMatches)
        {
            BestMbrMatches = bestMbrMatches;
            MaxQuantAnalysis = true;
            FlashLfqResults = null;
        }

        /// <summary>
        /// Creates a dictionary mapping modified peptide sequences (string) to the MBR spectral angle for peptides
        /// in each file (list of strings, with list length == number of files)
        /// </summary>
        private void PopulatePeptideScoreDict()
        {
            PeptideScoreDict = new();
            var flashLfqPeptides = FlashLfqResults.PeptideModifiedSequences.
                Select(p => p.Key).Distinct();
            foreach (string peptide in flashLfqPeptides) PeptideScoreDict.TryAdd(peptide, new());
            PeptideScoreDict.TryAdd("fileName", new());
            foreach (SpectraFileInfo spectraFile in FlashLfqResults.SpectraFiles)
            {
                PeptideScoreDict["fileName"].Add(spectraFile.FilenameWithoutExtension);
                var peaksForFile = BestMbrMatches.Select(kvp => kvp.Key).
                    Where(p => p.SpectraFileInfo == spectraFile && p.Identifications.Any()).
                    ToDictionary(p => p.Identifications.First().ModifiedSequence, p => p);
                IEnumerable<string> peptidesNotInFile = flashLfqPeptides.Except(peaksForFile.Keys);

                foreach (var peptide in peaksForFile.Keys)
                {
                    PeptideScoreDict.TryGetValue(peptide, out List<string> scoreList);
                    BestMbrMatches.TryGetValue(peaksForFile[peptide], out var mbrSpectralMatch);
                    string score = mbrSpectralMatch.spectralLibraryMatch != null && mbrSpectralMatch.spectralLibraryMatch.SpectralAngle > -1
                        ? mbrSpectralMatch.spectralLibraryMatch.SpectralAngle.ToString()
                        : "Spectrum Not Found";
                    scoreList.Add(score);
                }

                foreach (var peptide in peptidesNotInFile)
                {
                    if (PeptideScoreDict.TryGetValue(peptide, out List<string> scoreList))
                    {
                        scoreList.Add("NA");
                    }
                }
            }
        }

        public static (double?, double?) CalculateRtShift(ChromatographicPeak acceptorPeak)
        {
            if (acceptorPeak.RtHypothesis == null || acceptorPeak.Apex == null || acceptorPeak.Apex.IndexedPeak == null) return (null, null);
            double retentionTimeShift = acceptorPeak.Apex.IndexedPeak.RetentionTime - (double)acceptorPeak.RtHypothesis;

            // Interquartile range is approximately equal to 1.35 standard deviations
            double? rtStdDev = acceptorPeak.RtInterquartileRange != null
                ? acceptorPeak.RtInterquartileRange / 1.35
                : acceptorPeak.RtStdDev;
            if (rtStdDev != null) return (retentionTimeShift, Math.Abs(retentionTimeShift / (double)rtStdDev));
            else return (retentionTimeShift, null);
        }

        public static (string apexI, string envelopeI, string peakShape) CalculateEnvelopeStatistics(ChromatographicPeak peak, PpmTolerance tolerance)
        {
            List<IndexedMassSpectralPeak> apexPeaks = peak.IsotopicEnvelopes
                .Select(e => e.IndexedPeak)
                .Where(p => tolerance.Within(p.Mz, peak.Apex.IndexedPeak.Mz))
                .OrderBy(p => p.RetentionTime)
                .ToList();

            string apexIntensityIntegral = apexPeaks
                .Sum(p => p.Intensity)
                .ToString();
            string envelopeIntensityIntegral = peak.IsotopicEnvelopes
                .Sum(e => e.IndexedPeak.Intensity)
                .ToString();

            StringBuilder apexDescription = new StringBuilder();
            apexDescription.Append("[ ");
            foreach (IndexedMassSpectralPeak imsPeak in apexPeaks)
            {
                apexDescription.Append(imsPeak.RetentionTime);
                apexDescription.Append(",");
                apexDescription.Append(imsPeak.Intensity);
                apexDescription.Append("; ");
            }
            apexDescription.Length--; // Removes the last "; " 
            apexDescription.Append("]");
            string apexIntensityVsTime = apexDescription.ToString();

            return (apexIntensityIntegral, envelopeIntensityIntegral, apexIntensityVsTime);
        }

        /// <summary>
        /// A Tab Separated Header that is similar to a ChromatographicPeak header,
        /// but with a new column appended at the zero-indexed 16th position
        /// </summary>
        public string PeaksTabSeparatedHeader
        {
            get
            {
                string[] peakHeaderSplit = ChromatographicPeak.TabSeparatedHeader.Split('\t');
                StringBuilder sb = new();
                sb.Append(string.Join('\t', peakHeaderSplit[0..16]));
                sb.Append('\t');
                sb.Append("Spectral Contrast Angle");
                sb.Append('\t');
                sb.Append("Retention Time Shift (min)");
                sb.Append('\t');
                if (MaxQuantAnalysis)
                {
                    sb.Append("Ppm Error");
                }
                else
                {
                    sb.Append("Retention Time Z-Score");

                    sb.Append('\t');
                    sb.Append("Ms2 Retention Time-Recovered Spectrum");
                    sb.Append('\t');
                    sb.Append("Precursor Mass - Recovered Spectrum");
                    sb.Append('\t');
                    sb.Append("Intensity Integral (Most Abundant Peak)");
                    sb.Append('\t');
                    sb.Append("Intensity Integral (Whole Envelope)");
                    if (IdsToPwsms != null)
                    {
                        sb.Append('\t');
                        sb.Append("Precursor Isotopic Angle");
                    }
                    //if (Ms1TicDictionary != null)
                    //{
                    //    sb.Append('\t');
                    //    sb.Append("Apex Ms1 TIC");
                    //}
                    sb.Append('\t');
                    sb.Append("Most Abundant Peak [Retention Time, Intensity]");

                }
                sb.Append('\t');
                sb.Append("Recovered Spectrum Ms2 Scan Time");
                sb.Append('\t');
                sb.Append(string.Join('\t', peakHeaderSplit[16..]));
                string header = sb.ToString();
                return header.Trim();
            }
        }

        /// <summary>
        /// Writes the peaks quantified by FlashLFQ to a .tsv file. Adds a column giving the spectralContrastAngle for MBR peaks
        /// </summary>
        /// <param name="outputFolder"> Output folder path </param>
        /// <param name="fileName"> Name of the output file (without extension) </param>
        public void WritePeakQuantificationResultsToTsv(string outputFolder, string fileName)
        {
            var fullSeqPath = Path.Combine(outputFolder, fileName + ".tsv");
            List<ChromatographicPeak> orderedPeaks = new();

            if (FlashLfqResults != null)
            {
                orderedPeaks = FlashLfqResults.Peaks.
                    SelectMany(p => p.Value)
                    .OrderBy(p => p.SpectraFileInfo.FilenameWithoutExtension)
                    .ThenByDescending(p => p.Intensity).ToList();
            }
            else
            {
                orderedPeaks = BestMbrMatches.Select(kvp => kvp.Key)
                    .OrderBy(p => p.SpectraFileInfo.FilenameWithoutExtension)
                    .ThenByDescending(p => p.Intensity).ToList();
            }

            if (fullSeqPath != null)
            {
                using (StreamWriter output = new StreamWriter(fullSeqPath))
                {
                    output.WriteLine(PeaksTabSeparatedHeader);
                    // TODO: Change this so it reflects the actual ppm tolerance
                    PpmTolerance tolerance = new PpmTolerance(5);

                    foreach (var peak in orderedPeaks)
                    {
                        string spectralContrastAngle = "Spectrum Not Found";
                        string retentionTimeShift = "";
                        string ppmError = "";
                        string rtZScore = "";

                        string recoveredSpectrumRt = "";
                        string recoveredSpectrumPrecursorMass = "";

                        string apexIntensityIntegral = "";
                        string envelopeIntensityIntegral = "";
                        string apexIntensityVsTime = "";

                        string precursorIsotopicAngle = "";

                        if (!MaxQuantAnalysis)
                        {
                            if (BestMbrMatches.TryGetValue(peak, out var mbrSpectralMatch))
                            {
                                if (mbrSpectralMatch.spectralLibraryMatch != null)
                                {
                                    spectralContrastAngle =
                                        mbrSpectralMatch.spectralLibraryMatch.SpectralAngle > -1
                                            ? mbrSpectralMatch.spectralLibraryMatch.SpectralAngle.ToString()
                                            : "Spectrum Not Found";
                                    recoveredSpectrumRt =
                                        mbrSpectralMatch.spectralLibraryMatch.ScanRetentionTime.ToString();
                                    recoveredSpectrumPrecursorMass = mbrSpectralMatch.spectralLibraryMatch
                                        .ScanPrecursorMonoisotopicPeakMz.ToString();
                                }
                            }
                            

                            (double?, double?) rtInfo = CalculateRtShift(peak);
                            retentionTimeShift = rtInfo.Item1 != null ? rtInfo.Item1.ToString() : "";
                            rtZScore = rtInfo.Item2 != null ? rtInfo.Item2.ToString() : "";

                            var envelopeResults = CalculateEnvelopeStatistics(peak, tolerance);
                            apexIntensityIntegral = envelopeResults.apexI;
                            envelopeIntensityIntegral = envelopeResults.envelopeI;
                            apexIntensityVsTime = envelopeResults.peakShape;

                            if (IdsToPwsms != null && IdsToPwsms.TryGetValue(peak.Identifications.First(), out var peakId))
                            {
                                if(peakId != null) precursorIsotopicAngle = GetIsotopeCorrelation(peakId, peak).ToString();
                            }
                            
                        }
                        else if (peak is MaxQuantChromatographicPeak mqPeak)
                        {
                            if (mqPeak.SpectralContrastAngle != null) spectralContrastAngle = mqPeak.SpectralContrastAngle.ToString();
                            if (mqPeak.RtShift != null) retentionTimeShift = mqPeak.RtShift.ToString();
                            if (mqPeak.PpmError != null) ppmError = mqPeak.PpmError.ToString();
                        }

                        string[] peakStringSplit = peak.ToString().Split('\t');
                        // Peaks generated from MaxQuant data don't have an Apex, so retention time has to be pulledfrom IndexedPeak
                        if (MaxQuantAnalysis && peakStringSplit[10].Equals("-"))
                        {
                            peakStringSplit[10] = peak.IsotopicEnvelopes.First().IndexedPeak.RetentionTime.ToString();
                        }
                        StringBuilder sb = new();
                        sb.Append(string.Join('\t', peakStringSplit[0..16]));
                        sb.Append('\t');
                        sb.Append(spectralContrastAngle);
                        sb.Append('\t');
                        sb.Append(retentionTimeShift);
                        sb.Append('\t');
                        if (MaxQuantAnalysis)
                        {
                            sb.Append(ppmError);
                        }
                        else
                        {
                            sb.Append(rtZScore);
                            sb.Append('\t');
                            sb.Append(recoveredSpectrumRt);
                            sb.Append('\t');
                            sb.Append(recoveredSpectrumPrecursorMass);
                            sb.Append('\t');
                            sb.Append(apexIntensityIntegral);
                            sb.Append('\t');
                            sb.Append(envelopeIntensityIntegral);
                            if (IdsToPwsms != null)
                            {
                                sb.Append('\t');
                                sb.Append(precursorIsotopicAngle);
                            }

                            //if (Ms1TicDictionary != null)
                            //{
                            //    sb.Append('\t');
                            //    sb.Append(precursorIsotopicAngle);
                            //}
                            sb.Append('\t');
                            sb.Append(apexIntensityVsTime);
                            

                        }
                        sb.Append('\t');
                        sb.Append(string.Join('\t', peakStringSplit[16..]));
                        output.WriteLine(sb.ToString().Trim());
                    }
                }
            }
        }

        public static double GetIsotopeCorrelation(PeptideWithSetModifications selectedPeptide, ChromatographicPeak peak)
        {
            double isotopeAngle = 0;
            double fineResolution = 0.01;
            double minimumProbability = 0.005;
            ChemicalFormula peptideFormula = null;
            try
            {
                peptideFormula = selectedPeptide.FullChemicalFormula;
            }
            //XL data throws a nullReferenceException when you try to access the FullChemicalFormula
            catch (NullReferenceException e) { }
            if (peptideFormula == null || peptideFormula.AtomCount > 0)
            {
                // calculate averagine (used for isotopic distributions for unknown modifications)
                double averageC = 4.9384;
                double averageH = 7.7583;
                double averageO = 1.4773;
                double averageN = 1.3577;
                double averageS = 0.0417;

                double averagineMass =
                    PeriodicTable.GetElement("C").AverageMass * averageC +
                    PeriodicTable.GetElement("H").AverageMass * averageH +
                    PeriodicTable.GetElement("O").AverageMass * averageO +
                    PeriodicTable.GetElement("N").AverageMass * averageN +
                    PeriodicTable.GetElement("S").AverageMass * averageS;

                if (!String.IsNullOrEmpty(selectedPeptide.BaseSequence))
                {
                    Proteomics.AminoAcidPolymer.Peptide baseSequence = new Proteomics.AminoAcidPolymer.Peptide(selectedPeptide.BaseSequence);
                    peptideFormula = baseSequence.GetChemicalFormula();
                    // add averagine for any unknown mass difference (i.e., a modification)
                    double massDiff = selectedPeptide.MonoisotopicMass - baseSequence.MonoisotopicMass;

                    // Magic numbers are bad practice, but I pulled this directly from flashLfq
                    if (Math.Abs(massDiff) >= 20)
                    {
                        double averagines = massDiff / averagineMass;

                        peptideFormula.Add("C", (int)Math.Round(averagines * averageC, 0));
                        peptideFormula.Add("H", (int)Math.Round(averagines * averageH, 0));
                        peptideFormula.Add("O", (int)Math.Round(averagines * averageO, 0));
                        peptideFormula.Add("N", (int)Math.Round(averagines * averageN, 0));
                        peptideFormula.Add("S", (int)Math.Round(averagines * averageS, 0));
                    }
                }
            }

            List<IndexedMassSpectralPeak> imsPeaksOrderedByMz = peak.IsotopicEnvelopes
                .Select(e => e.IndexedPeak)
                .Where(p => Math.Abs(p.RetentionTime - peak.Apex.IndexedPeak.RetentionTime) < 0.0001)
                .OrderBy(p => p.Mz)
                .ToList();

            if (peptideFormula != null && imsPeaksOrderedByMz != null && imsPeaksOrderedByMz.Count > 0)
            {
                IsotopicDistribution peptideDistribution = IsotopicDistribution
                    .GetDistribution(peptideFormula, fineResolution, minimumProbability);
                SpectralSimilarity isotopeSimilarity = new(
                    imsPeaksOrderedByMz.Select(p => p.Mz).ToArray(),
                    imsPeaksOrderedByMz.Select(p => p.Intensity).ToArray(),
                    peptideDistribution.Masses.ToArray(),
                    peptideDistribution.Intensities.ToArray(),
                    SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum,
                    toleranceInPpm: 10.0,
                    allPeaks: true,
                    filterOutBelowThisMz: 200);
                isotopeAngle = (double?)isotopeSimilarity.SpectralContrastAngle() ?? double.NaN;
                return isotopeAngle;
            }

            return double.NaN;

        }


        /// <summary>
        /// Writes the peptides quantified by FlashLFQ to a .tsv file. Appends extra columns giving the spectral angle
        /// for MBR identified peptides for each file
        /// </summary>
        /// <param name="outputFolder"> Output folder path </param>
        /// <param name="fileName"> Name of the output file (without extension) </param>
        public void WritePeptideQuantificationResultsToTsv(string outputFolder, string fileName)
        {
            var peptidePath = Path.Combine(outputFolder, fileName + ".tsv");

            using (StreamWriter output = new StreamWriter(peptidePath))
            {
                output.WriteLine(PeptidesTabSeparatedHeader);

                var stringPeptidePairs = FlashLfqResults.PeptideModifiedSequences.
                    OrderBy(p => p.Key);
                foreach (var peptide in stringPeptidePairs)
                {
                    string scores = PeptideScoreDict.TryGetValue(peptide.Key, out var scoreList)
                        ? string.Join('\t', scoreList)
                        : new('\t', FlashLfqResults.SpectraFiles.Count-1);
                    output.WriteLine(peptide.Value.ToString(FlashLfqResults.SpectraFiles) + '\t' + scores.TrimEnd());
                }
            }
        }
    }
}
