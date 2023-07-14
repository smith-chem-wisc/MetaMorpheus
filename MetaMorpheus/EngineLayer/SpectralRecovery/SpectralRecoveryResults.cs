using FlashLFQ;
using EngineLayer.SpectralRecovery;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using EngineLayer;
using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry.MzSpectra;
using Proteomics.ProteolyticDigestion;
using System;
using System.Globalization;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.SpectralRecovery
{
    public class SpectralRecoveryResults
    {
        public readonly ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> BestMbrMatches;
        public readonly FlashLfqResults FlashLfqResults;
        public readonly IEnumerable<ChromatographicPeak> MbrPeaks; 
        private Dictionary<string, List<string>> PeptideScoreDict;

        /// <summary>
        /// A Tab Separated Header that is similar to a ChromatographicPeak header,
        /// but with a new column appended at the zero-indexed 16th position
        /// </summary>
        public static string PeaksTabSeparatedHeader
        {
            get
            {
                string[] peakHeaderSplit = ChromatographicPeak.TabSeparatedHeader.Split('\t');
                StringBuilder sb = new();
                sb.Append(string.Join('\t', peakHeaderSplit[0..16]));
                sb.Append('\t');
                sb.Append("Spectral Contrast Angle");
                sb.Append('\t');
                sb.Append(string.Join('\t', peakHeaderSplit[16..]));
                string header = sb.ToString();
                return header.Trim();
            }
        }

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

        public SpectralRecoveryResults(ConcurrentDictionary<ChromatographicPeak, SpectralRecoveryPSM> bestMbrMatches, FlashLfqResults flashLfqResults)
        {
            BestMbrMatches = bestMbrMatches;
            FlashLfqResults = flashLfqResults;
            MbrPeaks = FlashLfqResults != null
                ? FlashLfqResults.Peaks
                    .SelectMany(p => p.Value)
                    .OrderBy(p => p.SpectraFileInfo.FilenameWithoutExtension)
                    .ThenByDescending(p => p.Intensity)
                : bestMbrMatches.Keys
                    .OrderBy(p => p.SpectraFileInfo.FilenameWithoutExtension)
                    .ThenByDescending(p => p.Intensity);
            PopulatePeptideScoreDict();
        }

        /// <summary>
        /// Creates a dictionary mapping modified peptide sequences (string) to the MBR spectral angle for peptides
        /// in each file (list of strings, with list length == number of files)
        /// </summary>
        private void PopulatePeptideScoreDict()
        {
            if (FlashLfqResults == null) return; //TODO: Make this work with MaxQuantResults

            PeptideScoreDict = FlashLfqResults
                .PeptideModifiedSequences
                .Select(p => p.Key)
                .Distinct()
                .ToDictionary(peptide => peptide, x => new List<string>());

            // This is used for when writing the header for the AllQuantPeptidesFile
            PeptideScoreDict.Add("fileName", new());

            foreach (SpectraFileInfo spectraFile in FlashLfqResults.SpectraFiles)
            {
                // Adding filenames in this way ensures the column label & values are concordant
                PeptideScoreDict["fileName"].Add(spectraFile.FilenameWithoutExtension);
                var peaksForFile = BestMbrMatches
                    .Select(kvp => kvp.Key)
                    .Where(p => p.SpectraFileInfo.Equals(spectraFile) && p.Identifications.Any())
                    .ToDictionary(p => p.Identifications.First().ModifiedSequence, p => p);

                IEnumerable<string> peptidesNotInFile = PeptideScoreDict.Keys.Except(peaksForFile.Keys);

                foreach (var peptide in peaksForFile.Keys)
                {
                    if(!PeptideScoreDict.TryGetValue(peptide, out List<string> scoreList) 
                       | !BestMbrMatches.TryGetValue(peaksForFile[peptide], out var mbrSpectralMatch))
                    {
                        continue;
                    }
                    string score = mbrSpectralMatch != null && mbrSpectralMatch.SpectralAngle > -1
                        ? mbrSpectralMatch.SpectralAngle.ToString()
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

        /// <summary>
        /// Writes the peaks quantified by FlashLFQ to a .tsv file. Adds a column giving the spectralContrastAngle for MBR peaks
        /// </summary>
        /// <param name="outputFolder"> Output folder path </param>
        /// <param name="fileName"> Name of the output file (without extension) </param>
        public void WritePeakQuantificationResultsToTsv(string outputFolder, string fileName)
        {
            var fullSeqPath = Path.Combine(outputFolder, fileName + ".tsv");

            using StreamWriter output = new StreamWriter(fullSeqPath);
            output.WriteLine(PeaksTabSeparatedHeader);

            foreach (var peak in MbrPeaks)
            {
                string spectralContrastAngle = 
                    BestMbrMatches.TryGetValue(peak, out var mbrSpectralMatch) && mbrSpectralMatch != null && mbrSpectralMatch.SpectralAngle > -1
                    ? mbrSpectralMatch.SpectralAngle.ToString()
                    : "Spectrum Not Found";

                string[] peakStringSplit = peak.ToString().Split('\t');
                StringBuilder sb = new();
                sb.Append(string.Join('\t', peakStringSplit[0..16]));
                sb.Append('\t');
                sb.Append(spectralContrastAngle);
                sb.Append('\t');
                sb.Append(string.Join('\t', peakStringSplit[16..]));
                output.WriteLine(sb.ToString().Trim());
            }
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

                var stringPeptidePairs = FlashLfqResults
                    .PeptideModifiedSequences
                    .OrderBy(p => p.Key);
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
