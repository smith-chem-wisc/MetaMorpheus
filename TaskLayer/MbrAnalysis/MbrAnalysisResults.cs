using FlashLFQ;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer.MbrAnalysis
{
    public class MbrAnalysisResults
    {
        public readonly ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> BestMbrMatches;
        public readonly FlashLfqResults FlashLfqResults;
        private Dictionary<string, List<string>> PeptideScoreDict;

        /// <summary>
        /// A Tab Separated Header that is similar to a ChromatographicPeak header,
        /// but with a new column appended at the zero-indexed 16th position
        /// </summary>
        public static string TabSeparatedHeader
        {
            get
            {
                string[] peakHeaderSplit = ChromatographicPeak.TabSeparatedHeader.Split('\t');
                StringBuilder sb = new();
                sb.Append(string.Join('\t', peakHeaderSplit[0..15]));
                sb.Append("Spectral Contrast Angle\t");
                sb.Append(string.Join('\t', peakHeaderSplit[16..(peakHeaderSplit.Length-1)]));
                string header = sb.ToString();
                return header.Trim();
            }
        }

        public MbrAnalysisResults(ConcurrentDictionary<ChromatographicPeak, MbrSpectralMatch> bestMbrMatches, FlashLfqResults flashLfqResults)
        {
            BestMbrMatches = bestMbrMatches;
            FlashLfqResults = flashLfqResults;
            PopulatePeptideScoreDict();
        }

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
                //IEnumerable<string> peptidesInFile = peaksForFile.Where(p => p.Identifications.Any())
                //    .Select(p => p.Identifications.First().ModifiedSequence);
                IEnumerable<string> peptidesNotInFile = flashLfqPeptides.Except(peaksForFile.Keys);

                foreach (var peptide in peaksForFile.Keys)
                {
                    PeptideScoreDict.TryGetValue(peptide, out List<string> scoreList);
                    BestMbrMatches.TryGetValue(peaksForFile[peptide], out var mbrSpectralMatch);
                    scoreList.Add(mbrSpectralMatch.spectralLibraryMatch.SpectralAngle.ToString());
                }

                foreach (var peptide in peptidesNotInFile)
                {
                    PeptideScoreDict.TryGetValue(peptide, out List<string> scoreList);
                    scoreList.Add("NA");
                }
            }
        }

        /// <summary>
        /// Writes the peaks quantified by FlashLFQ to a .tsv file. Adds a column giving the spectralContrastAngle for MBR peaks
        /// </summary>
        /// <param name="outputFolder"> Output folder path </param>
        /// <param name="fileName"> Name of the output file </param>
        private void WritePeakQuantificationResultsToTsv(string outputFolder, string fileName)
        {
            var fullSeqPath = Path.Combine(outputFolder, fileName + ".tsv");

            IEnumerable<ChromatographicPeak> orderedPeaks = FlashLfqResults.Peaks.
                SelectMany(p => p.Value)
                .OrderBy(p => p.SpectraFileInfo.FilenameWithoutExtension)
                .ThenByDescending(p => p.Intensity);


            if (fullSeqPath != null)
            {
                using (StreamWriter output = new StreamWriter(fullSeqPath))
                {
                    output.WriteLine();

                    foreach (var peak in orderedPeaks)
                    {
                        string spectralContrastAngle = "NA";
                        if (BestMbrMatches.TryGetValue(peak, out var mbrSpectralMatch))
                        {
                            spectralContrastAngle = mbrSpectralMatch.spectralLibraryMatch.SpectralAngle.ToString();
                        }

                        string[] peakStringSplit = peak.ToString().Split('\t');
                        StringBuilder sb = new();
                        sb.Append(string.Join('\t', peakStringSplit[0..15]));
                        sb.Append(spectralContrastAngle);
                        sb.Append('t');
                        sb.Append(string.Join('\t', peakStringSplit[16..(peakStringSplit.Length - 1)]));
                        output.WriteLine(sb.ToString().Trim());
                    }
                }
            }
        }

        

        private void WritePeptideQuantificationResultsToTsv(string outputFolder, string fileName)
        {
            var peaksPath = Path.Combine(outputFolder, fileName + ".tsv");

            FlashLfqResults.PeptideModifiedSequences

            flashLFQResults.WriteResults(peaksPath, null, null, null, true);

            FinishedWritingFile(peaksPath, nestedIds);
        }
    }
}
