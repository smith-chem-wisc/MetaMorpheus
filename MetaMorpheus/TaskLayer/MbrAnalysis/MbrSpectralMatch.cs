using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using EngineLayer;
using FlashLFQ;

namespace TaskLayer.MbrAnalysis
{
    public class MbrSpectralMatch
    {
        public PeptideSpectralMatch spectralLibraryMatch { get; set; }

        public ChromatographicPeak acceptorPeak { get; set; }

        public bool originalMatchFound { get; private set; }

        public PeptideSpectralMatch originalSpectralMatch { get; private set; }

        public double RetentionTimeShift { get; private set; }

        public double? RetentionTimeZScore { get; private set; }

        public MbrSpectralMatch(PeptideSpectralMatch spectralLibraryMatch, ChromatographicPeak acceptorPeak)
        {
            this.spectralLibraryMatch = spectralLibraryMatch;
            this.acceptorPeak = acceptorPeak;
            SetRetentionTimeShift();
            originalMatchFound = false;
        }

        public void FindOriginalPsm(List<PeptideSpectralMatch> originalSearchPsms)
        {
            if (spectralLibraryMatch == null) return;
            int oneBasedScanNumber = spectralLibraryMatch.ScanNumber;
            string spectraFile = spectralLibraryMatch.FullFilePath;
            string fullSequence = spectralLibraryMatch.FullSequence;
            PeptideSpectralMatch originalPsm = null;


            originalPsm =
                originalSearchPsms.Where(
                    p => p.FullFilePath.Equals(spectraFile) &&
                         p.ScanNumber == oneBasedScanNumber &&
                         p.FullSequence != null &&
                         (p.FullSequence == fullSequence || p.FullSequence.Split('|').Any(m => m == fullSequence))
                ).FirstOrDefault();

            if (originalPsm != null)
            {
                originalMatchFound = true;
                originalSpectralMatch = originalPsm;
            }
        }

        /// <summary>
        /// Calculates the difference between the actual and expected retention times in minutes.
        /// Calculates the "Z-Score" by dividing this difference by the standard deviation of the alignment between
        /// peaks in the donor and acceptor files
        /// </summary>
        public void SetRetentionTimeShift()
        {
            if (acceptorPeak.RtHypothesis == null || acceptorPeak.Apex == null) return;
            RetentionTimeShift = acceptorPeak.Apex.IndexedPeak.RetentionTime - (double)acceptorPeak.RtHypothesis;

            // Interquartile range is approximately equal to 1.35 standard deviations
            double rtStdDev = acceptorPeak.RtStdDev ?? acceptorPeak.RtInterquartileRange / 1.35 ?? Double.NaN;
            if (rtStdDev == Double.NaN) return;
            RetentionTimeZScore = Math.Abs(RetentionTimeShift / rtStdDev);
        }

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append(PeptideSpectralMatch.GetTabSeparatedHeader());
                //sb.Append('\t');
                //sb.Append("Retention Time Shift");
                //sb.Append('\t');
                //sb.Append("Rt Shift Z-Score");
                sb.Append('\t');
                sb.Append("Initial Search Q-Value");
                sb.Append('\t');
                sb.Append("Initial Search PEP");
                sb.Append('\t');
                sb.Append("Initial Search PEP Q-Value");
                return sb.ToString();
            }
        }

        public override string ToString()
        {
            string retentionTimeShift = RetentionTimeShift.ToString();
            string retentionTimeZScore = RetentionTimeZScore == null ? "" : RetentionTimeZScore.ToString();
            string originalPsmQ = originalSpectralMatch == null ? "Spectrum Not Found" : originalSpectralMatch.FdrInfo.QValue.ToString();
            string originalPsmPEP = originalSpectralMatch == null ? "" : originalSpectralMatch.FdrInfo.PEP.ToString();
            string originalPsmPEPQ = originalSpectralMatch == null ? "" : originalSpectralMatch.FdrInfo.PEP_QValue.ToString();
            var sb = new StringBuilder();
            sb.Append(spectralLibraryMatch.ToString());
            sb.Append('\t');
            sb.Append(retentionTimeShift);
            sb.Append('\t');
            sb.Append(retentionTimeZScore);
            sb.Append('\t');
            sb.Append(originalPsmQ);
            sb.Append('\t');
            sb.Append(originalPsmPEP);
            sb.Append('\t');
            sb.Append(originalPsmPEPQ);
            return sb.ToString();
        }
    }
}
