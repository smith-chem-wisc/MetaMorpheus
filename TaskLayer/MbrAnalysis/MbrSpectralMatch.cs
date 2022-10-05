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

        public MbrSpectralMatch(PeptideSpectralMatch spectralLibraryMatch, ChromatographicPeak acceptorPeak)
        {
            this.spectralLibraryMatch = spectralLibraryMatch;
            this.acceptorPeak = acceptorPeak;
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

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append(PeptideSpectralMatch.GetTabSeparatedHeader());
                sb.Append('\t');
                sb.Append("Initial Search Q-Value");
                return sb.ToString();
            }
        }

        public override string ToString()
        {
            string originalPsmString = originalSpectralMatch == null ? "NA" : originalSpectralMatch.FdrInfo.QValue.ToString();
            var sb = new StringBuilder();
            sb.Append(spectralLibraryMatch.ToString());
            sb.Append('\t');
            sb.Append(originalPsmString);
            return sb.ToString();
        }
    }
}
