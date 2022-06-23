using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Globalization;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using FlashLFQ;
using EngineLayer.FdrAnalysis;

namespace EngineLayer
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
            this.originalMatchFound = false;
        }

        public void FindOriginalPsm(List<PeptideSpectralMatch> originalSearchPsms)
        {
            if (spectralLibraryMatch == null) return;
            int oneBasedScanNumber = spectralLibraryMatch.ScanNumber;
            string spectraFile = spectralLibraryMatch.FullFilePath;
            string fullSequence = spectralLibraryMatch.FullSequence;

            PeptideSpectralMatch originalPsm = originalSearchPsms
                .Where(p => p.FullFilePath.Equals(spectraFile))
                .Where(p => p.ScanNumber == oneBasedScanNumber)
                .Where(p => p.FullSequence == fullSequence)
                .FirstOrDefault();
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
                sb.Append(OriginalFdrHeader);
                sb.Append('\t');
                sb.Append(ChromatographicPeak.TabSeparatedHeader);
                return sb.ToString();
            }
        }

        public static string OriginalFdrHeader
        {
            get
            {
                string psmHeader = PeptideSpectralMatch.GetTabSeparatedHeader();
                string[] psmHeaderSplit = psmHeader.Split('\t');
                string[] newHeaderSplit = new string[psmHeaderSplit.Length];
                for (int i = 0; i < psmHeaderSplit.Length; i++) newHeaderSplit[i] = "Original " + psmHeaderSplit[i];
                string newHeader = string.Join('\t', newHeaderSplit); List<string> originalFdrHeader = new();
                return newHeader;
           }
        }

        public override string ToString()
        {
            string nullPsm = new string('\t', 55);
            string originalPsmString = originalSpectralMatch == null ? nullPsm : originalSpectralMatch.ToString();
            var sb = new StringBuilder();
            sb.Append(spectralLibraryMatch.ToString());
            sb.Append('\t');
            sb.Append(originalPsmString);
            sb.Append('\t');
            sb.Append(acceptorPeak.ToString());
            return sb.ToString();
        }
    }
}
