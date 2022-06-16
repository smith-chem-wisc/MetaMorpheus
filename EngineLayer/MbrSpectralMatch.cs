using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using FlashLFQ;
using EngineLayer.FdrAnalysis;

namespace EngineLayer
{
    public class MbrSpectralMatch
    {
        public PeptideSpectralMatch spectralLibraryMatch;

        public ChromatographicPeak acceptorPeak;

        public bool originalMatchFound;

        public FdrInfo originalFdrInfo;

        public MbrSpectralMatch(PeptideSpectralMatch spectralLibraryMatch, ChromatographicPeak acceptorPeak)
        {
            this.spectralLibraryMatch = spectralLibraryMatch;
            this.acceptorPeak = acceptorPeak;
            this.originalMatchFound = false;
        }

        public void FindOriginalPsm(List<PeptideSpectralMatch> originalSearchPsms)
        {
            int oneBasedScanNumber = spectralLibraryMatch.ScanNumber;
            string spectraFile = spectralLibraryMatch.FullFilePath;

            PeptideSpectralMatch originalPsm = originalSearchPsms.Where(p => p.FullFilePath.Equals(spectraFile)).Where(p => p.ScanNumber == oneBasedScanNumber).FirstOrDefault();
            if (originalPsm != null)
            {
                originalMatchFound = true;
                originalFdrInfo = originalPsm.FdrInfo;
            }
        }

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append(PeptideSpectralMatch.GetTabSeparatedHeader());
                sb.Append('\t');
                sb.Append(ChromatographicPeak.TabSeparatedHeader);
                return sb.ToString();
            }
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(spectralLibraryMatch.ToString());
            sb.Append('\t');
            sb.Append(acceptorPeak.ToString());
            return sb.ToString();
        }



    }
}
