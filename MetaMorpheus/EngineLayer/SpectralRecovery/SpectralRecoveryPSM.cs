using System.Collections.Generic;
using System.Linq;
using System.Text;
using FlashLFQ;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.SpectralRecovery
{
    public class SpectralRecoveryPSM : PeptideSpectralMatch
    {
        public ChromatographicPeak AcceptorPeak { get; set; }

        public PeptideSpectralMatch OriginalSpectralMatch { get; private set; }

        public SpectralRecoveryPSM(
            ChromatographicPeak acceptorPeak,
            PeptideWithSetModifications peptide,
            double score,
            RecoveredMs2ScanWithSpecificMass scan,
            CommonParameters commonParameters,
            List<MatchedFragmentIon> matchedFragmentIons,
            int notch = -1,
            int scanIndex = -1,
            double xcorr = 0) :
            base(peptide, notch, score, scanIndex, scan, commonParameters, matchedFragmentIons, xcorr)
        {
            AcceptorPeak = acceptorPeak;
        }

        public void FindOriginalPsm(List<PeptideSpectralMatch> originalSearchPsms)
        {
            OriginalSpectralMatch = originalSearchPsms.FirstOrDefault(
                    p => p.FullFilePath.Equals(FullFilePath) &&
                    p.MsDataScan.OneBasedScanNumber == MsDataScan.OneBasedScanNumber &&
                    (FullSequence.Equals(p.FullSequence) ||
                    p.FullSequence.Split('|').Any(seq => seq.Equals(FullSequence))));
        }

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append(GetTabSeparatedHeader());
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
            string originalPsmQ = OriginalSpectralMatch == null ? "Spectrum Not Found" : OriginalSpectralMatch.FdrInfo.QValue.ToString();
            string originalPsmPEP = OriginalSpectralMatch == null ? "" : OriginalSpectralMatch.FdrInfo.PEP.ToString();
            string originalPsmPEPQ = OriginalSpectralMatch == null ? "" : OriginalSpectralMatch.FdrInfo.PEP_QValue.ToString();
            var sb = new StringBuilder();
            sb.Append(base.ToString().Trim());
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
