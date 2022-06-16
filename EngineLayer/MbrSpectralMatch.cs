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

        public FdrInfo originalFdrInfo { get; private set; }

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
                List<string> originalFdrHeader = new();
                originalFdrHeader.Add("Original Cumulative Target");
                originalFdrHeader.Add("Original Cumulative Decoy");
                originalFdrHeader.Add("Original Cumulative Target Notch");
                originalFdrHeader.Add("Original Cumulative Decoy Notch");
                originalFdrHeader.Add("Original QValue");
                originalFdrHeader.Add("Original QValue Notch");
                originalFdrHeader.Add("Original PEP");
                originalFdrHeader.Add("Original PEP_QValue");
                return string.Join('\t', originalFdrHeader);
           }
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(spectralLibraryMatch.ToString());
            sb.Append('\t');
            sb.Append(OriginalFdrWriter());
            sb.Append('\t');
            sb.Append(acceptorPeak.ToString());
            return sb.ToString();
        }

        private string OriginalFdrWriter()
        {
            string cumulativeTarget = " ";
            string cumulativeDecoy = " ";
            string cumulativeTargetNotch = " ";
            string cumulativeDecoyNotch = " ";
            string qValue = " ";
            string qValueNotch = " ";
            string PEP = " ";
            string PEP_Qvalue = " ";

            if (originalMatchFound && originalFdrInfo != null)
            {
                cumulativeTarget = originalFdrInfo.CumulativeTarget.ToString(CultureInfo.InvariantCulture);
                cumulativeDecoy = originalFdrInfo.CumulativeDecoy.ToString(CultureInfo.InvariantCulture);
                cumulativeTargetNotch = originalFdrInfo.CumulativeTargetNotch.ToString(CultureInfo.InvariantCulture);
                qValue = originalFdrInfo.QValue.ToString("F6", CultureInfo.InvariantCulture);
                cumulativeDecoyNotch = originalFdrInfo.CumulativeDecoyNotch.ToString(CultureInfo.InvariantCulture);
                qValueNotch = originalFdrInfo.QValueNotch.ToString("F6", CultureInfo.InvariantCulture);
                PEP = originalFdrInfo.PEP.ToString();
                PEP_Qvalue = originalFdrInfo.PEP_QValue.ToString();
            }

            return string.Join('\t', new List<string> { cumulativeTarget, cumulativeDecoy,
                                                        cumulativeTargetNotch, cumulativeDecoyNotch,
                                                        qValue, qValueNotch, PEP, PEP_Qvalue});

        }



    }
}
