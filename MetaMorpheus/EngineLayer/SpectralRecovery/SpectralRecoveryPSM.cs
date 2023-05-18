using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.SpectralRecovery
{
    public class SpectralRecoveryPSM : PeptideSpectralMatch
    {
        public ChromatographicPeak AcceptorPeak { get; set; }

        public PeptideSpectralMatch OriginalSpectralMatch { get; private set; }
        /// <summary>
        /// In spectral recovery, we see situations where a mass was selected for fragmentation
        /// but could not be deconvoluted. As a result, base.ScanPrecursorMass == 0. However, the
        /// RecoveredMs2ScanWithSpecificMass contains information on the precursor peak that was closest
        /// to the theoretical precursor m/z. ClosestPrecursorPeak stores that info
        /// </summary>
        public MzPeak ClosestPrecursorPeak { get; }
        public PeptideWithSetModifications DonorPeptide { get; }

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
            DonorPeptide = peptide;
            ClosestPrecursorPeak = scan.PeakClosestToDonor;
        }

        public void FindOriginalPsm(List<PeptideSpectralMatch> originalSearchPsms)
        {
            OriginalSpectralMatch = originalSearchPsms.FirstOrDefault(
                    p => p.FullFilePath.Equals(FullFilePath) &&
                    p.MsDataScan.OneBasedScanNumber == MsDataScan.OneBasedScanNumber &&
                    (FullSequence.Equals(p.FullSequence) ||
                    p.FullSequence.Split('|').Any(seq => seq.Equals(FullSequence))));
        }

        public static string TabSeparatedHeader => string.Join('\t', SpectralRecoveryDataDictionary(null, null).Keys);

        public override string ToString()
        {
            return ToString(new Dictionary<string, int>());
        }

        public override string ToString(IReadOnlyDictionary<string, int> modsToWritePruned)
        {
            return string.Join('\t', SpectralRecoveryDataDictionary(this, modsToWritePruned).Values);
        }

        public static Dictionary<string, string> SpectralRecoveryDataDictionary(SpectralRecoveryPSM srPsm,
            IReadOnlyDictionary<string, int> modsToWritePruned)
        {
            // Get information from base PSM class
            Dictionary<string, string> psmDictionary = DataDictionary(srPsm, modsToWritePruned);

            // Populate fields specific to RecoveredPSMs
            psmDictionary[PsmTsvHeader_SpectralRecovery.IsolationWindowCenter] = 
                srPsm == null ? "" : srPsm.MsDataScan.IsolationMz.NullableToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.PrecursorOffset] =
                srPsm == null ? "" : CalculatePrecursorOffset(srPsm).NullableToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.IsolationWindowWidth] = 
                srPsm == null ? "" : srPsm.MsDataScan.IsolationMz.NullableToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.OriginalPsmQ] = 
                srPsm?.OriginalSpectralMatch == null ? "" : srPsm.OriginalSpectralMatch.FdrInfo.QValue.ToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.OriginalPsmPEP] =
                srPsm?.OriginalSpectralMatch == null ? "" : srPsm.OriginalSpectralMatch.FdrInfo.PEP.ToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.OriginalPsmPEP_QValue] =
                srPsm?.OriginalSpectralMatch == null ? "" : srPsm.OriginalSpectralMatch.FdrInfo.PEP_QValue.ToString(CultureInfo.InvariantCulture);
            return psmDictionary;
        }

        private static double? CalculatePrecursorOffset(SpectralRecoveryPSM srPsm)
        {
            if (srPsm?.MsDataScan.IsolationMz == null) return null;
            double precursorMz = Math.Abs(srPsm.ScanPrecursorMonoisotopicPeakMz) > 0.1
                ? srPsm.ScanPrecursorMonoisotopicPeakMz
                : srPsm.ClosestPrecursorPeak.Mz;
            return precursorMz - srPsm.MsDataScan.IsolationMz;
        }

    }
}
