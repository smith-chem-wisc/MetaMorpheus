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
        public bool DeconvolutablePrecursor { get; }
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
            DeconvolutablePrecursor = true;

            //In cases where the precursor wasn't deconvoluted, the precursor mz and mass are set to 0 and 0 - z*1.008, respectively.
            // This updates the precursor mz and mass using the ClosestPrecursorPeak mz
            if (Math.Abs(ScanPrecursorMonoisotopicPeakMz) < 0.1)
            {
                DeconvolutablePrecursor = false;

                ScanPrecursorMonoisotopicPeakMz = ClosestPrecursorPeak.Mz;
                ScanPrecursorMass = ScanPrecursorMonoisotopicPeakMz.ToMass(AcceptorPeak.Apex.ChargeState);
            }
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
            Dictionary<string, string> psmDictionary = new();

            PsmTsvWriter.AddBasicMatchData(psmDictionary, srPsm);
            AddSpectralRecoveryData(psmDictionary, srPsm);
            PsmTsvWriter.AddPeptideSequenceData(psmDictionary, srPsm, modsToWritePruned);
            PsmTsvWriter.AddMatchedIonsData(psmDictionary, srPsm?.MatchedFragmentIons);
            PsmTsvWriter.AddMatchScoreData(psmDictionary, srPsm);

            return psmDictionary;
        }

        /// <summary>
        /// Populate fields specific to RecoveredPSMs
        /// </summary>
        /// <param name="psmDictionary"></param>
        /// <param name="srPsm"></param>
        /// <returns></returns>
        private static void AddSpectralRecoveryData(Dictionary<string, string> psmDictionary,
            SpectralRecoveryPSM srPsm)
        {
            // Chromatographic Peak Info
            psmDictionary[PsmTsvHeader_SpectralRecovery.PeakApexRt] = srPsm?.AcceptorPeak?.Apex == null 
                ? " " 
                : srPsm.AcceptorPeak.Apex.IndexedPeak.RetentionTime.ToString(CultureInfo.InvariantCulture);

            // Scan Isolation Window info
            psmDictionary[PsmTsvHeader_SpectralRecovery.PrecursorDeconvoluted] = srPsm == null 
                ? " " 
                : srPsm.DeconvolutablePrecursor ? "Y" : "N";
            psmDictionary[PsmTsvHeader_SpectralRecovery.IsolationWindowCenter] =
                srPsm == null ? " " : srPsm.MsDataScan.IsolationMz.NullableToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.PrecursorOffset] =
                srPsm == null ? " " : CalculatePrecursorOffset(srPsm).NullableToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.IsolationWindowWidth] =
                srPsm == null ? " " : srPsm.MsDataScan.IsolationWidth.NullableToString(CultureInfo.InvariantCulture);

            // Original Psm Info
            psmDictionary[PsmTsvHeader_SpectralRecovery.OriginalPsmQ] =
                srPsm?.OriginalSpectralMatch == null ? " " : srPsm.OriginalSpectralMatch.FdrInfo.QValue.ToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.OriginalPsmPEP] =
                srPsm?.OriginalSpectralMatch == null ? " " : srPsm.OriginalSpectralMatch.FdrInfo.PEP.ToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.OriginalPsmPEP_QValue] =
                srPsm?.OriginalSpectralMatch == null ? " " : srPsm.OriginalSpectralMatch.FdrInfo.PEP_QValue.ToString(CultureInfo.InvariantCulture);
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
