using System.Collections.Generic;
using System.Linq;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;

namespace EngineLayer
{
    public class PeptideSpectralMatch : SpectralMatch
    {
        public PeptideSpectralMatch(IBioPolymerWithSetMods peptide, int notch, double score, int scanIndex,
            Ms2ScanWithSpecificMass scan, CommonParameters commonParameters,
            List<MatchedFragmentIon> matchedFragmentIons, double xcorr = 0) : base(peptide, notch, score, scanIndex,
            scan, commonParameters, matchedFragmentIons, xcorr)
        {
            _BestMatchingPeptides = new List<(int, PeptideWithSetModifications)>();
            ScanIndex = scanIndex;
            FullFilePath = scan.FullFilePath;
            ScanNumber = scan.OneBasedScanNumber;
            PrecursorScanNumber = scan.OneBasedPrecursorScanNumber;
            ScanRetentionTime = scan.RetentionTime;
            ScanExperimentalPeaks = scan.NumPeaks;
            TotalIonCurrent = scan.TotalIonCurrent;
            ScanPrecursorCharge = scan.PrecursorCharge;
            ScanPrecursorMonoisotopicPeakMz = scan.PrecursorMonoisotopicPeakMz;
            ScanPrecursorMass = scan.PrecursorMass;
            DigestionParams = commonParameters.DigestionParams;
            PeptidesToMatchingFragments = new Dictionary<PeptideWithSetModifications, List<MatchedFragmentIon>>();
            Xcorr = xcorr;
            NativeId = scan.NativeId;
            RunnerUpScore = commonParameters.ScoreCutoff;
            MsDataScan = scan.TheScan;
            SpectralAngle = -1;

        }
        #region Silac
            

        /// <summary>
        /// This method changes the base and full sequences to reflect heavy silac labels
        /// translates SILAC sequence into the proper peptide sequence ("PEPTIDEa" into "PEPTIDEK(+8.014)")
        /// </summary>
        public void ResolveHeavySilacLabel(List<SilacLabel> labels, IReadOnlyDictionary<string, int> modsToWritePruned)
        {
            //FullSequence
            FullSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.FullSequence)).ResolvedString; //string, not value
            FullSequence = SilacConversions.GetAmbiguousLightSequence(FullSequence, labels, false);

            //BaseSequence
            BaseSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.BaseSequence)).ResolvedString; //string, not value
            BaseSequence = SilacConversions.GetAmbiguousLightSequence(BaseSequence, labels, true);

            //EssentialSequence
            EssentialSequence = PsmTsvWriter.Resolve(_BestMatchingBioPolymersWithSetMods.Select(b => b.Pwsm.EssentialSequence(modsToWritePruned))).ResolvedString; //string, not value
            EssentialSequence = SilacConversions.GetAmbiguousLightSequence(EssentialSequence, labels, false);
        }

        /// <summary>
        /// This method is used by SILAC quantification to add heavy/light psms
        /// Don't have access to the scans at that point, so a new contructor is needed
        /// </summary>
        public PeptideSpectralMatch Clone(List<(int Notch, IBioPolymerWithSetMods Peptide)> bestMatchingPeptides) => new PeptideSpectralMatch(this, bestMatchingPeptides);
        
        protected PeptideSpectralMatch(SpectralMatch psm, List<(int Notch, IBioPolymerWithSetMods Peptide)> bestMatchingPeptides) 
            : base(psm, bestMatchingPeptides)
        {
        }
        #endregion


    }
}
